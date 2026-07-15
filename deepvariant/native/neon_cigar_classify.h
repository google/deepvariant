// neon_cigar_classify.h — NEON byte-level classifier for the M-block
// inner loop of AlleleCounter::Add (A2.2).
//
// For an ALIGNMENT_MATCH / SEQUENCE_MATCH / SEQUENCE_MISMATCH CIGAR
// element of length `n`, upstream's per-base inner loop
// (`deepvariant/allelecounter.cc:902-942`) does, for each base offset
// i in [0, n):
//
//   1. Quality check (`CanBasesBeUsed`):
//        canonical = read[i] ∈ {A,C,G,T}
//        if legacy:  use_base = canonical && qual[i] >= min_quality
//        else:       use_base = canonical
//                    is_low_quality_i = (qual[i] < min_quality)
//   2. Type:           is_ref = (ref[i] == read[i])
//   3. If `IsValidRefOffset && use_base`, emit a ReadAllele.
//
// Steps 1-2 are pure byte-level comparisons over three contiguous
// arrays (read[], ref[], qual[]) — perfect for NEON. The actual
// emit (step 3) is upstream scalar code that can read the per-base
// bitmasks produced here.
//
// This kernel produces four uint8 output arrays of length n:
//
//   use_base[i]        — 1 if base passes canonical+qual gates
//   is_low_quality[i]  — 1 if base passes canonical but is low-qual
//                        (only meaningful when !legacy)
//   is_ref[i]          — 1 if read[i] == ref[i]
//   canonical[i]       — 1 if read[i] ∈ {A,C,G,T} (debug)
//
// The kernel does NOT touch `to_add`, `methylation`, or any
// upstream bookkeeping. It is a *pre-classification* pass that lets
// the outer scalar code skip per-base function calls.
//
// Production wiring: this kernel is live in the allele-counting path —
// `deepvariant/allelecounter.cc` includes this header and calls
// `ClassifyMBlockNeon` from the M-block inner loop of AlleleCounter::Add.

#pragma once

#include <cstddef>
#include <cstdint>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#  include <arm_neon.h>
#  define DV_NEON_CIGAR_AVAILABLE 1
#else
#  define DV_NEON_CIGAR_AVAILABLE 0
#endif

namespace deepvariant {
namespace neon_cigar {

struct ClassifyMasks {
  uint8_t* use_base;        // length n: 1 or 0
  uint8_t* is_low_quality;  // length n: 1 or 0  (only set when !legacy)
  uint8_t* is_ref;          // length n: 1 or 0
  uint8_t* canonical;       // length n: 1 or 0
};

// Reference scalar implementation. Bit-exactly matches upstream's
// `CanBasesBeUsed(len=1)` semantics, used by both `legacy` and
// non-legacy modes.
inline void ClassifyMBlockScalar(const char* read, const char* ref,
                                 const uint8_t* qual, size_t n,
                                 uint8_t min_quality, bool legacy,
                                 const ClassifyMasks& out) {
  for (size_t i = 0; i < n; ++i) {
    const uint8_t b = static_cast<uint8_t>(read[i]);
    const uint8_t r = static_cast<uint8_t>(ref[i]);
    const uint8_t q = qual[i];

    // canonical = b ∈ {A,C,G,T} (uppercase only — matches
    // CanonicalBases::ACGT default).
    const uint8_t can =
        (b == 'A' || b == 'C' || b == 'G' || b == 'T') ? 1u : 0u;
    out.canonical[i] = can;

    if (!can) {
      out.use_base[i] = 0;
      out.is_low_quality[i] = 0;
      out.is_ref[i] = 0;
      continue;
    }

    if (legacy) {
      out.use_base[i] = (q >= min_quality) ? 1u : 0u;
      out.is_low_quality[i] = 0;
    } else {
      out.use_base[i] = 1u;
      out.is_low_quality[i] = (q < min_quality) ? 1u : 0u;
    }
    out.is_ref[i] = (r == b) ? 1u : 0u;
  }
}

// NEON 16-byte chunk-fill path. Falls through to scalar tail.
//
// vceqq_u8/vcgeq_u8 produce 0xFF/0x00 masks; we right-shift by 7 to
// turn them into 0x01/0x00 so downstream code can OR/AND them as
// natural 0/1 booleans (matching the scalar reference layout).
inline void ClassifyMBlockNeon(const char* read, const char* ref,
                               const uint8_t* qual, size_t n,
                               uint8_t min_quality, bool legacy,
                               const ClassifyMasks& out) {
#if DV_NEON_CIGAR_AVAILABLE
  size_t i = 0;
  if (n >= 16) {
    const uint8x16_t v_a = vdupq_n_u8('A');
    const uint8x16_t v_c = vdupq_n_u8('C');
    const uint8x16_t v_g = vdupq_n_u8('G');
    const uint8x16_t v_t = vdupq_n_u8('T');
    const uint8x16_t v_minq = vdupq_n_u8(min_quality);
    const uint8x16_t v_one = vdupq_n_u8(1);

    for (; i + 16 <= n; i += 16) {
      uint8x16_t b = vld1q_u8(reinterpret_cast<const uint8_t*>(read + i));
      uint8x16_t r = vld1q_u8(reinterpret_cast<const uint8_t*>(ref + i));
      uint8x16_t q = vld1q_u8(qual + i);

      // canonical = b == any of {A,C,G,T}
      uint8x16_t is_a = vceqq_u8(b, v_a);
      uint8x16_t is_c = vceqq_u8(b, v_c);
      uint8x16_t is_g = vceqq_u8(b, v_g);
      uint8x16_t is_t = vceqq_u8(b, v_t);
      uint8x16_t can_mask = vorrq_u8(vorrq_u8(is_a, is_c),
                                     vorrq_u8(is_g, is_t));
      // canonical → 0/1
      uint8x16_t can = vandq_u8(can_mask, v_one);
      vst1q_u8(out.canonical + i, can);

      // is_ref = (r == b) AND canonical (so non-canonical → 0).
      uint8x16_t eq_mask = vceqq_u8(b, r);
      uint8x16_t is_ref_mask = vandq_u8(eq_mask, can_mask);
      vst1q_u8(out.is_ref + i, vandq_u8(is_ref_mask, v_one));

      uint8x16_t qual_ok_mask = vcgeq_u8(q, v_minq);     // 0xFF if qual >= min
      uint8x16_t qual_low_mask = vmvnq_u8(qual_ok_mask); // 0xFF if qual < min

      uint8x16_t use_mask;
      uint8x16_t low_mask;
      if (legacy) {
        // legacy: emit only if canonical AND qual ok
        use_mask = vandq_u8(can_mask, qual_ok_mask);
        low_mask = vdupq_n_u8(0);
      } else {
        // non-legacy: emit if canonical (always); low-quality flag
        // tracks the slow path.
        use_mask = can_mask;
        low_mask = vandq_u8(can_mask, qual_low_mask);
      }
      vst1q_u8(out.use_base + i, vandq_u8(use_mask, v_one));
      vst1q_u8(out.is_low_quality + i, vandq_u8(low_mask, v_one));
    }
  }

  // Scalar tail.
  if (i < n) {
    ClassifyMasks tail{
        out.use_base + i,
        out.is_low_quality + i,
        out.is_ref + i,
        out.canonical + i,
    };
    ClassifyMBlockScalar(read + i, ref + i, qual + i, n - i, min_quality,
                         legacy, tail);
  }
#else
  ClassifyMBlockScalar(read, ref, qual, n, min_quality, legacy, out);
#endif
}

}  // namespace neon_cigar
}  // namespace deepvariant
