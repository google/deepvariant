// global_align_simd-inl.h — SIMD GlobalAlign DP  // NOLINT(build/header_guard)
/*
 * Copyright 2026 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

// This file is included multiple times (once per SIMD target) via
// foreach_target.h. It implements the Farrar striped Smith-Waterman-like
// DP used by FastPassAligner::GlobalAlign.
//
// The key difference from the prototype benchmark: this version stores
// ALL columns in striped layout (not just a rolling 2-column buffer),
// enabling BackTrackBestAlignmentStriped to read directly from the
// striped layout without costly O(n×m) conversion.

// clang-format off
#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "hwy/aligned_allocator.h"
#include "hwy/highway.h"
// clang-format on

HWY_BEFORE_NAMESPACE();
namespace learning {
namespace genomics {
namespace deepvariant {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

// PopulateDpMatrixSimd — Farrar striped SIMD DP (int16)
//
// Stores full (m+1) columns of M, E, F in striped layout.
// Each matrix occupies (m+1) * stripe_size int16_t values.
// Column j data is at offset j * stripe_size.
//
// Args:
//   query: query sequence (0-indexed, length n)
//   n: query length
//   target: target sequence (0-indexed, length m)
//   m: target length
//   match_score: positive match score (e.g., 4)
//   mismatch_pen: negative mismatch penalty (e.g., -6)
//   gap_open: negative gap opening penalty (e.g., -8)
//   gap_extend: negative gap extension penalty (e.g., -1)
//   edge_range: edge range for edge penalties
//   edge_indel_penalty: negative edge indel penalty
//   M_striped, E_striped, F_striped: output arrays, each (m+1)*stripe_size
//   M_row0, E_row0, F_row0: row 0 values, each (m+1) int16_t
//   stripe_size: S * W (precomputed)
//   S: ceil(n/W) (precomputed)
inline void PopulateDpMatrixSimd(const char* query, int n, const char* target,
                                 int m, int match_score, int mismatch_pen,
                                 int gap_open, int gap_extend, int edge_range,
                                 int edge_indel_penalty, int16_t* M_striped,
                                 int16_t* E_striped, int16_t* F_striped,
                                 int16_t* M_row0, int16_t* E_row0,
                                 int16_t* F_row0, size_t stripe_size, int S) {
  const hn::ScalableTag<int16_t> d;
  const size_t W = hn::Lanes(d);

  const int16_t kNegInf = -30000;

  // --- Initialize row 0 ---
  for (int j = 0; j <= m; ++j) {
    M_row0[j] = 0;
    int gap_val = gap_open + j * gap_extend;
    E_row0[j] = static_cast<int16_t>(
        std::max(gap_val, static_cast<int>(kNegInf)));
    F_row0[j] = static_cast<int16_t>(
        std::max(gap_val, static_cast<int>(kNegInf)));
  }

  // --- Helper: is_edge ---
  auto is_edge = [n, edge_range](int i) -> bool {
    return i > 0 && (i <= edge_range || i > n - edge_range);
  };

  // --- Precompute edge penalties in striped layout ---
  auto a_edge_open = hwy::AllocateAligned<int16_t>(stripe_size);
  auto a_edge_extend = hwy::AllocateAligned<int16_t>(stripe_size);
  int16_t* edge_open = a_edge_open.get();
  int16_t* edge_extend = a_edge_extend.get();

  // Score lookup tables: one per DNA base, in striped layout.
  auto a_lookupA = hwy::AllocateAligned<int16_t>(stripe_size);
  auto a_lookupC = hwy::AllocateAligned<int16_t>(stripe_size);
  auto a_lookupG = hwy::AllocateAligned<int16_t>(stripe_size);
  auto a_lookupT = hwy::AllocateAligned<int16_t>(stripe_size);
  int16_t* lookupA = a_lookupA.get();
  int16_t* lookupC = a_lookupC.get();
  int16_t* lookupG = a_lookupG.get();
  int16_t* lookupT = a_lookupT.get();

  for (int s = 0; s < S; ++s) {
    for (size_t l = 0; l < W; ++l) {
      int p = s + static_cast<int>(l) * S;  // 0-indexed query position
      int qi = p + 1;                       // 1-indexed query position
      size_t off = static_cast<size_t>(s) * W + l;
      if (p < n) {
        edge_open[off] = static_cast<int16_t>(
            is_edge(qi) ? edge_indel_penalty : 0);
        edge_extend[off] = static_cast<int16_t>(
            (is_edge(qi) && !is_edge(qi - 1)) ? edge_indel_penalty : 0);
        char qc = query[p];
        lookupA[off] =
            static_cast<int16_t>(qc == 'A' ? match_score : mismatch_pen);
        lookupC[off] =
            static_cast<int16_t>(qc == 'C' ? match_score : mismatch_pen);
        lookupG[off] =
            static_cast<int16_t>(qc == 'G' ? match_score : mismatch_pen);
        lookupT[off] =
            static_cast<int16_t>(qc == 'T' ? match_score : mismatch_pen);
      } else {
        // Padding beyond query length — use neutral values.
        edge_open[off] = 0;
        edge_extend[off] = 0;
        lookupA[off] = static_cast<int16_t>(mismatch_pen);
        lookupC[off] = static_cast<int16_t>(mismatch_pen);
        lookupG[off] = static_cast<int16_t>(mismatch_pen);
        lookupT[off] = static_cast<int16_t>(mismatch_pen);
      }
    }
  }

  // --- Initialize column 0 of striped storage ---
  // Column 0: M[i,0] = -kInf for i>0, E[i,0] = -kInf, F[i,0] = -kInf
  for (size_t k = 0; k < stripe_size; ++k) {
    M_striped[k] = kNegInf;
    E_striped[k] = kNegInf;
    F_striped[k] = kNegInf;
  }

  // --- Previous column buffer (starts as column 0 state) ---
  auto a_prevM = hwy::AllocateAligned<int16_t>(stripe_size);
  auto a_prevE = hwy::AllocateAligned<int16_t>(stripe_size);
  int16_t* prevM = a_prevM.get();
  int16_t* prevE = a_prevE.get();

  for (size_t k = 0; k < stripe_size; ++k) {
    prevM[k] = kNegInf;
    prevE[k] = kNegInf;
  }

  const auto vGapOpenExtend =
      hn::Set(d, static_cast<int16_t>(gap_open + gap_extend));
  const auto vGapExtend = hn::Set(d, static_cast<int16_t>(gap_extend));
  const auto vNegInf = hn::Set(d, kNegInf);
  const auto vMask0 = hn::FirstN(d, 1);

  // Dynamic lookup buffer for non-ACGT target characters.
  // Allocated once, reused per column when needed.
  auto a_lookupDyn = hwy::AllocateAligned<int16_t>(stripe_size);
  int16_t* lookupDyn = a_lookupDyn.get();

  // --- Process each column j = 1..m ---
  for (int j = 1; j <= m; ++j) {
    const size_t col_off = static_cast<size_t>(j) * stripe_size;

    // Select the score profile for the current target base.
    // For non-ACGT characters, build a dynamic profile matching the scalar
    // behavior (character equality: query[i-1] == target[j-1]).
    const int16_t* profile;
    char target_char = target[j - 1];
    switch (target_char) {
      case 'A': profile = lookupA; break;
      case 'C': profile = lookupC; break;
      case 'G': profile = lookupG; break;
      case 'T': profile = lookupT; break;
      default:
        // Non-ACGT: build dynamic lookup using character equality.
        for (int s = 0; s < S; ++s) {
          for (size_t l = 0; l < W; ++l) {
            int p = s + static_cast<int>(l) * S;
            size_t off = static_cast<size_t>(s) * W + l;
            if (p < n) {
              lookupDyn[off] = static_cast<int16_t>(
                  query[p] == target_char ? match_score : mismatch_pen);
            } else {
              lookupDyn[off] = static_cast<int16_t>(mismatch_pen);
            }
          }
        }
        profile = lookupDyn;
        break;
    }

    const int16_t M_row0_j = 0;
    const int16_t F_row0_j = static_cast<int16_t>(
        std::max(static_cast<int>(kNegInf), gap_open + j * gap_extend));
    const int16_t M_row0_prev = 0;

    // Pre-compute diagonal for segment 0.
    auto vLastPrevM = hn::Load(d, prevM + static_cast<size_t>(S - 1) * W);
    auto vDiag = hn::SlideUpLanes(d, vLastPrevM, 1);
    vDiag = hn::IfThenElse(vMask0, hn::Set(d, M_row0_prev), vDiag);

    // Initial "above" for F at segment 0.
    auto vM_above = hn::IfThenElse(vMask0, hn::Set(d, M_row0_j), vNegInf);
    auto vF_above = hn::IfThenElse(vMask0, hn::Set(d, F_row0_j), vNegInf);

    // ========== SINGLE FUSED PASS: E + H + F per segment ==========
    for (int s = 0; s < S; ++s) {
      const size_t off = static_cast<size_t>(s) * W;

      auto vPrevM = hn::Load(d, prevM + off);
      auto vPrevE = hn::Load(d, prevE + off);
      auto vEdgeOpen = hn::Load(d, edge_open + off);
      auto vEdgeExtend = hn::Load(d, edge_extend + off);
      auto vScore = hn::Load(d, profile + off);

      // E[i,j] = max(M[i,j-1] + gap_open + gap_extend + edge_penalty_open,
      //              E[i,j-1] + gap_extend + edge_penalty_extend)
      auto vE = hn::Max(hn::Add(hn::Add(vPrevM, vGapOpenExtend), vEdgeOpen),
                        hn::Add(hn::Add(vPrevE, vGapExtend), vEdgeExtend));

      // H[i,j] = max(M[i-1,j-1] + score, E[i,j])
      auto vH = hn::Max(hn::Add(vDiag, vScore), vE);

      // F[i,j] = max(M[i-1,j] + gap_open + gap_extend + edge_penalty_open,
      //              F[i-1,j] + gap_extend + edge_penalty_extend)
      auto vF = hn::Max(hn::Add(hn::Add(vM_above, vGapOpenExtend), vEdgeOpen),
                        hn::Add(hn::Add(vF_above, vGapExtend), vEdgeExtend));

      // M[i,j] = max(H, F)
      auto vM = hn::Max(vH, vF);

      // Store to column j in striped layout
      hn::Store(vE, d, E_striped + col_off + off);
      hn::Store(vF, d, F_striped + col_off + off);
      hn::Store(vM, d, M_striped + col_off + off);

      // Carry forward for next segment.
      vDiag = vPrevM;
      vF_above = vF;
      vM_above = vM;
    }

    // ========== LAZY-F CORRECTION ==========
    // The F matrix has a wraparound dependency: segment 0 depends on
    // segment S-1 from the same column. We resolve this iteratively.
    for (int iter = 0; iter < static_cast<int>(W); ++iter) {
      auto vLastM = hn::Load(d, M_striped + col_off +
                                    static_cast<size_t>(S - 1) * W);
      auto vLastF = hn::Load(d, F_striped + col_off +
                                    static_cast<size_t>(S - 1) * W);

      auto vM_sh = hn::SlideUpLanes(d, vLastM, 1);
      vM_sh = hn::IfThenElse(vMask0, hn::Set(d, M_row0_j), vM_sh);
      auto vF_sh = hn::SlideUpLanes(d, vLastF, 1);
      vF_sh = hn::IfThenElse(vMask0, hn::Set(d, F_row0_j), vF_sh);

      auto vF_new =
          hn::Max(hn::Add(hn::Add(vM_sh, vGapOpenExtend),
                         hn::Load(d, edge_open)),
                  hn::Add(hn::Add(vF_sh, vGapExtend),
                         hn::Load(d, edge_extend)));

      auto vF_old = hn::Load(d, F_striped + col_off);
      auto vF_upd = hn::Max(vF_old, vF_new);
      if (hn::AllFalse(d, hn::Ne(vF_upd, vF_old))) break;

      hn::Store(vF_upd, d, F_striped + col_off);
      hn::Store(hn::Max(hn::Load(d, M_striped + col_off), vF_upd), d,
                M_striped + col_off);

      bool any_changed = true;
      for (int s = 1; s < S && any_changed; ++s) {
        const size_t off = static_cast<size_t>(s) * W;
        const size_t poff = off - W;
        auto vFsn = hn::Max(
            hn::Add(hn::Add(hn::Load(d, M_striped + col_off + poff),
                            vGapOpenExtend),
                    hn::Load(d, edge_open + off)),
            hn::Add(hn::Add(hn::Load(d, F_striped + col_off + poff),
                            vGapExtend),
                    hn::Load(d, edge_extend + off)));
        auto vFso = hn::Load(d, F_striped + col_off + off);
        auto vFs = hn::Max(vFso, vFsn);
        any_changed = !hn::AllFalse(d, hn::Ne(vFs, vFso));
        hn::Store(vFs, d, F_striped + col_off + off);
        hn::Store(hn::Max(hn::Load(d, M_striped + col_off + off), vFs), d,
                  M_striped + col_off + off);
      }
    }

    // Update prev buffers for next column (read from stored striped data).
    for (size_t k = 0; k < stripe_size; ++k) {
      prevM[k] = M_striped[col_off + k];
      prevE[k] = E_striped[col_off + k];
    }
  }  // end column loop
}

// Returns the SIMD width (number of int16 lanes) for the best available target.
inline size_t SimdWidth() {
  const hn::ScalableTag<int16_t> d;
  return hn::Lanes(d);
}

}  // namespace HWY_NAMESPACE
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
HWY_AFTER_NAMESPACE();
