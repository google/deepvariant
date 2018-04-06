/*
 * Copyright 2017 Google Inc.
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

#include "deepvariant/pileup_image_native.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "tensorflow/core/platform/logging.h"

using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::CigarUnit;
using std::vector;

using learning::genomics::deepvariant::DeepVariantCall;

namespace learning {
namespace genomics {
namespace deepvariant {

using tensorflow::uint8;

namespace {

// Does this read support one of the alternative alleles?
inline bool ReadSupportsAlt(const DeepVariantCall& dv_call,
                            const Read& read,
                            const std::vector<string>& alt_alleles) {
  string key = (read.fragment_name() + "/" +
                std::to_string(read.read_number()));
  for (const string& alt_allele : alt_alleles) {
    const auto& allele_support = dv_call.allele_support();
    const bool alt_allele_present_in_call =
        allele_support.find(alt_allele) != allele_support.cend();
    if (alt_allele_present_in_call) {
      const auto& supp_read_names = allele_support.at(alt_allele).read_names();
      for (const string& read_name : supp_read_names) {
        if (read_name == key) return true;
      }
    }
  }
  return false;
}

}  // namespace


// The maximum value a pixel can have as a float. We use the 254.0
// value as originally set in DeepVariant v1. This means our pixel
// values can go from 0 to 254. Which, when converted to an int,
// gives us 255 or 256 possible pixel values.
const float kMaxPixelValueAsFloat = 254.0;


ImageRow::ImageRow(int width)
    : base(width, 0),
      base_quality(width, 0),
      mapping_quality(width, 0),
      on_positive_strand(width, 0),
      supports_alt(width, 0),
      matches_ref(width, 0)
{}

int ImageRow::Width() const {
  CHECK(base.size() == base_quality.size() &&
        base.size() == mapping_quality.size() &&
        base.size() == on_positive_strand.size() &&
        base.size() == supports_alt.size() &&
        base.size() == matches_ref.size());
  return base.size();
}

PileupImageEncoderNative::PileupImageEncoderNative(
    const PileupImageOptions& options)
    : options_(options) {
    CHECK((options_.width() % 2 == 1) && options_.width() >= 3)
        << "Width must be odd; found " << options_.width();
}

// Gets the pixel color (int) for a base.
int PileupImageEncoderNative::BaseColor(char base) const {
  switch (base) {
    case 'A': return (options_.base_color_offset_a_and_g() +
                      options_.base_color_stride() * 3);
    case 'G': return (options_.base_color_offset_a_and_g() +
                      options_.base_color_stride() * 2);
    case 'T': return (options_.base_color_offset_t_and_c() +
                      options_.base_color_stride() * 1);
    case 'C': return (options_.base_color_offset_t_and_c() +
                      options_.base_color_stride() * 0);
    default: return 0;
  }
}


int PileupImageEncoderNative::BaseColor(const string& base) const {
  CHECK_EQ(base.size(), 1) << "'base' string should be a single character";
  return BaseColor(base[0]);
}

int PileupImageEncoderNative::MatchesRefColor(bool base_matches_ref) const {
  float alpha = (base_matches_ref ?
                 options_.reference_matching_read_alpha() :
                 options_.reference_mismatching_read_alpha());
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

int PileupImageEncoderNative::SupportsAltColor(bool read_supports_alt) const {
  float alpha = (read_supports_alt ?
                 options_.allele_supporting_read_alpha() :
                 options_.allele_unsupporting_read_alpha());
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

int PileupImageEncoderNative::BaseQualityColor(int base_qual) const {
  float capped = static_cast<float>(
      std::min(options_.base_quality_cap(), base_qual));
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (capped / options_.base_quality_cap()));
}

int PileupImageEncoderNative::MappingQualityColor(int mapping_qual) const {
  float capped = static_cast<float>(
      std::min(options_.mapping_quality_cap(), mapping_qual));
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (capped / options_.mapping_quality_cap()));
}

int PileupImageEncoderNative::StrandColor(bool on_positive_strand) const {
  return (on_positive_strand ?
          options_.positive_strand_color() :
          options_.negative_strand_color());
}

std::unique_ptr<ImageRow>
PileupImageEncoderNative::EncodeRead(const DeepVariantCall& dv_call,
                                     const string& ref_bases,
                                     const Read& read,
                                     int image_start_pos,
                                     const vector<string>& alt_alleles) {
  ImageRow img_row(ref_bases.size());
  const bool supports_alt = ReadSupportsAlt(dv_call, read, alt_alleles);
  const int mapping_quality = read.alignment().mapping_quality();
  const bool is_forward_strand = !read.alignment().position().reverse_strand();
  const uint8 alt_color = SupportsAltColor(supports_alt);
  const uint8 mapping_color = MappingQualityColor(mapping_quality);
  const uint8 strand_color = StrandColor(is_forward_strand);
  const int min_base_quality = options_.read_requirements().min_base_quality();

  // Handler for each component of the CIGAR string, as subdivided
  // according the rules below.
  // Side effect: draws in img_row
  // Return value: true on normal exit; false if we determine that we
  // have a low quality base at the call position (in which case we
  // should return null) from EncodeRead.
  std::function<bool(int, int, const CigarUnit::Operation&)>
  action_per_cigar_unit = [&](int ref_i,
                              int read_i,
                              const CigarUnit::Operation& cigar_op) {
    char read_base = 0;
    if (cigar_op == CigarUnit::INSERT) {
      // redacted
      read_base = options_.indel_anchoring_base_char()[0];
    } else if (cigar_op == CigarUnit::DELETE) {
      ref_i -= 1;  // Adjust anchor base on reference
      read_base = options_.indel_anchoring_base_char()[0];
    } else if (cigar_op == CigarUnit::ALIGNMENT_MATCH ||
               cigar_op == CigarUnit::SEQUENCE_MATCH  ||
               cigar_op == CigarUnit::SEQUENCE_MISMATCH) {
      read_base = read.aligned_sequence()[read_i];
    }

    size_t col = ref_i - image_start_pos;
    if (read_base && 0 <= col && col < ref_bases.size()) {
      int base_quality = read.aligned_quality(read_i);
      int qual = std::min(base_quality, mapping_quality);
      if (ref_i == dv_call.variant().start() && qual < min_base_quality) {
        return false;
      }
      bool matches_ref = (read_base == ref_bases[col]);

      // Draw the pixel
      img_row.base[col]               = BaseColor(read_base);
      img_row.base_quality[col]       = BaseQualityColor(base_quality);
      img_row.mapping_quality[col]    = mapping_color;
      img_row.on_positive_strand[col] = strand_color;
      img_row.supports_alt[col]       = alt_color;
      img_row.matches_ref[col]        = MatchesRefColor(matches_ref);
    }
    return true;
  };

  // In the following, we iterate over alignment information for each
  // base of read, invoking action_per_cigar_unit for every segment of
  // the alignment.
  //
  // The handling of each cigar element type is given below, assuming
  // it has length n.
  //
  // ALIGNMENT_MATCH, SEQUENCE_MATCH, SEQUENCE_MISMATCH:
  //   Provide a segment ref_i, read_i for each of the n bases in the
  //   operator, where ref_i is the position on the genome where this
  //   base aligns.
  //
  // INSERT, CLIP_SOFT:
  //   Provides a single ref_i, read_i segment regardless of n. ref_i
  //   is set to the preceding base of the insertion; i.e., the anchor
  //   base. Beware that ref_i could be -1 if the insertion is aligned
  //   to the first base of a contig.  read_i points to the first base
  //   of the insertion. So if our cigar is 1M2I1M for a read starting
  //   at S, we'd see first (S, 0, '1M'), followed by one (S, 1,
  //   '2I'), and then (S + 1, 3, '1M').
  //
  // DELETE, SKIP:
  //   Provides a single ref_i, read_i segment regardless of n. ref_i
  //   is set to the first base of the deletion, just like in an
  //   ALIGNMENT_MATCH. read_i points to the previous base in the
  //   read, as there's no actual read sequence associated with a
  //   deletion. Beware that read_i could be -1 if the deletion is the
  //   first cigar of the read. So if our cigar is 1M2D1M for a read
  //   starting at S, we'd see first (S, 0, '1M'), followed by one (S
  //   + 1, 0, '2D'), and then (S + 3, 1, '1M').
  //
  // CLIP_HARD, PAD:
  //   These operators are ignored by as they don't impact the
  //   alignment of the read w.r.t. the reference.
  //
  // Any other CIGAR op:
  //   Fatal error, at present; later we should fail with a status encoding.
  int ref_i = read.alignment().position().position();
  int read_i = 0;
  bool ok = true;

  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();

    switch (op) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::SEQUENCE_MISMATCH:
        // Alignment op.
        for (int i = 0; i < op_len; i++) {
          ok = ok && action_per_cigar_unit(ref_i, read_i, op);
          ref_i++;
          read_i++;
        }
        break;
      case CigarUnit::INSERT:
      case CigarUnit::CLIP_SOFT:
        // Insert op.
        ok = action_per_cigar_unit(ref_i - 1, read_i, op);
        read_i += op_len;
        break;
      case CigarUnit::DELETE:
      case CigarUnit::SKIP:
        // Delete op.
        ok = action_per_cigar_unit(ref_i, read_i - 1, op);
        ref_i += op_len;
        break;
      case CigarUnit::CLIP_HARD:
      case CigarUnit::PAD:
        // Ignored ops.  Do nothing.
        break;
      default:
        LOG(FATAL) << "Unrecognized CIGAR op";
    }

    // Bail out if we found this read had a low-quality base at the
    // call site.
    if (!ok) {
      return nullptr;
    }
  }

  return std::unique_ptr<ImageRow>(new ImageRow(img_row));
}


std::unique_ptr<ImageRow>
PileupImageEncoderNative::EncodeReference(const string& ref_bases) {
  int ref_qual = options_.reference_base_quality();
  uint8 base_quality_color = BaseQualityColor(ref_qual);
  uint8 mapping_quality_color = MappingQualityColor(ref_qual);
  // We use "+" strand color for the reference.
  uint8 strand_color = StrandColor(true);
  uint8 alt_color = SupportsAltColor(false);
  uint8 ref_color = MatchesRefColor(true);

  ImageRow img_row(ref_bases.size());
  for (size_t i = 0; i < ref_bases.size(); ++i) {
    img_row.base[i]               = BaseColor(ref_bases[i]);
    img_row.base_quality[i]       = base_quality_color;
    img_row.mapping_quality[i]    = mapping_quality_color;
    img_row.on_positive_strand[i] = strand_color;
    img_row.supports_alt[i]       = alt_color;
    img_row.matches_ref[i]        = ref_color;
  }

  return std::unique_ptr<ImageRow>(new ImageRow(img_row));
}


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
