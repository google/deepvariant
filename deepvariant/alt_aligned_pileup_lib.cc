/*
 * Copyright 2023 Google LLC.
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

#include "deepvariant/alt_aligned_pileup_lib.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/realigner/fast_pass_aligner.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace {

// Margin added to the reference sequence for the aligner module.
constexpr int64_t kRefAlignMargin = 20;

using ::nucleus::genomics::v1::CigarUnit;
using ::nucleus::genomics::v1::Read;
using ::nucleus::genomics::v1::Range;
using ::nucleus::genomics::v1::Variant;

bool IsOperationRefAdvancing(const CigarUnit& unit) {
  static const absl::flat_hash_set<int>* kRefAdvancingOps =
      new absl::flat_hash_set<int>(
          {nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MATCH,
           nucleus::genomics::v1::CigarUnit::DELETE,
           nucleus::genomics::v1::CigarUnit::SKIP,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MISMATCH});
  return kRefAdvancingOps->contains(unit.operation());
}
bool IsOperationReadAdvancing(const CigarUnit& unit) {
  static const absl::flat_hash_set<int>* kReadAdvancingOps =
      new absl::flat_hash_set<int>(
          {nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MATCH,
           nucleus::genomics::v1::CigarUnit::INSERT,
           nucleus::genomics::v1::CigarUnit::CLIP_SOFT,
           nucleus::genomics::v1::CigarUnit::SEQUENCE_MISMATCH});
  return kReadAdvancingOps->contains(unit.operation());
}
}  // namespace

void TrimCigar(const ::google::protobuf::RepeatedPtrField<CigarUnit>& cigar,
               int64_t ref_start, int64_t ref_length,
               ::google::protobuf::RepeatedPtrField<CigarUnit>* new_cigar,
               int64_t& read_start, int64_t& new_read_length) {
  // First consume the ref until the trim is covered.
  int64_t trim_remaining = ref_start;
  // Then consume the ref until the ref_length is covered.
  int64_t ref_to_cover_remaining = ref_length;
  read_start = 0;
  new_read_length = 0;
  for (const auto& cigar_unit : cigar) {
    int64_t c_operation_length = cigar_unit.operation_length();
    // Each operation moves forward in the ref, the read, or both.
    bool advances_ref = IsOperationRefAdvancing(cigar_unit);
    bool advances_read = IsOperationReadAdvancing(cigar_unit);
    int64_t ref_step = advances_ref ? c_operation_length : 0L;
    // First, use up each operation until the trimmed area is covered.
    if (trim_remaining > 0) {
      if (ref_step <= trim_remaining) {
        // Fully apply to the trim.
        trim_remaining -= ref_step;
        read_start += (advances_read ? c_operation_length : 0L);
        continue;
      } else {
        // Partially apply to finish the trim.
        ref_step -= trim_remaining;
        read_start += (advances_read ? trim_remaining : 0L);
        // If trim finishes here, the rest of the ref_step can apply to the
        // next stage and count towards covering the given ref window.
        c_operation_length = ref_step;
        trim_remaining = 0;
      }
    }
    // Once the trim is done, start applying cigar entries to covering the ref
    // window.
    if (trim_remaining == 0) {
      if (ref_step <= ref_to_cover_remaining) {
        // Fully apply to the window.
        // CigarUnit new_cigar_unit;
        auto new_cigar_unit = new_cigar->Add();
        new_cigar_unit->set_operation(cigar_unit.operation());
        new_cigar_unit->set_operation_length(c_operation_length);
        ref_to_cover_remaining -= ref_step;
        new_read_length += (advances_read ? c_operation_length : 0L);
      } else {
        // Partially apply to finish the window.
        c_operation_length = ref_to_cover_remaining;
        auto new_cigar_unit = new_cigar->Add();
        new_cigar_unit->set_operation(cigar_unit.operation());
        new_cigar_unit->set_operation_length(c_operation_length);
        new_read_length += (advances_read ? c_operation_length : 0L);
        ref_to_cover_remaining = 0;
        break;
      }
    }
  }
}

Read TrimRead(const Read& read, const Range& region) {
  int64_t read_start = read.alignment().position().position();
  // Ref position where trimmed read should start.
  int64_t trim_left = std::max(region.start() - read_start, 0L);
  // Ref length of the trimmed read.
  int64_t ref_length = region.end() - std::max(region.start(), read_start);
  CHECK_GT(ref_length, 0);

  // Calculated read relative position after trimming.
  int64_t read_trim = 0;
  // Calculated read length after trimming.
  int64_t new_read_length = 0;
  // Manually copy the read proto to avoid copying some large fields that we
  // don't  need.
  Read new_read;
  auto new_cigar = new_read.mutable_alignment()->mutable_cigar();
  TrimCigar(read.alignment().cigar(), trim_left, ref_length, new_cigar,
            read_trim, new_read_length);
  new_read.set_fragment_name(read.fragment_name());
  new_read.set_id(read.id());
  new_read.set_read_group_id(read.read_group_id());
  new_read.set_read_group_set_id(read.read_group_set_id());
  new_read.set_read_number(read.read_number());
  new_read.set_fragment_length(read.fragment_length());
  new_read.set_number_reads(read.number_reads());
  for (const auto& [each_info_key, value] : read.info()) {
    (*new_read.mutable_info())[each_info_key] = value;
  }
  *new_read.mutable_alignment()->mutable_position() =
      read.alignment().position();
  new_read.mutable_alignment()->set_mapping_quality(
      read.alignment().mapping_quality());
  // Following fields are not needed but we copy them for consistency:
  if (read.has_next_mate_position()) {
    *new_read.mutable_next_mate_position() = read.next_mate_position();
  }
  new_read.set_proper_placement(read.proper_placement());
  new_read.set_duplicate_fragment(read.duplicate_fragment());
  new_read.set_failed_vendor_quality_checks(
      read.failed_vendor_quality_checks());
  new_read.set_secondary_alignment(read.secondary_alignment());
  new_read.set_supplementary_alignment(read.supplementary_alignment());

  // Set new read's alignment position.
  if (trim_left != 0) {
    new_read.mutable_alignment()->mutable_position()->set_position(
        region.start());
  }
  // Set aligned_sequence. Aligned string is trimmed starting from read_trim
  // and ending at read_trim + new_read_length.
  CHECK(read_trim >= 0 &&
        read_trim + new_read_length <= read.aligned_sequence().size());
  new_read.set_aligned_sequence(
      read.aligned_sequence().substr(read_trim, new_read_length));
  // Set aligned_quality, a repeated integer:
  CHECK(read_trim >= 0 &&
        read_trim + new_read_length <= read.aligned_quality().size());
  int64_t pos = 0;
  for (const auto& quality : read.aligned_quality()) {
    if (pos >= read_trim) {
      new_read.mutable_aligned_quality()->Add(quality);
    }
    if (pos == read_trim + new_read_length - 1) break;
    pos++;
  }
  return new_read;
}

Range CalculateAlignmentRegion(const Variant& variant, int half_width,
                               const nucleus::GenomeReference& ref_reader) {
  // Alignment region has the same width as pileup image.
  Range alignment_region;
  int64_t ref_start = variant.start();
  int64_t n_ref_bases = variant.reference_bases().size();
  int64_t ref_end = ref_start + n_ref_bases;
  alignment_region.set_reference_name(variant.reference_name());
  alignment_region.set_start(std::max(variant.start() - half_width, 0L));
  alignment_region.set_end(std::min(
      ref_reader.Contig(variant.reference_name()).ValueOrDie()->n_bases(),
      ref_end + half_width));
  return alignment_region;
}

// Helper function to trim a vector of reads. For the general use case we don't
// want to keep reads that are less than 15 bases long after trimming.
// original_alignment_positions is a vector where read alignment starting
// positions before trimming are stored.
std::vector<Read> TrimReads(const std::vector<const Read*>& reads,
                            const Range& region,
                            std::vector<int64_t>& original_alignment_positions,
                            int min_overlap) {
  std::vector<Read> ret;
  for (const Read* read : reads) {
    Read trimmed_read;
    trimmed_read = TrimRead(*read, region);
    if (trimmed_read.aligned_sequence().size() >= min_overlap) {
      original_alignment_positions.push_back(
          read->alignment().position().position());
      ret.push_back(trimmed_read);
    }
  }
  return ret;
}

namespace {
Range MakeRange(const std::string& ref_name, int64_t start, int64_t end) {
  Range range;
  range.set_reference_name(ref_name);
  range.set_start(start);
  range.set_end(end);
  return range;
}
}  // namespace

std::vector<Read> RealignReadsToHaplotype(
    absl::string_view haplotype, const std::vector<Read>& reads,
    absl::string_view contig, int64_t ref_start, int64_t ref_end,
    const nucleus::GenomeReference& ref_reader,
    const MakeExamplesOptions& options) {
  FastPassAligner realigner;
  AlignerOptions aln_config = options.realigner_options().aln_config();
  if (!reads.empty()) {
    aln_config.set_read_size(reads[0].aligned_sequence().size());
  }
  aln_config.set_force_alignment(true);
  realigner.set_options(aln_config);
  // Both reference and haplotype are padded with typically 20 bases from the
  // reference.
  int64_t ref_start_ext = std::max(0L, ref_start - kRefAlignMargin);
  int64_t ref_end_ext =
      std::min(ref_reader.Contig(std::string(contig)).ValueOrDie()->n_bases(),
               ref_end + kRefAlignMargin);

  std::string ref_prefix =
      ref_reader.GetBases(MakeRange(std::string(contig),
                                    ref_start_ext, ref_start)).ValueOrDie();
  std::string ref_suffix = ref_reader.GetBases(
      MakeRange(std::string(contig), ref_end, ref_end_ext)).ValueOrDie();
  realigner.set_reference(absl::StrCat(ref_prefix, haplotype, ref_suffix));
  realigner.set_ref_start(std::string(contig), ref_start_ext);
  realigner.set_ref_prefix_len(ref_start - ref_start_ext);
  realigner.set_ref_suffix_len(ref_end_ext - ref_end);
  realigner.set_haplotypes({absl::StrCat(ref_prefix, haplotype, ref_suffix)});
  return *realigner.AlignReads(reads);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
