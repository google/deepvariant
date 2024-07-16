/*
 * Copyright 2017 Google LLC.
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

// Implementation of allelecounter.h.
#include "deepvariant/allelecounter.h"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include "absl/log/log.h"
#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Separator string that will appear between the fragment name and read number
// the string key constructed from a Read with ReadKey().
static constexpr char kFragmentNameReadNumberSeparator[] = "/";

using absl::string_view;

using absl::StrCat;
using nucleus::GenomeReference;
using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;

// TODO Consolidate SumAlleleCounts functions into one since the
// functionality is identical.
std::vector<Allele> SumAlleleCounts(const AlleleCount& allele_count,
                                    bool include_low_quality) {
  std::map<std::pair<string_view, AlleleType>, int> allele_sums;
  for (const auto& entry : allele_count.read_alleles()) {
    if (include_low_quality || !entry.second.is_low_quality()) {
      ++allele_sums[{entry.second.bases(), entry.second.type()}];
    }
  }

  std::vector<Allele> to_return;
  to_return.reserve(allele_sums.size());
  for (const auto& entry : allele_sums) {
    to_return.push_back(
        MakeAllele(entry.first.first, entry.first.second, entry.second));
  }

  // TODO SumAlleleCounts is only used in one place in variant_calling.cc
  // where ref alleles are filtered out. The code below is redundant.
  // Verify that there are no other usages of ref alleles and remove this code.
  //
  // Creates a synthetic reference Allele if we saw any reference containing
  // alleles, whose count is tracked (for performance reasons) as an integer
  // in the AlleleCount.ref_supporting_read_count field of the proto. This
  // synthetic allele allows us to provide the same API from this function: a
  // vector of the Alleles observed in allele_count without having to track the
  // read names for reference containing reads, which is very memory-intensive.
  if (allele_count.ref_supporting_read_count() > 0 &&
      !allele_count.track_ref_reads()) {
    to_return.push_back(MakeAllele(allele_count.ref_base(),
                                   AlleleType::REFERENCE,
                                   allele_count.ref_supporting_read_count()));
  }

  return to_return;
}

std::vector<Allele> SumAlleleCounts(absl::Span<const AlleleCount> allele_counts,
                                    bool include_low_quality) {
  std::map<std::pair<string_view, AlleleType>, int> allele_sums;
  for (const AlleleCount& allele_count : allele_counts) {
    for (const auto& entry : allele_count.read_alleles()) {
      if (include_low_quality || !entry.second.is_low_quality()) {
        ++allele_sums[{entry.second.bases(), entry.second.type()}];
      }
    }
  }

  std::vector<Allele> to_return;
  to_return.reserve(allele_sums.size());
  for (const auto& entry : allele_sums) {
    to_return.push_back(
        MakeAllele(entry.first.first, entry.first.second, entry.second));
  }

  // TODO SumAlleleCounts is only used in one place in variant_calling.cc
  // where ref alleles are filtered out. The code below is redundant.
  // Verify that there are no other usages of ref alleles and remove this code.
  //
  // Creates a synthetic reference Allele if we saw any reference containing
  // alleles, whose count is tracked (for performance reasons) as an integer
  // in the AlleleCount.ref_supporting_read_count field of the proto. This
  // synthetic allele allows us to provide the same API from this function: a
  // vector of the Alleles observed in allele_count without having to track the
  // read names for reference containing reads, which is very memory-intensive.
  int ref_support_for_all_samples = 0;
  for (const AlleleCount& allele_count : allele_counts) {
    ref_support_for_all_samples += allele_count.ref_supporting_read_count();
  }
  if (ref_support_for_all_samples > 0 && !allele_counts.empty() &&
      !allele_counts[0].track_ref_reads()) {
    to_return.push_back(MakeAllele(allele_counts[0].ref_base(),
                                   AlleleType::REFERENCE,
                                   ref_support_for_all_samples));
  }

  return to_return;
}

// TODO Consolidate TotalAlleleCounts functions into one since the
// functionality is identical.
// Allele counter tracks reads supporting alt alleles. Simple counter is used
// for ref supporting reads. If track_ref_reads flag is set then ref supporting
// reads are tracked as well but only for positions marked as potential
// candidates.
int TotalAlleleCounts(const AlleleCount& allele_count,
                      bool include_low_quality) {
  int total_allele_counts = std::count_if(
      allele_count.read_alleles().begin(), allele_count.read_alleles().end(),
      [include_low_quality](google::protobuf::Map<string, Allele>::value_type e) {
        return (!e.second.is_low_quality() || include_low_quality) &&
               e.second.type() != AlleleType::REFERENCE;
      });
  total_allele_counts += allele_count.ref_supporting_read_count();
  return total_allele_counts;
}

// Allele counter tracks reads supporting alt alleles. Simple counter is used
// for ref supporting reads. If track_ref_reads flag is set then ref supporting
// reads are tracked as well but only for positions marked as potential
// candidates.
int TotalAlleleCounts(const std::vector<AlleleCount>& allele_counts,
                      bool include_low_quality) {
  int total_allele_count = 0;
  for (const AlleleCount& allele_count : allele_counts) {
    total_allele_count += std::count_if(
        allele_count.read_alleles().begin(), allele_count.read_alleles().end(),
        [include_low_quality](google::protobuf::Map<string, Allele>::value_type e) {
          return (!e.second.is_low_quality() || include_low_quality) &&
                 e.second.type() != AlleleType::REFERENCE;
        });
    total_allele_count += allele_count.ref_supporting_read_count();
  }
  return total_allele_count;
}

// Returns false if any of the bases from offset to offset+len are canonical
// bases.
// If `keep_legacy_behavior` is set to true, this function will also return
// false when any of the bases in read from offset to offset + len is below
// the quality threshold.
//
// There is a separate bool output `is_low_quality`, which will be set to
// true if all the bases in read from offset to offset+len is lower than
// the quality threshold to be used for generating alleles for our counts.
bool CanBasesBeUsed(const nucleus::genomics::v1::Read& read, int offset,
                    int len, const AlleleCounterOptions& options,
                    bool& is_low_quality) {
  CHECK_LE(offset + len, read.aligned_quality_size());

  const int min_base_quality = options.read_requirements().min_base_quality();
  int indel_base_quality = 0;
  for (int i = 0; i < len; i++) {
    indel_base_quality += read.aligned_quality(offset + i);
    if (read.aligned_quality(offset + i) < min_base_quality &&
        options.keep_legacy_behavior()) {
      return false;
    }
    if (!nucleus::IsCanonicalBase(read.aligned_sequence()[offset + i])) {
      return false;
    }
  }
  is_low_quality = false;
  if (!options.keep_legacy_behavior()) {
    if (indel_base_quality < min_base_quality * len) {
      is_low_quality = true;
    }
  }
  return true;
}

bool allele_pos_cmp(const AlleleCount& allele_count, int64_t pos) {
  return allele_count.position().position() < pos;
}

// Return the allele index by base position in allele_counts vector.
int AlleleIndex(absl::Span<const AlleleCount> allele_counts, int64_t pos) {
  auto idx = std::lower_bound(allele_counts.begin(), allele_counts.end(), pos,
                              allele_pos_cmp);
  if (idx == allele_counts.end() || idx->position().position() != pos) {
    return -1;
  }
  return std::distance(allele_counts.begin(), idx);
}

void AlleleCounter::Init() {
  // Initialize our counts vector of AlleleCounts with proper position and
  // reference base information. Initially the alleles repeated field is empty.
  const int64_t len = IntervalLength();
  counts_.reserve(len);
  // Set candidate positions relative to the interval.
  for (auto& candidate_position : candidate_positions_) {
    candidate_position -= interval_.start();
  }
  auto full_interval_offset = interval_.start() - reads_interval_.start();
  for (int i = 0; i < len; ++i) {
    AlleleCount allele_count;
    const int64_t pos = interval_.start() + i;
    *(allele_count.mutable_position()) =
        nucleus::MakePosition(interval_.reference_name(), pos);
    allele_count.set_ref_base(ref_bases_.substr(i + full_interval_offset, 1));
    allele_count.set_track_ref_reads(options_.track_ref_reads());
    counts_.push_back(allele_count);
  }
}

// AlleleCounter objects are passed to Python by pointers. We need to return
// a raw pointer here in order to test a Python specific API.
AlleleCounter* AlleleCounter::InitFromAlleleCounts(
    absl::Span<const AlleleCount> allele_counts) {
  auto allele_counter = new AlleleCounter();
  allele_counter->counts_.assign(allele_counts.begin(), allele_counts.end());
  return allele_counter;
}

// This constructor is only used for unit testing, therefore it is defined as
// private.
AlleleCounter::AlleleCounter() : ref_(nullptr) {}

AlleleCounter::AlleleCounter(const GenomeReference* const ref,
                             const Range& range,
                             const nucleus::genomics::v1::Range& full_range,
                             const std::vector<int>& candidate_positions,
                             const AlleleCounterOptions& options)
    : ref_(ref),
      interval_(range),
      reads_interval_(full_range),
      candidate_positions_(candidate_positions),
      options_(options),
      ref_bases_(ref_->GetBases(full_range).ValueOrDie()) {
  Init();
}

AlleleCounter::AlleleCounter(const GenomeReference* const ref,
                             const Range& range,
                             const std::vector<int>& candidate_positions,
                             const AlleleCounterOptions& options)
    : ref_(ref),
      interval_(range),
      reads_interval_(range),
      candidate_positions_(candidate_positions),
      options_(options),
      ref_bases_(ref_->GetBases(range).ValueOrDie()) {
  Init();
}

string AlleleCounter::RefBases(const int64_t rel_start, const int64_t len) {
  CHECK_GT(len, 0) << "Length must be >= 1";

  // If our region isn't valid (e.g., it is off the end of the chromosome),
  // return an empty string, otherwise get the actual bases from reference.
  const int abs_start = reads_interval_.start() + rel_start;
  const Range region = nucleus::MakeRange(reads_interval_.reference_name(),
                                          abs_start, abs_start + len);
  if (!ref_->IsValidInterval(region)) {
    return "";
  } else {
    return ref_->GetBases(region).ValueOrDie();
  }
}

string AlleleCounter::GetPrevBase(const Read& read, const int read_offset,
                                  const int interval_offset) {
  CHECK_GE(read_offset, 0) << "read_offset should be 0 or greater";
  if (read_offset == 0) {
    // The read_offset case is here to handle the case where the insertion/
    // deletion/soft_clip is the first cigar element of the read, and there's no
    // previous base in the read, and so we take our previous base from the
    // reference genome instead.
    return RefBases(interval_offset - 1, 1);
  } else {
    // In all other cases we actually take our previous base from the read
    // itself.
    return read.aligned_sequence().substr(read_offset - 1, 1);
  }
}

ReadAllele AlleleCounter::MakeIndelReadAllele(const Read& read,
                                              const int interval_offset,
                                              const int ref_offset,
                                              const int read_offset,
                                              const CigarUnit& cigar) {
  const int op_len = cigar.operation_length();
  const string prev_base = GetPrevBase(read, read_offset, ref_offset);
  bool is_low_quality_read_allele = false;

  if (prev_base.empty() || !nucleus::AreCanonicalBases(prev_base) ||
      (cigar.operation() != CigarUnit::DELETE &&
       !CanBasesBeUsed(read, read_offset, op_len, options_,
                       is_low_quality_read_allele))) {
    // There is no prev_base (we are at the start of the contig), or the bases
    // are unusable, so don't actually add the indel allele.
    return ReadAllele();
  }

  AlleleType type;
  string bases;
  switch (cigar.operation()) {
    case CigarUnit::DELETE:
      type = AlleleType::DELETION;
      bases = RefBases(ref_offset, op_len);
      if (bases.empty()) {
        // We couldn't get the ref bases for the deletion (which can happen if
        // the deletion spans off the end of the contig), so abort now without
        // considering this read any longer. It's rare but such things happen in
        // genomes but they do occur in practice, such as when: (1) the reads
        // spans off the chromosome, but because there's more sequence there
        // (the chromosome isn't complete), which means the read can have
        // whatever CIGAR it likes, which may include a deletion; (2) the
        // chromosome is actually circular, and the aligner is clever enough to
        // know that, and the read's cigar reflect true differences of the read
        // to the alignment at the start of the contig.  Nasty, I know.
        VLOG(2) << "Deletion spans off the chromosome for read: "
                << read.ShortDebugString() << " at cigar "
                << cigar.ShortDebugString() << " with interval "
                << Interval().ShortDebugString() << " with interval_offset "
                << interval_offset << " and read_offset " << read_offset;
        return ReadAllele();
      }

      if (!nucleus::AreCanonicalBases(bases)) {
        // The reference genome has non-canonical bases that are being deleted.
        // We don't add deletions with non-canonical bases so we return an empty
        // ReadAllele().
        return ReadAllele();
      }

      break;
    case CigarUnit::INSERT:
      type = AlleleType::INSERTION;
      bases = read.aligned_sequence().substr(read_offset, op_len);
      break;
    case CigarUnit::CLIP_SOFT:
      type = AlleleType::SOFT_CLIP;
      bases = read.aligned_sequence().substr(read_offset, op_len);
      break;
    default:
      LOG(FATAL) << "Unexpected cigar operation: " << cigar.DebugString();
  }

  return ReadAllele(interval_offset - 1, StrCat(prev_base, bases), type,
                    is_low_quality_read_allele);
}

void AlleleCounter::AddReadAlleles(const Read& read, absl::string_view sample,
                                   const std::vector<ReadAllele>& to_add) {
  for (size_t i = 0; i < to_add.size(); ++i) {
    const ReadAllele& to_add_i = to_add[i];

    // The read can span beyond and after the interval, so don't add counts
    // outside our interval boundaries.
    if (to_add_i.skip() || !IsValidIntervalOffset(to_add_i.position())) {
      continue;
    }

    // If sequential alleles have the same position, skip the first one. This
    // occurs, for example, when we observe a base at position p on the genome
    // which is enqueued as the ith element of our to_add vector. But the next
    // allele is an indel allele which, because of VCF convention, occurs at
    // position p, is enqueued at i+1 and supersedes the previous base
    // substitution. Resolving these conflicts here allows us to keep the
    // Read => ReadAllele algorithm logic simple.
    if (i + 1 < to_add.size() &&
        to_add_i.position() == to_add[i + 1].position()) {
      continue;
    }

    AlleleCount& allele_count = counts_[to_add_i.position()];

    if (to_add_i.type() == AlleleType::REFERENCE) {
      if (!to_add_i.is_low_quality()) {
        const int prev_count = allele_count.ref_supporting_read_count();
        allele_count.set_ref_supporting_read_count(prev_count + 1);
      }
    }

    // Always create non reference alleles.
    // Reference alleles are created only when the track_ref_reads flag is set
    // and we know that this position contains a potential candidate.
    if (to_add_i.type() != AlleleType::REFERENCE ||
        (options_.track_ref_reads() &&
         std::binary_search(candidate_positions_.begin(),
                            candidate_positions_.end(), to_add_i.position()))) {
      auto* read_alleles = allele_count.mutable_read_alleles();
      auto* sample_alleles = allele_count.mutable_sample_alleles();
      const string key = ReadKey(read);
      const Allele allele = MakeAllele(to_add_i.bases(), to_add_i.type(), 1,
                                       to_add_i.is_low_quality());

      // Naively, there should never be multiple counts for the same read key.
      // We detect such a situation here but only write out a warning. It would
      // be better to have a stronger response (FATAL), but unfortunately we see
      // data in the wild that we need to process that has duplicates.
      if (read_alleles->count(key)) {
        // Not thread safe.
        static int counter = 0;
        if (counter++ < 1) {
          VLOG(2) << "Found duplicate read: " << key << " at "
                  << allele_count.position().ShortDebugString();
        }
      }

      (*read_alleles)[key] = allele;
      // Update sample to allele map. This may allows us to determine set of
      // samples that support each allele.
      Allele* new_allele = (*sample_alleles)[std::string(sample)].add_alleles();
      *new_allele = allele;
    }
  }
}

// Convenience function to check if operations is match. Note, that we treat
// SEQUENCE_MISMATCH as ALIGNMENT_MATCH.
bool IsOperationMatch(const nucleus::genomics::v1::CigarUnit& op) {
  return (op.operation() == CigarUnit::ALIGNMENT_MATCH ||
          op.operation() == CigarUnit::SEQUENCE_MATCH ||
          op.operation() == CigarUnit::SEQUENCE_MISMATCH);
}

// Merge two operations. If operations are the same type then first operation's
// length is icreased and second operation length's is set to zero. If
// operations are different types then M operation of length MIN(op1, op2) is
// added instead of op1. Op2 is converted to the type of a larger operation
// and length is set to Max(op) - Min(op).
// Function returns true if operations are merged.
bool MergeOperations(nucleus::genomics::v1::CigarUnit& op1,
                     nucleus::genomics::v1::CigarUnit& op2) {
  // Simple merge if operations are of the same type. There are three different
  // types of "match" operation therefore it is not enough to just compare
  // operation types.
  if (op1.operation() == op2.operation() ||
      (IsOperationMatch(op1) && IsOperationMatch(op2))) {
    op1.set_operation_length(op1.operation_length() + op2.operation_length());
    op2.set_operation_length(0);
  } else if ((op1.operation() == CigarUnit::DELETE ||
              op1.operation() == CigarUnit::INSERT) &&
             (op2.operation() == CigarUnit::DELETE ||
              op2.operation() == CigarUnit::INSERT)) {
    auto min_indel_len =
        std::min(op1.operation_length(), op2.operation_length());
    auto new_indel_len =
        std::max(op1.operation_length(), op2.operation_length()) -
        min_indel_len;
    if (op1.operation_length() > op2.operation_length()) {
      op2.set_operation(op1.operation());
    }
    op1.set_operation(CigarUnit::ALIGNMENT_MATCH);
    op1.set_operation_length(min_indel_len);
    op2.set_operation_length(new_indel_len);
  } else {
    return false;
  }
  return true;
}

// Advance reference and read pointers depending on the cigar operation and
// its length.
void AdvanceReadReferencePointers(
    int increment, const nucleus::genomics::v1::CigarUnit_Operation& operation,
    int& read_offset, int& ref_offset) {
  switch (operation) {
    case CigarUnit::ALIGNMENT_MATCH:
    case CigarUnit::SEQUENCE_MATCH:
    case CigarUnit::SEQUENCE_MISMATCH:
      read_offset += increment;
      ref_offset += increment;
      break;
    case CigarUnit::CLIP_SOFT:
    case CigarUnit::INSERT:
      read_offset += increment;
      // No interval offset change, since an insertion doesn't move us on ref.
      break;
    case CigarUnit::DELETE:
    case CigarUnit::PAD:
    case CigarUnit::SKIP:
      // No read offset change, since a del/pad/skip don't consume read bases.
      ref_offset += increment;
      break;
    case CigarUnit::CLIP_HARD:
      break;
    default:
      // Lots of misc. enumerated values from proto that aren't useful such as
      // enumeration values INT_MIN_SENTINEL_DO_NOT_USE_ and
      // OPERATION_UNSPECIFIED.
      break;
  }
}

// Handle the case when INDEL is at the head of a cigar.
// DEL is removed and alignment position is shifted to the right.
// INS is converted into a REF and alignment position is shifted to the left.
int HandleHeadingIndel(
    std::vector<nucleus::genomics::v1::CigarUnit>::iterator it,
    std::vector<nucleus::genomics::v1::CigarUnit>& norm_cigar) {
  int read_alignment_shift = 0;
  // it must be a first operation or the first op following soft clip.
  CHECK(it == norm_cigar.begin() ||
        (!norm_cigar.empty() &&
         norm_cigar.begin()->operation() == CigarUnit::CLIP_SOFT &&
         it == norm_cigar.begin() + 1));
  if (it->operation() == CigarUnit::DELETE) {
    read_alignment_shift = it->operation_length();
    norm_cigar.erase(it);
  } else if (it->operation() == CigarUnit::INSERT) {
    read_alignment_shift = -it->operation_length();
    it->set_operation(CigarUnit::ALIGNMENT_MATCH);
  }
  return read_alignment_shift;
}
// Shift cigar operation according to the shift parameter. Only INDELs are
// shifted and only to the left. It is expected that operation to the left is
// REF or SOFT_CLIP or there is no operation. Operation to the left is decreased
// in length, operation to the right is increased in length. If there is no
// operation to the right  REF is created.
int ShiftOperation(int shift,
                   std::vector<nucleus::genomics::v1::CigarUnit>::iterator it,
                   std::vector<nucleus::genomics::v1::CigarUnit>& norm_cigar) {
  // If previous operation is ref or soft clip then it is reduced in length.
  // If it is the first operation then read alignment is shifted. In this case
  // it is removed if it is del or turned into ref if it is ins.
  int read_alignment_shift = 0;
  // If indel is first operation or second after a soft clip then indel is
  // treated specially.
  if (it == norm_cigar.begin() ||
      (!norm_cigar.empty() && it == norm_cigar.begin() + 1 &&
       norm_cigar.begin()->operation() == CigarUnit::CLIP_SOFT)) {
    return HandleHeadingIndel(it, norm_cigar);
  } else {
    auto prev_op = it - 1;
    // Previous operation should not be a soft clip. Soft clip in the middle of
    // cigar is an error in alignment. It should not happen.
    CHECK(prev_op->operation() != CigarUnit::CLIP_SOFT);
    if (IsOperationMatch(*prev_op)) {
      CHECK(shift <= prev_op->operation_length());
      prev_op->set_operation_length(prev_op->operation_length() - shift);
    } else {
      // Do nothing if prev operation is not REF.
      return read_alignment_shift;
    }
  }

  // Expand existing ref following it or add a new one if it is the last element
  auto post_op = it + 1;
  if (post_op == norm_cigar.end()) {
    nucleus::genomics::v1::CigarUnit post_ref;
    post_ref.set_operation_length(shift);
    post_ref.set_operation(CigarUnit::ALIGNMENT_MATCH);
    norm_cigar.insert(it + 1, post_ref);
  } else {
    post_op->set_operation_length(post_op->operation_length() + shift);
  }
  return read_alignment_shift;
}

// Iterate cigar operations and attempt merging adjacent operations of the same
// type or indels. Returns true of merging was done.
bool FindAndMergeOperations(
    std::vector<nucleus::genomics::v1::CigarUnit>& norm_cigar) {
  for (auto op = norm_cigar.begin(); op != norm_cigar.end(); op++) {
    auto next_op = op + 1;
    if (next_op != norm_cigar.end() && MergeOperations(*op, *next_op)) {
        return true;
    }
  }
  return false;
}

// Remove all zero length operations and merge operations that can be merged.
// Operations of the same type are merged by adding their lengths. If DEL and
// INS has to be merged then their overlapping part is converted to REF and
// non overlapping part is preserved.
// For example. 3D5I (3 del and 5 ins) is merged into 3M2I (3 ref and 2 ins).
// Return true if any change was made to the cigar.
bool SwipeAndMerge(std::vector<nucleus::genomics::v1::CigarUnit>& norm_cigar) {
  // Repeat the loop until nothing is merged.
  bool merged = true;
  bool is_modified = false;
  while (merged) {
    merged = false;
    // First remove all operations of length zero
    auto before_size = norm_cigar.size();
    norm_cigar.erase(
        std::remove_if(norm_cigar.begin(), norm_cigar.end(),
                       [](const nucleus::genomics::v1::CigarUnit& op) {
                         return (op.operation_length() == 0);
                       }),
        norm_cigar.end());
    if (norm_cigar.size() < before_size) {
      is_modified = true;
    }
    // Then merge operations that are mergable.
    merged  = FindAndMergeOperations(norm_cigar);
    if (merged) {
      is_modified = true;
    }
  }  // while(merged)
  return is_modified;
}

bool AlleleCounter::CanDelBeShifted(
    const absl::string_view read_seq,
    std::vector<nucleus::genomics::v1::CigarUnit>::const_iterator cigar_elt,
    int read_offset,
    int interval_offset,
    int op_len) const {
  if (cigar_elt->operation() != CigarUnit::DELETE) {
    return false;
  }
  if (read_offset <= 0) {
    return false;
  }
  if (interval_offset + op_len - 1 >= ref_bases_.size()) {
    return false;
  }

  return  read_seq[read_offset - 1] ==
       ref_bases_[interval_offset + op_len - 1];
}

bool AlleleCounter::CanInsBeShifted(
    const absl::string_view read_seq,
    std::vector<nucleus::genomics::v1::CigarUnit>::const_iterator cigar_elt,
    int read_offset,
    int interval_offset,
    int op_len) const {
  if (cigar_elt->operation() != CigarUnit::INSERT) {
    return false;
  }
  if (interval_offset <= 0) {
    return false;
  }
  if (read_offset + op_len - 1 >= read_seq.size()) {
    return false;
  }

  return read_seq[read_offset + op_len - 1] ==
                     ref_bases_[interval_offset - 1];
}

// Normalize cigar of a given read following
// https://genome.sph.umich.edu/wiki/Variant_Normalization. As a result
// alignment position may need to be adjusted, in this case read_shift parameter
// is set to a non-zero value. As a result of shifting indel operations
// sometimes merging of adjecent indels may be performed.
bool AlleleCounter::NormalizeCigar(
    const absl::string_view read_seq, int interval_offset,
    std::vector<nucleus::genomics::v1::CigarUnit>& norm_cigar,
    int& read_shift) const {
  bool is_modified = false;
  read_shift = 0;
  if (norm_cigar.empty()) {
    return is_modified;
  }
  int iteration = 0;  // while loop will run up to 10 times.
  while (iteration++ < 10) {
    int read_offset = 0;
    int cur_interval_offset = interval_offset + read_shift;
    // Iterate cigar operations and shift indels if possible
    // If shift occurred break the loop and recalculate
    int prev_op_len = norm_cigar.front().operation_length();
    bool is_shifted = false;
    for (auto cigar_elt = norm_cigar.begin(); cigar_elt != norm_cigar.end();
         cigar_elt++) {
      const int op_len = cigar_elt->operation_length();
      int shift = 0;
      if (cigar_elt->operation() == CigarUnit::INSERT ||
          cigar_elt->operation() == CigarUnit::DELETE) {
        // Move INS/DEL to the left until last base of INS/DEL is the same
        // as REF base at the start of the operation. Moving left can only
        // be done if cigar has a REF operation preceding the INS/DEL. By
        // moving INS>DEL to the left we consume the length of a preceding
        // REF operation.
        while (prev_op_len > 0 &&
               ( CanDelBeShifted(read_seq, cigar_elt, read_offset,
                                 cur_interval_offset, op_len) ||
                CanInsBeShifted(read_seq, cigar_elt, read_offset,
                                 cur_interval_offset, op_len))) {
          cur_interval_offset--;
          prev_op_len--;
          read_offset--;
          shift++;
        }
        if (shift > 0) {
          read_shift += ShiftOperation(shift, cigar_elt, norm_cigar);
          is_modified = true;
          is_shifted = true;
          break;
        }
      }
      prev_op_len = cigar_elt->operation_length();
      AdvanceReadReferencePointers(op_len, cigar_elt->operation(), read_offset,
                                   cur_interval_offset);
    }  // for

    // Input BAM may contain non-normalized records, therefore we try to
    // SwipeAndMerge even if there was no shift.
    bool is_merged = SwipeAndMerge(norm_cigar);
    if (is_merged) {
      is_modified = true;
    }
    if (!is_shifted && !is_merged) {
      break;
    }
  }  // while (iteration < 10)
  // Call shift to deal with an indel at the beginning of cigar.
  read_shift += HandleHeadingIndel(norm_cigar.begin(), norm_cigar);
  return is_modified;
}

void AlleleCounter::NormalizeAndAdd(
    const nucleus::genomics::v1::Read& read, absl::string_view sample,
    std::unique_ptr<std::vector<nucleus::genomics::v1::CigarUnit>>& norm_cigar,
    int& read_shift) {
  // Make sure our incoming read has a mapping quality above our min. threshold.
  if (read.alignment().mapping_quality() <
      options_.read_requirements().min_mapping_quality()) {
    return;
  }

  const LinearAlignment& aln = read.alignment();
  std::vector<ReadAllele> to_add;
  to_add.reserve(read.aligned_quality_size());
  int interval_offset = aln.position().position() - ReadsInterval().start();
  const string_view read_seq(read.aligned_sequence());
  // Copy input cigar into the local variable since it can be modified.
  std::vector<CigarUnit> input_output_cigar(aln.cigar().begin(),
                                            aln.cigar().end());
  bool is_modified =
      NormalizeCigar(read_seq, interval_offset, input_output_cigar, read_shift);
  if (is_modified) {
    norm_cigar->assign(input_output_cigar.begin(), input_output_cigar.end());
  }
  Add(read, sample, &input_output_cigar, read_shift);
}

void AlleleCounter::Add(const nucleus::genomics::v1::Read& read,
                        absl::string_view sample,
                        const std::vector<CigarUnit>* cigar_to_use,
                        int read_shift) {
  // Make sure our incoming read has a mapping quality above our min. threshold.
  if (read.alignment().mapping_quality() <
      options_.read_requirements().min_mapping_quality()) {
    return;
  }

  const LinearAlignment& aln = read.alignment();
  std::vector<ReadAllele> to_add;
  to_add.reserve(read.aligned_quality_size());
  int read_offset = 0;
  int ref_interval_offset =
      aln.position().position() + read_shift - ReadsInterval().start();
  int interval_offset =
      aln.position().position() + read_shift - Interval().start();
  const string_view read_seq(read.aligned_sequence());
  std::vector<CigarUnit> cigar;
  if (cigar_to_use != nullptr) {
    cigar.assign(cigar_to_use->begin(), cigar_to_use->end());
  } else {
    cigar.assign(aln.cigar().begin(), aln.cigar().end());
  }

  for (const auto& cigar_elt : cigar) {
    const int op_len = cigar_elt.operation_length();
    switch (cigar_elt.operation()) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::SEQUENCE_MISMATCH:
        for (int i = 0; i < op_len; ++i) {
          const int ref_offset = ref_interval_offset + i;
          const int base_offset = read_offset + i;
          bool is_low_quality_read_allele = false;
          if (IsValidRefOffset(ref_offset) &&
              CanBasesBeUsed(read, base_offset, 1, options_,
                             is_low_quality_read_allele)) {
            const AlleleType type =
                ref_bases_[ref_offset] == read_seq[base_offset]
                    ? AlleleType::REFERENCE
                    : AlleleType::SUBSTITUTION;
            to_add.emplace_back(interval_offset + i,
                                string(read_seq.substr(base_offset, 1)), type,
                                is_low_quality_read_allele);
          }
        }
        read_offset += op_len;
        ref_interval_offset += op_len;
        interval_offset += op_len;
        break;
      case CigarUnit::CLIP_SOFT:
      case CigarUnit::INSERT:
        // Note, by convention VCF insertion/deletion are at the preceding base.
        to_add.push_back(MakeIndelReadAllele(read, interval_offset,
                                             ref_interval_offset, read_offset,
                                             cigar_elt));
        read_offset += op_len;
        // No interval offset change, since an insertion doesn't move us on ref.
        break;
      case CigarUnit::DELETE:
        // By convention VCF insertion/deletion are at the preceding base.
        to_add.push_back(MakeIndelReadAllele(read, interval_offset,
                                             ref_interval_offset, read_offset,
                                             cigar_elt));
        // No read offset change, since a deletion doesn't consume read bases.
        ref_interval_offset += op_len;
        interval_offset += op_len;
        break;
      case CigarUnit::PAD:
      case CigarUnit::SKIP:
        // No read offset change, since a pad/skip don't consume read bases.
        ref_interval_offset += op_len;
        interval_offset += op_len;
        break;
      case CigarUnit::CLIP_HARD:
        break;
      default:
        // Lots of misc. enumerated values from proto that aren't useful such as
        // enumeration values INT_MIN_SENTINEL_DO_NOT_USE_ and
        // OPERATION_UNSPECIFIED.
        break;
    }
  }

  AddReadAlleles(read, sample, to_add);
  ++n_reads_counted_;
}

string AlleleCounter::ReadKey(const Read& read) {
  return StrCat(read.fragment_name(), kFragmentNameReadNumberSeparator,
                read.read_number());
}

std::vector<AlleleCountSummary> AlleleCounter::SummaryCounts(
    int left_padding, int right_padding) const {
  std::vector<AlleleCountSummary> summaries;
  CHECK_GE(left_padding, 0);
  CHECK_GE(right_padding, 0);
  CHECK_LT(left_padding + right_padding, counts_.size());
  summaries.reserve(counts_.size() - left_padding - right_padding);
  for (int i = left_padding; i < counts_.size() - right_padding; i++) {
    const AlleleCount& allele_count = counts_[i];
    AlleleCountSummary summary;
    summary.set_reference_name(allele_count.position().reference_name());
    summary.set_position(allele_count.position().position());
    summary.set_ref_base(allele_count.ref_base());
    summary.set_ref_supporting_read_count(
        allele_count.ref_supporting_read_count());
    summary.set_total_read_count(TotalAlleleCounts(allele_count));
    summary.set_ref_nonconfident_read_count(
        allele_count.ref_nonconfident_read_count());
    summaries.push_back(summary);
  }
  return summaries;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
