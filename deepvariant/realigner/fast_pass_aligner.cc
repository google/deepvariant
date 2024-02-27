/*
 * Copyright 2018 Google LLC.
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

#include "deepvariant/realigner/fast_pass_aligner.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/memory/memory.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "re2/re2.h"

namespace learning {
namespace genomics {
namespace deepvariant {

void FastPassAligner::set_reference(const string& reference) {
  this->reference_ = reference;
}

void FastPassAligner::set_reads(const std::vector<string>& reads) {
  this->reads_ = reads;
}

void FastPassAligner::set_ref_start(const string& chromosome,
                                    uint64_t position) {
  this->region_chromosome_ = chromosome;
  this->region_position_in_chr_ = position;
}

void FastPassAligner::set_haplotypes(const std::vector<string>& haplotypes) {
  this->haplotypes_ = haplotypes;
}

void FastPassAligner::set_options(const AlignerOptions& options) {
  // There is no is_set method in proto so we assume that value is set if it is
  // not zero.
  if (options.kmer_size() > 0) {
    this->kmer_size_ = options.kmer_size();
  }
  if (options.read_size() > 0) {
    this->read_size_ = options.read_size();
  }
  if (options.max_num_of_mismatches() > 0) {
    this->max_num_of_mismatches_ = options.max_num_of_mismatches();
  }
  if (options.realignment_similarity_threshold() > 0.0) {
    this->similarity_threshold_ = options.realignment_similarity_threshold();
  }
  if (options.match() > 0) {
    this->match_score_ = options.match();
  }
  if (options.mismatch() > 0) {
    this->mismatch_penalty_ = options.mismatch();
  }
  if (options.gap_open() > 0) {
    this->gap_opening_penalty_ = options.gap_open();
  }
  if (options.gap_extend() > 0) {
    this->gap_extending_penalty_ = options.gap_extend();
  }
  this->force_alignment_ = options.force_alignment();

  CHECK(kmer_size_ >= 3 && kmer_size_ <= 32);
  CHECK_GE(similarity_threshold_, 0.0);
  CHECK_LE(similarity_threshold_, 1.0);
  CHECK_GE(max_num_of_mismatches_, 0);
  CHECK(kmer_size_ >= 3 && kmer_size_ <= 32);
}

void FastPassAligner::CalculateSswAlignmentScoreThreshold() {
  ssw_alignment_score_threshold_ = match_score_
      * read_size_
      * similarity_threshold_
          - mismatch_penalty_
      * read_size_
      * (1 - similarity_threshold_);
  if (ssw_alignment_score_threshold_ < 0) {
    ssw_alignment_score_threshold_ = 1;
  }
}

// Fast align reads to haplotypes using reads index.
// Align reads that could not be aligned in the first step using ssw aligner.
// Keep the best alignment for each read, or preserve an original one if read
// could not be realigned with a high enough score.
std::unique_ptr<std::vector<nucleus::genomics::v1::Read>>
FastPassAligner::AlignReads(
    const std::vector<nucleus::genomics::v1::Read>& reads_param) {

  // Copy reads
  for (const auto& read : reads_param) {
    reads_.push_back(absl::AsciiStrToUpper(read.aligned_sequence()));
  }

  CalculateSswAlignmentScoreThreshold();

  // Build index
  BuildIndex();

  // Align reads to haplotypes using reads index. This is O(n) operation per
  // read, where n = read size.
  FastAlignReadsToHaplotypes();

  // Initialize ssw library. Set reference.
  InitSswLib();

  // Align haplotypes to the reference.
  AlignHaplotypesToReference();

  // calculate position shifts.
  CalculatePositionMaps();

  // Align reads that couldn't be aligned in FastAlignReadsToHaplotypes using
  // ssw library.
  SswAlignReadsToHaplotypes(ssw_alignment_score_threshold_);

  // Sort haplotypes by number of supporting reads. First haplotype is the one
  // that has fewer supporting reads.
  std::sort(read_to_haplotype_alignments_.begin(),
            read_to_haplotype_alignments_.end());

  // Realign reads that we could successfully realign in previous steps back to
  // reference. From all read to haplotype alignments the best one is picked.
  // In the case where read alignments are equally good to ref haplotype and
  // non-ref haplotype, a non-ref haplotype is preferred.
  std::unique_ptr<std::vector<nucleus::genomics::v1::Read>> realigned_reads(
      new std::vector<nucleus::genomics::v1::Read>());
  RealignReadsToReference(reads_param, &realigned_reads);

  return realigned_reads;
}

void FastPassAligner::InitSswLib() {
  // Initialize ssw library. Set reference.
  Filter filter;
  ssw_aligner_ =
      std::make_unique<Aligner>(match_score_, mismatch_penalty_,
                                gap_opening_penalty_, gap_extending_penalty_);
}

void FastPassAligner::SswSetReference(const string& reference) {
  CHECK(ssw_aligner_);
  ssw_aligner_->SetReferenceSequence(reference);
}

Alignment FastPassAligner::SswAlign(const string& target) const {
  CHECK(ssw_aligner_);
  Filter filter;
  Alignment alignment;
  if (ssw_aligner_->Align(target, filter, &alignment) == 0) {
    return alignment;
  } else {
    LOG(WARNING) << "SSW alignment failed for query: '" << target << "'";
    return Alignment();
  }
}

// For each haplotype try to find all reads that can be aligned using index.
void FastPassAligner::FastAlignReadsToHaplotypes() {
  std::vector<ReadAlignment> read_alignment_scores(reads_.size());
  for (int i = 0; i < haplotypes_.size(); i++) {
    const auto& haplotype = haplotypes_[i];
    int haplotype_score = 0;
    for (auto& readAlignment : read_alignment_scores) {
      readAlignment.reset();
    }
    FastAlignReadsToHaplotype(haplotype,
                              &haplotype_score,
                              &read_alignment_scores);

    // haplotype_score == 0 means we found a problem with this haplotype. In
    // this case we need to discard of this haplotype.
    if (haplotype_score == 0) {
      for (auto& readAlignment : read_alignment_scores) {
        readAlignment.reset();
      }
    }

    read_to_haplotype_alignments_.push_back(
        HaplotypeReadsAlignment(i, haplotype_score, read_alignment_scores));
  }
}

void FastPassAligner::FastAlignReadsToHaplotype(
    absl::string_view haplotype, int* haplotype_score,
    std::vector<ReadAlignment>* haplotype_read_alignment_scores) {
  CHECK(haplotype_score != nullptr);
  CHECK(haplotype_read_alignment_scores != nullptr);

  bool is_ref = (haplotype == reference_);
  std::vector<int> coverage(haplotype.size(), 0);
  // In the loop we try to align reads for each position in haplotype up to
  // lastPos.
  const auto& lastPos = haplotype.length() - kmer_size_;
  for (int i = 0; i <= lastPos; i++) {
    // get all reads that are aligned against i-th position
    auto index_it = kmer_index_.find(haplotype.substr(i, kmer_size_));
    if (index_it == kmer_index_.end()) {
      continue;
    }
    // Iterate through all the reads that are found in the index for the current
    // kmer.
    for (const auto& it : index_it->second) {
      uint64_t read_id_index = static_cast<uint64_t>(it.read_id);
      CHECK(read_id_index < reads_.size() && it.read_id.is_set);
      size_t target_start_pos = std::max(
          static_cast<int64_t>(0),
          static_cast<int64_t>(i) - static_cast<int64_t>(it.read_pos.pos));
      size_t cur_read_size = reads_[read_id_index].size();
      size_t span = cur_read_size;
      if (target_start_pos + cur_read_size > haplotype.length()) {
        continue;
      }
      auto& read_alignment =
          (*haplotype_read_alignment_scores)[read_id_index];

      // This read is already aligned, skip it.
      if (read_alignment.position != ReadAlignment::kNotAligned &&
          read_alignment.position == target_start_pos) {
        continue;
      }
      CHECK(target_start_pos + span <= haplotype.size());
      int num_of_mismatches = 0;
      int new_read_alignment_score = FastAlignStrings(
          haplotype.substr(target_start_pos, span),
          reads_[read_id_index],
          max_num_of_mismatches_ + 1, &num_of_mismatches);

      if (num_of_mismatches <= max_num_of_mismatches_) {
        CHECK(it.read_id.is_set &&
            read_id_index < haplotype_read_alignment_scores->size());
        int oldScore = read_alignment.score;

        for (auto pos = target_start_pos; pos < target_start_pos + span;
             pos++) {
          coverage[pos]++;
        }

        if (oldScore < new_read_alignment_score) {
          read_alignment.score = new_read_alignment_score;
          *haplotype_score -= oldScore;
          *haplotype_score += read_alignment.score;
          read_alignment.position = target_start_pos;
          read_alignment.cigar = std::to_string(cur_read_size) + "=";
        }
      }
    }  // for (matching reads)
    // We want to discard haplotypes that don't have good read support.
    // At the same time we don't want to discard reference haplotype because,
    // there might be cases were not all haplotypes were generated and we can
    // get away with that by aligning reads to reference haplotype.
    if (coverage[i] == 0 && i >= ref_prefix_len_ &&
        i < haplotype.size() - ref_suffix_len_ && !is_ref) {
      *haplotype_score = 0;
      return;
    }
  }    // for (all k-mer positions)
}

// Align 2 same length strings by comparing each character.
int FastPassAligner::FastAlignStrings(absl::string_view s1,
                                      absl::string_view s2,
                                      int max_mismatches,
                                      int* num_of_mismatches) const {
  int num_of_matches = 0;
  *num_of_mismatches = 0;
  CHECK(s1.size() == s2.size());
  for (int i = 0; i < s1.size(); i++) {
    const auto& c1 = s1[i];
    const auto& c2 = s2[i];
    if (c1 != c2 && (c1 != 'N' && c2 != 'N')) {
      if (c1 != c2) {
        (*num_of_mismatches)++;
      }
      if (*num_of_mismatches == max_mismatches) {
        return 0;
      }
    } else {
      num_of_matches++;
    }
  }
  return num_of_matches * match_score_ - *num_of_mismatches * mismatch_penalty_;
}

CigarUnit::Operation CigarOperationFromChar(char op) {
  switch (op) {
    case '=':
    case 'X':
      return nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH;
    case 'S':
      return nucleus::genomics::v1::CigarUnit::CLIP_SOFT;
    case 'D':
      return nucleus::genomics::v1::CigarUnit::DELETE;
    case 'I':
      return nucleus::genomics::v1::CigarUnit::INSERT;
    default:
      return nucleus::genomics::v1::CigarUnit::OPERATION_UNSPECIFIED;
  }
}

std::list<CigarOp> CigarStringToVector(absl::string_view cigar) {
  std::list<CigarOp> cigarOps;
  absl::string_view input(cigar);
  RE2 pattern("(\\d+)([XIDS=])");
  int opLen;
  string opType;
  while (RE2::Consume(&input, pattern, &opLen, &opType)) {
    CHECK_EQ(opType.length(), 1);
    CigarUnit::Operation op = CigarOperationFromChar(opType[0]);
    cigarOps.push_back(CigarOp(op, opLen));
  }
  return cigarOps;
}

inline bool AlignmentIsRef(absl::string_view cigar, size_t target_len) {
  return cigar == absl::StrCat(target_len, "=");
}

// Align haplotypes to reference using ssw library.
void FastPassAligner::AlignHaplotypesToReference() {
  SswSetReference(reference_);

  // Initialize read_to_haplotype_alignments_ if it is not initialized yet.
  if (read_to_haplotype_alignments_.empty()) {
    for (int i = 0; i < haplotypes_.size(); i++) {
      read_to_haplotype_alignments_.push_back(HaplotypeReadsAlignment(
          i, -1, std::vector<ReadAlignment>(reads_.size())));
    }
  }

  for (auto& haplotype_alignment : read_to_haplotype_alignments_) {
    Filter filter;
    CHECK(haplotype_alignment.haplotype_index < haplotypes_.size());
    auto hap_len = haplotypes_[haplotype_alignment.haplotype_index].size();
    // Most of the time, one haplotype will perfectly match the reference.
    if (haplotypes_[haplotype_alignment.haplotype_index] == reference_) {
      haplotype_alignment.is_reference = true;
      haplotype_alignment.cigar = absl::StrCat(hap_len, "=");
      haplotype_alignment.cigar_ops =
          CigarStringToVector(haplotype_alignment.cigar);
      haplotype_alignment.ref_pos = 0;
    } else {
      Alignment alignment =
          SswAlign(haplotypes_[haplotype_alignment.haplotype_index]);
      if (alignment.sw_score > 0) {
        // In rare cases, the ref haplotype will be a substring of the ref, and
        // therefore not caught by the string equality check above.
        haplotype_alignment.is_reference =
            AlignmentIsRef(alignment.cigar_string, hap_len);

        haplotype_alignment.cigar = alignment.cigar_string;
        haplotype_alignment.cigar_ops =
            CigarStringToVector(haplotype_alignment.cigar);
        haplotype_alignment.ref_pos = alignment.ref_begin;
      }
    }
  }
}

void FastPassAligner::SswAlignReadsToHaplotypes(uint16_t score_threshold) {
  // For each read
  for (int i = 0; i < reads_.size(); i++) {
    bool has_at_least_one_alignment = false;
    // Check if this read is aligned to at least one haplotype
    for (const auto& hap_alignment : read_to_haplotype_alignments_) {
      if (hap_alignment.read_alignment_scores[i].score > 0) {
        has_at_least_one_alignment = true;
        break;
      }
    }
    // If this read is not aligned to any of the haplotypes we try SSW.
    if (!has_at_least_one_alignment) {
      for (auto& hap_alignment : read_to_haplotype_alignments_) {
        // Skip haplotypes with no read support (score=0), except if
        // force_alignment, then compute an alignment against the reference no
        // matter what.
        if (hap_alignment.haplotype_score == 0 &&
            !(force_alignment_ && hap_alignment.is_reference)) {
          continue;
        }
        CHECK(hap_alignment.haplotype_index < haplotypes_.size());
        SswSetReference(haplotypes_[hap_alignment.haplotype_index]);
        Alignment alignment = SswAlign(reads_[i]);
        if (alignment.sw_score > 0) {
          // TODO Remove score_threshold condition. It is effectively
          // not used.
          if (alignment.sw_score >= score_threshold ||
              (force_alignment_ && hap_alignment.is_reference)) {
            hap_alignment.read_alignment_scores[i].score = alignment.sw_score;
            hap_alignment.read_alignment_scores[i].cigar =
                alignment.cigar_string;
            hap_alignment.read_alignment_scores[i].position =
                alignment.ref_begin;
          }
        } else if (force_alignment_ && hap_alignment.is_reference) {
        }
      }
    }
  }  // for all reads
}

// Each operation except MATCH is checked in the input cigar.
// If reference base preceding operation is the same as the last base of the
// operation alt bases then alignment is not normalized.
// Example:
// Reference: GATCGAGAGAGAGATC
// Read.    : GATCGA--GAGAGATC
// This alignment is not normalized because last base in deletion 'GA' is the
// same as the base preceding the deletion in the read.
// The reason we sometimes see these non-normalized alignments is because reads
// are aligned to the haplotype first. Continuing the example above, the
// haplotype could be the following:
// Reference: GATCGAGAGAGAGATC
// Haplotype: GATCGT--GAGAGATC
// Read.    : GATCGA--GAGAGATC
// In this case haplotype to the reference alignment is normalized (cannot be
// shifted left), read to haplotype alignment is also normalized, but read to
// reference alignment is not normalized.
// Normalizing the alignment could be done but it is not easy. If we shift
// cigar operation to the left we may need to deal with merging it with other
// operations that may exist.
bool FastPassAligner::IsAlignmentNormalized(
    const std::list<CigarOp>& cigar,
    int ref_offset,
    const absl::string_view read_sequence) const {
  if (ref_offset < 0) {
    return true;
  }
  int cur_ref_offset = ref_offset;
  int cur_read_offset = 0;
  for (const auto& op : cigar) {
    if (op.operation == nucleus::genomics::v1::CigarUnit::CLIP_SOFT) {
      cur_read_offset += op.length;
      continue;
    }
    if (op.operation != nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH) {
      std::string op_sequence;
      if  (op.operation == nucleus::genomics::v1::CigarUnit::DELETE) {
        if (cur_ref_offset + op.length > reference_.size()) {
          return false;
        }
        CHECK(cur_ref_offset + op.length <= reference_.size());
        op_sequence = reference_.substr(cur_ref_offset, op.length);
      } else {
        CHECK(cur_read_offset + op.length <= read_sequence.size());
        op_sequence =  read_sequence.substr(cur_read_offset, op.length);
      }
      if ((op.operation == nucleus::genomics::v1::CigarUnit::INSERT
          && op_sequence.back() == reference_[cur_ref_offset - 1]) ||
         (cur_read_offset > 0 &&
          op.operation == nucleus::genomics::v1::CigarUnit::DELETE
          && op_sequence.back() == read_sequence[cur_read_offset - 1])) {
        return false;
      }
    }
    if (op.operation != nucleus::genomics::v1::CigarUnit::INSERT) {
      cur_ref_offset += op.length;
    }
    if (op.operation != nucleus::genomics::v1::CigarUnit::DELETE) {
      cur_read_offset += op.length;
    }
  }
  return true;
}

void FastPassAligner::RealignReadsToReference(
    const std::vector<nucleus::genomics::v1::Read>& reads,
    std::unique_ptr<std::vector<nucleus::genomics::v1::Read>>*
        realigned_reads) {
  // Loop through all reads
  for (size_t read_index = 0; read_index < reads.size(); read_index++) {
    const nucleus::genomics::v1::Read& read = reads[read_index];
    nucleus::genomics::v1::Read realigned_read;
    realigned_read.MergeFrom(read);
    int best_hap_index = -1;
    // See if we have a better alignment
    if (GetBestReadAlignment(read_index, &best_hap_index)) {
      const HaplotypeReadsAlignment& bestHaplotypeAlignments =
          read_to_haplotype_alignments_[best_hap_index];
      std::unique_ptr<LinearAlignment> new_alignment(new LinearAlignment());
      new_alignment->MergeFrom(read.alignment());
      new_alignment->clear_cigar();
      // Calculate new alignment position.
      std::unique_ptr<nucleus::genomics::v1::Position> new_position(
          new nucleus::genomics::v1::Position());
      new_position->MergeFrom(read.alignment().position());
      auto read_to_hap_pos = bestHaplotypeAlignments
          .read_alignment_scores[read_index]
          .position;
      CHECK(read_to_hap_pos < bestHaplotypeAlignments
          .hap_to_ref_positions_map.size());
      int hap_to_ref_position = bestHaplotypeAlignments
          .hap_to_ref_positions_map[read_to_hap_pos];
      // We only change position of original read alignment and don't change
      // chromosome, it shouldn't change anyway!
      new_position->set_position(
          region_position_in_chr_
              + bestHaplotypeAlignments.ref_pos
              + read_to_hap_pos
              + hap_to_ref_position);
      new_alignment->set_allocated_position(new_position.release());
      std::list<CigarOp> readToRefCigarOps;
      // Calculate new cigar by merging read to haplotype and haplotype to ref
      // alignments.
      CalculateReadToRefAlignment(
          read_index, bestHaplotypeAlignments.read_alignment_scores[read_index],
          bestHaplotypeAlignments.cigar_ops, &readToRefCigarOps);

      // The following block is only executed if normalize_reads flag is not
      // set. This is because if --normalize_reads is true, they will be
      // normalize later on.
      if (!normalize_reads_) {
        // If read is not normalized (any cigar operation can be shifter left)
        // we discard the alignment.
        if (!IsAlignmentNormalized(
            readToRefCigarOps,
            bestHaplotypeAlignments.ref_pos
                + read_to_hap_pos
                + hap_to_ref_position,
            reads_[read_index])) {
            readToRefCigarOps.clear();
        }
      }

      for (auto& op : readToRefCigarOps) {
        CigarUnit* cu = new_alignment->add_cigar();
        cu->set_operation(op.operation);
        cu->set_operation_length(op.length);
      }

      if (!readToRefCigarOps.empty()) {
        realigned_read.set_allocated_alignment(new_alignment.release());
      } else if (force_alignment_) {
      }
      (*realigned_reads)->push_back(realigned_read);
    } else {  // Could not find a new alignment.
      if (force_alignment_) {
      } else {
        // Keeping original alignment (force_alignment is off).
        (*realigned_reads)->push_back(realigned_read);
      }
    }
  }  // for
}

void FastPassAligner::AddKmerToIndex(absl::string_view kmer,
                                     ReadId read_id, KmerOffset pos) {
  kmer_index_[kmer].push_back(KmerOccurrence(read_id, pos));
}

void FastPassAligner::AddReadToIndex(absl::string_view read, ReadId read_id) {
  // Ignoring reads that are too short for a kmer size. Those reads will still
  // be realigned with SSW.
  if (read.length() <= kmer_size_) {
    return;
  }
  auto last_pos = read.length() - kmer_size_;
  for (int i = 0; i <= last_pos; i++) {
    AddKmerToIndex(read.substr(i, kmer_size_), read_id, KmerOffset(i));
  }
}

void FastPassAligner::BuildIndex() {
  size_t read_id = 0;
  for (const auto& read : reads_) {
    AddReadToIndex(read, ReadId(read_id++));
  }
}

void SetPositionsMap(size_t haplotype_size,
                     HaplotypeReadsAlignment* hyplotype_alignment) {
  std::vector<int>& positions_map =
      hyplotype_alignment->hap_to_ref_positions_map;
  positions_map.resize(haplotype_size);
  RE2 pattern("(\\d+)([XIDS=])");  // matches cigar operation
  absl::string_view input(hyplotype_alignment->cigar);
  int cur_shift = 0;
  int haplotype_pos = 0;
  int last_pos = 0;
  int operation_len;
  string operation_type;
  while (RE2::Consume(&input, pattern, &operation_len, &operation_type)) {
    CHECK_EQ(operation_type.length(), 1);

    char op = operation_type[0];
    switch (op) {
      case '=':
      case 'X':
        last_pos = haplotype_pos + operation_len;
        while (haplotype_pos != last_pos) {
          positions_map[haplotype_pos] = cur_shift;
          haplotype_pos++;
        }
        break;
      case 'S':
        last_pos = haplotype_pos + operation_len;
        cur_shift -= operation_len;
        while (haplotype_pos != last_pos) {
          positions_map[haplotype_pos] = cur_shift;
          haplotype_pos++;
        }
        break;
      case 'D':
        cur_shift += operation_len;
        break;
      case 'I':
        last_pos = haplotype_pos + operation_len;
        while (haplotype_pos != last_pos) {
          positions_map[haplotype_pos] = cur_shift;
          cur_shift--;
          haplotype_pos++;
        }
        break;
    }
  }
}

void FastPassAligner::CalculatePositionMaps() {
  for (auto& hyplotype_alignment : read_to_haplotype_alignments_) {
    SetPositionsMap(haplotypes_[hyplotype_alignment.haplotype_index].size(),
                    &hyplotype_alignment);
  }
}

bool FastPassAligner::GetBestReadAlignment(
    size_t readId,
    int* best_hap_index) const {
  int best_score = 0;
  bool best_haplotype_found = false;
  for (int hap_index = 0; hap_index < haplotypes_.size(); hap_index++) {
    int hap_score = read_to_haplotype_alignments_[hap_index]
                        .read_alignment_scores[readId]
                        .score;
    // If compared scores are equal, preference is given to a read alignment
    // against a non-reference haplotype.
    if (hap_score > best_score ||
        (best_score > 0 && hap_score == best_score &&
         !read_to_haplotype_alignments_[hap_index].is_reference)) {
      best_score = read_to_haplotype_alignments_[hap_index]
                       .read_alignment_scores[readId]
                       .score;
      *best_hap_index = hap_index;
      best_haplotype_found = true;
    }
  }
  return best_haplotype_found;
}

// Calculate aligned length from cigar.
int AlignedLength(const std::list<CigarOp>& cigar) {
  int len = 0;
  for (auto& op : cigar) {
    if (op.operation != nucleus::genomics::v1::CigarUnit::DELETE) {
      len += op.length;
    }
  }
  return len;
}

// Merge cigar op to the end of the output cigar.
// <op> input parameter is always a one base operation to be merged.
// <cigar> input/output parameter is our merged output cigar.
// - If op.operation is the same as the last one in the cigar then the length
//   of the last operation is increased by op.length.
// - If op.operation is not the same then new operation is added to the cigar.
// For all operations except DELETE we need to make sure that aligned length
// does not go over the read length. This can happen for example when we merge
// a large INS (larger than a read itself)
void MergeCigarOp(const CigarOp& op, int read_len, std::list<CigarOp>* cigar) {
  const auto& last_cigar_op =
      cigar->empty() ? nucleus::genomics::v1::CigarUnit::OPERATION_UNSPECIFIED
                     : cigar->back().operation;
  int aligned_length_before_merge = AlignedLength(*cigar);
  int new_op_length = 0;
  if (op.operation != nucleus::genomics::v1::CigarUnit::DELETE) {
    new_op_length = std::min(op.length, read_len - aligned_length_before_merge);
  } else {
    new_op_length = op.length;
  }

  // Nothing is merged if we already aligned all positions of the read.
  if (new_op_length <= 0 || aligned_length_before_merge == read_len) {
    return;
  }

  // Special processing is used when merging INS and DEL. Going one base at a
  // time INS and DEL cancel each other. In addition we need to add new
  // match operations in place of removed INS and DEL. New Match operation is
  // added right before the trailing INDEL operation.
  if ((op.operation == nucleus::genomics::v1::CigarUnit::INSERT
      && last_cigar_op == nucleus::genomics::v1::CigarUnit::DELETE)
      ||
     (op.operation == nucleus::genomics::v1::CigarUnit::DELETE
     && last_cigar_op == nucleus::genomics::v1::CigarUnit::INSERT)) {
    // Determine the last element in the current cigar list.
    std::list<CigarOp>::iterator last_element = cigar->end();
    if (!cigar->empty()) {
      last_element--;
    }

    // Determine one before last element in cigar list.
    std::list<CigarOp>::iterator one_before_last = cigar->end();
    if (cigar->size() > 1) {
      one_before_last = std::prev(one_before_last, 2);
    } else {
      one_before_last = std::prev(one_before_last, 1);
    }
    if (one_before_last->operation !=
        nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH) {
      cigar->insert(last_element,
          CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 1));
    } else {
      one_before_last->length += 1;
    }

    // Reduce the size of the last operation. If last operation's length was 1
    // then it is removed.
    if (cigar->back().length == 1) {
      cigar->pop_back();
    } else {
      cigar->back().length -= 1;
    }
  } else if (op.operation == last_cigar_op) {
    cigar->back().length += new_op_length;

    // If the op we are adding is not the same as the last, set the proper
    // length and add it to the list.
  }  else {
    cigar->push_back(CigarOp(op.operation, new_op_length));
  }
}

// Following functions are for internal use only.
namespace {

std::list<CigarOp> LeftTrimHaplotypeToRefAlignment(
    const std::list<CigarOp>& haplotype_to_ref_cigar_ops_input,
    int read_to_haplotype_pos) {
  int cur_pos = 0;
  std::list<CigarOp> haplotype_to_ref_cigar_ops(
      haplotype_to_ref_cigar_ops_input);
  while (cur_pos != read_to_haplotype_pos) {
    CHECK(!haplotype_to_ref_cigar_ops.empty());
    CigarOp cur_hap_op = haplotype_to_ref_cigar_ops.front();
    haplotype_to_ref_cigar_ops.pop_front();
    if (cur_hap_op.operation ==
            nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH ||
        cur_hap_op.operation == nucleus::genomics::v1::CigarUnit::CLIP_HARD ||
        cur_hap_op.operation == nucleus::genomics::v1::CigarUnit::CLIP_SOFT ||
        cur_hap_op.operation == nucleus::genomics::v1::CigarUnit::INSERT) {
      if (cur_hap_op.length + cur_pos > read_to_haplotype_pos) {
        haplotype_to_ref_cigar_ops.push_front(
            CigarOp(cur_hap_op.operation,
                    cur_hap_op.length - (read_to_haplotype_pos - cur_pos)));
      }
      cur_pos = std::min(cur_hap_op.length + cur_pos, read_to_haplotype_pos);
    }
  }

  // If after trimming the first operation is DEL we need to remove it,
  // because read alignment cannot start with DEL.
  if (haplotype_to_ref_cigar_ops.front().operation ==
      nucleus::genomics::v1::CigarUnit::DELETE) {
    haplotype_to_ref_cigar_ops.pop_front();
  }

  return haplotype_to_ref_cigar_ops;
}

// Merging one base operations.
// This function handles all possible combinations of one base merges except
// INS+DEL and DEL+INS. Below is the list of all possible combinations of
// operations and how they are resolved. Note, that all mergings below are
// symmetrical (DEL + MATCH is the same as MATCH + DEL).
// DEL + MATCH = DEL
// INS + MATCH = INS
// DEL + DEL = DEL
// INS + INS = INS
// CLIP_SOFT + MATCH = CLIP_SOFT
// CLIP_SOFT + CLIP_SOFT = CLIP_SOFT
// MATCH + MATCH = MATCH
// INS + DEL = Exception!
// DEL + INS = Exception!
inline void MergeOneBaseOperations(
    const CigarOp& cur_read_to_hap_op,
    const CigarOp& cur_hap_to_ref_op,
    int read_len,
    std::list<CigarOp>* read_to_ref_cigar_ops) {
    // Assert that operations are not INS and DEL
    std::vector<CigarUnit::Operation> input_operations{
      cur_read_to_hap_op.operation, cur_hap_to_ref_op.operation};
    CHECK((input_operations != std::vector<CigarUnit::Operation>{
      nucleus::genomics::v1::CigarUnit::DELETE,
      nucleus::genomics::v1::CigarUnit::INSERT}) &&
         (input_operations != std::vector<CigarUnit::Operation>{
      nucleus::genomics::v1::CigarUnit::INSERT,
      nucleus::genomics::v1::CigarUnit::DELETE}));

    std::vector<CigarUnit::Operation> operations{
        nucleus::genomics::v1::CigarUnit::CLIP_SOFT,
        nucleus::genomics::v1::CigarUnit::DELETE,
        nucleus::genomics::v1::CigarUnit::INSERT,
        nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH
        };
    for (auto op : operations) {
      if (cur_read_to_hap_op.operation == op
          || cur_hap_to_ref_op.operation == op) {
        MergeCigarOp(CigarOp(op, 1), read_len, read_to_ref_cigar_ops);
        break;
      }
    }
}

}  // namespace

void FastPassAligner::CalculateReadToRefAlignment(
    size_t read_index,
    const ReadAlignment& read_to_haplotype_alignment,
    const std::list<CigarOp>& haplotype_to_ref_cigar_ops_input,
    std::list<CigarOp>* read_to_ref_cigar_ops) const {
  CHECK(read_index < reads_.size());
  int read_len = reads_[read_index].length();
  int read_to_haplotype_pos = read_to_haplotype_alignment.position;
  std::list<CigarOp> read_to_haplotype_cigar_ops =
      CigarStringToVector(read_to_haplotype_alignment.cigar);

  // Left trim haplotype to reference cigar to match read to haplotype
  // alignment position.
  std::list<CigarOp> haplotype_to_ref_cigar_ops =
      LeftTrimHaplotypeToRefAlignment(haplotype_to_ref_cigar_ops_input,
                                      read_to_haplotype_pos);

  // Sanity check. By design haplotype is built from reads. Therefore it should
  // be impossible that read does not overlap with haplotype.
  CHECK(!haplotype_to_ref_cigar_ops.empty());

  // Skip heading soft clips.
  if (!read_to_haplotype_cigar_ops.empty() &&
      read_to_haplotype_cigar_ops.front().operation ==
          nucleus::genomics::v1::CigarUnit::CLIP_SOFT) {
    MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::CLIP_SOFT,
                         read_to_haplotype_cigar_ops.front().length),
                 read_len, read_to_ref_cigar_ops);
    read_to_haplotype_cigar_ops.pop_front();
  }

  CigarOp cur_read_to_hap_op(
      nucleus::genomics::v1::CigarUnit::OPERATION_UNSPECIFIED, 0);
  CigarOp cur_hap_to_ref_op(
      nucleus::genomics::v1::CigarUnit::OPERATION_UNSPECIFIED, 0);

  // Build read to reference cigar by iterating CigarOp overlaps.
  while ((!read_to_haplotype_cigar_ops.empty() ||
          !haplotype_to_ref_cigar_ops.empty()) &&
         AlignedLength(*read_to_ref_cigar_ops) < read_len) {
    // TODO Need to verify this logic.
    // This can happen if read was aligned to hyplotype partially. In this case
    // The tail (or head) of read to haplotype alignment would be soft-clipped.
    if (!read_to_haplotype_cigar_ops.empty() &&
        haplotype_to_ref_cigar_ops.empty() && cur_hap_to_ref_op.length == 0) {
      MergeCigarOp(read_to_haplotype_cigar_ops.front(), read_len,
                   read_to_ref_cigar_ops);
      read_to_haplotype_cigar_ops.pop_front();
      continue;
    }

    // Read is aligned completely, we are done.
    if (read_to_haplotype_cigar_ops.empty() && cur_read_to_hap_op.length == 0 &&
        !haplotype_to_ref_cigar_ops.empty()) {
      break;
    }

    // Assign current Cigar Ops for each alignment.
    if (cur_read_to_hap_op.length == 0) {
      cur_read_to_hap_op = read_to_haplotype_cigar_ops.front();
      read_to_haplotype_cigar_ops.pop_front();
    }
    if (cur_hap_to_ref_op.length == 0) {
      cur_hap_to_ref_op = haplotype_to_ref_cigar_ops.front();
      haplotype_to_ref_cigar_ops.pop_front();
    }

    while (cur_read_to_hap_op.length > 0 && cur_hap_to_ref_op.length > 0) {
      if ((cur_read_to_hap_op.operation
              == nucleus::genomics::v1::CigarUnit::DELETE &&
          cur_hap_to_ref_op.operation
              == nucleus::genomics::v1::CigarUnit::INSERT)
           ||
          (cur_read_to_hap_op.operation
              == nucleus::genomics::v1::CigarUnit::INSERT &&
          cur_hap_to_ref_op.operation
              == nucleus::genomics::v1::CigarUnit::DELETE)) {
        cur_hap_to_ref_op.length--;
        cur_read_to_hap_op.length--;
        // When hap_to_ref deletion is consumed by read_to_hap insertion we
        // need to convert the deletion to match.
        if (cur_hap_to_ref_op.operation ==
            nucleus::genomics::v1::CigarUnit::DELETE) {
          haplotype_to_ref_cigar_ops.push_front(
              CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 1));
          read_to_haplotype_cigar_ops.push_front(
              CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 1));
        }
        continue;
      }

      MergeOneBaseOperations(cur_read_to_hap_op, cur_hap_to_ref_op,
                             read_len, read_to_ref_cigar_ops);

      // If read_to_hap is INS or DEL then we consume these operations first.
      if (cur_read_to_hap_op.operation ==
              nucleus::genomics::v1::CigarUnit::INSERT) {
        cur_read_to_hap_op.length--;
      } else if (cur_hap_to_ref_op.operation ==
              nucleus::genomics::v1::CigarUnit::DELETE) {
        cur_hap_to_ref_op.length--;
      } else {
        cur_hap_to_ref_op.length--;
        cur_read_to_hap_op.length--;
      }
    }  // while
  }  // while
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
