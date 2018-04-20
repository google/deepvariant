/*
 * Copyright 2018 Google Inc.
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

#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include "re2/re2.h"

#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/platform/logging.h"

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

void FastPassAligner::set_kmer_size(int kmer_size) {
  // Having kmer size less than 32 allow to make a better index in the future.
  // Otherwise there is no need for this restriction.
  // Kmer size should be less than read size.
  CHECK(kmer_size >= 3 && kmer_size <= 32);
  kmer_size_ = kmer_size;
}

void FastPassAligner::set_read_size(int read_size) {
  CHECK_GT(read_size, 0);
  read_size_ = read_size;
}

void FastPassAligner::set_max_num_of_mismatches(int max_num_of_mismatches) {
  CHECK_GE(max_num_of_mismatches, 0);
  this->max_num_of_mismatches_ = max_num_of_mismatches;
}

void FastPassAligner::set_score_schema(uint8_t match_score,
                                       uint8_t mismatch_penalty,
                                       uint8_t gap_opening_penalty,
                                       uint8_t gap_extending_penalty) {
  this->match_score_ = match_score;
  this->mismatch_penalty_ = mismatch_penalty;
  this->gap_opening_penalty_ = gap_opening_penalty;
  this->gap_extending_penalty_ = gap_extending_penalty;
}

// Align reads to all haplotypes using reads index. Based on a number of aligned
// reads choose x[polidy] number of best haplotypes. Having x[ploidy] haplotypes
// align reads that could not be aligned in the first step using ssw aligner.
// Keep the best alignment for each read.
std::unique_ptr<ReadsVectorType> FastPassAligner::AlignReads(
    const std::vector<nucleus::genomics::v1::Read>& reads_param) {
  // redacted
  return nullptr;
}

void FastPassAligner::InitSswLib() {
  // Initialize ssw library. Set reference.
  // redacted
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
    read_to_haplotype_alignments_.push_back(
        HaplotypeReadsAlignment(i, haplotype_score, read_alignment_scores));
  }
}

void FastPassAligner::FastAlignReadsToHaplotype(
    const string& haplotype, int* haplotype_score,
    std::vector<ReadAlignment>* haplotype_read_alignment_scores) {
  CHECK(haplotype_score != nullptr);
  CHECK(haplotype_read_alignment_scores != nullptr);
  tensorflow::StringPiece bases_view(haplotype);

  // In the loop we try to align reads for each position in haplotype up to
  // lastPos.
  const auto& lastPos = haplotype.length() - kmer_size_;
  for (int i = 0; i <= lastPos; i++) {
    // get all reads that are aligned against i-th position
    auto index_it = kmer_index_.find(bases_view.substr(i, kmer_size_));
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
      size_t readStartPos = 0;
      size_t cur_read_size = reads_[read_id_index].size();
      size_t span = cur_read_size;
      if (target_start_pos + cur_read_size > haplotype.length()) {
        continue;
      }
      CHECK(target_start_pos + span <= bases_view.size());
      int num_of_mismatches = 0;
      int new_read_alignment_score = FastAlignStrings(
          bases_view.substr(target_start_pos, span),
          reads_[read_id_index].substr(readStartPos, span),
          max_num_of_mismatches_ + 1, &num_of_mismatches);

      if (num_of_mismatches <= max_num_of_mismatches_) {
        CHECK(it.read_id.is_set &&
            read_id_index < haplotype_read_alignment_scores->size());
        auto& read_alignment =
            (*haplotype_read_alignment_scores)[read_id_index];
        int oldScore = read_alignment.score;
        if (oldScore < new_read_alignment_score) {
          read_alignment.score = new_read_alignment_score;
          *haplotype_score -= oldScore;
          *haplotype_score += read_alignment.score;
          read_alignment.position = target_start_pos;
          read_alignment.cigar = std::to_string(cur_read_size) + "=";
        }
      }
    }  // for (matching reads)
  }    // for (all k-mer positions)
}

// Align 2 same length strings by comparing each character.
int FastPassAligner::FastAlignStrings(const tensorflow::StringPiece& s1,
                                      const tensorflow::StringPiece& s2,
                                      int max_mismatches,
                                      int* num_of_mismatches) const {
  int num_of_matches = 0;
  *num_of_mismatches = 0;
  CHECK(s1.size() == s2.size());
  for (int i = 0; i < s1.size(); i++) {
    if (s1[i] != s2[i] && (s1[i] != 'N' && s2[i] != 'N')) {
      if (s1[i] != s2[i]) {
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

// Align haplotypes to reference using ssw library and update bestHaplotypes
// vector.
void FastPassAligner::AlignHaplotypesToReference() {
  // redacted
}

void FastPassAligner::AlignReadsToHaplotypes(uint16_t score_threshold) {
  // redacted
}

void FastPassAligner::RealignReadsToReference(
    const std::vector<nucleus::genomics::v1::Read>& reads,
    std::unique_ptr<ReadsVectorType> realigned_reads) {
  // redacted
}

void FastPassAligner::AddKmerToIndex(tensorflow::StringPiece kmer,
                                     ReadId read_id,
                                     KmerOffset pos) {
  kmer_index_[kmer].push_back(KmerOccurrence(read_id, pos));
}

void FastPassAligner::AddReadToIndex(const string& read, ReadId read_id) {
  CHECK(read.length() > kmer_size_);
  auto last_pos = read.length() - kmer_size_;
  tensorflow::StringPiece bases_view(read);
  for (int i = 0; i <= last_pos; i++) {
    AddKmerToIndex(bases_view.substr(i, kmer_size_), read_id, KmerOffset(i));
  }
}

void FastPassAligner::BuildIndex() {
  size_t read_id = 0;
  for (const auto& read : reads_) {
    AddReadToIndex(read, ReadId(read_id++));
  }
}

CigarUnit::Operation CigarOperationFromChar(char op) {
  // redacted
  return nucleus::genomics::v1::CigarUnit_Operation_OPERATION_UNSPECIFIED;
}

void FastPassAligner::CigarStringToVector(const string& cigar,
                                          std::list<CigarOp>* cigar_ops) {
  // redacted
}

void FastPassAligner::CalculatePositionMaps() {
  // redacted
}

// Return true if cigar has at least one structural variation
bool hasIndels(const std::list<CigarOp>& cigar) {
  // redacted
  return false;
}

bool FastPassAligner::GetBestReadAlignment(
    size_t readId,
    int* best_hap_index) const {
  // redacted
  return false;
}

// Calculate aligned length from cigar.
int alignedLength(const std::list<CigarOp>& cigar) {
  // redacted
  return -1;
}

// Merge op with the last CigarOp if they are the same, otherwise adds a new
// CigarOp at the end.
void mergeCigarOp(const CigarOp& op, std::list<CigarOp>* cigar, int read_size,
                  bool debug = false) {
  // redacted
}

void FastPassAligner::CalculateCigarForRead(
    size_t read_index,
    std::list<CigarOp>* read_to_ref_cigar) const {
  // redacted
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
