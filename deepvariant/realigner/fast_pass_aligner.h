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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_FAST_PASS_ALIGNER_H_
#define LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_FAST_PASS_ALIGNER_H_

#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "deepvariant/realigner/ssw.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/hash/hash.h"
#include "tensorflow/core/platform/types.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Position;

// redacted

struct Kmer {
  string sequence;
};

struct ReadId {
  ReadId() : is_set(false), id(0) {}
  explicit ReadId(size_t id) : is_set(true), id(id) {}

  explicit operator int64_t() const { return id; }
  explicit operator uint64_t() const { return id; }

  bool operator<(const ReadId& that) const { return id < that.id; }
  bool operator==(const ReadId& that) const { return id == that.id; }

  bool is_set;
  size_t id;
};

struct KmerOffset {
  KmerOffset() : is_set(false), pos(0) {}
  explicit KmerOffset(size_t pos) : is_set(true), pos(pos) {}

  bool operator==(const KmerOffset& that) const {
    return pos == that.pos && is_set == that.is_set;
  }

  bool is_set;
  size_t pos;
};

// Used in Reads kmer index
struct KmerOccurrence {
  KmerOccurrence() {}
  KmerOccurrence(ReadId read_id, KmerOffset pos)
      : read_id(read_id), read_pos(pos) {}

  bool operator==(const KmerOccurrence& that) const {
    return read_id == that.read_id && read_pos == that.read_pos;
  }

  ReadId read_id;       // ID if the read this Kmer appears in.
  KmerOffset read_pos;  // Position in the read.
};

// Store read alignment to haplotype
// redacted
// renamed. That would be better because default position = 0 is a valid
// position.
struct ReadAlignment {
  ReadAlignment() : position(0), cigar(""), score(0) {}

  ReadAlignment(uint16_t position_param, const string& cigar_param,
                int score_param)
      : position(position_param), cigar(cigar_param), score(score_param) {}

  bool operator==(const ReadAlignment& that) const {
    return score == that.score && position == that.position &&
           cigar == that.cigar;
  }

  void reset() {
    score = 0;
    position = 0;
    cigar = "";
  }

  uint16_t position;
  string cigar;
  int score;
};

// Cigar operation is defined by operation type and it's length.
struct CigarOp {
  CigarOp()
      : operation(
      nucleus::genomics::v1::CigarUnit_Operation_OPERATION_UNSPECIFIED),
        length(0) {}
  CigarOp(CigarUnit::Operation op, int len) : operation(op), length(len) {}

  CigarUnit::Operation operation;
  int length;
};

// Stores alignment scores for all reads to one haplotype.
// haplotype_score = sum of all aligned read scores, where read score
// is a number of matched bases. This way we can order haplotypes by number of
// supporting reads.
struct HaplotypeReadsAlignment {
  HaplotypeReadsAlignment() : haplotype_index(0), haplotype_score(0) {}
  HaplotypeReadsAlignment(
      size_t haplotype_index,
      int score,
      const std::vector<ReadAlignment>& read_alignment_scores)
        : haplotype_index(haplotype_index), haplotype_score(score) {
    this->read_alignment_scores.assign(read_alignment_scores.begin(),
                                       read_alignment_scores.end());
  }

  // halpotype index in haplotypes member of FastPassAligner class
  size_t haplotype_index;

  // sum of all aligned read scores. Each read's
  // haplotype_score = number of matched bases.
  int haplotype_score;

  // Simple alignment haplotype_score for each read.
  std::vector<ReadAlignment> read_alignment_scores;

  // Cigar string for haplotype against reference alignment.
  string cigar;

  // Cigar as a list of CigarOp
  std::list<CigarOp> cigar_ops;

  // Hypolotype to reference position.
  uint64_t ref_pos;

  // Map of shifts that allow to easily calculate read to reference position
  // from read to haplotype position. read pos = read_to_haplotype_pos
  //    + haplotype_to_ref_position
  //    + hap_to_ref_positions_map[read_to_haplotype_pos]
  std::vector<int> hap_to_ref_positions_map;
};

using KmerIndexType =
std::unordered_map<tensorflow::StringPiece, std::vector<KmerOccurrence>,
                   tensorflow::StringPieceHasher>;

using ReadsVectorType = std::vector<nucleus::genomics::v1::Read>;

// Align a set of reads to a target sequence.
// This class is intended for realigning reads to haplotypes (graph paths)
// generated by DeBrujn graph. Since graph's paths are constructed from the same
// reads all "correct" paths should have almost perfect alignment with exception
// of reads with sequencing errors. Alignment is done by using reads kmer index.
// Using the index we can quickly (in O(n)) find all reads overlapping all kmers
// in a given haplotype. Haplotype is quickly discarded If we cannot get a 100%
// coverage using read kmer index. Details of 1st step: For each kmer in
// haplotype:
//  1. Find all reads overlapping given kmer
//  2. Compare each read with target sequence. Keep the read if we can have (x)M
//  alignment, allowing up to 3 mismatches, where
//      x is the length of the read. This operation can be very fast if we
//      encode each base with 2 bits. We may use special instructions that count
//      set bits in one cycle.
//  3. Keep only reads that have perfect alignment.
//  4. Due to sequencing errors some reads may be left unaligned. For those read
//  we can use a SW
//      alignment. Number of such reads should be very low (less than error rate
//      because only structural errors may contribute)
// Next step we need to align haplotype to reference.
// Next step using haplotype to reference alignment we can correct all reads
// alignment for the final result.
class FastPassAligner {
 public:
  void set_reference(const string& reference);
  void set_reads(const std::vector<string>& reads);
  std::vector<string> get_reads() const { return reads_; }
  void set_ref_start(const string& chromosome, uint64_t position);
  void set_haplotypes(const std::vector<string>& haplotypes);
  void set_score_schema(uint8_t match_score, uint8_t mismatch_penalty,
                        uint8_t gap_opening_penalty,
                        uint8_t gap_extending_penalty);
  uint8_t get_match_score() const { return match_score_; }
  uint8_t get_mismatch_penalty() const {return mismatch_penalty_;}
  void set_kmer_size(int kmer_size);
  void set_read_size(int read_size);
  void set_max_num_of_mismatches(int max_num_of_mismatches);
  void set_is_debug(bool is_debug) { debug_out_ = is_debug; }
  void set_debug_read_id(int read_id) { debug_read_id_ = read_id; }

  // Align reads to the reference by first aligning reads to haplotypes and
  // then by merging haplotype to reference cigars and reads to haplotype
  // cigars.
  // This function is an entry point for FastPassAligner.
  std::unique_ptr<ReadsVectorType> AlignReads(
      const ReadsVectorType& reads_param);

  // Build K-mer index for all reads.
  void BuildIndex();

  KmerIndexType GetKmerIndex() const { return kmer_index_; }

  // Align all reads to a haplotype using fast pass alignment.
  void FastAlignReadsToHaplotype(
      const string& haplotype, int* haplotype_score,
      std::vector<ReadAlignment>* haplotype_read_alignment_scores);

 private:
  // Reference sequence for the window
  string reference_;

  // Chromosome name for the window start position
  string region_chromosome_;

  // Position in chromosome for the window start position
  uint64_t region_position_in_chr_;

  // Vector of haploptypes sequences
  std::vector<string> haplotypes_;

  std::vector<HaplotypeReadsAlignment> read_to_haplotype_alignments_;

  // index of reads. Allows to find all reads that contain a given k-mer
  // and their align position.
  KmerIndexType kmer_index_;

  // Vector of reads that need to be realigned
  std::vector<string> reads_;

  // K-mer size that is used for indexing input reads
  int kmer_size_;

  // Expected read length. It is needed for sanity checking. Actual reads may
  // be different sizes. Although, actual read sizes should be close to
  // read_size.
  int read_size_;

  int max_num_of_mismatches_ = 2;

  // Four parameters below specify alignment scoring schema.
  uint8_t match_score_ = 4;
  uint8_t mismatch_penalty_ = 6;
  uint8_t gap_opening_penalty_ = 8;
  uint8_t gap_extending_penalty_ = 1;

  // SSW aligner instance
  std::unique_ptr<Aligner> ssw_aligner_;

  // These attributes allow debug output.
  bool debug_out_ = false;
  int debug_read_id_ = 0;

  void InitSswLib();

  void FastAlignReadsToHaplotypes();

  void AlignReadsToHaplotypes(uint16_t score_threshold);

  void AlignHaplotypesToReference();

  // Changes alignment for each read that we could realign.
  // realigned_reads are eventually passed to Python wrapped in unique_ptr
  // as per clif requirements.
  void RealignReadsToReference(
      const std::vector<nucleus::genomics::v1::Read>& reads,
      std::unique_ptr<std::vector<nucleus::genomics::v1::Read>>
        realigned_reads);

  void AddReadToIndex(const string& read, ReadId read_id);

  void AddKmerToIndex(tensorflow::StringPiece kmer,
                      ReadId read_id,
                      KmerOffset pos);

  int FastAlignStrings(const tensorflow::StringPiece& s1,
                       const tensorflow::StringPiece& s2,
                       int max_mismatches, int* num_of_mismatches) const;

  void UpdateBestHaplotypes(
      size_t haplotype_index, int haplotype_score,
      const std::vector<ReadAlignment>& current_read_scores);

  void CalculatePositionMaps();

  bool GetBestReadAlignment(size_t readId, int* best_hap_index) const;

  void CalculateCigarForRead(
      size_t read_index,
      std::list<CigarOp>* read_to_ref_cigar) const;

  static void CigarStringToVector(const string& cigar,
                                  std::list<CigarOp>* cigar_ops);
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_FAST_PASS_ALIGNER_H_
