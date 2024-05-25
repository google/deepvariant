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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_FAST_PASS_ALIGNER_H_
#define LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_FAST_PASS_ALIGNER_H_

#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include "deepvariant/realigner/ssw.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/node_hash_map.h"
#include "absl/memory/memory.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "re2/re2.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Position;

// TODO Add ChromosomePosition type instead of uint64_t

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
// TODO position could be of type KmerOffset, but then it needs to be
// renamed. That would be better because default position = 0 is a valid
// position.
struct ReadAlignment {
  static constexpr uint16_t kNotAligned = std::numeric_limits<uint16_t>::max();
  ReadAlignment() : position(kNotAligned), cigar(""), score(0) {}

  ReadAlignment(uint16_t position_param, absl::string_view cigar_param,
                int score_param)
      : position(position_param), cigar(cigar_param), score(score_param) {}

  bool operator==(const ReadAlignment& that) const {
    return score == that.score && position == that.position &&
           cigar == that.cigar;
  }

  void reset() {
    score = 0;
    position = kNotAligned;
    cigar = "";
  }

  uint16_t position;
  string cigar;
  int score;
};

// Cigar operation is defined by operation type and its length.
struct CigarOp {
  CigarOp()
      : operation(nucleus::genomics::v1::CigarUnit::OPERATION_UNSPECIFIED),
        length(0) {}
  CigarOp(CigarUnit::Operation op, int len) : operation(op), length(len) {}

  bool operator==(const CigarOp& that) const {
    return operation == that.operation && length == that.length;
  }

  CigarUnit::Operation operation;
  int length;
};

// TODO HaplotypeReadsAlignment is not a good name for this structure.
// It contains not only read alignments to haplotype but also haplotype to
// reference alignment.
// Stores alignment scores for all reads to one haplotype.
// haplotype_score = sum of all aligned read scores, where read score
// is a number of matched bases. This way we can order haplotypes by number of
// supporting reads.
struct HaplotypeReadsAlignment {
  HaplotypeReadsAlignment() : haplotype_index(0), haplotype_score(0) {}
  HaplotypeReadsAlignment(size_t haplotype_index, int score,
                          absl::Span<const ReadAlignment> read_alignment_scores)
      : haplotype_index(haplotype_index), haplotype_score(score) {
    this->read_alignment_scores.assign(read_alignment_scores.begin(),
                                       read_alignment_scores.end());
  }

  bool operator==(const HaplotypeReadsAlignment& that) const {
    return haplotype_index == that.haplotype_index &&
           haplotype_score == that.haplotype_score &&
           read_alignment_scores == that.read_alignment_scores &&
           cigar == that.cigar && cigar_ops == that.cigar_ops &&
           is_reference == that.is_reference &&
           hap_to_ref_positions_map == that.hap_to_ref_positions_map;
  }

  bool operator<(const HaplotypeReadsAlignment& that) const {
    return haplotype_score < that.haplotype_score;
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

  // If true the haplotype is a reference.
  bool is_reference;
};

// Calculate a shift for each position of a haplotype from a haplotype to
// reference alignment.
void SetPositionsMap(size_t haplotype_size,
                     HaplotypeReadsAlignment* hyplotype_alignment);

void MergeCigarOp(const CigarOp& op, int read_len, std::list<CigarOp>* cigar);

using KmerIndexType =
    absl::flat_hash_map<absl::string_view, std::vector<KmerOccurrence>>;

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
  void set_normalize_reads(bool normalize_reads) {
    normalize_reads_ = normalize_reads;
  }
  uint8_t get_match_score() const { return match_score_; }
  uint8_t get_mismatch_penalty() const { return mismatch_penalty_; }
  void set_options(const AlignerOptions& options);
  void set_is_debug(bool is_debug) { debug_out_ = is_debug; }
  void set_debug_read_id(int read_id) { debug_read_id_ = read_id; }
  void set_ref_prefix_len(int ref_prefix_len) {
    ref_prefix_len_ = ref_prefix_len;
  }
  void set_ref_suffix_len(int ref_suffix_len) {
    ref_suffix_len_ = ref_suffix_len;
  }
  int16_t get_ssw_alignment_score_threshold() const {
    return ssw_alignment_score_threshold_;
  }

  // Align reads to the reference by first aligning reads to haplotypes and
  // then by merging haplotype to reference cigars and reads to haplotype
  // cigars.
  // This function is an entry point for FastPassAligner.
  std::unique_ptr<std::vector<nucleus::genomics::v1::Read>> AlignReads(
      const std::vector<nucleus::genomics::v1::Read>& reads_param);

  // Build K-mer index for all reads.
  void BuildIndex();

  KmerIndexType GetKmerIndex() const { return kmer_index_; }

  // Align all reads to a haplotype using fast pass alignment.
  void FastAlignReadsToHaplotype(
      absl::string_view haplotype, int* haplotype_score,
      std::vector<ReadAlignment>* haplotype_read_alignment_scores);

  // Align reads to haplotypes using SSW library. Only reads that could not
  // be aligned with FastAlignReadsToHaplotype are aligned here.
  // Only alignment with better than score_threshold score are kept.
  void SswAlignReadsToHaplotypes(uint16_t score_threshold);

  // Initialize SSW library.
  void InitSswLib();

  void SswSetReference(const string& reference);

  // Align target sequence to the reference using SSW library. InitSswLib must
  // be called prior to calling SswAlign.
  Alignment SswAlign(const string& target) const;

  // Align all haplotypes to the reference using SSW library.
  void AlignHaplotypesToReference();

  const std::vector<HaplotypeReadsAlignment>& GetReadToHaplotypeAlignments()
      const {
    return read_to_haplotype_alignments_;
  }

  // Finds the best alignment by iterating all haplotype alignments.
  // TODO Can be optimized if we store the best alignment with the read.
  bool GetBestReadAlignment(size_t readId, int* best_hap_index) const;

  // Calculate a read alignment by merging read to haplotype and haplotype to
  // reference alignments.
  // 1. Extracts a portion of haplotype to reference cigar that overlaps
  //    positions of read to haplotype alignment.
  // 2. Iterates through each cigar operation for both alignments, merging 2
  //    operations of length 1 at a time.
  // 3. Different logic is implemented for each type of merges: =:=, DEL:=,
  //    =:DEL, INS:=, =:INS, DEL:DEL, INS:INS, DEL:INS
  void CalculateReadToRefAlignment(
      size_t read_index,
      const ReadAlignment& read_to_haplotype_alignment,
      const std::list<CigarOp>& haplotype_to_ref_cigar_ops_input,
      std::list<CigarOp>* read_to_ref_cigar_ops) const;

  bool IsAlignmentNormalized(
    const std::list<CigarOp>& cigar,
    int ref_offset,
    absl::string_view read_sequence) const;

  // Replace each read alignment with the better one (if available), otherwise
  // an original alignment is preserved.
  // realigned_reads param is eventually passed to Python wrapped in unique_ptr
  // as per clif requirements.
  void RealignReadsToReference(
      absl::Span<const nucleus::genomics::v1::Read> reads,
      std::unique_ptr<std::vector<nucleus::genomics::v1::Read>>*
          realigned_reads);

  void CalculateSswAlignmentScoreThreshold();

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
  int kmer_size_ = 32;

  // Expected read length. It is needed for sanity checking. Actual reads may
  // be different sizes. Although, actual read sizes should be close to
  // read_size.
  int read_size_;

  int max_num_of_mismatches_ = 2;

  // Only read alignments with the score greater than
  // ssw_alignment_score_threshold_ are kept. ssw_alignment_score_threshold_
  // is calculated from similarity_threshold_ and read_size.
  // Most of the reads should almost perfectly align to haplotypes. Read may
  // not align to haplotype perfectly if:
  //  - there is a sequencing error;
  //  - read does not really come from this region;
  //  - haplotype does not represent a real target genome sequence.
  int16_t ssw_alignment_score_threshold_;

  // Four parameters below specify alignment scoring schema.
  uint8_t match_score_ = 4;
  uint8_t mismatch_penalty_ = 6;
  uint8_t gap_opening_penalty_ = 8;
  uint8_t gap_extending_penalty_ = 1;
  bool force_alignment_ = false;

  // Threshold is calculated from this flag using the following formula.
  // score_threshold = match_score_ * read_size_ * <similarity_threshold_>
  //    - mismatch_penalty_ * read_size_ * (1 - <similarity_threshold_ >);
  // This threshold is used for alignments performed with SSW.
  double similarity_threshold_ = 0.85;

  // SSW aligner instance
  std::unique_ptr<Aligner> ssw_aligner_;

  // These attributes allow debug output.
  bool debug_out_ = false;
  int debug_read_id_ = 0;

  int ref_prefix_len_;
  int ref_suffix_len_;

  // Set to the same value of --normalize_reads flag from make_examples.
  bool normalize_reads_ = false;

  // Alingn reads to haplotypes by simply comparing strings. This way we will
  // be able align all the reads that are aligned to haplotypes w/o indels.
  void FastAlignReadsToHaplotypes();

  void AddReadToIndex(absl::string_view read, ReadId read_id);

  void AddKmerToIndex(absl::string_view kmer, ReadId read_id,
                      KmerOffset pos);

  int FastAlignStrings(absl::string_view s1, absl::string_view s2,
                       int max_mismatches, int* num_of_mismatches) const;

  void UpdateBestHaplotypes(
      size_t haplotype_index, int haplotype_score,
      const std::vector<ReadAlignment>& current_read_scores);

  // Update position map for each haplotype. Position map stores shifts for
  // each position of a haplotype in respect to haplotype to reference
  // alignment. This map helps to quickly calculate reference position for
  // a read from a read to haplotype alignment.
  void CalculatePositionMaps();
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_FAST_PASS_ALIGNER_H_
