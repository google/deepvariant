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

// Compute AlleleCounts over an interval on the genome given the reads that
// align somewhere within that interval.
//
#ifndef LEARNING_GENOMICS_DEEPVARIANT_ALLELECOUNTER_H_
#define LEARNING_GENOMICS_DEEPVARIANT_ALLELECOUNTER_H_

#ifndef FRIEND_TEST
#define FRIEND_TEST(test_case_name, test_name)\
friend class test_case_name##_##test_name##_Test
#endif

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using std::string;

// Summarizes the counts of all of the distinct alleles present in allele_count.
//
// Computes and returns a vector of Allele objects, one for each distinct Allele
// across all reads in allele_count. Each allele in the vector has its count
// field set to the number of reads that carried that allele.
std::vector<Allele> SumAlleleCounts(const AlleleCount& allele_count,
                                    bool include_low_quality = false);

// Summarizes the counts of all of the distinct alleles present in allele_count
// for one position combined for all DeepTrio samples. Effectively this function
// merges allele_count from all DeepTrio samples.
// This function is similar to SumAlleleCounts(const AlleleCount& allele_count)
std::vector<Allele> SumAlleleCounts(absl::Span<const AlleleCount> allele_counts,
                                    bool include_low_quality = false);

// Gets the total count of observed alleles in this allele_count, which is the
// sum of the observed non-reference alleles in read_alleles + the total number
// of reference supporting reads.
int TotalAlleleCounts(const AlleleCount& allele_count,
                      bool include_low_quality = false);

// Gets the total count of observed alleles in allele_count from all DeepTrio
// samples, which is the sum of the observed non-reference alleles in
// read_alleles + the total number of reference supporting reads.
int TotalAlleleCounts(const std::vector<AlleleCount>& allele_counts,
                      bool include_low_quality = false);

// Binary search for allele index by position.
int AlleleIndex(absl::Span<const AlleleCount> allele_counts, int64_t pos);

// Represents an Allele observed in a read at a specific position in our
// interval. Supports the concept that the site should be skipped but still
// needs to be represented in a data processing chain. ReadAlleles marked as
// skip() won't be added to the AlleleCount at position but since they are
// explicitly represented can be reasoned about relative to other ReadAlleles.
// The pos() of a ReadAllele is an offset in our interval where this Allele
// should be added. For example, if pos() == 1, the Allele implied by this
// ReadAllele would be added to the second AlleleCount of this AlleleCounter.
// The position of the ReadAllele can be < 0 and beyond the length of the
// interval, indicating that some read that overlaps our interval carries an
// Allele but it shouldn't be added to our AlleleCounts.
class ReadAllele {
 public:
  ReadAllele() = default;

  // Creates a ReadAllele with position, bases, and type.
  ReadAllele(int position, absl::string_view bases, const AlleleType& type,
             bool is_low_quality = false, uint8_t mapping_quality = 0,
             uint8_t avg_base_quality = 0)
      : position_(position),
        bases_(bases),
        type_(type),
        low_quality_allele_(is_low_quality),
        mapping_quality_(mapping_quality),
        avg_base_quality_(avg_base_quality) {}

  // Gets the position of this ReadAllele. Can be < 0 or >= IntervalLength(),
  // indicating that the ReadAllele refers to a position outside of the
  // interval.
  int position() const { return position_; }

  // Returns true if this ReadAllele should be skipped.
  bool skip() const { return position_ == kInvalidPosition; }

  // Gets the bases of the allele observed at this position.
  const string& bases() const { return bases_; }

  // Gets the type of the allele observed at this position.
  const AlleleType& type() const { return type_; }

  bool is_low_quality() const { return low_quality_allele_; }

  uint8_t mapping_quality() const { return mapping_quality_; }

  uint8_t avg_base_quality() const { return avg_base_quality_; }

 private:
  static constexpr int kInvalidPosition = -1;

  const int position_ = kInvalidPosition;
  const string bases_ = "";
  const AlleleType type_ = AlleleType::UNSPECIFIED;
  bool low_quality_allele_ = false;
  uint8_t mapping_quality_ = 0;
  uint8_t avg_base_quality_ = 0;
};

// Workhorse class to compute AlleleCounts over an interval on the genome.
//
// AlleleCounter works roughly as follows:
//
// AlleleCounter counter = AlleleCounter(ref, interval)
//
// // Add reads that overlap our interval
// for ( read : reads_that_overlap_interval ) {
//   counter.Add(read)
//
// // Tell the counter that we are done adding reads
// counter.Finalize()
// vector<AlleleCount> counter.Counts()
//
// The returned vector of AlleleCounts has an AlleleCount for every position in
// the interval.  Each AlleleCount contains a position, the reference base at
// that position, and a repeated field of Allele protos which contain the bases,
// type, and counts of the alleles. Each observed allele comes from the
// alignment of a read at that position in the genome.  For example, if we have:
//
//   pos:  0123
//   ref:  ACGT
//   read:  CG
//
// We could produce AlleleCounts with 0 observed alleles at position 0, one C
// at position 1, one G at position 2, and no alleles at position 3.
//
// This becomes more complex when insertion and deletion alleles occur, so that
// you have:
//
//   pos:  0123
//   ref:  ACGT
//   read:  C-T
//
// Would actually produce an allele CG deletion at position 1, indicating that
// the reference base was C and that a G base was deleted. No allele count would
// occur at G (after all it was never observed), and there would be 1 T at
// position 3.  The AlleleCount proto and associated other protos provide
// additional details.
//
// The AlleleCount's generated by adding each read simply sum up independently
// with multiple reads, which is a very attractive property of the AlleleCount
// representation.
//
// Note that this code can diverge from the left-alignment requirement of VCF /
// variant protos when the input read cigars are themselves not  left aligned.
// For example, suppose we have:
//
// ref:    TAAAC
// sample: TAAC
//
// Ideally an aligner would place the DELETION cigar elements for any read
// covering this site to the left-most position:
//
// ref:    TAAAC
// read1:  T-AAC [cigar = 1M1D3M, ideal]
// read2:  TA-AC [cigar = 2M1D2M, pretty far from ideal, but equivalent]
//
// This code doesn't try to do anything clever by left-aligning CIGAR elements
// in order to fix this problem. This is largely ok because (1) the standard
// aligner (BWA) does in fact do consistent left alignment and (2) we anticipate
// feeding this AlleleCounter reads that have been assembled into a consistent
// alignment across all of the reads in a way that is impossible for a read-by-
// read aligner to do. So beware with the aligner you use, unless you've cleaned
// up the reads in some way (left aligning the cigars or just doing assembly).
//
// It is ok to send reads that only partially overlap the interval; the
// AlleleCounter will only add counts from the part that overlap the interval.
//
// This code assumes that the reference genome and the reads have only upper
// case bases. By construction our GenomeReference will not have lower case
// bases. Lower case bases are allowed in read sequences (see
// https://samtools.github.io/hts-specs/SAMv1.pdf),
// so it may be worth upper casing the read sequence at some point to make this
// code more robust.
//
// AlleleCounter performs some specialized logic to deal with non-canonical
// bases that might occur in the reference and/or the reads. By canonical we
// mean bases that are one of {A,C,G,T}. The basic constraint is that we will
// not add Allele's that contain non-canonical bases. So if we have a read that
// has (e.g.) an N, we will not produce Allele objects that contain the N. So
// if it doesn't match the reference at a site, it will be skipped as a
// SUBSTITUTION allele, and if it were part of an insertion, that entire
// insertion will be skipped. We do support AlleleCount objects that have a N
// base as its reference base, and those bases get alleles just as a normal base
// would. Downstream processing should look at the reference base of the
// AlleleCounts produced by this code for non-canonical reference bases and
// handle if special handling is needed. Finally, if a read's deletion CIGAR
// spans across part of the reference genome that has an N, the corresponding
// DELETION allele will be dropped.
class AlleleCounter {
 public:
  // Creates an AlleleCounter.
  //
  // The counter will track AlleleCounts over the interval range, using ref to
  // get reference data in that interval, using the options to determine which
  // reads and bases are counted.
  // Candidate_positions parameter is optional. When it is present allele
  // counter saves ref alleles for positions found in this vector.
  //
  // The GenomeReference must be available throughout the lifetime of this
  // AlleleCounter object.
  AlleleCounter(const nucleus::GenomeReference* ref,
                const nucleus::genomics::v1::Range& range,
                const std::vector<int>& candidate_positions,
                const AlleleCounterOptions& options);

  // An alternative constructor that allows to use a wider reference region for
  // allele counter. This is needed for read normalization for those reads that
  // only partially overlap allele counter region.
  AlleleCounter(const nucleus::GenomeReference* ref,
                const nucleus::genomics::v1::Range& range,
                const nucleus::genomics::v1::Range& full_range,
                const std::vector<int>& candidate_positions,
                const AlleleCounterOptions& options);

  // This Init is used by unit tests only.
  static AlleleCounter* InitFromAlleleCounts(
      absl::Span<const AlleleCount> allele_counts);

  // Adds the alleles from read to our AlleleCounts. This method is also called
  // by NormalizeAndAdd. In that case allele counts are created using a
  // normalized cigar and update read alignment position passed as an optional
  // parameter.
  void Add(const nucleus::genomics::v1::Read& read, absl::string_view sample,
           const std::vector<nucleus::genomics::v1::CigarUnit>* cigar_to_use =
               nullptr,
           int read_shift = 0);

  // Wrapper around Add() that normalize the input read first and then calls
  // Add().
  void NormalizeAndAdd(
      const nucleus::genomics::v1::Read& read, absl::string_view sample,
      std::unique_ptr<std::vector<nucleus::genomics::v1::CigarUnit>>&
          norm_cigar,
      int& read_shift);

  // Python wrapper around NormalizeAndAdd. It allows to avoid serialization of
  // protos when calling from Python.
  std::unique_ptr<std::vector<nucleus::genomics::v1::CigarUnit>>
  NormalizeAndAddPython(const nucleus::ConstProtoPtr<
                            const nucleus::genomics::v1::Read>& wrapped,
                        const string& sample, int* read_shift) {
    auto norm_cigar =
        std::make_unique<std::vector<nucleus::genomics::v1::CigarUnit>>(
            std::vector<nucleus::genomics::v1::CigarUnit>());
    NormalizeAndAdd(*(wrapped.p_), sample, norm_cigar, *read_shift);
    return norm_cigar;
  }

  // Simple wrapper around Add() that allows us to efficiently pass large
  // protobufs in from Python. Simply unwraps the ConstProtoPtr objects and
  // calls Add(read).
  void AddPython(const nucleus::ConstProtoPtr<
                     const nucleus::genomics::v1::Read>& wrapped,
                 const string& sample) {
    Add(*(wrapped.p_), sample, nullptr);
  }

  // Gets the options in use by this AlleleCounter
  const AlleleCounterOptions& Options() const { return options_; }

  // Gets the interval we are counting alleles over.
  const nucleus::genomics::v1::Range& Interval() const { return interval_; }

  // Gets the interval overlapping all the reads.
  const nucleus::genomics::v1::Range& ReadsInterval() const {
    return reads_interval_;
  }

  // Returns the number of basepairs in our interval.
  int64_t IntervalLength() const { return interval_.end() - interval_.start(); }

  // Return the number of basepairs of the interval overlapping all the reads.
  int64_t ReadsIntervalLength() const {
    return reads_interval_.end() - reads_interval_.start();
  }

  // Gets the completed AlleleCounts over this counter's interval.
  //
  // This function returns a vector of AlleleCount objects, one for each
  // basepair in our interval, filled in according to the reads that have been
  // added via calls to Add*() routines.
  //
  // Calling this routine finalizes the AlleleCounts(), so that subsequent calls
  // to Add*() will fail.
  const std::vector<AlleleCount>& Counts() const { return counts_; }

  // Similar to Counts() function but returns a lighter-weight summary proto.
  //
  // This function has all of the behavior of calling Counts() but instead of
  // returning the heavy-weight AlleleCount proto this returns a simpler proto.
  // See the proto description for more information about the proto fields.
  std::vector<AlleleCountSummary> SummaryCounts(int left_padding = 0,
                                                int right_padding = 0) const;

  // How many reads have been added to this counter?
  int NCountedReads() const { return n_reads_counted_; }

  // Constructs a unique string key for this read. The key is the concatenation
  // of fragment_name, "/", and read_number.
  string ReadKey(const nucleus::genomics::v1::Read& read);

 private:
  // This constructor is used for unit testing only.
  AlleleCounter();

  // Initialize allele counter.
  void Init();

  // Helper function to get the reference bases between offsets rel_start
  // (inclusive) and rel_end (exclusive). The offsets are both relative to our
  // interval, so rel_start = 0 means the first base in our interval.  Because
  // of their relative nature, it's meaningful to ask for offsets that are
  // negative (gets bases before the start of the interval) or offsets that are
  // longer than the interval. Will return the empty string if the actual
  // genomic coordinates implied by the offsets aren't all on the chromosome.
  string RefBases(int64_t rel_start, int64_t len);

  // Returns True if ref_offset (where 0 indicates the first position in the
  // interval, which could be base 1234 in genomic coordinates, for example), is
  // within our interval. This means that ref_offset >= 0 and ref_offset <
  // IntervalLength().
  bool IsValidRefOffset(int ref_offset) {
    return ref_offset >= 0 && ref_offset < ReadsIntervalLength();
  }

  // Returns True if interval_offset is within our allele counter interval.
  bool IsValidIntervalOffset(int interval_offset) {
    return interval_offset >= 0 && interval_offset < IntervalLength();
  }

  // Gets the base before read_offset in read, or if that would be before the
  // start of the read (i.e., read_offset == 0) then return the previous base on
  // the reference genome (at interval_offset - 1).
  string GetPrevBase(const nucleus::genomics::v1::Read& read, int read_offset,
                     int interval_offset);

  // Creates a ReadAllele for an indel (type based on cigar) from read starting
  // at read_offset position in the read to the AlleleCount at interval_offset.
  // Does all of the necessary quality checks to ensure we only add good bases
  // to the our AlleleCounts, as well as manages the complexity of determining
  // the correct allele to add. May return a ReadAllele marked as skip() if the
  // implied allele isn't valid for some reason (e.g., bases are too low
  // quality).
  ReadAllele MakeIndelReadAllele(
      const nucleus::genomics::v1::Read& read, int interval_offset,
      int ref_offset, int read_offset,
      const nucleus::genomics::v1::CigarUnit& cigar);

  // Adds the ReadAlleles in to_add to our AlleleCounts.
  void AddReadAlleles(const nucleus::genomics::v1::Read& read,
                      absl::string_view sample,
                      absl::Span<const ReadAllele> to_add);

  // Normalize cigar by shifting INDELs in the middle of a repeat all the way
  // to the left. As a result of shifting two INDELs may become merged. Merged
  // INDEL may become non-normalized so the process is repeated up to 10 times.
  // If INDEL is shifted all the way to the beginning of the read then this
  // INDEL is removed and read alignment position has to be shifted.
  bool NormalizeCigar(const absl::string_view read_seq, int interval_offset,
                      std::vector<nucleus::genomics::v1::CigarUnit>& cigar,
                      int& read_shift) const;

  // Helper function used in NormalizeCigar function. Returns true if deletion
  // operation can be shifted left preserving the alignment.
  bool CanDelBeShifted(
      const absl::string_view read_seq,
      std::vector<nucleus::genomics::v1::CigarUnit>::const_iterator cigar_elt,
      int read_offset,
      int interval_offset,
      int op_len) const;

  // Helper function used in NormalizeCigar function. Returns true if insertion
  // operation can be shifted left preserving the alignment.
  bool CanInsBeShifted(
      const absl::string_view read_seq,
      std::vector<nucleus::genomics::v1::CigarUnit>::const_iterator cigar_elt,
      int read_offset,
      int interval_offset,
      int op_len) const;

  // Our GenomeReference, which we use to get information about the reference
  // bases in our interval.
  const nucleus::GenomeReference* const ref_;

  // The interval chr from start (0-based, inclusive) to end (0-based,
  // exclusive) describing where we are counting on the genome. We will produce
  // one AlleleCount for each base in interval, from start to end (exclusive).
  const nucleus::genomics::v1::Range interval_;

  // The interval chr from start (0-based, inclusive) to end (0-based,
  // exclusive) of the available ref bases. By default this interval is equal to
  // to the interval_. If read normalization is needed then reads_interval_ may
  // be an etentsion of interval_. reads_interval_ spans from the first position
  // of the very first read to the last position of the last read.
  const nucleus::genomics::v1::Range reads_interval_;

  // Vector of potential candidate positions. Ref alleles are stored only for
  // positions found in this vector. This functionality is optional and
  // governed by track_ref_reads flag.
  std::vector<int> candidate_positions_;

  // The options that are controlling how we count reads.
  const AlleleCounterOptions options_;

  // The number of reads we've added to this interval.
  int n_reads_counted_ = 0;

  // Our AlleleCount objects, one for each base in our interval, in order.
  std::vector<AlleleCount> counts_;

  // The reference bases covering our interval;
  const string ref_bases_;

  // Following tests call protected method NormalizeCigar.
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarDel);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarIns);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarInsDel);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarInsertAtTheEnd);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarTwoDelsMerged);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarDelInsMerged);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarInsShiftedToEdge);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarInsShiftedAllTheWayToSoftClip);
  FRIEND_TEST(AlleleCounterTest, NormalizeCigarDelInsMergedNoShift);
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_ALLELECOUNTER_H_
