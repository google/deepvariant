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
 *
 */

// Get basic information about a reference genome as well as make it cpu and
// memory efficient and scalable to get the reference bases for an interval on
// the genome.
//
// The GenomeReference provides the core functionality needed to use a reference
// genome for data processing and analyses tools:
//
//   -- Get information about the contigs (aka chromosomes) present the FASTA,
//      such as its name, description, and number of basepairs.
//   -- Efficiently lookup the sequence of bases in an interval in the reference
//      genome. For example, GetBases("chr1", 0, 10) gets the basepair sequence
//      from the first base to the ninth base on chr1. This function call has
//      cost roughly proportional to the size of the query interval, regardless
//      of its position in the original FASTA file.
//
// The code here makes some strong assumptions about what a client could want.
// It doesn't record the position of bases in the original FASTA, and it doesn't
// track line breaks, comments, and other features of the FASTA. It uppercases
// the basepair sequences, so complexity or other information encoded in the
// case of the bases is lost. It also ensures that all of the bases in the
// reference are either {A,C,G,T,N} by refusing to import reference sequences
// containing other characters. The code assumes that random accesses of
// reasonably small chunks of sequence is important at the expense of low-cost
// (but still reasonably efficient) iteration of all sequences in the FASTA.
#ifndef THIRD_PARTY_NUCLEUS_IO_REFERENCE_H_
#define THIRD_PARTY_NUCLEUS_IO_REFERENCE_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "absl/types/optional.h"
#include "htslib/faidx.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/text_reader.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/fasta.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/vendor/statusor.h"

namespace nucleus {

constexpr int INDEXED_FASTA_READER_DEFAULT_CACHE_SIZE = 64 * 1024;

// Alias for the abstract base class for FASTA record iterables, which
// corresponds to (name, sequence) pairs.
using GenomeReferenceRecord = std::pair<string, string>;
using GenomeReferenceRecordIterable = Iterable<GenomeReferenceRecord>;

class GenomeReference : public Reader {
 public:
  GenomeReference(const GenomeReference&) = delete;
  GenomeReference& operator=(const GenomeReference&) = delete;

  virtual ~GenomeReference() {}

  virtual const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const = 0;

  // Returns true if reference has a contig named contig_name.
  virtual bool HasContig(const string& contig_name) const;

  // Gets a vector of the contig names in this reference, in the order that they
  // appeared in the original FASTA file. For example, would return:
  //
  //    {"chrM", "chr1", "chr2"}
  //
  // for a reference with three contigs named chrM, chr1, and chr2.
  virtual std::vector<string> ContigNames() const;

  // Gets the metadata about a contig in the fasta, such as its name,
  // description, length, etc. If contig_name isn't found in this reference,
  // returns a value whose status is not ok().
  virtual StatusOr<const nucleus::genomics::v1::ContigInfo*> Contig(
      const string& contig_name) const;

  // Gets the basepairs in the FASTA file from Range range.
  //
  // This follows the Range convention of getting bases from start
  // (inclusive) to end (exclusive), both 0-based. That is,
  // GetBases(Range("chr1", 2, 4)) gets a string starting with the *3rd* base on
  // chr1 and extending through the *4th* base (excluding 4 as an offset).
  // If chr isn't present in this reference, start is invalid or end is beyond
  // the length of chr, returns a value whose status is not ok().
  virtual StatusOr<string> GetBases(
      const nucleus::genomics::v1::Range& range) const = 0;

  // Gets all of the FASTA records in this file in order.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction. Returns an OK status if the iterable can be
  // constructed, or not OK otherwise.
  virtual StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>> Iterate()
      const = 0;

  // Returns true iff the Range chr:start-end is a valid interval on chr and chr
  // is a known contig in this reference.
  bool IsValidInterval(const nucleus::genomics::v1::Range& range) const;

  // Close the underlying resource descriptors.
  virtual ::nucleus::Status Close() { return ::nucleus::Status(); }

  // This no-op function is needed only for Python context manager support.  Do
  // not use it!
  virtual ::nucleus::Status PythonEnter() const { return ::nucleus::Status(); }

 protected:
  // Default constructor for subclasses.
  GenomeReference() {}
};

// A FASTA reader backed by a htslib FAI index.
//
// FASTA files store information about DNA/RNA/Amino Acid sequences:
//
// https://en.wikipedia.org/wiki/FASTA_format
//
// This reader is specialized for the FASTA variant used in NGS analyses, which
// has a FAI index created by samtools that allows efficient query() operations
// to obtain the subsequence of the FASTA on a specific contig between a start
// and end offsets:
//
// http://www.htslib.org/doc/faidx.html
// http://www.htslib.org/doc/samtools.html [faidx section]
//
// The FASTA file can be optionally block-gzipped compressed.
//
// This class provides methods to iterate through a the FASTA records and to
// also query() for the bases spanning a specific region on the genome.
//
// Uses the htslib C API for reading the FASTA and FAI. For details of
// the API, see:
//
// https://github.com/samtools/htslib/tree/develop/htslib
//
// The objects returned by iterate() or query() are strings containing the
// bases, all upper-cased.
class IndexedFastaReader : public GenomeReference {
 public:
  // Creates a new GenomeReference backed by the FASTA file fasta_path.
  //
  // Returns this newly allocated IndexedFastaReader object, passing ownership
  // to the caller via a unique_ptr.
  //
  // htslib currently assumes that the FAI file is named fasta_path + '.fai',
  // so that file must exist and be readable by htslib.
  //
  // We maintain a single entry cache of the bases from the last FASTA fetch, to
  // reduce the number of file reads, which can be quite costly for remote
  // filesystems.  64K is the default block size for htslib faidx fetches, so
  // there is no penalty to rounding up all small access sizes to 64K.  The
  // cache can be disabled using `cache_size=0`.
  static StatusOr<std::unique_ptr<IndexedFastaReader>> FromFile(
      const string& fasta_path, const string& fai_path,
      const nucleus::genomics::v1::FastaReaderOptions& options,
      int cache_size_bases = INDEXED_FASTA_READER_DEFAULT_CACHE_SIZE);
  static StatusOr<std::unique_ptr<IndexedFastaReader>> FromFile(
      const string& fasta_path, const string& fai_path,
      int cache_size_bases = INDEXED_FASTA_READER_DEFAULT_CACHE_SIZE);

  ~IndexedFastaReader();

  // Disable copy and assignment operations
  IndexedFastaReader(const IndexedFastaReader& other) = delete;
  IndexedFastaReader& operator=(const IndexedFastaReader&) = delete;

  const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const override {
    return contigs_;
  }

  StatusOr<string> GetBases(
      const nucleus::genomics::v1::Range& range) const override;

  // Get the options controlling the behavior of this FastaReader.
  const nucleus::genomics::v1::FastaReaderOptions& Options() const {
    return options_;
  }

  StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>> Iterate()
      const override;

  // Close the underlying resource descriptors.
  ::nucleus::Status Close() override;

 private:
  // Allow iteration to access the underlying reader.
  friend class IndexedFastaReaderIterable;

  // Must use one of the static factory methods.
  IndexedFastaReader(const string& fasta_path, faidx_t* faidx,
                     const nucleus::genomics::v1::FastaReaderOptions& options,
                     int cache_size_bases);

  // Path to the FASTA file containing our genomic bases.
  const string fasta_path_;

  // Initialized htslib faidx_t index for our fasta. Effectively constant during
  // the life of this object, but htslib API doesn't allow us to mark this as
  // const.
  faidx_t* faidx_;

  // The options controlling the behavior of this FastaReader.
  const nucleus::genomics::v1::FastaReaderOptions options_;

  // A list of ContigInfo, each of which contains the information about the
  // contigs used by this BAM file.
  const std::vector<nucleus::genomics::v1::ContigInfo> contigs_;

  // Size, in bases, of the read cache.
  const int cache_size_bases_;

  // Cache of the last "small" read from the FASTA (<= kFastaCacheSize).
  mutable string small_read_cache_;

  // The range that is held in the cache, or "empty" if there is no range cached
  // yet.  Range must be <= kFastaCacheSize in length.
  mutable absl::optional<nucleus::genomics::v1::Range> cached_range_;
};

// A FASTA reader that is not backed by a htslib FAI index.
//
// FASTA files store information about DNA/RNA/Amino Acid sequences:
//
// https://en.wikipedia.org/wiki/FASTA_format
//
// This reader is for FASTA files that contain many small records and are
// explicitly not indexed. The FASTA files can be optionally block-gzipped
// compressed.
//
// This class provides methods to iterate through a the FASTA records but
// doesn't support query() for the bases spanning a specific region on the
// genome.
//
// The (name, bases) tuple returned by iterate() are strings containing the
// bases in uppercase.
class UnindexedFastaReader : public GenomeReference {
 public:
  // Creates a new GenomeReference backed by the FASTA file fasta_path.
  //
  // Returns this newly allocated UnindexedFastaReader object, passing ownership
  // to the caller via a unique_ptr.
  static StatusOr<std::unique_ptr<UnindexedFastaReader>> FromFile(
      const string& fasta_path);

  ~UnindexedFastaReader();

  // Disable copy and assignment operations
  UnindexedFastaReader(const UnindexedFastaReader& other) = delete;
  UnindexedFastaReader& operator=(const UnindexedFastaReader&) = delete;

  const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const override;

  StatusOr<string> GetBases(
      const nucleus::genomics::v1::Range& range) const override;

  StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>> Iterate()
      const override;

  // Close the underlying resource descriptors.
  ::nucleus::Status Close() override;

 private:
  // Allow iteration to access the underlying reader.
  friend class UnindexedFastaReaderIterable;

  // Must use one of the static factory methods.
  UnindexedFastaReader(std::unique_ptr<TextReader> text_reader);

  const std::vector<nucleus::genomics::v1::ContigInfo> contigs_;

  // Underlying file reader.
  std::unique_ptr<TextReader> text_reader_;
};

// An FASTA reader backed by in-memory ReferenceSequence protos.
//
// FASTA files store information about DNA/RNA/Amino Acid sequences:
//
// https://en.wikipedia.org/wiki/FASTA_format
//
//
// An InMemoryFastaReader provides the same API as GenomeReferenceFAI but
// doesn't fetch its data from an on-disk FASTA file but rather fetches the
// bases from an in-memory cache containing ReferenceSequence protos.
//
// In particular the GetBases(Range(chrom, start, end)) operation fetches bases
// from the tuple where chrom == chromosome, and then from the bases where the
// first base of bases starts at start. If start > 0, then the bases string is
// assumed to contain bases starting from that position in the region. For
// example, the record ('1', 10, 'ACGT') implies that
// GetBases(ranges.make_range('1', 11, 12)) will return the base 'C', as the 'A'
// base is at position 10. This makes it straightforward to cache a small region
// of a full chromosome without having to store the entire chromosome sequence
// in memory (potentially big!).
class InMemoryFastaReader : public GenomeReference {
 public:
  // Creates a new InMemoryFastaReader backed by ReferenceSequence protos.
  static StatusOr<std::unique_ptr<InMemoryFastaReader>> Create(
      const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
      const std::vector<nucleus::genomics::v1::ReferenceSequence>& seqs);

  // Disable copy and assignment operations
  InMemoryFastaReader(const InMemoryFastaReader& other) = delete;
  InMemoryFastaReader& operator=(const InMemoryFastaReader&) = delete;

  const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const override {
    return contigs_;
  }

  const std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence>&
  ReferenceSequences() const {
    return seqs_;
  }

  StatusOr<string> GetBases(
      const nucleus::genomics::v1::Range& range) const override;

  StatusOr<std::shared_ptr<GenomeReferenceRecordIterable>> Iterate()
      const override;

 private:
  // Allow iteration to access the underlying reader.
  friend class FastaFullFileIterable;

  // Must use one of the static factory methods.
  InMemoryFastaReader(
      const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
      std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence> seqs)
      : contigs_(contigs), seqs_(seqs) {}

  const std::vector<nucleus::genomics::v1::ContigInfo> contigs_;
  const std::unordered_map<string, nucleus::genomics::v1::ReferenceSequence>
      seqs_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_REFERENCE_H_
