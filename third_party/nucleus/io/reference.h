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

#include <vector>

#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/vendor/statusor.h"
#include "tensorflow/core/platform/types.h"

using tensorflow::string;
using tensorflow::int64;

namespace nucleus {

class GenomeReference {
 public:
  GenomeReference(const GenomeReference&) = delete;
  GenomeReference& operator=(const GenomeReference&) = delete;

  virtual ~GenomeReference() {}

  // Gets the path to the fasta file used to create this index.
  virtual const string& FastaPath() const = 0;

  // Gets a human-readable string describing the properties of this
  // GenomeReference.
  virtual string Info() const = 0;

  virtual const std::vector<nucleus::genomics::v1::ContigInfo>& Contigs()
      const = 0;

  // Gets the number of contigs in this reference.
  virtual int NContigs() const;

  // Returns true if reference has a contig named contig_name.
  virtual bool HasContig(const string& contig_name) const;

  // Gets the total number of basepairs across all contigs.
  virtual int64 NTotalBasepairs() const;

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

  // Returns true iff the Range chr:start-end is a valid interval on chr and chr
  // is a known contig in this reference.
  bool IsValidInterval(const nucleus::genomics::v1::Range& range) const;

  // Close the underlying resource descriptors.
  virtual tensorflow::Status Close() { return tensorflow::Status::OK(); }

  // This no-op function is needed only for Python context manager support.  Do
  // not use it!
  virtual tensorflow::Status PythonEnter() const {
    return tensorflow::Status::OK();
  }

 protected:
  // Default constructor for subclasses.
  GenomeReference() {}
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_REFERENCE_H_
