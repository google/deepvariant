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

#ifndef THIRD_PARTY_NUCLEUS_IO_SAM_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_SAM_READER_H_

#include <memory>
#include <string>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/util/samplers.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/vendor/statusor.h"

namespace nucleus {

// Alias for the abstract base class for SAM record iterables.
using SamIterable = Iterable<nucleus::genomics::v1::Read>;

// A SAM/BAM/CRAM reader.
//
// SAM/BAM/CRAM files store information about next-generation DNA sequencing
// info:
//
// https://samtools.github.io/hts-specs/SAMv1.pdf
// https://samtools.github.io/hts-specs/CRAMv3.pdf
//
// These files are block-gzipped series of records. When aligned they are
// frequently sorted and indexed:
//
// http://www.htslib.org/doc/samtools.html
//
// This class provides methods to iterate through a BAM file or, if indexed,
// to also query() for only read overlapping a specific region on the genome.
//
// Uses the htslib C API for reading NGS reads (BAM, SAM, CRAM etc). For details
// of the API, see:
//
// https://github.com/samtools/htslib/tree/develop/htslib
//
// The objects returned by iterate() or query() are nucleus.genomics.v1.Read
// objects parsed from the SAM/BAM/CRAM records in the file. Currently all
// fields except the extended key/value maps in each BAM fields are parsed.
//
class SamReader : public Reader {
 public:
  // Creates a new SamReader reading from the SAM/BAM/CRAM file reads_path.
  //
  // reads_path must point to an existing SAM/BAM/CRAM formatted file (text SAM,
  // compressed or uncompressed BAM file, CRAM files in all sorts of flavors).
  //
  // ref_path can be "", in which case the argument is ignored, or must point
  // to an existing FASTA file. If not "" and the reads_path points to a CRAM
  // file, the CRAM_OPT_REFERNECE field will be set to this path so that the
  // CRAM decoder uses ref_path to decode the reference-compressed read
  // sequences in the CRAM file. Because many low-level IO routines (e.g. stat)
  // are currently directly used in the CRAM implementation in htslib, ref_path
  // must be on a local (e.g., POSIX accessible) mount point. File system access
  // provided by htslib plugins (e.g., S3) won't work.
  //
  // If the filetype is BAM/CRAM, this constructor will attempt to load a BAI or
  // CRAI index from file reads_path + '.bai' or reads_path (without the .bam
  // extension) + '.bai'; if the index is not found, attempts to Query will
  // fail.
  //
  // Returns a StatusOr that is OK if the SamReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<SamReader>> FromFile(
      const string& reads_path, const string& ref_path,
      const nucleus::genomics::v1::SamReaderOptions& options);

  static StatusOr<std::unique_ptr<SamReader>> FromFile(
      const string& reads_path,
      const nucleus::genomics::v1::SamReaderOptions& options) {
    return FromFile(reads_path, "", options);
  }

  ~SamReader();

  // Disable assignment/copy operations
  SamReader(const SamReader& other) = delete;
  SamReader& operator=(const SamReader&) = delete;

  // Gets all of the reads in this file in order.
  //
  // This function allows one to iterate through all of the reads in this
  // SAM/BAM/CRAM file in order.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction. Returns an OK status if the iterable can be
  // constructed, or not OK otherwise.
  StatusOr<std::shared_ptr<SamIterable>> Iterate() const;

  // Gets all of the reads that overlap any bases in range.
  //
  // This function allows one to iterate through all of the reads in this
  // SAM/BAM/CRAM file in order that overlap a specific interval on the genome.
  // The query operation is efficient in that the cost is O(n) for n elements
  // that overlap range, and not O(N) for N elements in the entire file.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction.
  //
  // If no index was loaded by the constructor a non-OK status value will be
  // returned.
  //
  // If range isn't a valid interval in this BAM file a non-OK status value will
  // be returned.
  StatusOr<std::shared_ptr<SamIterable>> Query(
      const nucleus::genomics::v1::Range& region) const;

  // Returns True if this SamReader loaded an index file.
  bool HasIndex() const { return idx_ != nullptr; }

  // Close the underlying resource descriptors. Returns a Status to indicate if
  // everything went OK with the close.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support.  Do
  // not use it! Returns a Status indicating whether the enter was successful.
  ::nucleus::Status PythonEnter() const { return ::nucleus::Status(); }

  bool KeepRead(const nucleus::genomics::v1::Read& read) const;

  const nucleus::genomics::v1::SamReaderOptions& options() const {
    return options_;
  }

  // Returns a SamHeader message representing the structured header information.
  const nucleus::genomics::v1::SamHeader& Header() const { return sam_header_; }

 private:
  // Private constructor; use FromFile to safely create a SamReader from a
  // file.
  SamReader(const string& reads_path,
            const nucleus::genomics::v1::SamReaderOptions& options, htsFile* fp,
            bam_hdr_t* header, hts_idx_t* idx);

  // Our options that control the behavior of this class.
  const nucleus::genomics::v1::SamReaderOptions options_;

  // A pointer to the htslib file used to access the SAM/BAM data.
  htsFile* fp_;

  // A htslib header data structure obtained by parsing the header of this BAM.
  bam_hdr_t* header_;

  // The htslib index data structure for our indexed BAM file. May be NULL if no
  // index was loaded.
  hts_idx_t* idx_;

  // The sam.proto SamHeader message representing the structured header
  // information.
  nucleus::genomics::v1::SamHeader sam_header_;

  // For downsampling reads.
  mutable FractionalSampler sampler_;
};

namespace sam_reader_internal {

// Returns false if Read does not satisfy all of the ReadRequirements.
bool ReadSatisfiesRequirements(
    const nucleus::genomics::v1::Read& read,
    const nucleus::genomics::v1::ReadRequirements& requirements);

}  // namespace sam_reader_internal

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_SAM_READER_H_
