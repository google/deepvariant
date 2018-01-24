/*
 * Copyright 2017 Google Inc.
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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_CORE_SAM_READER_H_
#define LEARNING_GENOMICS_DEEPVARIANT_CORE_SAM_READER_H_

#include "deepvariant/core/protos/core.pb.h"
#include "deepvariant/core/reader_base.h"
#include "deepvariant/core/samplers.h"
#include "deepvariant/vendor/statusor.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "deepvariant/core/genomics/range.pb.h"
#include "deepvariant/core/genomics/reads.pb.h"
#include "tensorflow/core/platform/types.h"

namespace learning {
namespace genomics {
namespace core {

using tensorflow::string;

// Alias for the abstract base class for SAM record iterables.
using SamIterable = Iterable<nucleus::genomics::v1::Read>;

// A SAM/BAM reader.
//
// SAM/BAM files store information about next-generation DNA sequencing info:
//
// https://samtools.github.io/hts-specs/SAMv1.pdf
//
// These files are block-gzipped series of records. When aligned they are
// frequently sorted and indexed:
//
// http://www.htslib.org/doc/samtools.html
//
// This class provides methods to iterate through a BAM file or, if indexed,
// to also query() for only read overlapping a specific region on the genome.
//
// Uses the htslib C API for reading NGS reads (BAM, SAM, etc). For details of
// the API, see:
//
// https://github.com/samtools/htslib/tree/develop/htslib
//
// The objects returned by iterate() or query() are nucleus.genomics.v1.Read
// objects parsed from the SAM/BAM records in the file. Currently all fields
// except the extended key/value maps in each BAM fields are parsed.
//
class SamReader : public Reader {
 public:
  // Creates a new SamReader reading reads from the SAM/BAM file reads_path.
  //
  // reads_path must point to an existing SAM/BAM formatted file (text SAM or
  // compressed or uncompressed BAM file).
  //
  // If options.index_mode indicates we should load an index, this constructor
  // will attempt to load a BAI index from file reads_path + '.bai'.
  //
  // Returns a StatusOr that is OK if the SamReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<SamReader>> FromFile(
      const string& reads_path, const SamReaderOptions& options);

  ~SamReader();

  // Disable assignment/copy operations
  SamReader(const SamReader& other) = delete;
  SamReader& operator=(const SamReader&) = delete;

  // Gets all of the reads in this file in order.
  //
  // This function allows one to iterate through all of the reads in this
  // SAM/BAM file in order.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction. Returns an OK status if the iterable can be
  // constructed, or not OK otherwise.
  StatusOr<std::shared_ptr<SamIterable>> Iterate() const;

  // Gets all of the reads that overlap any bases in range.
  //
  // This function allows one to iterate through all of the reads in this
  // SAM/BAM file in order that overlap a specific interval on the genome. The
  // query operation is efficient in that the cost is O(n) for n elements that
  // overlap range, and not O(N) for N elements in the entire file.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction.
  //
  // This function is only available if an index was loaded. If no index was
  // loaded a non-OK status value will be returned.
  //
  // If range isn't a valid interval in this BAM file a non-OK status value will
  // be returned.
  StatusOr<std::shared_ptr<SamIterable>> Query(
      const nucleus::genomics::v1::Range& region) const;

  // Returns True if this SamReader loaded an index file.
  bool HasIndex() const { return idx_ != nullptr; }

  // Gets a list of the contigs used by this SAM file.
  const std::vector<ContigInfo>& Contigs() const { return contigs_; }

  // Close the underlying resource descriptors. Returns a Status to indicate if
  // everything went OK with the close.
  tensorflow::Status Close();

  // This no-op function is needed only for Python context manager support.  Do
  // not use it! Returns a Status indicating whether the enter was successful.
  tensorflow::Status PythonEnter() const { return tensorflow::Status::OK(); }

  bool KeepRead(const nucleus::genomics::v1::Read& read) const;

  const SamReaderOptions& options() const { return options_; }

  // Returns all unique sample names from the read groups of the SAM/BAM.
  const std::set<string>& Samples() const { return samples_; }

 private:
  // Private constructor; use FromFile to safely create a SamReader from a
  // file.
  SamReader(const string& reads_path, const SamReaderOptions& options,
            htsFile* fp, bam_hdr_t* header, hts_idx_t* idx);

  // Our options that control the behavior of this class.
  const SamReaderOptions options_;

  // A pointer to the htslib file used to access the SAM/BAM data.
  htsFile * fp_;

  // A htslib header data structure obtained by parsing the header of this BAM.
  bam_hdr_t * header_;

  // The htslib index data structure for our indexed BAM file. May be NULL if no
  // index was loaded.
  hts_idx_t* idx_;

  // A list of ContigInfo, each of which contains the information about the
  // contigs used by this BAM file.
  std::vector<ContigInfo> contigs_;

  // A set of sample names for which this SAM/BAM file contains reads.
  std::set<string> samples_;

  void ParseSamplesFromHeader();

  // For downsampling reads.
  mutable core::PhiloxFractionalSampler sampler_;
};

}  // namespace core
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_CORE_SAM_READER_H_
