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

#ifndef THIRD_PARTY_NUCLEUS_IO_VCF_READER_H_
#define THIRD_PARTY_NUCLEUS_IO_VCF_READER_H_

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/io/vcf_conversion.h"
#include "third_party/nucleus/protos/index.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/vendor/statusor.h"
#include "tensorflow/core/platform/types.h"

namespace nucleus {

using tensorflow::string;

// Alias for the abstract base class for VCF record iterables.
using VariantIterable = Iterable<nucleus::genomics::v1::Variant>;

// A VCF reader that provides access to Tabix indexed VCF files.
//
// VCF files store information about genetic variation:
//
// https://samtools.github.io/hts-specs/VCFv4.2.pdf
//
// These files are commonly block-gzipped and indexed with Tabix:
//
// https://academic.oup.com/bioinformatics/article/27/5/718/262743/Tabix-fast-retrieval-of-sequence-features-from
// https://github.com/samtools/htslib
// http://www.htslib.org/doc/tabix.html
//
// This class provides methods to iterate through a VCF file or, if indexed
// with Tabix, to also query() for only variants overlapping a specific regions
// on the genome.
//
// The objects returned by iterate() or query() are nucleus.genomics.v1.Variant
// objects parsed from the VCF records in the file. Currently all fields except
// the INFO key/value maps in the VCF variant and genotype fields are parsed.
//
class VcfReader : public Reader {
 public:
  // Creates a new VcfReader reading variants from the VCF file variantsPath.
  //
  // variantsPath must point to an existing VCF formatted file (text or
  // bgzip compressed VCF file).
  //
  // If options.index_mode indicates we should load an index, this constructor
  // will attempt to load an Tabix index from file variantsPath + '.tbi'.
  //
  // Returns a StatusOr that is OK if the VcfReader could be successfully
  // created or an error code indicating the error that occurred.
  static StatusOr<std::unique_ptr<VcfReader>> FromFile(
      const string& variants_path,
      const nucleus::genomics::v1::VcfReaderOptions& options);

  ~VcfReader();


  // Disable copy or assignment
  VcfReader(const VcfReader& other) = delete;
  VcfReader& operator=(const VcfReader&) = delete;

  // Gets all of the variants in this file in order.
  //
  // This function allows one to iterate through all of the variants in this
  // VCF file in order.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction. Returns an OK status if the iterable can be
  // constructed, or not OK otherwise.
  StatusOr<std::shared_ptr<VariantIterable>> Iterate();

  // Gets all of the variants that overlap any bases in range.
  //
  // This function allows one to iterate through all of the variants in this
  // VCF file in order that overlap a specific iterval on the genome. The query
  // operation is efficient in that the cost is O(n) for n elements that overlap
  // range, and not O(N) for N elements in the entire file.
  //
  // The specific parsing, filtering, etc behavior is determined by the options
  // provided during construction.
  //
  // This function is only available if an index was loaded. If no index was
  // loaded a non-OK status value will be returned.
  //
  // If range isn't a valid interval in this VCF file a non-OK status value will
  // be returned.
  StatusOr<std::shared_ptr<VariantIterable>> Query(
      const nucleus::genomics::v1::Range& region);

  // Returns True if this VcfReader loaded an index file.
  bool HasIndex() const { return idx_ != nullptr; }

  // Returns the VCF header associated with this reader.
  const nucleus::genomics::v1::VcfHeader Header() const { return vcf_header_; }

  // Get the options controlling the behavior of this VcfReader.
  const nucleus::genomics::v1::VcfReaderOptions& Options() const {
    return options_;
  }

  // Close the underlying resource descriptors. Returns a Status to indicate if
  // everything went OK with the close.
  tensorflow::Status Close();

  // This no-op function is needed only for Python context manager support.  Do
  // not use it! Returns a Status indicating whether the enter was successful.
  tensorflow::Status PythonEnter() const { return tensorflow::Status::OK(); }

  // Access to the record converter.
  const VcfRecordConverter& RecordConverter() const {
    return record_converter_;
  }

 private:
  VcfReader(const string& variants_path,
            const nucleus::genomics::v1::VcfReaderOptions& options, htsFile* fp,
            bcf_hdr_t* header, tbx_t* idx);

  // The options controlling the behavior of this VcfReader.
  const nucleus::genomics::v1::VcfReaderOptions options_;

  // A pointer to the htslib file used to access the VCF data.
  htsFile * fp_;

  // A htslib header data structure obtained by parsing the header of this VCF.
  bcf_hdr_t * header_;

  // The htslib tbx_t data structure for tabix indexed files. May be NULL if no
  // index was loaded.
  tbx_t* idx_;

  // The VcfHeader data structure that represents the information in the header
  // of the VCF.
  nucleus::genomics::v1::VcfHeader vcf_header_;

  // Object for converting VCF records to to Variant proto.
  VcfRecordConverter record_converter_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_VCF_READER_H_
