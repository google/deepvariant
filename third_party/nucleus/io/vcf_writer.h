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

#ifndef THIRD_PARTY_NUCLEUS_IO_VCF_WRITER_H_
#define THIRD_PARTY_NUCLEUS_IO_VCF_WRITER_H_

#include <memory>
#include <string>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "third_party/nucleus/io/vcf_conversion.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace nucleus {

// A VCF writer, allowing us to write VCF files.
class VcfWriter {
 public:
  // Creates a new VcfWriter writing to the file at variants_path, which is
  // opened and created if needed. Returns either a unique_ptr to the VcfWriter
  // or a Status indicating why an error occurred.
  static StatusOr<std::unique_ptr<VcfWriter>> ToFile(
      const string& variants_path,
      const nucleus::genomics::v1::VcfHeader& header,
      const nucleus::genomics::v1::VcfWriterOptions& options);
  ~VcfWriter();

  // Disable copy or assignment
  VcfWriter(const VcfWriter& other) = delete;
  VcfWriter& operator=(const VcfWriter&) = delete;

  // Write a variant record to the VCF.
  // Note that variant calls must be provided in the same order as samples
  // listed in the options. Returns Status::OK() if the write was successful;
  // otherwise the status provides information about what error occurred.
  ::nucleus::Status Write(
      const nucleus::genomics::v1::Variant& variant_message);
  ::nucleus::Status WritePython(
      const ConstProtoPtr<const nucleus::genomics::v1::Variant>& wrapped) {
    return Write(*(wrapped.p_));
  }
  // These API extends Write() method to enable somatic processing.
  ::nucleus::Status WriteSomatic(
      const nucleus::genomics::v1::Variant& variant_message);
  ::nucleus::Status WriteSomaticPython(
      const ConstProtoPtr<const nucleus::genomics::v1::Variant>& wrapped) {
    return WriteSomatic(*(wrapped.p_));
  }

  // Close the underlying resource descriptors. Returns Status::OK() if the
  // close was successful; otherwise the status provides information about what
  // error occurred.
  ::nucleus::Status Close();

  // This no-op function is needed only for Python context manager support.  Do
  // not use it!
  void PythonEnter() const {}

  // Access to the record converter.
  const VcfRecordConverter& RecordConverter() const {
    return record_converter_;
  }

  // Infers htsFile open mode from the given file path.
  // Returns one of the following.
  //  "wb"  for compressed BCF if path ends with .bcf.gz;
  //  "wbu" for uncompressed BCF if  path ends with .bcf;
  //  "wz"  for compressed VCF if file path ends with .gz and not .bcf.gz;
  //  "w"   for all other cases default to uncompressed VCF.
  // This method is public for testing.
  static const char* GetOpenMode(const string& file_path);

 private:
  VcfWriter(const nucleus::genomics::v1::VcfHeader& header,
            const nucleus::genomics::v1::VcfWriterOptions& options,
            htsFile* fp);

  ::nucleus::Status WriteHeader();

  // A pointer to the htslib file used to write the VCF data.
  htsFile* fp_;

  // The options controlling the behavior of this VcfWriter.
  const nucleus::genomics::v1::VcfWriterOptions options_;

  // The VcfHeader proto representation of the VCF header.
  const nucleus::genomics::v1::VcfHeader vcf_header_;

  // A pointer to the VCF header object.
  bcf_hdr_t* header_;

  // VCF record interconverter.
  VcfRecordConverter record_converter_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_VCF_WRITER_H_
