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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_CORE_VCF_WRITER_H_
#define LEARNING_GENOMICS_DEEPVARIANT_CORE_VCF_WRITER_H_

#include "deepvariant/core/genomics/range.pb.h"
#include "deepvariant/core/genomics/variants.pb.h"
#include "deepvariant/core/protos/core.pb.h"
#include "deepvariant/vendor/statusor.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/platform/types.h"

namespace learning {
namespace genomics {
namespace core {

using tensorflow::string;

// A VCF writer, allowing us to write VCF files.
class VcfWriter {
 public:
  // Creates a new VcfWriter writing to the file at variants_path, which is
  // opened and created if needed. Returns either a unique_ptr to the VcfWriter
  // or a Status indicating why an error occurred.
  static StatusOr<std::unique_ptr<VcfWriter>> ToFile(
      const string& variants_path, const VcfWriterOptions& options);
  ~VcfWriter();

  // Disable copy or assignment
  VcfWriter(const VcfWriter& other) = delete;
  VcfWriter& operator=(const VcfWriter&) = delete;

  // Write a variant record to the VCF.
  // Note that variant calls must be provided in the same order as samples
  // listed in the options. Returns Status::OK() if the write was successful;
  // otherwise the status provides information about what error occurred.
  tensorflow::Status Write(
      const learning::genomics::v1::Variant& variant_message);

  // Close the underlying resource descriptors. Returns Status::OK() if the
  // close was successful; otherwise the status provides information about what
  // error occurred.
  tensorflow::Status Close();

  // This no-op function is needed only for Python context manager support.  Do
  // not use it!
  void PythonEnter() const {}

 private:
  VcfWriter(const VcfWriterOptions& options, htsFile* fp);

  tensorflow::Status WriteHeader();

  // A pointer to the htslib file used to write the VCF data.
  htsFile* fp_;

  // The options controlling the behavior of this VcfReader.
  const VcfWriterOptions options_;

  // A pointer to the VCF header object.
  bcf_hdr_t* header_;
};

}  // namespace core
}  // namespace genomics
}  // namespace learning


#endif  // LEARNING_GENOMICS_DEEPVARIANT_CORE_VCF_WRITER_H_
