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

// Implementation of vcf_writer.h
#include "third_party/nucleus/io/vcf_writer.h"

#include <cmath>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/strings/substitute.h"
#include "absl/log/check.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/io/vcf_conversion.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/status.h"
#include "google/protobuf/map.h"
#include "google/protobuf/repeated_field.h"

namespace nucleus {

using nucleus::genomics::v1::Variant;

namespace {

// Modes used to write a htsFile. See hts_open in hts.h for more info.
constexpr char kBCFOpenModeCompressed[] = "wb";
constexpr char kBCFOpenModeUncompressed[] = "wbu";
constexpr char kOpenModeCompressed[] = "wz";
constexpr char kOpenModeUncompressed[] = "w";

// RAII wrapper on top of bcf1_t* to always perform cleanup.
class BCFRecord {
 public:
  BCFRecord() : bcf1_(bcf_init()) {}
  ~BCFRecord() {
    if (bcf1_) {
      bcf_destroy(bcf1_);
    }
  }

  // Disable copy or assignment
  BCFRecord(const BCFRecord& other) = delete;
  BCFRecord& operator=(const BCFRecord&) = delete;

  bcf1_t* get_bcf1() { return bcf1_; }

 private:
  bcf1_t* bcf1_;
};

}  // namespace

StatusOr<std::unique_ptr<VcfWriter>> VcfWriter::ToFile(
    const string& variants_path, const nucleus::genomics::v1::VcfHeader& header,
    const nucleus::genomics::v1::VcfWriterOptions& options) {
  const char* const open_mode = GetOpenMode(variants_path);
  htsFile* fp = hts_open_x(variants_path, open_mode);
  if (fp == nullptr) {
    return ::nucleus::Unknown(
        absl ::StrCat("Could not open variants_path: ", variants_path));
  }

  auto writer = absl::WrapUnique(new VcfWriter(header, options, fp));
  NUCLEUS_RETURN_IF_ERROR(writer->WriteHeader());
  return std::move(writer);
}

VcfWriter::VcfWriter(const nucleus::genomics::v1::VcfHeader& header,
                     const nucleus::genomics::v1::VcfWriterOptions& options,
                     htsFile* fp)
    : fp_(fp),
      options_(options),
      vcf_header_(header),
      record_converter_(
          vcf_header_,
          std::vector<string>(options_.excluded_info_fields().begin(),
                              options_.excluded_info_fields().end()),
          std::vector<string>(options_.excluded_format_fields().begin(),
                              options_.excluded_format_fields().end()),
          options_.retrieve_gl_and_pl_from_info_map()) {
  CHECK(fp != nullptr);

  NUCLEUS_CHECK_OK(VcfHeaderConverter::ConvertFromPb(vcf_header_, &header_));
}

::nucleus::Status VcfWriter::WriteHeader() {
  if (options_.exclude_header()) {
    return ::nucleus::Status();
  }
  if (bcf_hdr_write(fp_, header_) < 0) {
    return ::nucleus::Unknown("Failed to write header");
  }
  return ::nucleus::Status();
}

VcfWriter::~VcfWriter() {
  if (fp_) {
    // There's nothing we can do but assert fail if there's an error during
    // the Close() call here.
    NUCLEUS_CHECK_OK(Close());
  }
}

::nucleus::Status VcfWriter::Write(const Variant& variant_message) {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition("Cannot write to closed VCF stream.");
  BCFRecord v;
  if (v.get_bcf1() == nullptr) {
    return ::nucleus::Unknown("bcf_init call failed");
  }
  NUCLEUS_RETURN_IF_ERROR(
      RecordConverter().ConvertFromPb(variant_message, *header_, v.get_bcf1()));
  if (options_.round_qual_values() &&
      !bcf_float_is_missing(v.get_bcf1()->qual)) {
    // Round quality value printed out to one digit past the decimal point.
    double rounded_quality = floor(variant_message.quality() * 10 + 0.5) / 10;
    v.get_bcf1()->qual = rounded_quality;
  }
  if (bcf_write(fp_, header_, v.get_bcf1()) != 0) {
    return ::nucleus::Unknown("bcf_write call failed");
  }
  return ::nucleus::Status();
}

::nucleus::Status VcfWriter::Close() {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition(
        "Cannot close an already closed VcfWriter");
  if (hts_close(fp_) < 0) return ::nucleus::Unknown("hts_close call failed");
  fp_ = nullptr;
  bcf_hdr_destroy(header_);
  header_ = nullptr;
  return ::nucleus::Status();
}

// static
const char* VcfWriter::GetOpenMode(const string& file_path) {
  if (EndsWith(file_path, ".bcf")) {
    return kBCFOpenModeUncompressed;
  } else if (EndsWith(file_path, ".bcf.gz")) {
    return kBCFOpenModeCompressed;
  } else if (EndsWith(file_path, ".gz")) {
    return kOpenModeCompressed;
  }
  return kOpenModeUncompressed;
}

}  // namespace nucleus
