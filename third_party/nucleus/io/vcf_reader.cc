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

// Implementation of vcf_reader.h
#include "third_party/nucleus/io/vcf_reader.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "absl/memory/memory.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/io/vcf_conversion.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/math.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"
#include "google/protobuf/map.h"
#include "google/protobuf/repeated_field.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Variant;
using std::vector;

namespace {

bool FileTypeIsIndexable(htsFormat format) {
  return format.format == vcf && format.compression == bgzf;
}

}  // namespace

// Iterable class for traversing VCF records found in a query window.
class VcfQueryIterable : public VariantIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::Variant* out) override;

  // Constructor will be invoked via VcfReader::Query.
  VcfQueryIterable(const VcfReader* reader, htsFile* fp, bcf_hdr_t* header,
                   tbx_t* idx, hts_itr_t* iter);

  ~VcfQueryIterable() override;

 private:
  htsFile* fp_;
  bcf_hdr_t* header_;
  bcf1_t* bcf1_;
  tbx_t* idx_;
  hts_itr_t* iter_;
  kstring_t str_;
};

// Iterable class for traversing all VCF records in the file.
class VcfFullFileIterable : public VariantIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::Variant* out) override;

  // Constructor will be invoked via VcfReader::Iterate.
  VcfFullFileIterable(const VcfReader* reader, htsFile* fp, bcf_hdr_t* header);

  ~VcfFullFileIterable() override;

 private:
  htsFile* fp_;
  bcf_hdr_t* header_;
  bcf1_t* bcf1_;
};

StatusOr<std::unique_ptr<VcfReader>> VcfReader::FromFile(
    const string& vcf_filepath,
    const nucleus::genomics::v1::VcfReaderOptions& options) {
  return FromFileHelper(vcf_filepath, options, nullptr);
}

StatusOr<std::unique_ptr<VcfReader>> VcfReader::FromFileWithHeader(
    const string& vcf_filepath,
    const nucleus::genomics::v1::VcfReaderOptions& options,
    const nucleus::genomics::v1::VcfHeader& header) {
  bcf_hdr_t* h = nullptr;
  NUCLEUS_RETURN_IF_ERROR(VcfHeaderConverter::ConvertFromPb(header, &h));
  return FromFileHelper(vcf_filepath, options, h);
}

StatusOr<std::unique_ptr<VcfReader>> VcfReader::FromFileHelper(
    const string& vcf_filepath,
    const nucleus::genomics::v1::VcfReaderOptions& options, bcf_hdr_t* h) {
  htsFile* fp = hts_open_x(vcf_filepath, "r");
  if (fp == nullptr) {
    return ::nucleus::NotFound(absl::StrCat("Could not open ", vcf_filepath));
  }

  if (h == nullptr) {
    h = bcf_hdr_read(fp);
    if (h == nullptr) {
      return ::nucleus::Unknown(
          absl::StrCat("Couldn't parse header for ", fp->fn));
    }
  } else {
    // Call bcf_hdr_read to verify that this is a headerless VCF file.
    bcf_hdr_t* null_h = bcf_hdr_read(fp);
    if (null_h != nullptr) {
      // TODO: We need RAII wrappers for these raw htslib pointers.
      hts_close(fp);
      bcf_hdr_destroy(null_h);
      bcf_hdr_destroy(h);
      return ::nucleus::Unknown(
          absl::StrCat("Unexpected header in", vcf_filepath));
    }
    // Without the header, htslib fails to parse the format. Default to vcf.
    // TODO: support bcf files with no header.
    fp->format.format = htsExactFormat::vcf;
  }

  // Try to load the Tabix index if requested.
  tbx_t* idx = nullptr;
  if (FileTypeIsIndexable(fp->format)) {
    idx = tbx_index_load(fp->fn);
    // idx may be null; only an error if we try to Query later.
  }

  return absl::WrapUnique<VcfReader>(
      new VcfReader(vcf_filepath, options, fp, h, idx));
}

void VcfReader::NativeHeaderUpdated() {
  VcfHeaderConverter::ConvertToPb(header_, &vcf_header_);
  vector<string> infos_to_exclude(options_.excluded_info_fields().begin(),
                                  options_.excluded_info_fields().end());
  vector<string> formats_to_exclude(options_.excluded_format_fields().begin(),
                                    options_.excluded_format_fields().end());
  record_converter_ =
      VcfRecordConverter(vcf_header_, infos_to_exclude, formats_to_exclude,
                         options_.store_gl_and_pl_in_info_map());
}

VcfReader::VcfReader(const string& vcf_filepath,
                     const nucleus::genomics::v1::VcfReaderOptions& options,
                     htsFile* fp, bcf_hdr_t* header, tbx_t* idx)
    : vcf_filepath_(vcf_filepath),
      options_(options),
      fp_(fp),
      header_(header),
      idx_(idx),
      bcf1_(bcf_init()) {
  NativeHeaderUpdated();
}

VcfReader::~VcfReader() {
  bcf_destroy(bcf1_);
  if (fp_) {
    // We cannot return a value from the destructor, so the best we can do is
    // CHECK-fail if the Close() wasn't successful.
    NUCLEUS_CHECK_OK(Close());
  }
}

StatusOr<std::shared_ptr<VariantIterable>> VcfReader::Iterate() {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition("Cannot Iterate a closed VcfReader.");
  return StatusOr<std::shared_ptr<VariantIterable>>(
      MakeIterable<VcfFullFileIterable>(this, fp_, header_));
}

StatusOr<std::shared_ptr<VariantIterable>> VcfReader::Query(
    const Range& region) {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition("Cannot Query a closed VcfReader.");
  if (!HasIndex()) {
    return ::nucleus::FailedPrecondition("Cannot query without an index");
  }

  const char* reference_name = region.reference_name().c_str();
  if (bcf_hdr_name2id(header_, reference_name) < 0) {
    return ::nucleus::NotFound(
        absl::StrCat("Unknown reference_name '", region.reference_name(), "'"));
  }
  if (region.start() < 0 || region.start() >= region.end())
    return ::nucleus::InvalidArgument(
        absl::StrCat("Malformed region '", region.ShortDebugString(), "'"));

  // Get the tid (index of reference_name in our tabix index),
  const int tid = tbx_name2id(idx_, reference_name);
  hts_itr_t* iter = nullptr;
  if (tid >= 0) {
    // Note that query is 0-based inclusive on start and exclusive on end,
    // matching exactly the logic of our Range.
    iter = tbx_itr_queryi(idx_, tid, region.start(), region.end());
    if (iter == nullptr) {
      return ::nucleus::NotFound(
          absl::StrCat("region '", region.ShortDebugString(),
                       "' returned an invalid hts_itr_queryi result"));
    }
  }  // implicit else case:
  // The chromosome isn't reflected in the tabix index (meaning, no
  // variant records) => return an *empty* iterable by leaving iter empty.
  return StatusOr<std::shared_ptr<VariantIterable>>(
      MakeIterable<VcfQueryIterable>(this, fp_, header_, idx_, iter));
}

::nucleus::Status VcfReader::FromString(const absl::string_view& vcf_line,
                                        nucleus::genomics::v1::Variant* v) {
  size_t len = vcf_line.length();
  std::unique_ptr<char[]> cstr{new char[len + 1]};
  std::strncpy(cstr.get(), vcf_line.data(), len);
  *(cstr.get() + len) = '\0';
  kstring_t str = {.l = len + 1, .m = len + 1, .s = cstr.get()};

  // vcf_parse1 returns -1 on critical errors and 0 otherwise. BCF_ERR_CTG_UNDEF
  // and BCF_ERR_TAG_UNDEF indicate missing header definitions, and are
  // non-critical errors. Ignore these missing header definitions because they
  // are common in the wild.
  if (vcf_parse1(&str, header_, bcf1_) < 0) {
    return ::nucleus::DataLoss(
        absl::StrCat("Failed to parse VCF record: ", cstr.get()));
  }
  if (bcf1_->errcode == BCF_ERR_CTG_UNDEF ||
      bcf1_->errcode == BCF_ERR_TAG_UNDEF) {
    bcf1_->errcode = 0;
    NativeHeaderUpdated();
  }

  if (bcf1_->errcode != 0) {
    return ::nucleus::DataLoss(absl::StrCat(
        "Failed to parse VCF record with errcode: ", bcf1_->errcode));
  }

  NUCLEUS_RETURN_IF_ERROR(RecordConverter().ConvertToPb(header_, bcf1_, v));
  return ::nucleus::Status();
}

StatusOr<bool> VcfReader::FromStringPython(const absl::string_view& vcf_line,
                                           nucleus::genomics::v1::Variant* v) {
  ::nucleus::Status s = FromString(vcf_line, v);
  if (!s.ok()) {
    return s;
  }
  return true;
}

::nucleus::Status VcfReader::Close() {
  if (fp_ == nullptr)
    return ::nucleus::FailedPrecondition("VcfReader already closed");
  if (HasIndex()) {
    tbx_destroy(idx_);
    idx_ = nullptr;
  }
  bcf_hdr_destroy(header_);
  header_ = nullptr;
  int retval = hts_close(fp_);
  fp_ = nullptr;
  if (retval < 0) {
    return ::nucleus::Internal("hts_close() failed");
  } else {
    return ::nucleus::Status();
  }
}

// Iterable class definitions.

StatusOr<bool> VcfQueryIterable::Next(Variant* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  if (tbx_itr_next(fp_, idx_, iter_, &str_) < 0) return false;
  if (vcf_parse1(&str_, header_, bcf1_) < 0) {
    return ::nucleus::DataLoss(
        absl::StrCat("Failed to parse VCF record: ", str_.s));
  }
  const VcfReader* reader = static_cast<const VcfReader*>(reader_);
  NUCLEUS_RETURN_IF_ERROR(
      reader->RecordConverter().ConvertToPb(header_, bcf1_, out));
  return true;
}

VcfQueryIterable::~VcfQueryIterable() {
  hts_itr_destroy(iter_);
  bcf_destroy(bcf1_);
  if (str_.s != nullptr) {
    free(str_.s);
  }
}

VcfQueryIterable::VcfQueryIterable(const VcfReader* reader, htsFile* fp,
                                   bcf_hdr_t* header, tbx_t* idx,
                                   hts_itr_t* iter)
    : Iterable(reader),
      fp_(fp),
      header_(header),
      bcf1_(bcf_init()),
      idx_(idx),
      iter_(iter),
      str_({0, 0, nullptr}) {}

StatusOr<bool> VcfFullFileIterable::Next(Variant* out) {
  NUCLEUS_RETURN_IF_ERROR(CheckIsAlive());
  if (bcf_read(fp_, header_, bcf1_) < 0) {
    if (bcf1_->errcode) {
      return ::nucleus::DataLoss("Failed to parse VCF record");
    } else {
      return false;
    }
  }
  const VcfReader* reader = static_cast<const VcfReader*>(reader_);
  NUCLEUS_RETURN_IF_ERROR(
      reader->RecordConverter().ConvertToPb(header_, bcf1_, out));
  return true;
}

VcfFullFileIterable::~VcfFullFileIterable() { bcf_destroy(bcf1_); }

VcfFullFileIterable::VcfFullFileIterable(const VcfReader* reader, htsFile* fp,
                                         bcf_hdr_t* header)
    : Iterable(reader), fp_(fp), header_(header), bcf1_(bcf_init()) {}

}  // namespace nucleus
