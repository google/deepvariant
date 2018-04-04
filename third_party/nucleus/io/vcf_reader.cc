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

// Implementation of vcf_reader.h
#include "third_party/nucleus/io/vcf_reader.h"

#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/io/vcf_conversion.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/math.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

namespace tf = tensorflow;

using std::vector;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using tensorflow::strings::StrCat;

namespace {

// Returns the hrec that contains information or nullptr if none does.
const bcf_hrec_t* GetPopulatedHrec(const bcf_idpair_t& idPair) {
  for (int i = 0; i < 3; i++) {
    const bcf_hrec_t* hrec = idPair.val->hrec[i];
    if (hrec != nullptr) {
      return hrec;
    }
  }
  LOG(ERROR) << "No populated hrec in idPair. Error in htslib.";
  return nullptr;
}

// Adds Contig information from the idPair to the ContigInfo object.
void AddContigInfo(const bcf_idpair_t& idPair,
                   nucleus::genomics::v1::ContigInfo* contig,
                   int pos_in_fasta) {
  // ID and length are special-cased in the idPair.
  contig->set_name(idPair.key);
  contig->set_n_bases(idPair.val->info[0]);
  contig->set_pos_in_fasta(pos_in_fasta);
  const bcf_hrec_t* hrec0 = GetPopulatedHrec(idPair);
  if (hrec0 != nullptr) {
    for (int j = 0; j < hrec0->nkeys; j++) {
      // Add any non-ID and non-length info to the structured map of additional
      // information. "IDX" is an htslib-internal key that should also be
      // ignored.
      if (string(hrec0->keys[j]) != "ID" &&
          string(hrec0->keys[j]) != "length" &&
          string(hrec0->keys[j]) != "IDX") {
        (*contig->mutable_extra())[hrec0->keys[j]] =
            Unquote(hrec0->vals[j]).ToString();
      }
    }
  }
}

// Adds FILTER information from the bcf_hrec_t to the VcfFilterInfo object.
void AddFilterInfo(const bcf_hrec_t* hrec,
                   nucleus::genomics::v1::VcfFilterInfo* filter) {
  if (hrec->nkeys >= 2 && string(hrec->keys[0]) == "ID" &&
      string(hrec->keys[1]) == "Description") {
    filter->set_id(hrec->vals[0]);
    // "Unquote" the description identifier.
    filter->set_description(Unquote(hrec->vals[1]).ToString());
  } else {
    LOG(WARNING) << "Malformed FILTER field detected in header, leaving this "
                    "filter empty";
  }
}

// Adds INFO information from the bcf_hrec_t to the VcfInfo object.
void AddInfo(const bcf_hrec_t* hrec, nucleus::genomics::v1::VcfInfo* info) {
  if (hrec->nkeys >= 4 && string(hrec->keys[0]) == "ID" &&
      string(hrec->keys[1]) == "Number" && string(hrec->keys[2]) == "Type" &&
      string(hrec->keys[3]) == "Description") {
    info->set_id(hrec->vals[0]);
    info->set_number(hrec->vals[1]);
    info->set_type(hrec->vals[2]);
    info->set_description(Unquote(hrec->vals[3]).ToString());
    for (int i = 4; i < hrec->nkeys; i++) {
      if (string(hrec->keys[i]) == "Source") {
        info->set_source(Unquote(hrec->vals[i]).ToString());
      } else if (string(hrec->keys[i]) == "Version") {
        info->set_version(Unquote(hrec->vals[i]).ToString());
      }
    }
  } else {
    LOG(WARNING) << "Malformed INFO field detected in header, leaving this "
                    "info empty";
  }
}

// Adds FORMAT information from the bcf_hrec_t to the VcfFormatInfo object.
void AddFormatInfo(const bcf_hrec_t* hrec,
                   nucleus::genomics::v1::VcfFormatInfo* format) {
  if (hrec->nkeys >= 4 && string(hrec->keys[0]) == "ID" &&
      string(hrec->keys[1]) == "Number" && string(hrec->keys[2]) == "Type" &&
      string(hrec->keys[3]) == "Description") {
    format->set_id(hrec->vals[0]);
    format->set_number(hrec->vals[1]);
    format->set_type(hrec->vals[2]);
    format->set_description(Unquote(hrec->vals[3]).ToString());
  } else {
    LOG(WARNING) << "Malformed FORMAT field detected in header, leaving this "
                    "format empty";
  }
}

// Adds structured information from the bcf_hrec_t to the VcfStructuredExtra.
void AddStructuredExtra(const bcf_hrec_t* hrec,
                        nucleus::genomics::v1::VcfStructuredExtra* extra) {
  extra->set_key(hrec->key);
  for (int i = 0; i < hrec->nkeys; i++) {
    nucleus::genomics::v1::VcfExtra& toAdd = *extra->mutable_fields()->Add();
    toAdd.set_key(hrec->keys[i]);
    toAdd.set_value(Unquote(hrec->vals[i]).ToString());
  }
}

// Adds unstructured information from the bcf_hrec_t to the VcfExtra object.
void AddExtra(const bcf_hrec_t* hrec, nucleus::genomics::v1::VcfExtra* extra) {
  extra->set_key(hrec->key);
  extra->set_value(hrec->value);
}

}  // namespace

// Iterable class for traversing VCF records found in a query window.
class VcfQueryIterable : public VariantIterable {
 public:
  // Advance to the next record.
  StatusOr<bool> Next(nucleus::genomics::v1::Variant* out) override;

  // Constructor will be invoked via VcfReader::Query.
  VcfQueryIterable(const VcfReader* reader,
                   htsFile* fp,
                   bcf_hdr_t* header,
                   tbx_t* idx,
                   hts_itr_t* iter);

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
  VcfFullFileIterable(const VcfReader* reader,
                      htsFile* fp,
                      bcf_hdr_t* header);

  ~VcfFullFileIterable() override;

 private:
  htsFile* fp_;
  bcf_hdr_t* header_;
  bcf1_t* bcf1_;
};

StatusOr<std::unique_ptr<VcfReader>> VcfReader::FromFile(
    const string& variants_path,
    const nucleus::genomics::v1::VcfReaderOptions& options) {
  htsFile* fp = hts_open_x(variants_path.c_str(), "r");
  if (fp == nullptr) {
    return tf::errors::NotFound(StrCat("Could not open ", variants_path));
  }

  bcf_hdr_t* header = bcf_hdr_read(fp);
  if (header == nullptr)
    return tf::errors::Unknown(StrCat("Couldn't parse header for ", fp->fn));

  // Try to load the Tabix index if requested.
  tbx_t* idx = nullptr;
  if (options.index_mode() ==
      nucleus::genomics::v1::IndexHandlingMode::INDEX_BASED_ON_FILENAME) {
    // tbx_index_load() computes the expected index filename from the
    // filename fn of the htsFile fp.
    idx = tbx_index_load(fp->fn);
    if (idx == nullptr)
      return tf::errors::NotFound(StrCat("No index found for ", fp->fn));
  }

  return std::unique_ptr<VcfReader>(
      new VcfReader(variants_path, options, fp, header, idx));
}

VcfReader::VcfReader(const string& variants_path,
                     const nucleus::genomics::v1::VcfReaderOptions& options,
                     htsFile* fp, bcf_hdr_t* header, tbx_t* idx)
    : options_(options), fp_(fp), header_(header), idx_(idx) {
  if (header_->nhrec < 1) {
    LOG(WARNING) << "Empty header, not a valid VCF.";
    return;
  }
  if (string(header_->hrec[0]->key) != "fileformat") {
    LOG(WARNING) << "Not a valid VCF, fileformat needed: " << variants_path;
  }
  vcf_header_.set_fileformat(header_->hrec[0]->value);

  // Fill in the contig info for each contig in the VCF header. Directly
  // accesses the low-level C struct because there are no indirection
  // macros/functions by htslib API.
  // BCF_DT_CTG: offset for contig (CTG) information in BCF dictionary (DT).
  const int n_contigs = header_->n[BCF_DT_CTG];
  for (int i = 0; i < n_contigs; ++i) {
    const bcf_idpair_t& idPair = header->id[BCF_DT_CTG][i];
    AddContigInfo(idPair, vcf_header_.add_contigs(), i);
  }

  // Iterate through all hrecs (except the first, which was 'fileformat') to
  // populate the rest of the headers.
  for (int i = 1; i < header_->nhrec; i++) {
    const bcf_hrec_t* hrec0 = header_->hrec[i];
    switch (hrec0->type) {
      case BCF_HL_CTG:
        // Contigs are populated above, since they store length in the
        // bcf_idinfo_t* structure.
        break;
      case BCF_HL_FLT:
        AddFilterInfo(hrec0, vcf_header_.add_filters());
        break;
      case BCF_HL_INFO:
        AddInfo(hrec0, vcf_header_.add_infos());
        break;
      case BCF_HL_FMT:
        AddFormatInfo(hrec0, vcf_header_.add_formats());
        break;
      case BCF_HL_STR:
        AddStructuredExtra(hrec0, vcf_header_.add_structured_extras());
        break;
      case BCF_HL_GEN:
        AddExtra(hrec0, vcf_header_.add_extras());
        break;
      default:
        LOG(WARNING) << "Unknown hrec0->type: " << hrec0->type;
    }
  }

  // Populate samples info.
  int n_samples = bcf_hdr_nsamples(header_);
  for (int i = 0; i < n_samples; i++) {
    vcf_header_.add_sample_names(header_->samples[i]);
  }
  vector<string> infos_to_exclude(options.excluded_info_fields().begin(),
                                  options.excluded_info_fields().end());
  vector<string> formats_to_exclude(options.excluded_format_fields().begin(),
                                    options.excluded_format_fields().end());
  record_converter_ =
      VcfRecordConverter(vcf_header_, infos_to_exclude, formats_to_exclude);
}

VcfReader::~VcfReader() {
  if (fp_) {
    // We cannot return a value from the destructor, so the best we can do is
    // CHECK-fail if the Close() wasn't successful.
    TF_CHECK_OK(Close());
  }
}

StatusOr<std::shared_ptr<VariantIterable>> VcfReader::Iterate() {
  if (fp_ == nullptr)
    return tf::errors::FailedPrecondition("Cannot Iterate a closed VcfReader.");
  return StatusOr<std::shared_ptr<VariantIterable>>(
      MakeIterable<VcfFullFileIterable>(this, fp_, header_));
}

StatusOr<std::shared_ptr<VariantIterable>> VcfReader::Query(
    const Range& region) {
  if (fp_ == nullptr)
    return tf::errors::FailedPrecondition("Cannot Query a closed VcfReader.");
  if (!HasIndex()) {
    return tf::errors::FailedPrecondition("Cannot query without an index");
  }

  const char* reference_name = region.reference_name().c_str();
  if (bcf_hdr_name2id(header_, reference_name) < 0) {
    return tf::errors::NotFound(
        StrCat("Unknown reference_name '", region.reference_name(), "'"));
  }
  if (region.start() < 0 || region.start() >= region.end())
    return tf::errors::InvalidArgument(
        StrCat("Malformed region '", region.ShortDebugString(), "'"));

  // Get the tid (index of reference_name in our tabix index),
  const int tid = tbx_name2id(idx_, reference_name);
  hts_itr_t* iter = nullptr;
  if (tid >= 0) {
    // Note that query is 0-based inclusive on start and exclusive on end,
    // matching exactly the logic of our Range.
    iter = tbx_itr_queryi(idx_, tid, region.start(), region.end());
    if (iter == nullptr) {
      return tf::errors::NotFound(
          StrCat("region '", region.ShortDebugString(),
                 "' returned an invalid hts_itr_queryi result"));
    }
  }  // implicit else case:
  // The chromosome isn't reflected in the tabix index (meaning, no
  // variant records) => return an *empty* iterable by leaving iter empty.
  return StatusOr<std::shared_ptr<VariantIterable>>(
      MakeIterable<VcfQueryIterable>(this, fp_, header_, idx_, iter));
}

tf::Status VcfReader::Close() {
  if (fp_ == nullptr)
    return tf::errors::FailedPrecondition("VcfReader already closed");
  if (HasIndex()) {
    tbx_destroy(idx_);
    idx_ = nullptr;
  }
  bcf_hdr_destroy(header_);
  header_ = nullptr;
  int retval = hts_close(fp_);
  fp_ = nullptr;
  if (retval < 0) {
    return tf::errors::Internal("hts_close() failed");
  } else {
    return tf::Status::OK();
  }
}

// Iterable class definitions.

StatusOr<bool> VcfQueryIterable::Next(Variant* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  if (tbx_itr_next(fp_, idx_, iter_, &str_) < 0) return false;
  if (vcf_parse1(&str_, header_, bcf1_) < 0) {
    return tf::errors::DataLoss(StrCat("Failed to parse VCF record: ", str_.s));
  }
  const VcfReader* reader = static_cast<const VcfReader*>(reader_);
  TF_RETURN_IF_ERROR(
      reader->RecordConverter().ConvertToPb(header_, bcf1_, out));
  return true;
}

VcfQueryIterable::~VcfQueryIterable() {
  hts_itr_destroy(iter_);
  bcf_destroy(bcf1_);
  if (str_.s != nullptr) { free(str_.s); }
}

VcfQueryIterable::VcfQueryIterable(const VcfReader* reader,
                                   htsFile* fp,
                                   bcf_hdr_t* header,
                                   tbx_t* idx,
                                   hts_itr_t* iter)
    : Iterable(reader),
      fp_(fp),
      header_(header),
      bcf1_(bcf_init()),
      idx_(idx),
      iter_(iter),
      str_({0, 0, nullptr})
{}


StatusOr<bool> VcfFullFileIterable::Next(Variant* out) {
  TF_RETURN_IF_ERROR(CheckIsAlive());
  if (bcf_read(fp_, header_, bcf1_) < 0) {
    if (bcf1_->errcode) {
      return tf::errors::DataLoss("Failed to parse VCF record");
    } else {
      return false;
    }
  }
  const VcfReader* reader = static_cast<const VcfReader*>(reader_);
  TF_RETURN_IF_ERROR(
      reader->RecordConverter().ConvertToPb(header_, bcf1_, out));
  return true;
}

VcfFullFileIterable::~VcfFullFileIterable() {
  bcf_destroy(bcf1_);
}

VcfFullFileIterable::VcfFullFileIterable(const VcfReader* reader,
                                         htsFile* fp,
                                         bcf_hdr_t* header)
    : Iterable(reader),
      fp_(fp),
      header_(header),
      bcf1_(bcf_init())
{}

}  // namespace nucleus
