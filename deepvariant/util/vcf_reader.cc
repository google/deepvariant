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

// Implementation of vcf_reader.h
#include "deepvariant/util/vcf_reader.h"

#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "deepvariant/util/genomics/range.pb.h"
#include "deepvariant/util/genomics/reference.pb.h"
#include "deepvariant/util/genomics/variants.pb.h"
#include "deepvariant/util/hts_path.h"
#include "deepvariant/util/math.h"
#include "deepvariant/util/utils.h"
#include "deepvariant/util/vcf_conversion.h"
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
    const string& variants_path, const VcfReaderOptions& options) {
  htsFile* fp = hts_open_x(variants_path.c_str(), "r");
  if (fp == nullptr) {
    return tf::errors::NotFound(StrCat("Could not open ", variants_path));
  }

  bcf_hdr_t* header = bcf_hdr_read(fp);
  if (header == nullptr)
    return tf::errors::Unknown(StrCat("Couldn't parse header for ", fp->fn));

  // Try to load the Tabix index if requested.
  tbx_t* idx = nullptr;
  if (options.index_mode() == IndexHandlingMode::INDEX_BASED_ON_FILENAME) {
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
                     const VcfReaderOptions& options, htsFile* fp,
                     bcf_hdr_t* header, tbx_t* idx)
    : options_(options),
      fp_(fp),
      header_(header),
      idx_(idx) {
  // Fill in the contig info for each contig in the VCF header. Directly
  // accesses the low-level C struct because there are no indirection
  // macros/functions by htslib API.
  // BCF_DT_CTG: offset for contig (CTG) information in BCF dictionary (DT).
  const int n_contigs = header_->n[BCF_DT_CTG];
  for (int i = 0; i < n_contigs; ++i) {
    nucleus::genomics::v1::ContigInfo contig;
    contig.set_name(header_->id[BCF_DT_CTG][i].key);
    contig.set_n_bases(header_->id[BCF_DT_CTG][i].val->info[0]);
    contig.set_pos_in_fasta(i);
    contigs_.push_back(contig);
  }

  // Populate FILTER and FORMAT info.
  int n_id_headers = header_->n[BCF_DT_ID];
  for (int i = 0; i < n_id_headers; i++) {
    const bcf_idpair_t& idPair = header_->id[BCF_DT_ID][i];
    const bcf_hrec_t* hrec0 = idPair.val->hrec[0];
    if (hrec0 != nullptr && hrec0->type == BCF_HL_FLT) {
      nucleus::genomics::v1::VcfFilterInfo filter;
      if (hrec0->nkeys >= 2 && string(hrec0->keys[0]) == "ID" &&
          string(hrec0->keys[1]) == "Description") {
        filter.set_id(hrec0->vals[0]);
        // "Unquote" the description identifier
        filter.set_description(Unquote(hrec0->vals[1]).ToString());
        filters_.push_back(filter);
      } else {
        LOG(WARNING) << "Malformed filter field detected in header, aborting "
                        "header parsing";
        break;
      }
    } else if (hrec0 != nullptr && hrec0->type == BCF_HL_FMT) {
      nucleus::genomics::v1::VcfFormatInfo format;
      if (hrec0->nkeys >= 4 && string(hrec0->keys[0]) == "ID" &&
          string(hrec0->keys[1]) == "Number" &&
          string(hrec0->keys[2]) == "Type" &&
          string(hrec0->keys[3]) == "Description") {
        format.set_id(hrec0->keys[0]);
        format.set_number(hrec0->keys[1]);
        format.set_type(hrec0->keys[2]);
        format.set_description(Unquote(hrec0->keys[3]).ToString());
        formats_.push_back(format);
      } else {
        LOG(WARNING) << "Malformed format field detected in header, aborting "
                        "header parsing";
        break;
      }
    }
  }

  // Populate samples info.
  int n_samples = bcf_hdr_nsamples(header_);
  for (int i = 0; i < n_samples; i++) {
    samples_.push_back(header_->samples[i]);
  }
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
        StrCat("Unknown reference_name ", region.ShortDebugString()));
  }
  if (region.start() < 0 || region.start() >= region.end())
    return tf::errors::InvalidArgument(
        StrCat("Malformed region ", region.ShortDebugString()));

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
  TF_RETURN_IF_ERROR(ConvertToPb(
      header_, bcf1_, reader->Options().desired_format_entries(), out));
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
      return tf::errors::DataLoss(StrCat("Failed to parse VCF record"));
    } else {
      return false;
    }
  }
  const VcfReader* reader = static_cast<const VcfReader*>(reader_);
  TF_RETURN_IF_ERROR(ConvertToPb(
      header_, bcf1_, reader->Options().desired_format_entries(), out));
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
