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

// Implementation of vcf_writer.h
#include "deepvariant/util/vcf_writer.h"

#include <algorithm>
#include <cstring>

#include "deepvariant/util/genomics/reference.pb.h"
#include "deepvariant/util/genomics/variants.pb.h"
#include "deepvariant/util/hts_path.h"
#include "deepvariant/util/utils.h"
#include "deepvariant/util/vcf_conversion.h"

#include "absl/memory/memory.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/lib/strings/stringprintf.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

namespace tf = tensorflow;

using nucleus::genomics::v1::Variant;
using tensorflow::strings::StrCat;

namespace {

constexpr char kOpenModeCompressed[] = "wz";
constexpr char kOpenModeUncompressed[] = "w";

constexpr char kFormatHeaderGT[] =
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
constexpr char kFormatHeaderGQ[] =
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
constexpr char kFormatHeaderDP[] =
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth of all "
    "passing filters reads.\">";  // NOLINT
constexpr char kFormatHeaderMINDP[] =
    "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP "
    "observed within the GVCF block.\">";
constexpr char kFormatHeaderAD[] =
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth of all "
    "passing filters reads for each allele.\">";  // NOLINT
constexpr char kFormatHeaderVAF[] =
    "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele "
    "fractions.\">";  // NOLINT
constexpr char kFormatHeaderPL[] =
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods, "
    "Phred encoded\">";  // NOLINT
constexpr char kFormatHeaderGL[] =
    "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods, "
    "log10 encoded\">";  // NOLINT
constexpr char kInfoHeaderEND[] =
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the "
    "interval\">";  // NOLINT

constexpr char kContigHeaderFmt[] = "##contig=<ID=%s,length=%lld>";
constexpr char kFilterHeaderFmt[] = "##FILTER=<ID=%s,Description=\"%s\">";

}  // namespace

StatusOr<std::unique_ptr<VcfWriter>> VcfWriter::ToFile(
    const string& variants_path, const VcfWriterOptions& options) {
  const char* const openMode = EndsWith(variants_path, ".gz")
                                   ? kOpenModeCompressed
                                   : kOpenModeUncompressed;
  htsFile* fp = hts_open_x(variants_path.c_str(), openMode);
  if (fp == nullptr)
    return tf::errors::Unknown(
        StrCat("Could not open variants_path ", variants_path));

  auto writer = absl::WrapUnique(new VcfWriter(options, fp));
  TF_RETURN_IF_ERROR(writer->WriteHeader());
  return std::move(writer);
}

VcfWriter::VcfWriter(const VcfWriterOptions& options, htsFile* fp)
    : fp_(fp), options_(options) {
  CHECK(fp != nullptr);

  header_ = bcf_hdr_init("w");
  for (const nucleus::genomics::v1::VcfFilterInfo& filter : options.filters()) {
    string filterStr = tensorflow::strings::Printf(
        kFilterHeaderFmt, filter.id().c_str(), filter.description().c_str());
    bcf_hdr_append(header_, filterStr.c_str());
  }

  bcf_hdr_append(header_, kFormatHeaderGT);
  bcf_hdr_append(header_, kFormatHeaderGQ);
  bcf_hdr_append(header_, kFormatHeaderDP);
  bcf_hdr_append(header_, kFormatHeaderMINDP);
  bcf_hdr_append(header_, kFormatHeaderAD);
  bcf_hdr_append(header_, kFormatHeaderVAF);
  bcf_hdr_append(header_, kFormatHeaderGL);
  bcf_hdr_append(header_, kFormatHeaderPL);
  bcf_hdr_append(header_, kInfoHeaderEND);

  for (const nucleus::genomics::v1::ContigInfo& contig : options.contigs()) {
    string ctgStr =
        tensorflow::strings::Printf(kContigHeaderFmt, contig.name().c_str(),
                                    static_cast<int64>(contig.n_bases()));
    bcf_hdr_append(header_, ctgStr.c_str());
  }

  for (const string& sampleName : options.sample_names()) {
    bcf_hdr_add_sample(header_, sampleName.c_str());
  }
  bcf_hdr_add_sample(header_, nullptr);
}

tf::Status VcfWriter::WriteHeader() {
  if (bcf_hdr_write(fp_, header_) < 0)
    return tf::errors::Unknown("Failed to write header");
  else
    return tf::Status::OK();
}

VcfWriter::~VcfWriter() {
  if (fp_) {
    // There's nothing we can do but assert fail if there's an error during
    // the Close() call here.
    TF_CHECK_OK(Close());
  }
}

tf::Status VcfWriter::Write(const Variant& variant_message) {
  if (fp_ == nullptr)
    return tf::errors::FailedPrecondition("Cannot write to closed VCF stream.");
  bcf1_t* v = bcf_init();
  if (v == nullptr)
    return tf::errors::Unknown("bcf_init call failed");
  TF_RETURN_IF_ERROR(ConvertFromPb(variant_message, *header_, v));
  if (options_.round_qual_values() && !bcf_float_is_missing(v->qual)) {
    // Round quality value printed out to one digit past the decimal point.
    double rounded_quality = floor(variant_message.quality() * 10 + 0.5) / 10;
    v->qual = rounded_quality;
  }
  if (bcf_write(fp_, header_, v) != 0)
    return tf::errors::Unknown("bcf_write call failed");
  bcf_destroy(v);
  return tf::Status::OK();
}

tf::Status VcfWriter::Close() {
  if (fp_ == nullptr)
    return tf::errors::FailedPrecondition(
        "Cannot close an already closed VcfWriter");
  if (hts_close(fp_) < 0)
    return tf::errors::Unknown("hts_close call failed");
  fp_ = nullptr;
  bcf_hdr_destroy(header_);
  header_ = nullptr;
  return tf::Status::OK();
}

}  // namespace nucleus
