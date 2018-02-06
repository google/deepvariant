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

#include "deepvariant/util/genomics/struct.pb.h"
#include "deepvariant/util/genomics/variants.pb.h"
#include "deepvariant/util/hts_path.h"
#include "deepvariant/util/math.h"
#include "deepvariant/util/utils.h"
#include "deepvariant/util/vcf_conversion.h"

#include "absl/memory/memory.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/lib/strings/stringprintf.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

namespace tf = tensorflow;

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using tensorflow::strings::StrCat;

namespace {

static const std::map<string, int>* FIELD_TYPE = new std::map<string, int>({
    {"AD", BCF_HT_INT},
    {"DP", BCF_HT_INT},
    {"MIN_DP", BCF_HT_INT},
    {"GQ", BCF_HT_INT},
    {"VAF", BCF_HT_REAL},
});

const char kOpenModeCompressed[] = "wz";
const char kOpenModeUncompressed[] = "w";

const char kFormatHeaderGT[] =
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
const char kFormatHeaderGQ[] =
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
const char kFormatHeaderDP[] =
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth of all "
    "passing filters reads.\">";  // NOLINT
const char kFormatHeaderMINDP[] =
    "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP "
    "observed within the GVCF block.\">";
const char kFormatHeaderAD[] =
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth of all "
    "passing filters reads for each allele.\">";  // NOLINT
const char kFormatHeaderVAF[] =
    "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele "
    "fractions.\">";  // NOLINT
const char kFormatHeaderPL[] =
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods, "
    "Phred encoded\">";  // NOLINT
const char kFormatHeaderGL[] =
    "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods, "
    "log10 encoded\">";  // NOLINT
const char kInfoHeaderEND[] =
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the "
    "interval\">";  // NOLINT

const char kContigHeaderFmt[] = "##contig=<ID=%s,length=%lld>";
const char kFilterHeaderFmt[] = "##FILTER=<ID=%s,Description=\"%s\">";

// Translate the variant call allele encoding used in the protobuf
// message into the htslib VCF constant.
int32_t vcfEncodeAllele(int pbAllele, bool isPhased) {
  // pbAllele is -1 for missing; 0 for ref, 1+ for alt.
  CHECK_GE(pbAllele, -1);
  if (isPhased) {
    return bcf_gt_phased(pbAllele);
  } else {
    return bcf_gt_unphased(pbAllele);
  }
}

// Write out one of the format tags to a VCF variant line.  The input encodes
// the format field values for every sample.  The vector of vectors encodes
// the field values as follows:
//   - if vv is empty, it is taken to mean that there are no format values to be
//     written, and this function call is effectively a no-op;
//   - vv[i] represents the field values for sample at index i; if vv[i] is
//     an empty vector, it means the values are MISSING for this sample.
//   - the subvectors of vv should all be the same length, except for potential
//     empty subvectors
// (This the inverse of the ReadFormatValues utility used by VcfReader)
template <class ValueType>
tf::Status EncodeFormatValues(const std::vector<std::vector<ValueType>>& values,
                              const char* tag, const bcf_hdr_t* h, bcf1_t* v) {
  using VT = VcfType<ValueType>;

  if (values.empty()) {
    return tf::Status::OK();
  }

  if (values.size() != bcf_hdr_nsamples(h))
    return tf::errors::FailedPrecondition("Values.size() != nsamples");
  size_t n_samples = values.size();

  size_t values_per_sample = 0;
  for (size_t s = 0; s < n_samples; s++) {
    values_per_sample =
        std::max(values_per_sample, values[s].size());
  }

  std::vector<ValueType> flat_values;  // Flat-encoded with missing/vector_ends
  for (size_t s = 0; s < n_samples; s++) {
    if (values[s].empty()) {
      for (size_t j = 0; j < values_per_sample; j++) {
        flat_values.push_back(ValueType());
        if (j == 0) {
          VT::SetMissing(&flat_values.back());
        } else {
          VT::SetVectorEnd(&flat_values.back());
        }
      }
    } else {
      if (values[s].size() != values_per_sample)
        return tf::errors::FailedPrecondition(
            "values[s].size() != values_per_sample");
      for (size_t j = 0; j < values_per_sample; j++) {
        flat_values.push_back(values[s][j]);
      }
    }
  }
  if (flat_values.size() != n_samples * values_per_sample)
    return tf::errors::FailedPrecondition(
        "flat_values.size() != n_samples * values_per_sample");
  return VT::PutFormatValues(tag, flat_values.data(), flat_values.size(), h, v);
}

// Private helper class for encoding VariantCall.info values in VCF FORMAT
// field values.
//
// The standard was to interact with this class is:
//
// Create an adaptor for an .info map key "DP".
// FormatFieldAdapter adapter("DP");
//
// For each variant, we check if we need to encode values.
// if (adapter.IsPresentInAnyVariantCalls(variant)) {
//   // And if so, encode them into our htslib bcf_record passing in the
//   // required htslib header object as well.
//   EncodeValues(variant, header, bcf_record);
// }
//
// redacted
// code to other data types. This code isn't templated on purpose, because we
// don't in fact know the proper types in C++ at compile time. We need to look
// up the type in the VCF header, and handle the including datatype based on the
// info in the header. We should look up information in the header when these
// are constructed, and then go through our VariantCall fields to fetch the
// appropriate data and get the appropriate values from our info fields here.
class FormatFieldAdapter {
 public:
  // Creates a new adapter for a field name field_name.
  explicit FormatFieldAdapter(const string& field_name)
      : field_name_(field_name) {}

  // Returns true if the format field field_name occurs in any VariantCall.info
  // maps present in Variant.
  bool IsPresentInAnyVariantCalls(const Variant& variant) const {
    for (const VariantCall& vc : variant.calls()) {
      if (vc.info().find(field_name_) != vc.info().end()) return true;
    }
    return false;
  }

  // Adds the values for our field_name from variant's calls into our bcf1_t
  // record bcf_record.
  //
  // This function supports both single value (DP) and multi-value (AD) info
  // fields.
  //
  // This function only works if the values in the VariantCall info maps don't
  // need to be modified in any way before adding them to the bcf_record. For
  // example, if DP in the VariantCall.info["DP"] map is 10, then we will write
  // a 10 in the bcf_record for the DP field. An example of a field that isn't
  // support by this class is the repeated genotype_likelihood field in the
  // VariantCall proto. These aren't stored in the info field, but are inlined
  // directly in the proto, even though they are written to the FORMAT field
  // in VCF.
  //
  // WARNING: This code currently assumes that all field values are real
  // numbers for "VAF" field. For all other field, it assumes the field values
  // are integers.
  template <class T>
  tf::Status EncodeValues(const Variant& variant, const bcf_hdr_t* header,
                          uint32_t header_type, bcf1_t* bcf_record) const {
    const int id = bcf_hdr_id2int(header, BCF_HL_FLT, field_name_.c_str());
    if (id < 0)
      return tf::errors::FailedPrecondition(
          StrCat("Field ", field_name_, " not in the VCF header"));
    if (bcf_hdr_id2type(header, BCF_HL_FMT, id) != header_type)
      return tf::errors::FailedPrecondition(
          StrCat("Field ", field_name_,
                 " isn't the right type. Expected header_type: ", header_type));

    const int n_calls = variant.calls().size();
    std::vector<std::vector<T>> values(n_calls, std::vector<T>{});
    for (int i = 0; i < n_calls; ++i) {
      const VariantCall& vc = variant.calls(i);
      auto found = vc.info().find(field_name_);
      if (found != vc.info().end()) {
        for (auto& list_value : (*found).second.values()) {
          values[i].push_back(list_value.number_value());
        }
      }
      // Since we don't have a field_name_ key/value pair in this sample, we
      // just leave the values[i] empty.
    }

    // Encode the values from our vector into the htslib bcf_t record.
    return EncodeFormatValues(values, field_name_.c_str(), header, bcf_record);
  }

  tf::Status EncodeValues(const Variant& variant, const bcf_hdr_t* header,
                          bcf1_t* bcf_record) const {
    if (FIELD_TYPE->count(field_name_) <= 0)
      return tf::errors::FailedPrecondition(
          StrCat("Field ", field_name_, " not in the FIELD_TYPE map."));
    int header_type = FIELD_TYPE->find(field_name_)->second;
    if (header_type == BCF_HT_REAL) {
      return EncodeValues<float>(variant, header, header_type, bcf_record);
    } else if (header_type == BCF_HT_INT) {
      return EncodeValues<int>(variant, header, header_type, bcf_record);
    } else {
      return tf::errors::FailedPrecondition(
          StrCat("Unrecognized type for field ", field_name_));
    }
    return tf::Status::OK();
  }

  // The name of our field, such as "DP", "AD", or "VAF".
  const string field_name_;
};

// Parses the Variant protobuf variant_message into an htslib VCF record
tf::Status ConvertFromPb(const Variant& variant_message, const bcf_hdr_t& h,
                         bcf1_t* v) {
  CHECK(v != nullptr) << "bcf1_t record cannot be null";

  v->rid = bcf_hdr_name2id(&h, variant_message.reference_name().c_str());
  if (v->rid < 0)
    return tf::errors::NotFound(
        "Record's reference name is not available in VCF header.");

  v->pos = variant_message.start();
  v->rlen = variant_message.end() - variant_message.start();

  // vcf_format on its own is not properly outputting the END
  // descriptor for GVCF records.  This is arguably a bug in htslib.
  // See b/62297987 for discussion.  For now we use this workaround.

  // Workaround: if our ALT is just a single symbolic allele, this is
  // a gVCF record and we need to populate END, because htslib won't.
  if (variant_message.alternate_bases_size() == 1 &&
      !variant_message.alternate_bases()[0].empty() &&
      variant_message.alternate_bases()[0][0] == '<') {
    int end = v->pos + v->rlen;
    if (bcf_update_info_int32(&h, v, "END", &end, 1) != 0)
      return tf::errors::Unknown("Failure to write END to vcf record");
  }

  // Some variants don't have names; these will get the placeholder "." in the
  // ID column
  if (variant_message.names_size() > 0) {
    string combined_names =
        tensorflow::str_util::Join(variant_message.names(), ";");
    v->d.id = strdup(combined_names.c_str());
  }

  // QUAL
  v->qual = variant_message.quality();

  // Alleles
  int nAlleles = 1 + variant_message.alternate_bases_size();
  auto alleles = absl::make_unique<const char*[]>(nAlleles);
  alleles.get()[0] = variant_message.reference_bases().c_str();
  for (int i = 1; i < nAlleles; i++) {
    alleles.get()[i] = variant_message.alternate_bases(i - 1).c_str();
  }
  bcf_update_alleles(&h, v, alleles.get(), nAlleles);

  // FILTER
  int nFilters = variant_message.filter_size();
  if (nFilters > 0) {
    auto filterIds = absl::make_unique<int32_t[]>(nFilters);
    for (int i = 0; i < nFilters; i++) {
      const char* filterName = variant_message.filter(i).c_str();
      int32_t filterId = bcf_hdr_id2int(&h, BCF_DT_ID, filterName);
      if (filterId < 0) {
        return tf::errors::NotFound("Filter must be found in header.");
      }
      filterIds.get()[i] = filterId;
    }
    bcf_update_filter(&h, v, filterIds.get(), nFilters);
  }

  // Variant calls
  int nCalls = variant_message.calls().size();
  if (nCalls != bcf_hdr_nsamples(&h))
    return tf::errors::FailedPrecondition(
        "Variant call count must match number of samples.");

  // We need to determine the effective ploidy (as the max number of GT calls
  // among samples at this variant); any genotypes shorter than this ploidy
  // will be padded (consult the VCF spec).
  int ploidy = 0;
  for (int c = 0; c < nCalls; c++) {
    ploidy = std::max(ploidy, variant_message.calls(c).genotype_size());
  }

  if (nCalls > 0) {
    // Write genotypes.
    auto gts = absl::make_unique<int32_t[]>(nCalls * ploidy);
    for (int c = 0; c < nCalls; c++) {
      const VariantCall& vc = variant_message.calls(c);

      if (vc.genotype_size() > ploidy)
        return tf::errors::FailedPrecondition(
            "Too many genotypes given the ploidy");
      if (vc.call_set_name() != h.samples[c])
        return tf::errors::FailedPrecondition(
            "Out-of-order call set names, or unrecognized call set name, with "
            "respect to samples declared in VCF header.");

      const bool isPhased = !vc.phaseset().empty();
      int a = 0;
      for (; a < vc.genotype_size(); a++) {
        gts.get()[c * ploidy + a] = vcfEncodeAllele(vc.genotype(a), isPhased);
      }
      for (; a < ploidy; a++) {
        gts.get()[c * ploidy + a] = bcf_int32_vector_end;
      }
    }
    bcf_update_genotypes(&h, v, gts.get(), nCalls * ploidy);

    // Write remaining FORMAT fields
    bool has_ll = false;
    for (int c = 0; c < nCalls; c++) {
      const VariantCall& vc = variant_message.calls(c);
      if (vc.genotype_likelihood_size() > 0) {
        has_ll = true;
      }
    }

    // The order of adapter definitions here determines the order of the fields
    // in the VCF.
    std::vector<FormatFieldAdapter> format_field_adapters = {
        FormatFieldAdapter("GQ"),     FormatFieldAdapter("DP"),
        FormatFieldAdapter("MIN_DP"), FormatFieldAdapter("AD"),
        FormatFieldAdapter("VAF"),
    };
    for (const FormatFieldAdapter& field : format_field_adapters) {
      if (field.IsPresentInAnyVariantCalls(variant_message)) {
        TF_RETURN_IF_ERROR(field.EncodeValues(variant_message, &h, v));
      }
    }

    std::vector<std::vector<int>> ll_values_phred;
    for (int c = 0; c < nCalls; c++) {
      const VariantCall& vc = variant_message.calls(c);
      // loglikelihood
      if (has_ll) {
        bool ll_not_missing = vc.genotype_likelihood_size() > 0;
        if (ll_not_missing) {
          std::vector<double> lls_this_call;
          for (double ll_val : vc.genotype_likelihood()) {
            lls_this_call.push_back(ll_val);
          }
          // "Normalize" likelihoods...
          std::vector<double> lls_this_call_normalized =
              ZeroShiftLikelihoods(lls_this_call);

          // Phred-transform them...
          std::vector<int> phreds_this_call(lls_this_call.size());
          std::transform(lls_this_call_normalized.cbegin(),
                         lls_this_call_normalized.cend(),
                         phreds_this_call.begin(), Log10PErrorToPhred);
          ll_values_phred.push_back(phreds_this_call);
        } else {
          ll_values_phred.push_back({});
        }
      }
    }

    TF_RETURN_IF_ERROR(EncodeFormatValues(ll_values_phred, "PL", &h, v));
  }
  return tf::Status::OK();
}

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
  for (const VcfFilterInfo& filter : options.filters()) {
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

  for (const ContigInfo& contig : options.contigs()) {
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
  if (options_.round_qual_values()) {
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
