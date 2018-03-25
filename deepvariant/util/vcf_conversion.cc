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

#include "deepvariant/util/vcf_conversion.h"

#include "absl/memory/memory.h"
#include "deepvariant/util/math.h"
#include "deepvariant/util/utils.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/lib/strings/strcat.h"

namespace nucleus {

namespace {

// -----------------------------------------------------------------------------
// "Raw" low-level interface to encoding/decoding to FORMAT fields.  We use
// these directly for FORMAT fields that have special semantics and so cannot be
// handled by VcfFormatFieldAdapter.

// Read in one of the format tags for a variant line of a VCF file, and return
// a vector of vectors, where each subvector records values for a single sample.
// An empty subvector represents that the field is "missing" for this sample.
// If the tag is not represented in this variant, a zero-length vector is
// returned.
template <class ValueType>
std::vector<std::vector<ValueType>> ReadFormatValues(const bcf_hdr_t* h,
                                                     const bcf1_t* v,
                                                     const char* tag) {
  using VT = VcfType<ValueType>;

  if (bcf_get_fmt(h, const_cast<bcf1_t*>(v), tag) == nullptr) {
    return std::vector<std::vector<ValueType>>();
  }

  int n_values, n_dst = 0;
  ValueType* dst = nullptr;
  n_values = VT::GetFormatValues(h, v, tag, &dst, &n_dst);
  // redacted
  CHECK_GE(n_values, 0);
  CHECK(dst != nullptr);

  int n_values_per_sample = n_values / v->n_sample;
  std::vector<std::vector<ValueType>> values(v->n_sample);

  for (int i = 0; i < v->n_sample; i++) {
    for (int j = 0; j < n_values_per_sample; j++) {
      ValueType* p = dst + n_values_per_sample * i + j;
      if (VT::IsVectorEnd(*p)) break;
      // We only support fields that are entirely missing, so if we encounter a
      // single missing field, we clear this subvector and break.
      if (VT::IsMissing(*p)) {
        values[i].clear();
        break;
      }
      values[i].push_back(*p);
    }
  }

  free(dst);
  return values;
}

// Sentinel value used to set variant.quality if one was not specified.
constexpr double kQualUnset = -1;

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
// (This the inverse of ReadFormatValues)
template <class ValueType>
tensorflow::Status EncodeFormatValues(
    const std::vector<std::vector<ValueType>>& values, const char* tag,
    const bcf_hdr_t* h, bcf1_t* v) {
  using VT = VcfType<ValueType>;

  if (values.empty()) {
    return tensorflow::Status::OK();
  }

  if (values.size() != bcf_hdr_nsamples(h))
    return tensorflow::errors::FailedPrecondition("Values.size() != nsamples");
  size_t n_samples = values.size();

  size_t values_per_sample = 0;
  for (size_t s = 0; s < n_samples; s++) {
    values_per_sample = std::max(values_per_sample, values[s].size());
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
        return tensorflow::errors::FailedPrecondition(
            "values[s].size() != values_per_sample");
      for (size_t j = 0; j < values_per_sample; j++) {
        flat_values.push_back(values[s][j]);
      }
    }
  }
  if (flat_values.size() != n_samples * values_per_sample)
    return tensorflow::errors::FailedPrecondition(
        "flat_values.size() != n_samples * values_per_sample");
  return VT::PutFormatValues(tag, flat_values.data(), flat_values.size(), h, v);
}


}  // namespace

// -----------------------------------------------------------------------------
// VcfFormatFieldAdapter implemenation.

VcfFormatFieldAdapter::VcfFormatFieldAdapter(const string& field_name,
                                             int vcf_type)
    : field_name_(field_name), vcf_type_(vcf_type) {}


tensorflow::Status VcfFormatFieldAdapter::EncodeValues(
    const nucleus::genomics::v1::Variant& variant,
    const bcf_hdr_t* header,
    bcf1_t* bcf_record) const {

  if (vcf_type_ == BCF_HT_REAL) {
    return EncodeValues<float>(variant, header, bcf_record);
  } else if (vcf_type_ == BCF_HT_INT) {
    return EncodeValues<int>(variant, header, bcf_record);
  } else {
    // redacted
    return tensorflow::errors::FailedPrecondition(
        "Unrecognized type for field ", field_name_);
  }
  return tensorflow::Status::OK();
}

// redacted
// the intermediate vectors contain variant objects (Value).
template <class T> tensorflow::Status VcfFormatFieldAdapter::EncodeValues(
    const nucleus::genomics::v1::Variant& variant,
    const bcf_hdr_t* header,
    bcf1_t* bcf_record) const {

  const int n_calls = variant.calls().size();
  std::vector<std::vector<T>> values(n_calls, std::vector<T>{});
  for (int i = 0; i < n_calls; ++i) {
    const nucleus::genomics::v1::VariantCall& vc = variant.calls(i);
    auto found = vc.info().find(field_name_);
    if (found != vc.info().end()) {
      for (auto& list_value : (*found).second.values()) {
        if (std::is_integral<T>::value) {
          values[i].push_back(list_value.int_value());
        } else {
          values[i].push_back(list_value.number_value());
        }
      }
    }
    // Since we don't have a field_name_ key/value pair in this sample, we
    // just leave the values[i] empty.
  }

  // Encode the values from our vector into the htslib bcf_t record.
  return EncodeFormatValues(values, field_name_.c_str(), header, bcf_record);
}


tensorflow::Status VcfFormatFieldAdapter::DecodeValues(
    const bcf_hdr_t* header, const bcf1_t* bcf_record,
    nucleus::genomics::v1::Variant* variant) const {

  if (vcf_type_ == BCF_HT_REAL) {
    return DecodeValues<float>(header, bcf_record, variant);
  } else if (vcf_type_ == BCF_HT_INT) {
    return DecodeValues<int>(header, bcf_record, variant);
  } else {
    // redacted
    return tensorflow::errors::FailedPrecondition(
        "Unrecognized type for field ", field_name_);
  }
  return tensorflow::Status::OK();
}

template <class T> tensorflow::Status VcfFormatFieldAdapter::DecodeValues(
    const bcf_hdr_t *header, const bcf1_t *bcf_record,
    nucleus::genomics::v1::Variant *variant) const {

  if (bcf_record->n_sample > 0) {
    std::vector<std::vector<T>> values =
        ReadFormatValues<T>(header, bcf_record, field_name_.c_str());
    for (int i = 0; i < bcf_record->n_sample; i++) {
      // Is the format field present for this variant, *and* non-missing for
      // this sample?
      bool have_field = !values.empty() && !values[i].empty();
      if (have_field) {
        nucleus::genomics::v1::VariantCall* call = variant->mutable_calls(i);
        SetInfoField(field_name_, values[i], call);
      }
    }
  }
  return tensorflow::Status::OK();
}

// -----------------------------------------------------------------------------
// VcfRecordConverter implementation.

VcfRecordConverter::VcfRecordConverter(
    const nucleus::genomics::v1::VcfHeader& vcf_header,
    const OptionalVariantFieldsToParse& desired_format_entries)
    : desired_format_entries_(desired_format_entries)
{
  for (const auto& format_spec : vcf_header.formats()) {
    string tag = format_spec.id();
    string type = format_spec.type();

    // These fields are handled specially.
    if (tag == "GT" || tag == "GL" || tag == "PL" || tag == "PS") continue;

    // Check if configuration has disabled this FORMAT field.
    if ((tag == "GQ" && desired_format_entries_.exclude_genotype_quality()) ||
        (tag == "AD" && desired_format_entries_.exclude_allele_depth()) ||
        (tag == "DP" && desired_format_entries_.exclude_read_depth()) ||
        (tag == "VAF" &&
         desired_format_entries_.exclude_variant_allele_frequencies()) ||
        (tag == "MIN_DP" && desired_format_entries_.exclude_allele_depth())) {
      continue;
    }

    // redacted
    int vcf_type;
    if (type == "Integer") {
      vcf_type = BCF_HT_INT;
    } else if (type == "Float") {
      vcf_type = BCF_HT_REAL;
    } else {
      // redacted
      LOG(WARNING) << "Unhandled FORMAT field type: field " << tag
                   << " of type " << type;
      continue;
    }
    format_adapters_.emplace_back(tag, vcf_type);
  }
}


tensorflow::Status VcfRecordConverter::ConvertToPb(
    const bcf_hdr_t* h, bcf1_t* v,
    nucleus::genomics::v1::Variant* variant_message) const {

  CHECK(h != nullptr) << "BCF header cannot be null";
  CHECK(v != nullptr) << "bcf1_t record cannot be null";
  CHECK(variant_message != nullptr) << "variant_message record cannot be null";

  variant_message->Clear();

  // Tell htslib to parse out all of the fields of the VCF record v.
  bcf_unpack(v, BCF_UN_ALL);

  variant_message->set_reference_name(bcf_hdr_id2name(h, v->rid));
  variant_message->set_start(v->pos);
  variant_message->set_end(v->pos + v->rlen);

  // Parse the ID field of the Variant.
  if (v->d.id && strcmp(v->d.id, ".") != 0) {
    // Don't add the missing "." marker to the id field.
    variant_message->add_names(v->d.id);
  }

  // Parse out the ref and alt alleles.
  if (v->n_allele > 0) {
    variant_message->set_reference_bases(v->d.allele[0]);
    for (int i = 1; i < v->n_allele; ++i) {
      variant_message->add_alternate_bases(v->d.allele[i]);
    }
  }

  // Parse out the QUAL field. QUAL is the only field where the unset default
  // proto value is a valid value, so we have to explicitly populate the result
  // with a sentinel.
  if (bcf_float_is_missing(v->qual)) {
    variant_message->set_quality(kQualUnset);
  } else {
    variant_message->set_quality(v->qual);
  }

  // Parse out the FILTER field.
  for (int i = 0; i < v->d.n_flt; ++i) {
    variant_message->add_filter(h->id[BCF_DT_ID][v->d.flt[i]].key);
  }

  bool want_ll = !desired_format_entries_.exclude_genotype_likelihood();
  bool want_gt = !desired_format_entries_.exclude_genotype();

  // Parse the calls of the variant.
  if (v->n_sample > 0) {
    int* gt_arr = nullptr;
    int ploidy, n_gts = 0;
    if (bcf_get_genotypes(h, v, &gt_arr, &n_gts) < 0) {
      free(gt_arr);
      return tensorflow::errors::DataLoss("Couldn't parse genotypes");
    }
    ploidy = n_gts / v->n_sample;

    std::vector<std::vector<int>> ps_values = ReadFormatValues<int>(h, v, "PS");

    for (int i = 0; i < v->n_sample; i++) {
      bool have_ps = !ps_values.empty() && !ps_values[i].empty();
      nucleus::genomics::v1::VariantCall* call = variant_message->add_calls();
      call->set_call_set_name(h->samples[i]);
      // Get the GT calls, if requested and available.
      if (want_gt) {
        bool gt_is_phased = false;
        for (int j = 0; j < ploidy; j++) {
          int gt_idx = gt_arr[i * ploidy + j];
          int gt = bcf_gt_allele(gt_idx);
          gt_is_phased = gt_is_phased || bcf_gt_is_phased(gt_idx);
          call->add_genotype(gt);
        }
        // We need to set the phaseset, if available in field PS. If this is a
        // phased genotype and there is no PS field, phaseset is "*".
        if (have_ps) {
          call->set_phaseset(std::to_string(ps_values[i][0]));
        } else if (gt_is_phased) {
          call->set_phaseset("*");
        }
      }
    }
    free(gt_arr);

    // Parse "generic" FORMAT fields.
    for (const auto& adapter : format_adapters_) {
      TF_RETURN_IF_ERROR(adapter.DecodeValues(h, v, variant_message));
    }

    // Handle FORMAT fields requiring special logic.
    std::vector<std::vector<int>> pl_values = ReadFormatValues<int>(h, v, "PL");
    std::vector<std::vector<float>> gl_values =
        ReadFormatValues<float>(h, v, "GL");

    for (int i = 0; i < v->n_sample; i++) {
      // Each indicator here is true iff the format field is present for this
      // variant, *and* is non-missing for this sample.
      bool have_gl = !gl_values.empty() && !gl_values[i].empty();
      bool have_pl = !pl_values.empty() && !pl_values[i].empty();

      nucleus::genomics::v1::VariantCall* call =
          variant_message->mutable_calls(i);

      if (want_ll) {
        // If GL and PL are *both* present, we populate the genotype_likelihood
        // fields with the GL values per the variants.proto spec, since PLs are
        // a lower resolution version of the same information.
        if (have_gl) {
          for (int gl : gl_values[i]) {
            call->add_genotype_likelihood(gl);
          }
        } else if (have_pl) {
          for (int pl : pl_values[i]) {
            call->add_genotype_likelihood(PhredToLog10PError(pl));
          }
        }
      }
    }
  }
  return tensorflow::Status::OK();
}


tensorflow::Status VcfRecordConverter::ConvertFromPb(
    const nucleus::genomics::v1::Variant& variant_message, const bcf_hdr_t& h,
    bcf1_t* v) const {

  CHECK(v != nullptr) << "bcf1_t record cannot be null";

  v->rid = bcf_hdr_name2id(&h, variant_message.reference_name().c_str());
  if (v->rid < 0)
    return tensorflow::errors::NotFound(
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
      return tensorflow::errors::Unknown("Failure to write END to vcf record");
  }

  // Some variants don't have names; these will get the placeholder "." in the
  // ID column
  if (variant_message.names_size() > 0) {
    string combined_names =
        tensorflow::str_util::Join(variant_message.names(), ";");
    v->d.id = strdup(combined_names.c_str());
  }

  // QUAL
  if (variant_message.quality() == kQualUnset) {
    bcf_float_set_missing(v->qual);
  } else {
    v->qual = variant_message.quality();
  }

  // Alleles
  int nAlleles = 1 + variant_message.alternate_bases_size();
  auto alleles = absl::make_unique<const char* []>(nAlleles);
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
        return tensorflow::errors::NotFound("Filter must be found in header.");
      }
      filterIds.get()[i] = filterId;
    }
    bcf_update_filter(&h, v, filterIds.get(), nFilters);
  }

  // Variant calls
  int nCalls = variant_message.calls().size();
  if (nCalls != bcf_hdr_nsamples(&h))
    return tensorflow::errors::FailedPrecondition(
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
      const nucleus::genomics::v1::VariantCall& vc = variant_message.calls(c);

      if (vc.genotype_size() > ploidy)
        return tensorflow::errors::FailedPrecondition(
            "Too many genotypes given the ploidy");
      if (vc.call_set_name() != h.samples[c])
        return tensorflow::errors::FailedPrecondition(
          "Out-of-order call set names, or unrecognized call set name, "
          "with respect to samples declared in VCF header. Variant has ",
          vc.call_set_name(), " at position ", c,
          " while the VCF header expected a sample named ",
          h.samples[c], " at this position");

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
      const nucleus::genomics::v1::VariantCall& vc = variant_message.calls(c);
      if (vc.genotype_likelihood_size() > 0) {
        has_ll = true;
      }
    }

    for (const VcfFormatFieldAdapter& field : format_adapters_) {
      TF_RETURN_IF_ERROR(field.EncodeValues(variant_message, &h, v));
    }

    std::vector<std::vector<int>> ll_values_phred;
    for (int c = 0; c < nCalls; c++) {
      const nucleus::genomics::v1::VariantCall& vc = variant_message.calls(c);
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
  return tensorflow::Status::OK();
}

}  // namespace nucleus
