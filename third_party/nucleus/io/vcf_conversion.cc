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

#include "third_party/nucleus/io/vcf_conversion.h"

#include <ctype.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <memory>

#include "google/protobuf/map.h"
#include "google/protobuf/repeated_field.h"
#include "absl/memory/memory.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "absl/strings/substitute.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/util/math.h"
#include "third_party/nucleus/util/utils.h"

namespace nucleus {

namespace {

constexpr char kFilterHeaderFmt[] = "##FILTER=<ID=$0,Description=\"$1\">";
constexpr char kInfoHeaderFmt[] =
    "##INFO=<ID=$0,Number=$1,Type=$2,Description=\"$3\"$4>";
constexpr char kFormatHeaderFmt[] =
    "##FORMAT=<ID=$0,Number=$1,Type=$2,Description=\"$3\">";
constexpr char kContigHeaderFmt[] = "##contig=<ID=$0$1>";
constexpr char kStructuredExtraHeaderFmt[] = "##$0=<$1>";
constexpr char kExtraHeaderFmt[] = "##$0=$1";

// Adds a FILTER field to the bcf_hdr_t header based on the VcfFilterInfo
// object.
void AddFilterToHeader(const nucleus::genomics::v1::VcfFilterInfo& filter,
                       bcf_hdr_t* header) {
  string filterStr = absl::Substitute(
      kFilterHeaderFmt, filter.id().c_str(), filter.description().c_str());
  bcf_hdr_append(header, filterStr.c_str());
}

// Adds an INFO field to the bcf_hdr_t header based on the VcfInfo object.
void AddInfoToHeader(const nucleus::genomics::v1::VcfInfo& info,
                     bcf_hdr_t* header) {
  string extra;
  if (!info.source().empty()) {
    absl::StrAppend(&extra, ",Source=\"", info.source(), "\"");
  }
  if (!info.version().empty()) {
    absl::StrAppend(&extra, ",Version=\"", info.version(), "\"");
  }
  string infoStr = absl::Substitute(
      kInfoHeaderFmt, info.id().c_str(), info.number().c_str(),
      info.type().c_str(), info.description().c_str(), extra.c_str());
  bcf_hdr_append(header, infoStr.c_str());
}

// Adds a FORMAT field to the bcf_hdr_t header based on the VcfFormatInfo
// object.
void AddFormatToHeader(const nucleus::genomics::v1::VcfFormatInfo& format,
                       bcf_hdr_t* header) {
  string formatStr = absl::Substitute(
      kFormatHeaderFmt, format.id().c_str(), format.number().c_str(),
      format.type().c_str(), format.description().c_str());
  bcf_hdr_append(header, formatStr.c_str());
}

// Adds a structured extra field to the bcf_hdr_t header based on the
// VcfStructuredExtra object.
void AddStructuredExtraToHeader(
    const nucleus::genomics::v1::VcfStructuredExtra& sExtra,
    bcf_hdr_t* header) {
  string fieldStr;
  for (auto const& kv : sExtra.fields()) {
    absl::StrAppend(&fieldStr, kv.key(), "=\"", kv.value(), "\",");
  }
  if (!fieldStr.empty()) {
    // Cut off the dangling comma.
    fieldStr.pop_back();
  }
  string result = absl::Substitute(
      kStructuredExtraHeaderFmt, sExtra.key().c_str(), fieldStr.c_str());
  bcf_hdr_append(header, result.c_str());
}

// Adds an unstructured extra field to the bcf_hdr_t header based on the
// VcfExtra object.
void AddExtraToHeader(const nucleus::genomics::v1::VcfExtra& extra,
                      bcf_hdr_t* header) {
  string result = absl::Substitute(
      kExtraHeaderFmt, extra.key().c_str(), extra.value().c_str());
  bcf_hdr_append(header, result.c_str());
}

// Adds a contig field to the bcf_hdr_t header based on the ContigInfo object.
void AddContigToHeader(const nucleus::genomics::v1::ContigInfo& contig,
                       bcf_hdr_t* header) {
  string extra;
  if (contig.n_bases()) {
    absl::StrAppend(&extra, ",length=", contig.n_bases());
  }
  if (!contig.description().empty()) {
    absl::StrAppend(&extra, ",description=\"", contig.description(), "\"");
  }
  for (auto const& kv : contig.extra()) {
    absl::StrAppend(&extra, ",", kv.first, "=\"", kv.second, "\"");
  }
  string contigStr = absl::Substitute(
      kContigHeaderFmt, contig.name().c_str(), extra.c_str());
  bcf_hdr_append(header, contigStr.c_str());
}

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
  if (dst == nullptr) {
    LOG(WARNING) << "Error reading format values (dst == nullptr) for tag "
                 << tag;
    return std::vector<std::vector<ValueType>>();
  }
  if (n_values < 0) {
    LOG(WARNING) << "Error reading format values (n_values < 0) for tag "
                 << tag;
    free(dst);
    return std::vector<std::vector<ValueType>>();
  }

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


// Specialized instantiation for string fields, which require different
// memory management and semantics.
template<>
std::vector<std::vector<string>> ReadFormatValues(const bcf_hdr_t* h,
                                                  const bcf1_t* v,
                                                  const char* tag) {
  if (bcf_get_fmt(h, const_cast<bcf1_t*>(v), tag) == nullptr) {
    return {};
  }

  int n_dst = 0;
  char** dst = nullptr;
  std::vector<std::vector<string>> values(v->n_sample);
  if (bcf_get_format_string(h, const_cast<bcf1_t*>(v), tag, &dst, &n_dst) > 0) {
    for (int i = 0; i < bcf_hdr_nsamples(h); i++) {
      // TODO: (1) validate the length of this list is as declared in
      // the header, and figure out what to do when declared length is smaller
      // than the actual length of the list.

      // According to https://samtools.github.io/hts-specs/VCFv4.3.pdf
      // Section 6.3.3 strings in VCF cannot contain ',' (a field separator).
      values[i] = absl::StrSplit(dst[i], ',');
    }
    // As noted in bcf_get_format_string declaration in vcf.h, the format
    // function we are using here allocates two arrays and both must be cleaned
    // by the user.
    free(dst[0]);
    free(dst);
  }
  return values;
}

// Sentinel value used to set variant.quality if one was not specified.
constexpr double kQualUnset = -1;

// Translate the variant call allele encoding used in the protobuf
// message into the htslib VCF constant.
int32 vcfEncodeAllele(int pbAllele, bool isPhased) {
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

  if (static_cast<int>(values.size()) != bcf_hdr_nsamples(h)) {
    return tensorflow::errors::FailedPrecondition("Values.size() != nsamples");
  }
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


// Specialized instantiation for string.
template<>
tensorflow::Status EncodeFormatValues(
    const std::vector<std::vector<string>>& values, const char* tag,
    const bcf_hdr_t* h, bcf1_t* v) {

  if (values.empty()) {
    return tensorflow::Status::OK();
  }

  if (static_cast<int>(values.size()) != bcf_hdr_nsamples(h)) {
    return tensorflow::errors::FailedPrecondition("Values.size() != nsamples");
  }
  size_t n_samples = values.size();

  size_t values_per_sample = 0;
  for (size_t s = 0; s < n_samples; s++) {
    values_per_sample = std::max(values_per_sample, values[s].size());
  }
  if (values_per_sample > 1) {
    return tensorflow::errors::FailedPrecondition(
        "Can't currently handle > 1 string format entry per sample.");
  }
  auto c_values_up = absl::make_unique<const char*[]>(n_samples);
  const char** c_values = c_values_up.get();
  CHECK(c_values != nullptr);
  for (size_t i = 0; i < n_samples; ++i) {
    if (values[i].empty()) {
      c_values[i] = ".";
    } else {
      c_values[i] = values[i][0].c_str();
    }
  }
  int rc = bcf_update_format_string(h, v, tag, c_values,
                                    values_per_sample * n_samples);
  if (rc < 0) {
    return tensorflow::errors::Internal(
        "Failure to write VCF FORMAT field");
  }

  return tensorflow::Status::OK();
}


// -----------------------------------------------------------------------------
// "Raw" low-level interface to encoding/decoding to INFO fields. These
// functions parallel the "FORMAT" functions above.

template <class ValueType>
std::vector<ValueType> ReadInfoValue(const bcf_hdr_t* h,
                                     const bcf1_t* v,
                                     const char* tag) {
  using VT = VcfType<ValueType>;

  if (bcf_get_info(h, const_cast<bcf1_t*>(v), tag) == nullptr) {
    return {};
  }

  int n_values, n_dst = 0;
  ValueType* dst = nullptr;
  n_values = VT::GetInfoValues(h, v, tag, &dst, &n_dst);
  if (dst == nullptr) {
    LOG(WARNING) << "Error reading info (dst == nullptr) value " << tag;
    return {};
  }
  if (n_values < 0) {
    free(dst);
    LOG(WARNING) << "Error reading info (n_values < 0) value " << tag;
    return {};
  }

  std::vector<ValueType> value(dst, dst + n_values);

  free(dst);
  return value;
}

template <>
std::vector<string> ReadInfoValue(const bcf_hdr_t* h,
                                  const bcf1_t* v,
                                  const char* tag) {
  if (bcf_get_info(h, const_cast<bcf1_t*>(v), tag) == nullptr) {
    return {};
  }

  int n_dst = 0;
  char* dst = nullptr;
  if (bcf_get_info_string(h, const_cast<bcf1_t*>(v), tag, &dst, &n_dst) < 0) {
    // TODO: cleanup error handling.
    LOG(FATAL) << "Failure to get INFO string";
  }
  std::string string_value(dst);
  free(dst);
  return std::vector<string> { string_value };
}

template <>
std::vector<bool> ReadInfoValue(const bcf_hdr_t* h,
                                const bcf1_t* v,
                                const char* tag) {
  void* dst;
  int n_dst = 0;
  int rc = bcf_get_info_flag(h, const_cast<bcf1_t*>(v), tag, &dst, &n_dst);
  if (rc == 1) {
    return {true};
  } else if (rc == 0) {
    return {false};
  } else {
    // TODO: cleanup error handling.
    LOG(FATAL) << "Failure to get INFO flag.";
  }
}

template <class ValueType>
tensorflow::Status EncodeInfoValue(
    const std::vector<ValueType>& value, const char* tag,
    const bcf_hdr_t* h, bcf1_t* v) {
  using VT = VcfType<ValueType>;

  if (value.empty()) {
    return tensorflow::Status::OK();
  }
  return VT::PutInfoValues(tag, value.data(), value.size(), h, v);
}

template <>
tensorflow::Status EncodeInfoValue(
    const std::vector<string>& value, const char* tag,
    const bcf_hdr_t* h, bcf1_t* v) {
  if (value.empty()) {
    return tensorflow::Status::OK();
  }
  if (value.size() != 1) {
    return tensorflow::errors::FailedPrecondition(
        "VCF string INFO fields can only contain a single string.");
  }
  const char* string_value = value[0].c_str();
  int rc = bcf_update_info_string(h, v, tag, string_value);
  if (rc < 0) {
    return tensorflow::errors::Internal(
        "Failure to write VCF INFO field");
  }
  return tensorflow::Status::OK();
}

template <>
tensorflow::Status EncodeInfoValue(const std::vector<bool>& value,
                                   const char* tag, const bcf_hdr_t* h,
                                   bcf1_t* v) {
  bool flag_setting;
  if (value.size() == 1 && value[0] == false) {
    flag_setting = false;
  } else if (value.size() == 1 && value[0] == true) {
    flag_setting = true;
  } else {
    return tensorflow::errors::FailedPrecondition(
        "Illegal setting of INFO FLAG value in Variant message.");
  }
  int rc = bcf_update_info_flag(h, v, tag, "", flag_setting);
  if (rc < 0) {
    return tensorflow::errors::Internal(
        "Failure to write VCF INFO field");
  }
  return tensorflow::Status::OK();
}

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
        // TODO: remove string conversion.
        (*contig->mutable_extra())[hrec0->keys[j]] =
            string(Unquote(hrec0->vals[j]));
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
    // TODO: remove string conversion
    filter->set_description(string(Unquote(hrec->vals[1])));
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
    // TODO: remove string conversions below.
    info->set_description(string(Unquote(hrec->vals[3])));
    for (int i = 4; i < hrec->nkeys; i++) {
      if (string(hrec->keys[i]) == "Source") {
        info->set_source(string(Unquote(hrec->vals[i])));
      } else if (string(hrec->keys[i]) == "Version") {
        info->set_version(string(Unquote(hrec->vals[i])));
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
  std::string id, number, type, description;
  for (int i = 0; i < hrec->nkeys; i++) {
    if (string(hrec->keys[i]) == "ID") {
      id = hrec->vals[i];
    } else if (string(hrec->keys[i]) == "Number") {
      number = hrec->vals[i];
    } else if (string(hrec->keys[i]) == "Type") {
      type = hrec->vals[i];
    } else if (string(hrec->keys[i]) == "Description") {
      // TODO: remove string conversions below.
      description = string(Unquote(hrec->vals[i]));
    }
  }
  if (!id.empty() && !number.empty() && !type.empty()) {
    format->set_id(id);
    format->set_number(number);
    format->set_type(type);
    format->set_description(description);
  } else if (!id.empty()) {
    LOG(WARNING) << "Malformed FORMAT field " + id +
                        ": Number or Type definition is missing. Leaving this "
                        "format empty.";
  } else {
    LOG(WARNING) << "Found malformed FORMAT field, all FORMAT fields must "
                    "define ID, Number, and Type. ";
  }
}

// Adds structured information from the bcf_hrec_t to the VcfStructuredExtra.
void AddStructuredExtra(const bcf_hrec_t* hrec,
                        nucleus::genomics::v1::VcfStructuredExtra* extra) {
  extra->set_key(hrec->key);
  for (int i = 0; i < hrec->nkeys; i++) {
    nucleus::genomics::v1::VcfExtra& toAdd = *extra->mutable_fields()->Add();
    toAdd.set_key(hrec->keys[i]);
    // TODO: remove string conversion.
    toAdd.set_value(string(Unquote(hrec->vals[i])));
  }
}

// Adds unstructured information from the bcf_hrec_t to the VcfExtra object.
void AddExtra(const bcf_hrec_t* hrec, nucleus::genomics::v1::VcfExtra* extra) {
  extra->set_key(hrec->key);
  extra->set_value(hrec->value);
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
  } else if (vcf_type_ == BCF_HT_STR) {
    return EncodeValues<string>(variant, header, bcf_record);
  } else {
    return tensorflow::errors::FailedPrecondition(
        "Unrecognized type for field ", field_name_);
  }
  return tensorflow::Status::OK();
}

// TODO: consider eliminating this templated function by making
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
      values[i] = ListValues<T>((*found).second);
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
  } else if (vcf_type_ == BCF_HT_STR) {
    return DecodeValues<string>(header, bcf_record, variant);
  } else {
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
// VcfInfoFieldAdapter implementation.

VcfInfoFieldAdapter::VcfInfoFieldAdapter(const string& field_name,
                                         int vcf_type)
    : field_name_(field_name), vcf_type_(vcf_type) {}


tensorflow::Status VcfInfoFieldAdapter::EncodeValues(
    const nucleus::genomics::v1::Variant& variant, const bcf_hdr_t* header,
    bcf1_t* bcf_record) const {
  if (vcf_type_ == BCF_HT_REAL) {
    return EncodeValues<float>(variant, header, bcf_record);
  } else if (vcf_type_ == BCF_HT_INT) {
    return EncodeValues<int>(variant, header, bcf_record);
  } else if (vcf_type_ == BCF_HT_STR) {
    return EncodeValues<string>(variant, header, bcf_record);
  } else if (vcf_type_ == BCF_HT_FLAG) {
    return EncodeValues<bool>(variant, header, bcf_record);
  } else {
    return tensorflow::errors::FailedPrecondition(
        "Unrecognized type for field ", field_name_);
  }
  return tensorflow::Status::OK();
}

template <class T> tensorflow::Status VcfInfoFieldAdapter::EncodeValues(
    const nucleus::genomics::v1::Variant& variant,
    const bcf_hdr_t* header,
    bcf1_t* bcf_record) const {

  auto found = variant.info().find(field_name_);
  if (found != variant.info().end()) {
    std::vector<T> value = ListValues<T>((*found).second);
    return EncodeInfoValue(value, field_name_.c_str(), header, bcf_record);
  } else {
    return tensorflow::Status::OK();
  }
}


tensorflow::Status VcfInfoFieldAdapter::DecodeValues(
    const bcf_hdr_t* header, const bcf1_t* bcf_record,
    nucleus::genomics::v1::Variant* variant) const {
  if (vcf_type_ == BCF_HT_REAL) {
    return DecodeValues<float>(header, bcf_record, variant);
  } else if (vcf_type_ == BCF_HT_INT) {
    return DecodeValues<int>(header, bcf_record, variant);
  } else if (vcf_type_ == BCF_HT_STR) {
    return DecodeValues<string>(header, bcf_record, variant);
  } else if (vcf_type_ == BCF_HT_FLAG) {
    return DecodeValues<bool>(header, bcf_record, variant);
  } else {
    return tensorflow::errors::FailedPrecondition(
        "Unrecognized type for field ", field_name_);
  }
  return tensorflow::Status::OK();
}

template <class T> tensorflow::Status VcfInfoFieldAdapter::DecodeValues(
    const bcf_hdr_t *header, const bcf1_t *bcf_record,
    nucleus::genomics::v1::Variant *variant) const {
  std::vector<T> value =
      ReadInfoValue<T>(header, bcf_record, field_name_.c_str());
  SetInfoField(field_name_, value, variant);
  return tensorflow::Status::OK();
}



// -----------------------------------------------------------------------------
// VcfRecordConverter implementation.

VcfRecordConverter::VcfRecordConverter(
    const nucleus::genomics::v1::VcfHeader& vcf_header,
    const std::vector<string>& infos_to_exclude,
    const std::vector<string>& formats_to_exclude,
    const bool gl_and_pl_in_info_map) {
  // Install adapters for INFO fields.
  for (const auto& format_spec : vcf_header.infos()) {
    string tag = format_spec.id();
    string type = format_spec.type();

    // Skip fields that are handled specially.
    if (tag == "END") continue;

    // Check if configuration has disabled this INFO field.
    if (std::find(infos_to_exclude.begin(), infos_to_exclude.end(), tag) !=
        infos_to_exclude.end())
      continue;

    int vcf_type;
    if (type == "Integer") {
      vcf_type = BCF_HT_INT;
    } else if (type == "Float") {
      vcf_type = BCF_HT_REAL;
    } else if (type == "String" || type == "Character") {
      vcf_type = BCF_HT_STR;
    } else if (type == "Flag") {
      vcf_type = BCF_HT_FLAG;
    } else {
      LOG(WARNING) << "Unhandled INFO field type: field " << tag
                   << " of type " << type;
      continue;
    }
    info_adapters_.emplace_back(tag, vcf_type);
  }

  // Install adapters for FORMAT fields.
  want_pl_ = false;
  want_gl_ = false;
  gl_and_pl_in_info_map_ = gl_and_pl_in_info_map;
  for (const auto& format_spec : vcf_header.formats()) {
    string tag = format_spec.id();
    string type = format_spec.type();

    // Check if configuration has disabled this FORMAT field.
    if (std::find(formats_to_exclude.begin(), formats_to_exclude.end(), tag) !=
        formats_to_exclude.end())
      continue;

    // These fields are handled specially.
    if (tag == "GT") continue;

    if (tag == "GL") {
      want_gl_ = true;
      if (!gl_and_pl_in_info_map) continue;
    }
    if (tag == "PL") {
      want_pl_ = true;
      if (!gl_and_pl_in_info_map) continue;
    }

    // TODO: how do we really want to encode the type here?
    int vcf_type;
    if (type == "Integer") {
      vcf_type = BCF_HT_INT;
    } else if (type == "Float") {
      vcf_type = BCF_HT_REAL;
    } else if (type == "String" || type == "Character") {
      vcf_type = BCF_HT_STR;
    } else {
      LOG(WARNING) << "Unhandled FORMAT field type: field " << tag
                   << " of type " << type;
      continue;
    }
    format_adapters_.emplace_back(tag, vcf_type);
  }

  // Update special-cased variant fields.
  want_variant_end_ =
      std::find(infos_to_exclude.begin(), infos_to_exclude.end(), "END") ==
      infos_to_exclude.end();
  want_genotypes_ =
      std::find(formats_to_exclude.begin(), formats_to_exclude.end(), "GT") ==
      formats_to_exclude.end();
}

// static
void VcfHeaderConverter::ConvertToPb(const bcf_hdr_t* hdr,
                                     genomics::v1::VcfHeader* vcf_header) {
  vcf_header->Clear();
  if (hdr->nhrec < 1) {
    LOG(WARNING) << "Empty header, not a valid VCF.";
    return;
  }
  if (string(hdr->hrec[0]->key) != "fileformat") {
    LOG(WARNING) << "Not a valid VCF, fileformat needed.";
  }
  vcf_header->set_fileformat(hdr->hrec[0]->value);

  // Fill in the contig info for each contig in the VCF header. Directly
  // accesses the low-level C struct because there are no indirection
  // macros/functions by htslib API.
  // BCF_DT_CTG: offset for contig (CTG) information in BCF dictionary (DT).
  const int n_contigs = hdr->n[BCF_DT_CTG];
  for (int i = 0; i < n_contigs; ++i) {
    const bcf_idpair_t& idPair = hdr->id[BCF_DT_CTG][i];
    AddContigInfo(idPair, vcf_header->add_contigs(), i);
  }

  // Iterate through all hrecs (except the first, which was 'fileformat') to
  // populate the rest of the headers.
  for (int i = 1; i < hdr->nhrec; i++) {
    const bcf_hrec_t* hrec0 = hdr->hrec[i];
    switch (hrec0->type) {
      case BCF_HL_CTG:
        // Contigs are populated above, since they store length in the
        // bcf_idinfo_t* structure.
        break;
      case BCF_HL_FLT:
        AddFilterInfo(hrec0, vcf_header->add_filters());
        break;
      case BCF_HL_INFO:
        AddInfo(hrec0, vcf_header->add_infos());
        break;
      case BCF_HL_FMT:
        AddFormatInfo(hrec0, vcf_header->add_formats());
        break;
      case BCF_HL_STR:
        AddStructuredExtra(hrec0, vcf_header->add_structured_extras());
        break;
      case BCF_HL_GEN:
        AddExtra(hrec0, vcf_header->add_extras());
        break;
      default:
        LOG(WARNING) << "Unknown hrec0->type: " << hrec0->type;
    }
  }

  // Populate samples info.
  int n_samples = bcf_hdr_nsamples(hdr);
  for (int i = 0; i < n_samples; i++) {
    vcf_header->add_sample_names(hdr->samples[i]);
  }
}

tensorflow::Status VcfHeaderConverter::ConvertFromPb(
    const nucleus::genomics::v1::VcfHeader& vcf_header, bcf_hdr_t** h) {
  // Note: bcf_hdr_init writes the fileformat= and the FILTER=<ID=PASS,...>
  // filter automatically.
  *h = bcf_hdr_init("w");
  for (const nucleus::genomics::v1::VcfFilterInfo& filter :
       vcf_header.filters()) {
    if (filter.id() != "PASS") {
      AddFilterToHeader(filter, *h);
    }
  }
  for (const nucleus::genomics::v1::VcfInfo& info : vcf_header.infos()) {
    AddInfoToHeader(info, *h);
  }
  for (const nucleus::genomics::v1::VcfFormatInfo& format :
       vcf_header.formats()) {
    AddFormatToHeader(format, *h);
  }
  for (const nucleus::genomics::v1::VcfStructuredExtra& extra :
       vcf_header.structured_extras()) {
    AddStructuredExtraToHeader(extra, *h);
  }
  for (const nucleus::genomics::v1::VcfExtra& extra : vcf_header.extras()) {
    AddExtraToHeader(extra, *h);
  }
  for (const nucleus::genomics::v1::ContigInfo& contig : vcf_header.contigs()) {
    AddContigToHeader(contig, *h);
  }

  for (const string& sampleName : vcf_header.sample_names()) {
    bcf_hdr_add_sample(*h, sampleName.c_str());
  }
  bcf_hdr_add_sample(*h, nullptr);
  int ret = bcf_hdr_sync(*h);
  if (ret < 0) {
    return tensorflow::errors::DataLoss("Couldn't sync bcf header");
  }
  return tensorflow::Status::OK();
}

// Convert a C string to uppercase, in place, unless it starts with "<".
char *uppercase_allele(char *s) {
  char *r = s;
  if (*s == '<') return r;
  while (*s) {
    *s = toupper((unsigned char) *s);
    s++;
  }
  return r;
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
  // Don't add the missing "." marker to the id field.
  if (v->d.id && strcmp(v->d.id, ".") != 0) {
    std::vector<string> names = absl::StrSplit(v->d.id, ';');
    for (string &n : names) {
      variant_message->add_names(n);
    }
  }

  // Parse out the ref and alt alleles.
  if (v->n_allele > 0) {
    variant_message->set_reference_bases(uppercase_allele(v->d.allele[0]));
    for (int i = 1; i < v->n_allele; ++i) {
      variant_message->add_alternate_bases(uppercase_allele(v->d.allele[i]));
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

  // Parse the generic INFO fields.
  for (const auto& adapter : info_adapters_) {
    TF_RETURN_IF_ERROR(adapter.DecodeValues(h, v, variant_message));
  }

  // Parse the calls of the variant.
  if (v->n_sample > 0) {
    int* gt_arr = nullptr;
    int n_gts = 0;
    if (bcf_get_genotypes(h, v, &gt_arr, &n_gts) < 0) {
      free(gt_arr);
      return tensorflow::errors::DataLoss("Couldn't parse genotypes");
    }
    int max_ploidy = n_gts / v->n_sample;

    for (int i = 0; i < v->n_sample; i++) {
      nucleus::genomics::v1::VariantCall* call = variant_message->add_calls();
      call->set_call_set_name(h->samples[i]);
      // Get the GT calls, if requested and available.
      if (want_genotypes_) {
        bool gt_is_phased = false;
        for (int j = 0; j < max_ploidy; j++) {
          int gt_idx = gt_arr[i * max_ploidy + j];
          // Check whether this sample has smaller ploidy.
          if (gt_idx == bcf_int32_vector_end) break;

          int gt = bcf_gt_allele(gt_idx);
          gt_is_phased = gt_is_phased || bcf_gt_is_phased(gt_idx);
          call->add_genotype(gt);
        }
        call->set_is_phased(gt_is_phased);
      }
    }
    free(gt_arr);

    // Parse "generic" FORMAT fields.
    for (const auto& adapter : format_adapters_) {
      TF_RETURN_IF_ERROR(adapter.DecodeValues(h, v, variant_message));
    }

    // Handle FORMAT fields requiring special logic.
    if (!gl_and_pl_in_info_map_) {
      std::vector<std::vector<int>> pl_values =
          ReadFormatValues<int>(h, v, "PL");
      std::vector<std::vector<float>> gl_values =
          ReadFormatValues<float>(h, v, "GL");

      for (int i = 0; i < v->n_sample; i++) {
        // Each indicator here is true iff the format field is present for this
        // variant, *and* is non-missing for this sample.
        bool have_gl = !gl_values.empty() && !gl_values[i].empty();
        bool have_pl = !pl_values.empty() && !pl_values[i].empty();

        nucleus::genomics::v1::VariantCall* call =
            variant_message->mutable_calls(i);

        if (want_gl_ || want_pl_) {
          // If GL and PL are *both* present, we populate the
          // genotype_likelihood fields with the GL values per the
          // variants.proto spec, since PLs are a lower resolution version of
          // the same information.
          if (have_gl) {
            for (const auto& gl : gl_values[i]) {
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
  // See internal for discussion.  For now we use this workaround.

  // Workaround: if our ALT is just a single symbolic allele, this is
  // a gVCF record and we need to populate END, because htslib won't.
  if (want_variant_end_ && variant_message.alternate_bases_size() == 1 &&
      !variant_message.alternate_bases()[0].empty() &&
      variant_message.alternate_bases()[0][0] == '<') {
    int end = v->pos + v->rlen;
    if (bcf_update_info_int32(&h, v, "END", &end, 1) != 0)
      return tensorflow::errors::Unknown("Failure to write END to vcf record");
  }

  // Some variants don't have names; these will get the placeholder "." in the
  // ID column
  if (variant_message.names_size() > 0) {
    string combined_names = absl::StrJoin(variant_message.names(), ";");
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
    auto filterIds = absl::make_unique<int32[]>(nFilters);
    for (int i = 0; i < nFilters; i++) {
      const char* filterName = variant_message.filter(i).c_str();
      int32 filterId = bcf_hdr_id2int(&h, BCF_DT_ID, filterName);
      if (filterId < 0) {
        return tensorflow::errors::NotFound("Filter must be found in header.");
      }
      filterIds.get()[i] = filterId;
    }
    bcf_update_filter(&h, v, filterIds.get(), nFilters);
  }

  // Generic INFO fields
  for (const VcfInfoFieldAdapter& field : info_adapters_) {
    TF_RETURN_IF_ERROR(field.EncodeValues(variant_message, &h, v));
  }

  // Variant calls
  int nCalls = variant_message.calls().size();
  int nSamples = bcf_hdr_nsamples(&h);
  if (nCalls != nSamples)
    return tensorflow::errors::FailedPrecondition(
        "Variant call count ", nCalls, " must match number of samples ",
        nSamples, ".");

  // We need to determine the effective ploidy (as the max number of GT calls
  // among samples at this variant); any genotypes shorter than this ploidy
  // will be padded (consult the VCF spec).
  int ploidy = 0;
  for (int c = 0; c < nCalls; c++) {
    ploidy = std::max(ploidy, variant_message.calls(c).genotype_size());
  }

  if (nCalls > 0) {
    // Write genotypes.
    auto gts = absl::make_unique<int32[]>(nCalls * ploidy);
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

      const bool isPhased = vc.is_phased();
      int a = 0;
      for (; a < vc.genotype_size(); a++) {
        gts.get()[c * ploidy + a] = vcfEncodeAllele(vc.genotype(a), isPhased);
      }
      for (; a < ploidy; a++) {
        gts.get()[c * ploidy + a] = bcf_int32_vector_end;
      }
    }
    if (bcf_update_genotypes(&h, v, gts.get(), nCalls * ploidy) < 0) {
      return tensorflow::errors::Unknown(
          "Failure to write genotypes to VCF record");
    }

    // Write remaining FORMAT fields.
    // NOTE: This assumes that the Variant is well-formed with respect to GL and
    // PL; i.e. that either the map contains those fields or the
    // genotype_likelihood field is populated, but not both. This is guaranteed
    // when reading Variants from a VCF file.
    for (const VcfFormatFieldAdapter& field : format_adapters_) {
      TF_RETURN_IF_ERROR(field.EncodeValues(variant_message, &h, v));
    }

    if (!gl_and_pl_in_info_map_) {
      bool has_ll = false;
      for (int c = 0; c < nCalls; c++) {
        const nucleus::genomics::v1::VariantCall& vc = variant_message.calls(c);
        if (vc.genotype_likelihood_size() > 0) {
          has_ll = true;
          break;
        }
      }

      if (want_gl_ && has_ll) {
        std::vector<std::vector<float>> gl_values;
        for (int c = 0; c < nCalls; c++) {
          const nucleus::genomics::v1::VariantCall& vc =
              variant_message.calls(c);
          std::vector<float> gl_this_call;
          for (double ll_val : vc.genotype_likelihood()) {
            gl_this_call.push_back(ll_val);
          }
          // Do we need to zero-shift?
          gl_values.push_back(gl_this_call);
        }
        TF_RETURN_IF_ERROR(EncodeFormatValues(gl_values, "GL", &h, v));
      }

      if (want_pl_ && has_ll) {
        std::vector<std::vector<int>> ll_values_phred;
        for (int c = 0; c < nCalls; c++) {
          const nucleus::genomics::v1::VariantCall& vc =
              variant_message.calls(c);
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

        TF_RETURN_IF_ERROR(EncodeFormatValues(ll_values_phred, "PL", &h, v));
      }
    }
  }
  return tensorflow::Status::OK();
}

}  // namespace nucleus
