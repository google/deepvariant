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

#ifndef THIRD_PARTY_NUCLEUS_UTIL_VCF_CONVERSION_H_
#define THIRD_PARTY_NUCLEUS_UTIL_VCF_CONVERSION_H_

#include "htslib/vcf.h"
#include "deepvariant/util/genomics/variants.pb.h"
#include "deepvariant/util/genomics/vcf.pb.h"
#include "deepvariant/util/vendor/statusor.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

// -----------------------------------------------------------------------------
// VCF type encoding utilities
template<class T>
struct VcfType {
  // Predicates for checking missing and sentinel entries.  Use these, not ==.
  // Is argument the "missing" value?
  static bool IsMissing(T);
  // Is argument the vector end sentinel value?
  static bool IsVectorEnd(T);
  // Missing and sentinel value assignment.  Use these, not = assignment
  // Set argument to "missing" value
  static void SetMissing(T* v);
  // Set argument to vector end sentinel
  static void SetVectorEnd(T* v);

  // Format field extraction: wrapper for bcf_get_format_*()
  // Given a VCF record, this will grab format tag "tag" from the record, a
  // buffer *dst for the contents, and copy them there.  The return value will
  // be either negative (for error) or the number of values copied.  *ndst
  // should not be used.
  static int GetFormatValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                             const char *tag, T **dst, int *ndst);

  // Format field writing: wrapper for bcf_update_format_*()
  // Given a VCF record, this routine will populate the format field specified
  // by "tag" with the values src[0]...src[nsrc-1].
  static tensorflow::Status PutFormatValues(const char *tag, const T *src,
                                            int nsrc, const bcf_hdr_t *hdr,
                                            bcf1_t *line);
};

// See interface description comment above.
template<>
struct VcfType<int> {
  static bool IsMissing(int v)   { return (v == bcf_int32_missing); }
  static bool IsVectorEnd(int v) { return (v == bcf_int32_vector_end); }
  static void SetMissing(int* v)     { *v = bcf_int32_missing;  }
  static void SetVectorEnd(int* v)   { *v = bcf_int32_vector_end; }

  static int GetFormatValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                             const char *tag, int **dst, int *ndst) {
    return bcf_get_format_int32(hdr, const_cast<bcf1_t *>(line),
                                tag, dst, ndst);
  }

  static tensorflow::Status PutFormatValues(const char *tag, const int *src,
                                            int nsrc, const bcf_hdr_t *hdr,
                                            bcf1_t *line) {
    if (bcf_update_format_int32(hdr, line, tag, src, nsrc) != 0)
      return tensorflow::errors::Internal("bcf_update_format_int32 failed");
    else
      return tensorflow::Status::OK();
  }
};

// See interface description comment above.
template<>
struct VcfType<float> {
  static bool IsMissing(float v)   { return bcf_float_is_missing(v); }
  static bool IsVectorEnd(float v) { return bcf_float_is_vector_end(v); }
  static void SetMissing(float* v) { bcf_float_set(v, bcf_float_missing); }
  static void SetVectorEnd(float* v) { bcf_float_set(v, bcf_float_vector_end); }

  static int GetFormatValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                             const char *tag, float **dst, int *ndst) {
    return bcf_get_format_float(hdr, const_cast<bcf1_t *>(line),
                                tag, dst, ndst);
  }

  static tensorflow::Status PutFormatValues(const char *tag, const float *src,
                                            int nsrc, const bcf_hdr_t *hdr,
                                            bcf1_t *line) {
    if (bcf_update_format_float(hdr, line, tag, src, nsrc) != 0)
      return tensorflow::errors::Internal("bcf_update_format_float failed");
    else
      return tensorflow::Status::OK();
  }
};


// -----------------------------------------------------------------------------
// Helper class for encoding VariantCall.info values in VCF FORMAT field values.
//
// The standard way to interact with this class is:
//
// Create an adaptor for an .info map key "DP".
// FormatFieldAdapter adapter("DP");
//
// For each variant, we check if we need to encode values.
// if (adapter.IsPresentInAnyVariantCalls(variant)) {
//   // And if so, encode them into our htslib bcf_record passing in the
//   // required htslib header object as well.
//   adapter.EncodeValues(variant, header, bcf_record);
// }
class VcfFormatFieldAdapter {
 public:
  // Creates a new adapter for a field name field_name.
  explicit VcfFormatFieldAdapter(const string& field_name, int vcf_type);

  // Returns true if the format field field_name occurs in any VariantCall.info
  // maps present in Variant.
  bool IsPresentInAnyVariantCalls(
      const nucleus::genomics::v1::Variant& variant) const;

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
  // supported by this class is the repeated genotype_likelihood field in the
  // VariantCall proto. These aren't stored in the info field, but are inlined
  // directly in the proto, even though they are written to the FORMAT field
  // in VCF.
  //
  // WARNING: This code currently assumes that all field values are real
  // numbers for "VAF" field. For all other field, it assumes the field values
  // are integers.
  tensorflow::Status EncodeValues(const nucleus::genomics::v1::Variant& variant,
                                  const bcf_hdr_t* header,
                                  bcf1_t* bcf_record) const;

 private:  // Non-API methods
  template <class T>
  tensorflow::Status EncodeValues(const nucleus::genomics::v1::Variant& variant,
                                  const bcf_hdr_t* header,
                                  bcf1_t* bcf_record) const;

 private:  // Fields
  // The name of our field, such as "DP", "AD", or "VAF".
  string field_name_;
  // The htslib/VCF "type" of this field, such as BCF_HT_INT.
  int vcf_type_;
};


// Helper class for converting between Variant proto messages and VCF records.
class VcfRecordConverter {
 public:
  // Constructor.
  VcfRecordConverter(
      const OptionalVariantFieldsToParse &desired_format_entries);

  // Convert a VCF line parsed by htslib into a Variant protocol buffer.
  // The parsed line is passed in v, and the parsed header is in h.
  tensorflow::Status ConvertToPb(
      const bcf_hdr_t *h, bcf1_t *v,
      nucleus::genomics::v1::Variant *variant_message) const;

  // Convert a Variant protocol buffer into htslib's representation of a VCF
  // line.
  tensorflow::Status ConvertFromPb(
      const nucleus::genomics::v1::Variant &variant_message, const bcf_hdr_t &h,
      bcf1_t *v) const;

 private:
  // Lookup table for genotype field adapters by VCF tag name.
  // The order of adapter definitions here determines the order of the fields
  // in a written VCF.
  const std::vector<VcfFormatFieldAdapter> format_adapters_;
  // Configuration of the FORMAT entries should we parse in from, or write out
  // to, VCF.
  const OptionalVariantFieldsToParse desired_format_entries_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_UTIL_VCF_CONVERSION_H_
