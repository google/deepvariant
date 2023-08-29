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

#ifndef THIRD_PARTY_NUCLEUS_IO_VCF_CONVERSION_H_
#define THIRD_PARTY_NUCLEUS_IO_VCF_CONVERSION_H_

#include <string>
#include <vector>

#include "htslib/vcf.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

// -----------------------------------------------------------------------------
// VCF type encoding utilities
template <class T>
struct VcfType {
  // Predicates for checking missing and sentinel entries.  Use these, not ==.
  // Is argument the "missing" value?
  static bool IsMissing(T);
  // Is argument the vector end sentinel value?
  static bool IsVectorEnd(T);
  // Missing and sentinel value assignment.  Use these, not = assignment
  // Set argument to "missing" value
  static void SetMissing(T *v);
  // Set argument to vector end sentinel
  static void SetVectorEnd(T *v);

  // FORMAT field extraction: wrapper for bcf_get_format_*()
  // Given a VCF record, this will grab FORMAT tag "tag" from the record, a
  // buffer *dst for the contents, and copy them there.  The return value will
  // be either negative (for error) or the number of values copied.  *ndst
  // should not be used.
  static int GetFormatValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                             const char *tag, T **dst, int *ndst);

  // FORMAT field writing: wrapper for bcf_update_format_*()
  // Given a VCF record, this routine will populate the FORMAT field specified
  // by "tag" with the values src[0]...src[nsrc-1].
  static ::nucleus::Status PutFormatValues(const char *tag, const T *src,
                                           int nsrc, const bcf_hdr_t *hdr,
                                           bcf1_t *line);

  // INFO field extraction: wrapper for bcf_get_info_*()
  // Given a VCF record, this will grab INFO tag "tag" from the record, a
  // buffer *dst for the contents, and copy them there.  The return value will
  // be either negative (for error) or the number of values copied.  *ndst
  // should not be used.
  static int GetInfoValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                           const char *tag, T **dst, int *ndst);

  // INFO field writing: wrapper for bcf_update_info_*()
  // Given a VCF record, this routine will populate the INFO field specified
  // by "tag" with the values src[0]...src[nsrc-1].
  static ::nucleus::Status PutInfoValues(const char *tag, const T *src,
                                         int nsrc, const bcf_hdr_t *hdr,
                                         bcf1_t *line);
};

// See interface description comment above.
template <>
struct VcfType<int> {
  static bool IsMissing(int v) { return (v == bcf_int32_missing); }
  static bool IsVectorEnd(int v) { return (v == bcf_int32_vector_end); }
  static void SetMissing(int *v) { *v = bcf_int32_missing; }
  static void SetVectorEnd(int *v) { *v = bcf_int32_vector_end; }

  static int GetFormatValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                             const char *tag, int **dst, int *ndst) {
    return bcf_get_format_int32(hdr, const_cast<bcf1_t *>(line), tag, dst,
                                ndst);
  }

  static ::nucleus::Status PutFormatValues(const char *tag, const int *src,
                                           int nsrc, const bcf_hdr_t *hdr,
                                           bcf1_t *line) {
    if (bcf_update_format_int32(hdr, line, tag, src, nsrc) != 0)
      return ::nucleus::Internal("bcf_update_format_int32 failed");
    else
      return ::nucleus::Status();
  }

  static int GetInfoValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                           const char *tag, int **dst, int *ndst) {
    return bcf_get_info_int32(hdr, const_cast<bcf1_t *>(line), tag, dst, ndst);
  }

  static ::nucleus::Status PutInfoValues(const char *tag, const int *src,
                                         int nsrc, const bcf_hdr_t *hdr,
                                         bcf1_t *line) {
    if (bcf_update_info_int32(hdr, line, tag, src, nsrc) != 0)
      return ::nucleus::Internal("bcf_update_info_int32 failed");
    else
      return ::nucleus::Status();
  }
};

// See interface description comment above.
template <>
struct VcfType<float> {
  static bool IsMissing(float v) { return bcf_float_is_missing(v); }
  static bool IsVectorEnd(float v) { return bcf_float_is_vector_end(v); }
  static void SetMissing(float *v) { bcf_float_set(v, bcf_float_missing); }
  static void SetVectorEnd(float *v) { bcf_float_set(v, bcf_float_vector_end); }

  static int GetFormatValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                             const char *tag, float **dst, int *ndst) {
    return bcf_get_format_float(hdr, const_cast<bcf1_t *>(line), tag, dst,
                                ndst);
  }

  static ::nucleus::Status PutFormatValues(const char *tag, const float *src,
                                           int nsrc, const bcf_hdr_t *hdr,
                                           bcf1_t *line) {
    if (bcf_update_format_float(hdr, line, tag, src, nsrc) != 0)
      return ::nucleus::Internal("bcf_update_format_float failed");
    else
      return ::nucleus::Status();
  }

  static int GetInfoValues(const bcf_hdr_t *hdr, const bcf1_t *line,
                           const char *tag, float **dst, int *ndst) {
    return bcf_get_info_float(hdr, const_cast<bcf1_t *>(line), tag, dst, ndst);
  }

  static ::nucleus::Status PutInfoValues(const char *tag, const float *src,
                                         int nsrc, const bcf_hdr_t *hdr,
                                         bcf1_t *line) {
    if (bcf_update_info_float(hdr, line, tag, src, nsrc) != 0)
      return ::nucleus::Internal("bcf_update_info_float failed");
    else
      return ::nucleus::Status();
  }
};

// -----------------------------------------------------------------------------
// Helper class for encoding VariantCall.info values in VCF FORMAT field values.
// This class is only intended for use with FORMAT fields that can be directly
// mapped between a VCF record and the FORMAT info dictionary, without special
// logic.  Where special logic is needed (e.g. for GT, GL/PL, etc.), the lower
// level functions `ReadFormatValues` and `EncodeFormatValues` are called
// directly.
//
// The standard way to interact with this class is as follows.
//
// Create an adaptor for FORMAT field "DP" of integer type:
//   VcfFormatFieldAdapter adapter("DP", BCF_HT_INT32);
//
// For each variant, we encode this format field into the vcf record:
//   adapter.EncodeValues(variant, header, bcf_record);
//
class VcfFormatFieldAdapter {
 public:
  // Creates a new adapter for a field name field_name.
  VcfFormatFieldAdapter(const string &field_name, int vcf_type);

  // Adds the values for our field_name from variant's calls into our bcf1_t
  // record bcf_record.
  ::nucleus::Status EncodeValues(const nucleus::genomics::v1::Variant &variant,
                                 const bcf_hdr_t *header,
                                 bcf1_t *bcf_record) const;

  // Add the values for this genotype field in the bcf1_t `bcf_record` to the
  // VariantCall info maps within this Variant proto message `variant`.
  ::nucleus::Status DecodeValues(const bcf_hdr_t *header,
                                 const bcf1_t *bcf_record,
                                 nucleus::genomics::v1::Variant *variant) const;

 private:  // Non-API methods
  template <class T>
  ::nucleus::Status EncodeValues(const nucleus::genomics::v1::Variant &variant,
                                 const bcf_hdr_t *header,
                                 bcf1_t *bcf_record) const;

  template <class T>
  ::nucleus::Status DecodeValues(const bcf_hdr_t *header,
                                 const bcf1_t *bcf_record,
                                 nucleus::genomics::v1::Variant *variant) const;

 private:  // Fields
  // The name of our field, such as "DP", "AD", or "VAF".
  string field_name_;
  // The htslib/VCF "type" of this field, such as BCF_HT_INT.
  int vcf_type_;
};

// -----------------------------------------------------------------------------
// Helper class for encoding Variant.info values in VCF INFO field values.
//
// (Usage of this class is completely analogous to the VcfFormatFieldAdapter
// class.)
class VcfInfoFieldAdapter {
 public:
  // Creates a new adapter for a field name field_name.
  VcfInfoFieldAdapter(const string &field_name, int vcf_type);

  // Adds the values for our field_name from the Variant into our bcf1_t
  // record bcf_record.
  ::nucleus::Status EncodeValues(const nucleus::genomics::v1::Variant &variant,
                                 const bcf_hdr_t *header,
                                 bcf1_t *bcf_record) const;

  // Add the values for this INFO field in the bcf1_t `bcf_record` to the
  // Variant message info map.
  ::nucleus::Status DecodeValues(const bcf_hdr_t *header,
                                 const bcf1_t *bcf_record,
                                 nucleus::genomics::v1::Variant *variant) const;

 private:  // Non-API methods
  template <class T>
  ::nucleus::Status EncodeValues(const nucleus::genomics::v1::Variant &variant,
                                 const bcf_hdr_t *header,
                                 bcf1_t *bcf_record) const;

  template <class T>
  ::nucleus::Status DecodeValues(const bcf_hdr_t *header,
                                 const bcf1_t *bcf_record,
                                 nucleus::genomics::v1::Variant *variant) const;

 private:  // Fields
  // The name of our info field, such as "H2" or "END"
  string field_name_;
  // The htslib/VCF "type" of this field, such as BCF_HT_INT.
  int vcf_type_;
};

// Helper class for converting between VcfHeader proto messages and bcf_hdr_t
// structs.
class VcfHeaderConverter {
 public:
  static void ConvertToPb(const bcf_hdr_t *h,
                          nucleus::genomics::v1::VcfHeader *vcf_header);

  // Converts a proto VcfHeader to a bcf_hdr_t. Caller needs to take ownership
  // of |h|.
  static ::nucleus::Status ConvertFromPb(
      const nucleus::genomics::v1::VcfHeader &vcf_header, bcf_hdr_t **h);
};

// Helper class for converting between Variant proto messages and VCF records.
class VcfRecordConverter {
 public:
  // Primary constructor.
  VcfRecordConverter(const nucleus::genomics::v1::VcfHeader &vcf_header,
                     const std::vector<string> &infos_to_exclude,
                     const std::vector<string> &formats_to_exclude,
                     const bool gl_and_pl_in_info_map);

  // Not the constructor you want.
  VcfRecordConverter() = default;

  // Convert a VCF line parsed by htslib into a Variant protocol buffer.
  // The parsed line is passed in v, and the parsed header is in h.
  ::nucleus::Status ConvertToPb(
      const bcf_hdr_t *h, bcf1_t *v,
      nucleus::genomics::v1::Variant *variant_message) const;

  // Convert a Variant protocol buffer into htslib's representation of a VCF
  // line.
  ::nucleus::Status ConvertFromPb(
      const nucleus::genomics::v1::Variant &variant_message, const bcf_hdr_t &h,
      bcf1_t *v) const;

 private:
  // Lookup table for variant INFO fields adapters by VCF tag name.
  // The order of adapter definitions here determines the order of the fields
  // in a written VCF.
  std::vector<VcfInfoFieldAdapter> info_adapters_;
  // Lookup table for genotype FORMAT field adapters by VCF tag name.
  // The order of adapter definitions here determines the order of the fields
  // in a written VCF.
  std::vector<VcfFormatFieldAdapter> format_adapters_;

  // Individual special-cased INFO fields.
  bool want_variant_end_;
  // Individual special-cased FORMAT fields.
  bool want_genotypes_;
  bool want_gl_;
  bool want_pl_;

  // Set to true if the GL and PL fields should be stored to and retrieved from
  // the info map with other FORMAT fields, rather than being special-cased as
  // first-class members of the proto.
  bool gl_and_pl_in_info_map_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_IO_VCF_CONVERSION_H_
