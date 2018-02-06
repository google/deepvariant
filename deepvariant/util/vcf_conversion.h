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

#ifndef THIRD_PARTY_NUCLEUS_UTIL_VCF_CONVERSION_H_
#define THIRD_PARTY_NUCLEUS_UTIL_VCF_CONVERSION_H_

#include "htslib/vcf.h"
#include "deepvariant/util/vendor/statusor.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

// redacted
// independent testing of those functions.

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
      return tensorflow::errors::Internal("bcf_update_format_float failed");
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

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_UTIL_VCF_CONVERSION_H_
