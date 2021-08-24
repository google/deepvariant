/*
 * Copyright 2017 Google LLC.
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

#include "deepvariant/utils.h"

#include "absl/strings/str_cat.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

Allele MakeAllele(const string& bases, const AlleleType type, const int count,
                  const bool is_low_quality) {
  Allele allele;
  allele.set_bases(bases);
  allele.set_type(type);
  allele.set_count(count);
  allele.set_is_low_quality(is_low_quality);
  return allele;
}

string SimplifyRefAlt(const string& ref, const string& alt) {
  int shortest_allele_len = ref.length();
  if (alt.length() < shortest_allele_len) {
    shortest_allele_len = alt.length();
  }
  int common_suffix_len = 0;
  for (int suffix_idx = 1; suffix_idx < shortest_allele_len; ++suffix_idx) {
    if (ref.at(ref.length() - suffix_idx) !=
        alt.at(alt.length() - suffix_idx)) {
      break;
    }
    common_suffix_len = suffix_idx;
  }

  if (common_suffix_len == 0) {
    return absl::StrCat(ref, "->", alt);
  } else {
    return absl::StrCat(ref.substr(0, ref.length() - common_suffix_len), "->",
                        alt.substr(0, alt.length() - common_suffix_len));
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
