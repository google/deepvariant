/*
 * Copyright 2026 Google LLC.
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

#include "deepvariant/variant_calling_utils.h"

#include <cstdint>
#include <string>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace variant_calling_utils {

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;

const char* const kGVCFAltAllele = "<*>";
const char* const kSupportingUncalledAllele = "UNCALLED_ALLELE";
const char* const kDPFormatField = "DP";
const char* const kADFormatField = "AD";
const char* const kVAFFormatField = "VAF";
const char* const kMFFormatField = "MF";
const char* const kMDFormatField = "MD";

const char* const kNoAltAllele = ".";

int DeletionSize(const Allele& allele) {
  return allele.type() == AlleleType::DELETION ? allele.bases().length() : -1;
}

std::string MakeAltAllele(absl::string_view prefix,
                          absl::string_view variant_ref, const uint32_t from) {
  const auto postfix =
      from >= variant_ref.length() ? "" : variant_ref.substr(from);
  return absl::StrCat(prefix, postfix);
}

void AddGenotypes(const std::string& sample_name,
                  absl::Span<const int> genotypes, Variant* variant) {
  CHECK(variant != nullptr);

  VariantCall* call = variant->add_calls();
  call->set_call_set_name(sample_name);
  for (const auto genotype : genotypes) {
    call->add_genotype(genotype);
  }
}

}  // namespace variant_calling_utils
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
