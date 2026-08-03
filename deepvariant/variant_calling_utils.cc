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

#include <algorithm>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/utils.h"
#include "absl/container/btree_map.h"
#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"

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

std::string CalcRefBases(absl::string_view ref_bases,
                         absl::Span<const Allele> alt_alleles,
                         absl::Span<const Allele> rejected_alleles) {
  if (alt_alleles.empty()) {
    // We don't have any alternate alleles, so used the provided ref_bases.
    return std::string(ref_bases);
  }

  const auto max_element_main =
      std::max_element(alt_alleles.cbegin(), alt_alleles.cend(),
                       [](const Allele& allele1, const Allele& allele2) {
                         return DeletionSize(allele1) < DeletionSize(allele2);
                       });
  const auto max_element_rejected =
      std::max_element(rejected_alleles.cbegin(), rejected_alleles.cend(),
                       [](const Allele& allele1, const Allele& allele2) {
                         return DeletionSize(allele1) < DeletionSize(allele2);
                       });
  // rejected_alleles may be empty, so we need to check for that.
  auto max_element =
      (max_element_rejected == rejected_alleles.cend() ||
       DeletionSize(*max_element_main) > DeletionSize(*max_element_rejected))
          ? max_element_main
          : max_element_rejected;

  if (max_element->type() != AlleleType::DELETION) {
    return std::string(ref_bases);
  } else {
    // Deletion alleles may have an anchor base that is the reference or some
    // other base, but a Variant must have a reference sequence that starts with
    // the reference base. The index 1 skips the first base of the deletion,
    // which is the anchor base of the deletion.
    CHECK(max_element->bases().size() > 1)
        << "Saw invalid deletion allele with too few bases"
        << max_element->ShortDebugString();
    return absl::StrCat(ref_bases, max_element->bases().substr(1));
  }
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

void AddReadDepths(const AlleleCount& allele_count, const AlleleMap& allele_map,
                   absl::string_view allele_map_refbases, Variant* variant) {
  // Set the DP to the total good reads seen at this position.
  VariantCall* call = variant->mutable_calls(0);
  nucleus::SetInfoField(kDPFormatField, TotalAlleleCounts(allele_count), call);

  if (variant->alternate_bases_size() == 1 &&
      (variant->alternate_bases(0) == kNoAltAllele ||
       variant->alternate_bases(0) == kGVCFAltAllele)) {
    // Variant has no alts or is a a gVCF record so only DP is meaningful.
    return;
  } else {
    int dp = TotalAlleleCounts(allele_count);
    // Build up AD and VAF.
    std::vector<int> ad;
    std::vector<double> vaf;
    ad.reserve(variant->alternate_bases_size() + 1);
    vaf.reserve(variant->alternate_bases_size());
    ad.push_back(allele_count.ref_supporting_read_count());

    absl::btree_map<std::string, const Allele*, std::less<>> alt_to_alleles;
    for (const auto& [allele, alt_bases] : allele_map) {
      const std::string key =
          SimplifyRefAlt(allele_map_refbases, alt_bases);
      alt_to_alleles[key] = &allele;
    }
    CHECK(alt_to_alleles.size() == allele_map.size())
        << "Non-unique alternative alleles!";
    for (const std::string& alt : variant->alternate_bases()) {
      const std::string simplified_ref_alt =
          SimplifyRefAlt(variant->reference_bases(), alt);
      int count_of_allele = 0;
      auto found = alt_to_alleles.find(simplified_ref_alt);
      if (found != alt_to_alleles.end()) {
        count_of_allele = (*found->second).count();
      }
      double this_vaf = 0.0;
      if (dp > 0) {
        this_vaf = 1.0 * count_of_allele / dp;
      }
      ad.push_back(count_of_allele);
      vaf.push_back(this_vaf);
    }

    nucleus::SetInfoField(kADFormatField, ad, call);
    nucleus::SetInfoField(kVAFFormatField, vaf, call);
  }
}

}  // namespace variant_calling_utils
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
