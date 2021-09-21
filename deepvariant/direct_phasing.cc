/*
 * Copyright 2021 Google LLC.
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

#include "deepvariant/direct_phasing.h"

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

absl::StatusOr<std::vector<int>> DirectPhasing::PhaseReads(
    const std::vector<DeepVariantCall>& candidates,
    const std::vector<
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>& reads) {
  // Not Implemented.
  return absl::Status(absl::StatusCode::kUnimplemented, "Not Implemented");
}

// Helper functions.
AlleleType AlleleTypeFromCandidate(
    std::string_view bases,
    const DeepVariantCall& candidate) {
  if (bases.size() > candidate.variant().end() - candidate.variant().start()) {
    return AlleleType::INSERTION;
  }
  if (bases.size() < candidate.variant().end() - candidate.variant().start()) {
    return AlleleType::DELETION;
  }
  if (bases.size() == candidate.variant().end() - candidate.variant().start()) {
    return AlleleType::SUBSTITUTION;
  }
  return  AlleleType::UNSPECIFIED;
}

int NumOfSubstitutionAlleles(const DeepVariantCall& candidate) {
  return std::count_if(candidate.allele_support_ext().begin(),
         candidate.allele_support_ext().end(),
         [candidate](std::pair<std::string,
                                     DeepVariantCall_SupportingReadsExt> it) {
           return (it.first != kUncalledAllele &&
               AlleleTypeFromCandidate(it.first, candidate) ==
                   AlleleType::SUBSTITUTION);
         });
}

int NumOfIndelAlleles(const DeepVariantCall& candidate) {
  return std::count_if(candidate.allele_support_ext().begin(),
         candidate.allele_support_ext().end(),
         [candidate](std::pair<std::string,
                                     DeepVariantCall_SupportingReadsExt> it) {
           return (it.first != kUncalledAllele &&
               (AlleleTypeFromCandidate(it.first, candidate) ==
                   AlleleType::DELETION
                   || AlleleTypeFromCandidate(it.first, candidate) ==
                  AlleleType::INSERTION));
         });
}

int SubstitutionAllelesDepth(const DeepVariantCall& candidate) {
  int count = 0;
  for (const auto& allele_info_it : candidate.allele_support_ext()) {
    if (allele_info_it.first != kUncalledAllele
        && AlleleTypeFromCandidate(allele_info_it.first, candidate) ==
        AlleleType::SUBSTITUTION) {
      // redacted
      count += allele_info_it.second.read_infos_size();
    }
  }
  return count;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
