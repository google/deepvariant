/*
 * Copyright 2025 Google LLC.
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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_METHYLATION_AWARE_PHASING_H_
#define LEARNING_GENOMICS_DEEPVARIANT_METHYLATION_AWARE_PHASING_H_

#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Performs the Wilcoxon Rank-Sum Test (Mann-Whitney U Test) between two
// haplotypes’ methylation values. Returns the two-sided p-value.
// Returns -1.0 if either haplotype has no data.
double WilcoxonRankSumTest(absl::Span<const double> hap1_methyl,
                           absl::Span<const double> hap2_methyl);

// Votes on the haplotype assignment for an unphased read using methylation
// similarity.
//
// For each informative site, compares the read's methylation level to the
// average levels for hap1 and hap2, and tallies votes. The read is assigned
// to the haplotype with the majority of votes (minimum 3 votes required).
//
// Returns:
//   1 = Haplotype 1
//   2 = Haplotype 2
//   0 = Cannot assign
int HaplotypeVoteWithMethylation(
    const DeepVariantCall_ReadSupport& unphased_read,
    absl::Span<const DeepVariantCall> informative_calls,
    const std::vector<const DeepVariantCall_ReadSupport*>& hap1_reads,
    const std::vector<const DeepVariantCall_ReadSupport*>& hap2_reads);

// Returns methylation level as a float between 0 and 1, or -1 if no
// methylation level is present (methylation_level == 0).
double GetMethylationLevelAtSite(const DeepVariantCall_ReadSupport& read);

// Identifies methylated reference sites that show differential methylation
// between haplotype 1 and haplotype 2 read sets.
std::vector<DeepVariantCall> IdentifyInformativeSitesAndUpdatePValues(
    const std::vector<const DeepVariantCall_ReadSupport*>& hap1_reads,
    const std::vector<const DeepVariantCall_ReadSupport*>& hap2_reads,
    std::vector<DeepVariantCall>& methylated_calls);

// Returns a unique key identifying a read based on its fragment name and
// read number.
std::string ReadKeyForMethylationAwarePhasing(
    const nucleus::genomics::v1::Read& read);

// Extracts reads assigned to a specific haplotype phase (0, 1, or 2).
std::vector<const DeepVariantCall_ReadSupport*> ExtractReadsByPhase(
    absl::Span<const DeepVariantCall_ReadSupport> reads,
    absl::Span<const int> phases, int target_phase);

// Performs iterative methylation-aware phasing on a set of reads.
//
// This function assigns haplotype phases to unphased reads based on their
// methylation profiles at informative CpG sites. It uses already-phased reads
// from direct phasing to define the haplotype-specific methylation signatures,
// identifies differentially methylated sites, and votes for each unphased
// read's haplotype assignment by comparing its methylation levels to those
// signatures.
std::tuple<std::vector<int>, std::vector<double>>
PerformMethylationAwarePhasing(
    absl::Span<const nucleus::genomics::v1::Read> reads_to_phase,
    const std::vector<int>& initial_read_phases,
    std::vector<DeepVariantCall>& methylated_ref_sites, int max_iter);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_METHYLATION_AWARE_PHASING_H_
