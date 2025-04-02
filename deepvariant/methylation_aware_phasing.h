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

#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Returns true if methylation values between the two haplotypes are
// significantly different (p < 0.05).
// Returns false if p >= threshold or inputs are invalid (e.g. one haplotype
// has no methylation data).
bool IsDifferentiallyMethylated(
    const std::vector<double>& hap1_methyl,
    const std::vector<double>& hap2_methyl);

// Performs the Wilcoxon Rank-Sum Test (Mann-Whitney U Test) between two
// haplotypes’ methylation values. Returns the two-sided p-value.
// Returns -1.0 if either haplotype has no data.
double WilcoxonRankSumTest(
    const std::vector<double>& hap1_methyl,
    const std::vector<double>& hap2_methyl);


// Performs iterative methylation-aware phasing on a set of reads.
//
// This function assigns haplotype phases to unphased reads based on their
// methylation profiles at informative CpG sites. It uses already-phased reads
// from direct phasing to define the haplotype-specific methylation signatures,
// identifies differentially methylated sites, and votes for each unphased
// read's haplotype assignment by comparing its methylation levels to those
// signatures.
std::vector<int> PerformMethylationAwarePhasing(
  const std::vector<nucleus::genomics::v1::Read>&
    reads_to_phase,
  const std::vector<int>& initial_read_phases,
  const std::vector<DeepVariantCall>& methylated_ref_sites,
  int max_iter);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_METHYLATION_AWARE_PHASING_H_
