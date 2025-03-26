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

#include "deepvariant/methylation_aware_phasing.h"

#include <vector>

#include "tensorflow/core/platform/test.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Test case: Distinct distributions should be differentially methylated
TEST(MethylationAwarePhasingTest,
     DistinctDistributionsIsDifferentiallyMethylated) {
  std::vector<double> hap1_methyl = {0.10, 0.15, 0.20, 0.12, 0.18};
  std::vector<double> hap2_methyl = {0.75, 0.80, 0.85, 0.78, 0.82};

  EXPECT_TRUE(IsDifferentiallyMethylated(hap1_methyl, hap2_methyl));
}

// Test case: Identical distributions should not be differentially methylated
TEST(MethylationAwarePhasingTest,
     IdenticalDistributionsIsDifferentiallyMethylated) {
  std::vector<double> hap1_methyl = {0.35, 0.40, 0.45, 0.50, 0.42};
  std::vector<double> hap2_methyl = {0.35, 0.40, 0.45, 0.50, 0.42};

  EXPECT_FALSE(IsDifferentiallyMethylated(hap1_methyl, hap2_methyl));
}

// Test case: Empty haplotypes should not be differentially methylated
TEST(MethylationAwarePhasingTest,
     EmptyHaplotypesIsNotDifferentiallyMethylated) {
  std::vector<double> hap1_methyl = {};
  std::vector<double> hap2_methyl = {};

  // Case 1: Both haplotypes are empty
  EXPECT_FALSE(IsDifferentiallyMethylated(
      hap1_methyl, hap2_methyl));

  // Case 2: One haplotype is empty, one is non-empty
  hap2_methyl = {0.2, 0.4, 0.6};

  EXPECT_FALSE(IsDifferentiallyMethylated(
    hap1_methyl, hap2_methyl));
}

TEST(MethylationAwarePhasingTest, WilcoxonRankSumSortOrderMatters) {
  // hap1 has clearly higher values than hap2, so we expect a very small p-value
  std::vector<double> hap1_methyl = {0.9, 0.85, 0.88, 0.95, 0.92};
  std::vector<double> hap2_methyl = {0.1, 0.12, 0.15, 0.05, 0.09};

  double p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);

  // The difference is large, so p should be near 0
  EXPECT_GT(0.01, p_val);
}

TEST(MethylationAwarePhasingTest, WilcoxonRankSumGroupAssignmentMatters) {
  // Hap1 has low methylation, hap2 has high methylation
  std::vector<double> hap1_methyl = {0.1, 0.2, 0.3};
  std::vector<double> hap2_methyl = {0.8, 0.9, 1.0};

  double p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);

  // Ensure test fails if group mapping is reversed (should be significant)
  EXPECT_LT(p_val, 0.05);
}



}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
