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

#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Helper function for creating a DeepVariant Call with the given parameters.
DeepVariantCall MakeMethylCall(
    int pos,
    int n_hap1_reads,
    int n_hap2_reads,
    int hap1_methyl,
    int hap2_methyl,
    std::vector<const DeepVariantCall_ReadSupport*>& hap1_reads,
    std::vector<const DeepVariantCall_ReadSupport*>& hap2_reads) {
  DeepVariantCall call;
  call.mutable_variant()->set_start(pos);
  auto* ext = call.mutable_ref_support_ext();
  for (int i = 0; i < n_hap1_reads; ++i) {
    auto* r = ext->add_read_infos();
    r->set_read_name("hap1_" + std::to_string(i));
    r->set_methylation_level(hap1_methyl);
    hap1_reads.push_back(r);
  }
  for (int i = 0; i < n_hap2_reads; ++i) {
    auto* r = ext->add_read_infos();
    r->set_read_name("hap2_" + std::to_string(i));
    r->set_methylation_level(hap2_methyl);
    hap2_reads.push_back(r);
  }
  return call;
}

// Test case: Distinct distributions should be differentially methylated
TEST(MethylationAwarePhasingTest,
     DistinctDistributionsIsDifferentiallyMethylated) {
  std::vector<double> hap1_methyl = {0.10, 0.15, 0.20, 0.12, 0.18};
  std::vector<double> hap2_methyl = {0.75, 0.80, 0.85, 0.78, 0.82};

  double p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);
  EXPECT_LT(p_val, 0.05);
}

// Test case: Identical distributions should not be differentially methylated
TEST(MethylationAwarePhasingTest,
     IdenticalDistributionsIsDifferentiallyMethylated) {
  std::vector<double> hap1_methyl = {0.35, 0.40, 0.45, 0.50, 0.42};
  std::vector<double> hap2_methyl = {0.35, 0.40, 0.45, 0.50, 0.42};
  double p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);
  EXPECT_GT(p_val, 0.05);
}
// Test case: Empty haplotypes should not be differentially methylated
TEST(MethylationAwarePhasingTest,
     EmptyHaplotypesIsNotDifferentiallyMethylated) {
  std::vector<double> hap1_methyl = {};
  std::vector<double> hap2_methyl = {};
  double p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);
  // Case 1: Both haplotypes are empty; p-value should be -1.0
  EXPECT_EQ(p_val, -1.0);
  // Case 2: One haplotype is empty, one is non-empty; p-value should be -1.0
  hap2_methyl = {0.2, 0.4, 0.6};
  p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);
  EXPECT_EQ(p_val, -1.0);
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

TEST(MethylationAwarePhasingTest, GetMethylationLevelAtSiteReturnsNormalized) {
  DeepVariantCall_ReadSupport support;
  support.set_methylation_level(128);
  EXPECT_NEAR(GetMethylationLevelAtSite(support), 128.0 / 255.0, 1e-6);
}

TEST(MethylationAwarePhasingTest, GetMethylationLevelAtSiteReturnsMinusOne) {
  DeepVariantCall_ReadSupport support;
  support.set_methylation_level(0);
  EXPECT_EQ(GetMethylationLevelAtSite(support), -1.0);
}

TEST(MethylationAwarePhasingTest, ReadKeyMatchesExpectedFormat) {
  nucleus::genomics::v1::Read read;
  read.set_fragment_name("frag");
  read.set_read_number(1);
  EXPECT_EQ(ReadKeyForMethylationAwarePhasing(read), "frag/1");
}

TEST(MethylationAwarePhasingTest, ExtractReadsByPhaseReturnsCorrectSubset) {
  std::vector<DeepVariantCall_ReadSupport> reads(3);
  reads[0].set_read_name("a");
  reads[1].set_read_name("b");
  reads[2].set_read_name("c");
  std::vector<int> phases = {0, 2, 2};

  auto result = ExtractReadsByPhase(reads, phases, 2);
  ASSERT_EQ(result.size(), 2);
  EXPECT_EQ(result[0]->read_name(), "b");
  EXPECT_EQ(result[1]->read_name(), "c");
}
TEST(MethylationAwarePhasingTest, IdentifyInformativeSitesFiltersCorrectly) {
  std::vector<const DeepVariantCall_ReadSupport*> hap1_reads, hap2_reads;
  DeepVariantCall informative_call =
      MakeMethylCall(100, 3, 3, 25, 230, hap1_reads, hap2_reads);
  DeepVariantCall low_coverage_call =
      MakeMethylCall(200, 1, 1, 25, 230, hap1_reads, hap2_reads);
  DeepVariantCall low_total_read_support =
      MakeMethylCall(300, 2, 2, 25, 230, hap1_reads, hap2_reads);
  DeepVariantCall low_mean_methyl_diff_call =
      MakeMethylCall(200, 3, 3, 125, 130, hap1_reads, hap2_reads);
  DeepVariantCall high_stddev_call =
      MakeMethylCall(400, 3, 3, 10, 200, hap1_reads, hap2_reads);
  // Add one more read to high_stddev_call to make the stddev high on one side.
  auto* r = high_stddev_call.mutable_ref_support_ext()->add_read_infos();
  r->set_read_name("hap1_3");
  r->set_methylation_level(250);
  hap1_reads.push_back(r);

  std::vector<DeepVariantCall> calls = {
      informative_call, low_coverage_call, low_total_read_support,
      low_mean_methyl_diff_call, high_stddev_call};
  auto informative_sites =
      IdentifyInformativeSitesAndUpdatePValues(hap1_reads, hap2_reads, calls);
  ASSERT_EQ(informative_sites.size(), 1);
  EXPECT_EQ(informative_sites[0].variant().start(), 100);
  EXPECT_NEAR(informative_sites[0].methylation_p_value(), 0.049534,
              1e-6);  // informative call
}

TEST(MethylationAwarePhasingTest, HaplotypeVoteWithMethylationVotesCorrectly) {
  std::vector<DeepVariantCall> informative_calls;

  for (int i = 0; i < 3; ++i) {
    DeepVariantCall call;
    call.mutable_variant()->set_start(1000 + i);
    auto* ext = call.mutable_ref_support_ext();

    auto* hap1 = ext->add_read_infos();
    hap1->set_read_name("hap1");
    hap1->set_methylation_level(25);  // ~0.1

    auto* hap2 = ext->add_read_infos();
    hap2->set_read_name("hap2");
    hap2->set_methylation_level(230);  // ~0.9

    auto* unphased = ext->add_read_infos();
    unphased->set_read_name("unphased");
    unphased->set_methylation_level(240);  // closer to hap2

    informative_calls.push_back(call);
  }

  std::vector<const DeepVariantCall_ReadSupport*> hap1_reads;
  std::vector<const DeepVariantCall_ReadSupport*> hap2_reads;

  for (const auto& call : informative_calls) {
    for (const auto& r : call.ref_support_ext().read_infos()) {
      if (r.read_name() == "hap1") hap1_reads.push_back(&r);
      if (r.read_name() == "hap2") hap2_reads.push_back(&r);
    }
  }

  // Pick unphased read from the first call to test
  const auto& unphased_read =
      informative_calls[0].ref_support_ext().read_infos(2);

  int vote = HaplotypeVoteWithMethylation(
    unphased_read, informative_calls, hap1_reads, hap2_reads);
  EXPECT_EQ(vote, 2);  // Should assign to hap2
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

