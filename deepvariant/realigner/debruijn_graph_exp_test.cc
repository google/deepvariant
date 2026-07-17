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

#include "deepvariant/realigner/debruijn_graph_exp.h"

#include <cstddef>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using ::testing::ElementsAre;
using ::testing::UnorderedElementsAre;

class DeBruijnGraphTest : public ::testing::Test {
 protected:
  DeBruijnGraphOptions options() {
    DeBruijnGraphOptions options;
    options.set_min_k(10);
    options.set_max_k(100);
    options.set_step_k(1);
    options.set_min_mapq(20);
    options.set_min_base_quality(20);
    options.set_min_edge_weight(1);
    options.set_max_num_paths(10);
    return options;
  }

  nucleus::genomics::v1::Read MakeRead(
      const std::string& seq) {
    nucleus::genomics::v1::Read read;
    read.set_aligned_sequence(seq);
    read.set_aligned_quality(std::string(seq.size(), 30));
    read.mutable_alignment()->set_mapping_quality(30);
    return read;
  }
};

TEST_F(DeBruijnGraphTest, TestCollapse) {
  std::string ref = "GATTACA";
  // Use k=3.
  // GAT -> ATT -> TTA -> TAC -> ACA
  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);

  auto read = MakeRead("GATTACA");
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads{
    nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read),
    nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read)
  };

  auto dbg = DeBruijnGraphExp::Build(ref, reads, opts);
  ASSERT_NE(dbg, nullptr);

  EXPECT_THAT(dbg->CandidateHaplotypesRanked(0), ElementsAre("GATTACA"));

  // Before collapse, it should have 5 vertices.
  // Actually, I can't easily check vertex count from public API, but I can
  // check haplotypes.

  dbg->Collapse();

  // After collapse, it should still have the same haplotype.
  EXPECT_THAT(dbg->CandidateHaplotypesRanked(0), ElementsAre("GATTACA"));
}

TEST_F(DeBruijnGraphTest, TestCollapseWithBranch) {
  std::string ref = "GATTACA";
  // Create reads that branches and then rejoins.
  // GAT -> ATG -> TGA -> GAC -> ACA
  // GAT -> ATT -> TTA -> TAC -> ACA

  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);
  // We need to use 2 reads for each path to keep edeges from trimming since all
  // edges with weigh less than 2 are trimmed.
  // In addition we need an extra base at the begining to avoid edge trimming at
  // the start.
  auto read1 = MakeRead("AGATGACA");
  auto read2 = MakeRead("GATGACAA");
  auto read3 = MakeRead("AGATTACA");
  auto read4 = MakeRead("GATTACAA");

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads{
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(
          &read1),
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(
          &read2),
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(
          &read3),
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(
          &read4)
  };

  auto dbg = DeBruijnGraphExp::Build(ref, reads, opts);
  ASSERT_NE(dbg, nullptr);

  EXPECT_THAT(dbg->CandidateHaplotypesRanked(0),
              UnorderedElementsAre("AGATGACAA", "AGATTACAA"));

  dbg->Collapse();

  EXPECT_THAT(dbg->CandidateHaplotypesRanked(0),
              UnorderedElementsAre("AGATGACAA", "AGATTACAA"));
}

struct CandidatePathsTestCase {
  std::string ref;
  std::vector<std::string> reads;
  std::vector<std::string> expected_haplotypes;
};

class DeBruijnGraphParameterizedTest
    : public DeBruijnGraphTest,
      public ::testing::WithParamInterface<CandidatePathsTestCase> {};

TEST_P(DeBruijnGraphParameterizedTest, TestCandidatePathsRanked) {
  const CandidatePathsTestCase& param = GetParam();

  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);

  std::vector<nucleus::genomics::v1::Read> reads_storage;
  for (const auto& read_seq : param.reads) {
    reads_storage.push_back(MakeRead(read_seq));
  }

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads;
  for (const auto& read : reads_storage) {
    reads.push_back(
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read));
  }

  auto dbg = DeBruijnGraphExp::Build(param.ref, reads, opts);
  ASSERT_NE(dbg, nullptr);

  EXPECT_THAT(dbg->CandidateHaplotypesRanked(0),
              ::testing::UnorderedElementsAreArray(param.expected_haplotypes));
}

INSTANTIATE_TEST_SUITE_P(
    DeBruijnGraphExpTestCases, DeBruijnGraphParameterizedTest,
    ::testing::Values(
        CandidatePathsTestCase{"GATTACA", {}, {}},
        CandidatePathsTestCase{"GATTACA", {"GATGACA", "GATGACA"},
            {"GATGACA"}}));

TEST_F(DeBruijnGraphTest, TestDeterminismWithBranch) {
  std::string ref = "GATTACA";
  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);

  auto read1 = MakeRead("AGATGACA");
  auto read2 = MakeRead("GATGACAA");
  auto read3 = MakeRead("AGATTACA");
  auto read4 = MakeRead("GATTACAA");

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads{
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read1),
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read2),
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read3),
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read4)};

  // Keep all dust alive across iterations to prevent the allocator from
  // recycling freed addresses, ensuring each iteration genuinely sees
  // different heap layouts.
  std::vector<std::vector<std::unique_ptr<char[]>>> all_dust;
  std::vector<std::vector<std::string>> all_results;
  for (int i = 0; i < 100; ++i) {
    std::vector<std::unique_ptr<char[]>> dust;
    std::mt19937 rng(i);
    for (int j = 0; j < 50; ++j) {
      dust.push_back(std::make_unique<char[]>(rng() % 1000 + 1));
    }
    all_dust.push_back(std::move(dust));

    auto dbg = DeBruijnGraphExp::Build(ref, reads, opts);
    ASSERT_NE(dbg, nullptr);
    auto results = dbg->CandidateHaplotypesRanked(0);
    all_results.push_back(results);
  }

  // Check if all results are identical (including order).
  for (size_t i = 1; i < all_results.size(); ++i) {
    EXPECT_EQ(all_results[i], all_results[0]);
  }
}

TEST_F(DeBruijnGraphTest, TestDeterminismWithBranchAndCollapse) {
  std::string ref = "GATTACA";
  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);

  auto read1 = MakeRead("AGATGACA");
  auto read2 = MakeRead("GATGACAA");
  auto read3 = MakeRead("AGATTACA");
  auto read4 = MakeRead("GATTACAA");

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads{
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read1),
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read2),
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read3),
      nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&read4)};

  // Keep all dust alive across iterations to prevent the allocator from
  // recycling freed addresses, ensuring each iteration genuinely sees
  // different heap layouts.
  std::vector<std::vector<std::unique_ptr<char[]>>> all_dust;
  std::vector<std::vector<std::string>> all_results;
  for (int i = 0; i < 100; ++i) {
    std::vector<std::unique_ptr<char[]>> dust;
    std::mt19937 rng(i);
    for (int j = 0; j < 50; ++j) {
      dust.push_back(std::make_unique<char[]>(rng() % 1000 + 1));
    }
    all_dust.push_back(std::move(dust));

    auto dbg = DeBruijnGraphExp::Build(ref, reads, opts);
    ASSERT_NE(dbg, nullptr);
    dbg->Collapse();
    auto results = dbg->CandidateHaplotypesRanked(0);
    all_results.push_back(results);
  }

  // Check if all results are identical (including order).
  for (size_t i = 1; i < all_results.size(); ++i) {
    EXPECT_EQ(all_results[i], all_results[0]);
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
