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

#include <sys/types.h>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
// #include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::Variant;
using ::testing::ElementsAreArray;
using ::testing::UnorderedElementsAreArray;
using AlleleSupportMap =
    absl::flat_hash_map<std::string, std::vector<std::string>>;
using RefSupport = std::vector<std::string>;

DeepVariantCall MakeCandidate(
    int64_t start, int64_t end,
    const AlleleSupportMap& allele_support = AlleleSupportMap(),
    const RefSupport& ref_support = RefSupport()) {
  DeepVariantCall candidate;
  Variant* variant = candidate.mutable_variant();
  variant->set_start(start);
  variant->set_end(end);
  if (!allele_support.empty()) {
    auto allele_support_field = candidate.mutable_allele_support_ext();
    for (const auto& one_allele_support : allele_support) {
      auto& support_infos = (*allele_support_field)[one_allele_support.first];
      for (const auto& read_name : one_allele_support.second) {
        auto* read_info = support_infos.add_read_infos();
        read_info->set_read_name(read_name);
        read_info->set_is_low_quality(false);
      }
    }
  }
  if (!ref_support.empty()) {
    auto ref_support_field = candidate.mutable_ref_support_ext();
    for (const auto& ref_support_read : ref_support) {
      auto* read_info = ref_support_field->add_read_infos();
      read_info->set_read_name(ref_support_read);
      read_info->set_is_low_quality(false);
    }
  }
  return candidate;
}

TEST(DirectPhasingTest, TestAlleleTypeFromCandidateSubstitution) {
  EXPECT_EQ(AlleleType::SUBSTITUTION,
            AlleleTypeFromCandidate("CC", MakeCandidate(100, 102)));
}

TEST(DirectPhasingTest, TestAlleleTypeFromCandidateDeletion) {
  EXPECT_EQ(AlleleType::DELETION,
            AlleleTypeFromCandidate("C", MakeCandidate(100, 102)));
}

TEST(DirectPhasingTest, TestAlleleTypeFromCandidateInsertion) {
  EXPECT_EQ(AlleleType::INSERTION,
            AlleleTypeFromCandidate("CCC", MakeCandidate(100, 101)));
}

TEST(DirectPhasingTest, TestAlleleTypeFromCandidateOneBaseSubstitution) {
  EXPECT_EQ(AlleleType::SUBSTITUTION,
            AlleleTypeFromCandidate("A", MakeCandidate(100, 101)));
}

struct TestNumOfSubstitutionAllelesTestCase {
  DeepVariantCall candidate;
  int expected_allele_depth;
};

TEST(DirectPhasingTest, TestNumOfSubstitutionAllelesMultipleSubAlleles) {
  EXPECT_EQ(2,
    NumOfSubstitutionAlleles(MakeCandidate(100, 101, {
      {"A", {"read1", "read2", "read3"}},  // SUB allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST(DirectPhasingTest, TestNumOfSubstitutionAllelesUncalledAllelePresent) {
  EXPECT_EQ(1,
    NumOfSubstitutionAlleles(MakeCandidate(100, 101, {
      {"UNCALLED_ALLELE", {"read1", "read2", "read3"}},  // Uncalled allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST(DirectPhasingTest, TestNumOfIndelAlleles2Sub1Indel) {
  EXPECT_EQ(1,
    NumOfIndelAlleles(MakeCandidate(100, 101, {
      {"A", {"read1", "read2", "read3"}},  // SUB allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST(DirectPhasingTest, TestNumOfIndelAllelesUncalledPresent) {
  EXPECT_EQ(2,
    NumOfIndelAlleles(MakeCandidate(100, 103, {
      {"UNCALLED_ALLELE", {"read1", "read2", "read3"}},  // Uncalled allele
      {"C", {"read4", "read5"}},          // DEL allele
      {"CCCC", {"read6", "read7"}}}         // INS allele
      )));
}

TEST(DirectPhasingTest, TestSubstitutionAllelesDepth2SubAlleles) {
  EXPECT_EQ(5,
    SubstitutionAllelesDepth(MakeCandidate(100, 101, {
      {"A", {"read1", "read2", "read3"}},  // SUB allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST(DirectPhasingTest, TestSubstitutionAllelesDepth1UncalledAnd2Indels) {
  EXPECT_EQ(0,
    SubstitutionAllelesDepth(MakeCandidate(100, 103, {
      {"UNCALLED_ALLELE", {"read1", "read2", "read3"}},  // Uncalled allele
      {"C", {"read4", "read5"}},          // DEL allele
      {"CCCC", {"read6", "read7"}}}         // INS allele
      )));
}

void PopulateReadSupportInfo(
    const std::vector<std::pair<std::string, bool>>& read_supports,
    google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport>&
        read_support_proto) {
  for (const auto& read_support : read_supports) {
    learning::genomics::deepvariant::DeepVariantCall_ReadSupport*
        read_support_item = read_support_proto.Add();
    read_support_item->set_read_name(read_support.first);
    read_support_item->set_is_low_quality(read_support.second);
  }
}

TEST(DirectPhasingTest, ReadSupportFromProtoSimple) {
  std::vector<ReadSupportInfo> expected_read_support_infos{
      {1, false},
      {2, false},
  };
  DirectPhasing direct_phasing;
  google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport> read_support_proto;
  PopulateReadSupportInfo({{"read1", false}, {"read2", false}},
                          read_support_proto);

  direct_phasing.read_to_index_ = {
      {"read1", 1}, {"read2", 2}, {"read3", 3}, {"read4", 4}, {"read5", 5}};

  EXPECT_THAT(direct_phasing.ReadSupportFromProto(read_support_proto),
              UnorderedElementsAreArray(expected_read_support_infos));
}

TEST(DirectPhasingTest, ReadSupportFromProtoLQReads) {
  std::vector<ReadSupportInfo> expected_read_support_infos{
      {1, false},
      {2, false},
  };
  DirectPhasing direct_phasing;
  google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport> read_support_proto;
  PopulateReadSupportInfo({{"read1", false}, {"read2", false},
                          {"read3", true}},
                          read_support_proto);

  direct_phasing.read_to_index_ = {
      {"read1", 1}, {"read2", 2}, {"read3", 3}, {"read4", 4}, {"read5", 5}};

  EXPECT_THAT(direct_phasing.ReadSupportFromProto(read_support_proto),
              UnorderedElementsAreArray(expected_read_support_infos));
}

struct ReadFields {
  std::string read_name;
  std::string chr;
  int position;
  std::string bases;
  std::vector<std::string> cigar;
};

// Creates test reads.
// Only read names are used in tests. All other read fields do not affect the
// logic of tests.
std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
CreateTestReads(int num_of_reads) {
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
      reads_out;
  for (int i = 1; i <= num_of_reads; i++) {
    reads_out.push_back(new nucleus::genomics::v1::Read(
        nucleus::MakeRead("chr1", 89 + i, "ACGTTGACTTGC", {"12M"},
                          absl::StrCat("read", std::to_string(i)))));
  }
  return reads_out;
}

TEST(DirectPhasingTest, BuildGraphSimple) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {
                     {"C", {"read1/0" , "read2/0", "read3/0"}}},
                     {"read4/0", "read5/0", "read6/0"}
                    ),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(6);

  // Manually define an expected vertices.
  std::vector<AlleleInfo> expected_vertices_list = {
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 100,
                 .bases = "A",
                 .read_support = {{0, false}, {1, false}, {2, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 100,
                 .bases = "C",
                 .read_support = {{3, false}, {4, false}, {5, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 105,
                 .bases = "C",
                 .read_support = {{0, false}, {1, false}, {2, false}}},
      AlleleInfo{.type = AlleleType::REFERENCE,
                 .position = 105,
                 .bases = "REF",
                 .read_support = {{3, false}, {4, false}, {5, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 110,
                 .bases = "T",
                 .read_support = {{0, false}, {1, false}, {2, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 110,
                 .bases = "G",
                 .read_support = {{3, false}, {4, false}}}};

  // Manually define expected edges.
  std::vector<std::pair<AlleleInfo, AlleleInfo>> expected_edges_list = {
      {expected_vertices_list[0], expected_vertices_list[2]},
      {expected_vertices_list[1], expected_vertices_list[3]},
      {expected_vertices_list[2], expected_vertices_list[4]},
      {expected_vertices_list[3], expected_vertices_list[5]},
  };

  direct_phasing.Build(candidates, reads);

  // Populate a list of edges that can be used by test comparator.
  std::vector<std::pair<AlleleInfo, AlleleInfo>> graph_edges;
  DirectPhasing::EdgeIterator ei, eend;
  std::tie(ei, eend) = boost::edges(direct_phasing.graph_);
  for (; ei != eend; ei++) {
    graph_edges.push_back(
        std::pair(direct_phasing.graph_[ei->m_source].allele_info,
                  direct_phasing.graph_[ei->m_target].allele_info));
    // gtest comparator does not output per field differences. If test fails
    // it is easier to debug if edges are printed here.
    // LOG(WARNING) << "Edge: "
    //     << direct_phasing.graph_[ei->m_source].allele_info.position << " "
    //     << direct_phasing.graph_[ei->m_source].allele_info.bases << "-"
    //     << direct_phasing.graph_[ei->m_target].allele_info.position << " "
    //     << direct_phasing.graph_[ei->m_target].allele_info.bases;
  }
  std::vector<AlleleInfo> graph_vertices;
  DirectPhasing::VertexIterator vi, vend;
  std::tie(vi, vend) = boost::vertices(direct_phasing.graph_);
  for (; vi != vend; ++vi) {
    // gtest comparator does not output per field differences. If test fails
    // it is easier to debug if vertices are printed here.
    // std::ostringstream ss;
    // for (auto read_info :
    //          direct_phasing.graph_[*vi].allele_info.read_support) {
    //   ss << read_info.read_index << ",";
    // }
    // LOG(WARNING) << "Vertex: "
    //     << direct_phasing.graph_[*vi].allele_info.position << " "
    //     << direct_phasing.graph_[*vi].allele_info.bases << " "
    //     << ss.str();
    graph_vertices.push_back(direct_phasing.graph_[*vi].allele_info);
    // gtest comparator does not output per field differences. If test fails
    // it is easier to debug if edges are printed here.
  }

  EXPECT_THAT(graph_vertices, UnorderedElementsAreArray(
      expected_vertices_list));
  EXPECT_THAT(graph_edges, UnorderedElementsAreArray(expected_edges_list));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

DirectPhasing::VertexIterator FindVertex(const DirectPhasing::BoostGraph& graph,
                                         const AlleleInfo& ai) {
  DirectPhasing::VertexIterator vi, vi_end;
  for (std::tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
    if (graph[*vi].allele_info.position == ai.position &&
        graph[*vi].allele_info.bases == ai.bases)
      return vi;
  }
  return vi_end;
}

bool operator==(const DirectPhasing::Score& score1,
                const DirectPhasing::Score& score2) {
  return score1.score == score2.score &&
         std::equal(std::begin(score1.from), std::end(score1.from),
                    std::begin(score2.from)) &&
         std::equal(std::begin(score1.read_support),
                    std::end(score1.read_support),
                    std::begin(score2.read_support));
}

TEST(DirectPhasingTest, CalculateScoreFirstIteration) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read1/0", "read2/0", "read4/0", "read5/0"}}},
                    {"read6/0", "read7/0", "read8/0"}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(8);

  direct_phasing.Build(candidates, reads);
  DirectPhasing::Vertex v_100_a = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 100, "A", {}});
  DirectPhasing::Vertex v_100_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 100, "C", {}});
  DirectPhasing::Vertex v_105_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 105, "C", {}});
  direct_phasing.UpdateStartingScore({v_100_a, v_100_c});
  DirectPhasing::Edge edge1, edge2;
  bool found = false;
  tie(edge1, found) = boost::edge(v_100_a, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);
  tie(edge2, found) = boost::edge(v_100_c, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);

  DirectPhasing::Score calculated_score =
      direct_phasing.CalculateScore(edge1, edge2);
  EXPECT_EQ(calculated_score,
            (DirectPhasing::Score{.score = 5 + 4,
                                  .from = {v_100_a, v_100_c},
                                  .read_support = {{0, 1}, {3, 4}}}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, CalculateScoreWithPreviousScore) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read1/0", "read2/0", "read4/0", "read5/0"}}},
                    {"read6/0", "read7/0", "read8/0"}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(8);

  // Build graph
  direct_phasing.Build(candidates, reads);

  // Find all vertices.
  DirectPhasing::Vertex v_100_a = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 100, "A", {}});
  DirectPhasing::Vertex v_100_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 100, "C", {}});
  DirectPhasing::Vertex v_105_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 105, "C", {}});
  DirectPhasing::Vertex v_110_t = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 110, "T", {}});
  DirectPhasing::Vertex v_110_g = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 110, "G", {}});

  // Update starting score.
  direct_phasing.UpdateStartingScore({v_100_a, v_100_c});
  DirectPhasing::Edge edge1, edge2;
  bool found = false;

  // Update the score for {edge1, edge2}
  tie(edge1, found) = boost::edge(v_100_a, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);
  tie(edge2, found) = boost::edge(v_100_c, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);
  direct_phasing.scores_[{v_105_c, v_105_c}] = direct_phasing.CalculateScore(
      edge1, edge2);

  // Verify scores for all combinations of edge1 and edge2.
  tie(edge1, found) = boost::edge(v_105_c, v_110_t, direct_phasing.graph_);
  EXPECT_TRUE(found);
  tie(edge2, found) = boost::edge(v_105_c, v_110_g, direct_phasing.graph_);
  EXPECT_TRUE(found);

  EXPECT_EQ(direct_phasing.CalculateScore(edge1, edge1),
            (DirectPhasing::Score{.score = 5 + 4 + 2,
                                  .from = {v_105_c, v_105_c},
                                  .read_support = {{0, 1}, {}}}));
  EXPECT_EQ(direct_phasing.CalculateScore(edge2, edge2),
            (DirectPhasing::Score{.score = 5 + 4 + 2,
                                  .from = {v_105_c, v_105_c},
                                  .read_support = {{}, {3, 4}}}));
  EXPECT_EQ(direct_phasing.CalculateScore(edge1, edge2),
            (DirectPhasing::Score{.score = 5 + 4 + 4,
                                  .from = {v_105_c, v_105_c},
                                  .read_support = {{0, 1}, {3, 4}}}));
  EXPECT_EQ(direct_phasing.CalculateScore(edge2, edge1),
            (DirectPhasing::Score{.score = 5 + 4 + 0,
                                  .from = {v_105_c, v_105_c},
                                  .read_support = {{}, {}}}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadSimpleTest) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read1/0", "read2/0", "read4/0", "read5/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(5);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({1, 1, 1, 2, 2}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadWithErrorCorrection) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  // read3 supports phase 1 in the candidate at 100, but it also supports
  // phase 2 in the candidate at 110.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(
          105, 106,
          {{"C", {"read1/0", "read2/0", "read3/0", "read4/0", "read5/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0"}},             // SUB allele
                     {"G", {"read3/0", "read4/0", "read5/0"}}}  // SUB allele
                    ),
      MakeCandidate(120, 121,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(5);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({1, 1, 1, 2, 2}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadChangedOrderOfAlleles) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  // read3 supports phase 1 in the candidate at 100, but it also supports
  // phase 2 in the candidate at 110.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(
          105, 106,
          {{"C", {"read1/0", "read2/0", "read3/0", "read4/0", "read5/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read4/0", "read5/0"}},             // SUB allele
                     {"G", {"read1/0", "read2/0", "read3/0"}}}  // SUB allele
                    ),
      MakeCandidate(120, 121,
                    {
                        {"G", {"read4/0", "read5/0"}},
                        {"T", {"read1/0", "read2/0", "read3/0"}}  // SUB allele
                    }                                             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(5);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({1, 1, 1, 2, 2}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadUnphasedRead) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  // read 3 overlaps one allele phase1, one allele phase 2 and pne homozygous
  // allele. Phase of read 3 has to be unassigned.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(
          105, 106,
          {{"C", {"read1/0", "read2/0", "read3/0", "read4/0", "read5/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0"}},             // SUB allele
                     {"G", {"read4/0", "read5/0", "read3/0"}}}  // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(5);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({1, 1, 0, 2, 2}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadBrokenPath) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  // No edge between A at 100 and G at 105
  // 100     105     110
  // A       G ----- T   Phase 1
  //
  // C ----- C ----- G   Phase 2
  // In this example reads 1,2,3 can be assigned any phase 1 or no phase, but
  // algorithm favors assigning phase1.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(
          105, 106,
          {{"C", {"read4/0", "read5/0"}}, {"G", {"read6/0", "read7/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read6/0", "read7/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}  // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(7);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({0, 0, 0, 2, 2, 1, 1}));
  DirectPhasing::Vertex v_105_g = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 105, "G", {}});
  DirectPhasing::Vertex v_105_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 105, "C", {}});

  EXPECT_THAT(direct_phasing.graph_[v_105_g].allele_info.read_support,
              UnorderedElementsAreArray(
      {
        ReadSupportInfo{
          .read_index = 5,
          .is_low_quality = false,
          .is_first_allele = true
        },
        ReadSupportInfo{
          .read_index = 6,
          .is_low_quality = false,
          .is_first_allele = true
        }
      }));

  EXPECT_THAT(direct_phasing.graph_[v_105_c].allele_info.read_support,
              UnorderedElementsAreArray(
      {
        ReadSupportInfo{
          .read_index = 3,
          .is_low_quality = false,
          .is_first_allele = false
        },
        ReadSupportInfo{
          .read_index = 4,
          .is_low_quality = false,
          .is_first_allele = false
        }
      }));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadBrokenPathNoConnection) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  // No edge between C at 105 and G at 110, and not edges between G at 105 and
  // C at 110. There is no connection between subragphs, so we should restart
  // phasing from position 110.
  // 100     105   110    120
  // A ---- C      G ----- T   Phase 1
  //
  // C ---- G      C ----- G   Phase 2
  // In this example reads 1,2,3 can be assigned any phase 1 or no phase, but
  // algorithm favors assigning phase1.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},
                     {"C", {"read4/0", "read5/0"}}}),
      MakeCandidate(105, 106,
                    {{"C", {"read1/0", "read2/0", "read3/0"}},
                     {"G", {"read4/0", "read5/0"}}}),
      MakeCandidate(
          110, 111,
          {{"C", {"read6/0", "read7/0"}}, {"G", {"read8/0", "read9/0"}}}),
      MakeCandidate(
          120, 121,
          {{"T", {"read6/0", "read7/0"}}, {"G", {"read8/0", "read9/0"}}})};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(9);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(),
              ElementsAreArray({1, 1, 1, 2, 2, 1, 1, 2, 2}));
  DirectPhasing::Vertex v_110_g = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 110, "G", {}});
  DirectPhasing::Vertex v_110_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 110, "C", {}});
  DirectPhasing::Score score_110_C_G =
      DirectPhasingPeer::FindScore(direct_phasing, v_110_c, v_110_g);

  EXPECT_EQ(score_110_C_G.score, 4);

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, NotPhasablePosition) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  // At position 110 all possible partitions yield the same score. In this case
  // we cannot phase this position.
  // 100     105   110    120    125
  // A ---- C ---- G      T ---- A  Phase 1
  //          \ /
  //          / \
  // C ---- G ---- C      G ---- T   Phase 2
  // In this example reads 1,2,3 can be assigned any phase 1 or no phase, but
  // algorithm favors assigning phase1.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0", "read10/0"}},
                     {"C", {"read4/0", "read5/0"}}}),
      MakeCandidate(
          105, 106,
          {{"C", {"read1/0", "read2/0", "read3/0", "read10/0", "read11/0"}},
           {"G", {"read4/0", "read5/0", "read12/0", "read13/0"}}}),
      MakeCandidate(
          110, 111,
          {{"C", {"read10/0", "read13/0"}}, {"G", {"read11/0", "read12/0"}}}),
      MakeCandidate(
          120, 121,
          {{"T", {"read6/0", "read7/0"}}, {"G", {"read8/0", "read9/0"}}}),
      MakeCandidate(
          125, 126,
          {{"A", {"read6/0", "read7/0"}}, {"T", {"read8/0", "read9/0"}}})};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(13);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  // EXPECT_THAT(phases.ValueOrDie(),
  //             ElementsAreArray({1, 1, 1, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0}));
  DirectPhasing::Vertex v_110_g = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 110, "G", {}});
  DirectPhasing::Vertex v_110_c = *FindVertex(
      direct_phasing.graph_, {AlleleType::SUBSTITUTION, 110, "C", {}});

  DirectPhasing::Score score_110_C_G =
      DirectPhasingPeer::FindScore(direct_phasing, v_110_c, v_110_g);
  DirectPhasing::Score score_110_C_C =
      DirectPhasingPeer::FindScore(direct_phasing, v_110_c, v_110_c);
  DirectPhasing::Score score_110_G_G =
      DirectPhasingPeer::FindScore(direct_phasing, v_110_g, v_110_g);

  // All scores at 110 should be equal.
  EXPECT_EQ(score_110_C_G.score, score_110_C_C.score);
  EXPECT_EQ(score_110_C_C.score, score_110_G_G.score);
  // Candidate at 110 should be unphased.
  EXPECT_EQ(direct_phasing.graph_[v_110_g].allele_info.phase, 0);
  EXPECT_EQ(direct_phasing.graph_[v_110_c].allele_info.phase, 0);

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadFullyConnectedGraph) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read4/0", "read5/0", "read1/0"}},
                     {"G", {"read2/0", "read3/0", "read6/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(6);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({1, 1, 1, 2, 2, 2}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadUnorderedInputFail) {
  DirectPhasing direct_phasing;

  // Create test candidates that are not ordered by position.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(105, 106,
                    {{"C", {"read4/0", "read5/0", "read1/0"}},
                     {"G", {"read2/0", "read3/0", "read6/0"}}}),
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    ),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(6);

  EXPECT_DEATH(direct_phasing.PhaseReads(candidates, reads), "");

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, PhaseReadCandidateOutOfOrderInTheMiddle) {
  DirectPhasing direct_phasing;

  // Create test candidates that are not ordered by position.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read4/0", "read5/0", "read1/0"}},
                     {"G", {"read2/0", "read3/0", "read6/0"}}}),
      MakeCandidate(104, 105,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    ),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(6);

  EXPECT_DEATH(direct_phasing.PhaseReads(candidates, reads), "");

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

// Test verifies that candidate that has only one allele and less than 3 reads
// supporting reference is filtered out.
TEST(DirectPhasingTest, FilterOneAlleleCandidate) {
  DirectPhasing direct_phasing;

  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {
                     {"C", {"read4/0", "read5/0", "read6/0"}}},  // SUB allele
                     {"read7/0"}  // One read supporting REF
                    ),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
  )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(7);

  AlleleInfo v_100_c = {AlleleType::SUBSTITUTION, 100, "C", {}};

  direct_phasing.Build(candidates, reads);

  // We expect that vertex at position 100 is not created.
  DirectPhasing::VertexIterator vi = FindVertex(direct_phasing.graph_, v_100_c);
  EXPECT_EQ(vi, direct_phasing.graph_.m_vertices.cend());

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

// Test verifies that candidate containing INDEL is filtered out.
TEST(DirectPhasingTest, FilterCandidateWithIndel) {
  DirectPhasing direct_phasing;

  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 102,
                     {
                      {"CC", {"read4/0", "read5/0", "read6/0"}},  // SUB allele
                      {"A", {"read1/0", "read2/0"}},  // INDEL allele
                     },
                     {"read7/0"}  // One read supporting REF
                    ),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
  )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(7);

  AlleleInfo v_100_c = {AlleleType::SUBSTITUTION, 100, "C", {}};

  direct_phasing.Build(candidates, reads);

  // We expect that vertex at positopm 100 is not created.
  DirectPhasing::VertexIterator vi = FindVertex(direct_phasing.graph_, v_100_c);
  EXPECT_EQ(vi, direct_phasing.graph_.m_vertices.cend());

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, DirectPhasingReuseObject) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read1/0", "read2/0", "read4/0", "read5/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(5);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());
  EXPECT_THAT(phases.ValueOrDie(), ElementsAreArray({1, 1, 1, 2, 2}));

  std::vector<DeepVariantCall> candidates2 = {
      MakeCandidate(120, 121,
                    {{"G", {"read1/0", "read2/0", "read3/0"}},
                     {"A", {"read4/0", "read5/0"}}}
                    ),
      MakeCandidate(130, 131,
                    {
                      {"T", {"read1/0", "read2/0", "read3/0",
                             "read4/0", "read5/0"}}
                    })};

  nucleus::StatusOr<std::vector<int>> phases2 =
      direct_phasing.PhaseReads(candidates2, reads);
  EXPECT_TRUE(phases2.ok());
  EXPECT_THAT(phases2.ValueOrDie(), ElementsAreArray({0, 0, 0, 0, 0}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

bool operator==(const PhasedVariant& a, const PhasedVariant& b) {
  return a.position == b.position && a.phase_1_bases == b.phase_1_bases &&
         a.phase_2_bases == b.phase_2_bases;
}

TEST(DirectPhasingTest, GetPhasedVariantsSanity) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {{"C", {"read4/0", "read5/0", "read1/0"}},
                     {"G", {"read2/0", "read3/0", "read6/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0", "read6/0"}}}  // SUB allele
                    )};

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads(6);

  nucleus::StatusOr<std::vector<int>> phases =
      direct_phasing.PhaseReads(candidates, reads);
  EXPECT_TRUE(phases.ok());

  EXPECT_THAT(
      direct_phasing.GetPhasedVariants(),
      ElementsAreArray(std::vector<PhasedVariant>{
          {.position = 100, .phase_1_bases = "A", .phase_2_bases = "C"},
          {.position = 105, .phase_1_bases = "G", .phase_2_bases = "C"},
          {.position = 110, .phase_1_bases = "T", .phase_2_bases = "G"}}));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
