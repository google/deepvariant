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

#include <string_view>

#include "deepvariant/protos/deepvariant.pb.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::Variant;
using ::testing::UnorderedElementsAreArray;

using AlleleSupportMap =
    absl::flat_hash_map<std::string, std::vector<std::string>>;

DeepVariantCall MakeCandidate(
    int64_t start, int64_t end,
    AlleleSupportMap allele_support = AlleleSupportMap()) {
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

std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
CreateTestReads(const std::vector<ReadFields>& reads) {
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
      reads_out;
  for (const auto& read : reads) {
    reads_out.push_back(new nucleus::genomics::v1::Read(nucleus::MakeRead(
        read.chr, read.position, read.bases, read.cigar, read.read_name)));
  }
  return reads_out;
}

TEST(DirectPhasingTest, BuildGraphSimple) {
  DirectPhasing direct_phasing;

  // Create test candidates.
  std::vector<DeepVariantCall> candidates = {
      MakeCandidate(100, 101,
                    {{"A", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"C", {"read4/0", "read5/0"}}}             // SUB allele
                    ),
      MakeCandidate(105, 106,
                    {
                     {"C", {"read1/0" , "read2/0", "read4/0", "read5/0"}}}
                    ),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};


  // Create test reads.
  // NOTE: The read content "ACGTTGACTTGC" here isn't actually used in the logic
  // for this test, because the `candidates` are already created. You can ignore
  // "ACGTTGACTTGC" when reading this test.
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads({
          {"read1", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read2", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read3", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read4", "chr1", 109, "ACGTTGACTTGC", {"12M"}},
          {"read5", "chr1", 109, "ACGTTGACTTGC", {"12M"}},
      });

  // Manually define an expected vertices.
  std::vector<AlleleInfo> expected_vertices_list = {
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 100,
                 .bases = "A",
                 .read_support = {{0, false}, {1, false}, {2, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 110,
                 .bases = "T",
                 .read_support = {{0, false}, {1, false}, {2, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 105,
                 .bases = "C",
                 .read_support = {{0, false}, {1, false}, {3, false},
                                  {4, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 100,
                 .bases = "C",
                 .read_support = {{3, false}, {4, false}}},
      AlleleInfo{.type = AlleleType::SUBSTITUTION,
                 .position = 110,
                 .bases = "G",
                 .read_support = {{3, false}, {4, false}}}};

  // Manually define expected edges.
  std::vector<std::pair<AlleleInfo, AlleleInfo>> expected_edges_list = {
      {expected_vertices_list[0], expected_vertices_list[2]},
      {expected_vertices_list[3], expected_vertices_list[2]},
      {expected_vertices_list[2], expected_vertices_list[1]},
      {expected_vertices_list[2], expected_vertices_list[4]}
  };

  direct_phasing.Build(candidates, reads);

  // Populate a list of edges that can be used by test comparator.
  std::vector<std::pair<AlleleInfo, AlleleInfo>> graph_edges;
  EdgeIterator ei, eend;
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
  VertexIterator vi, vend;
  std::tie(vi, vend) = boost::vertices(direct_phasing.graph_);
  for (; vi != vend; ++vi) {
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

VertexIterator FindVertex(const BoostGraph& graph, const AlleleInfo& ai) {
  VertexIterator vi, vi_end;
  for (std::tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
    if (graph[*vi].allele_info.position == ai.position &&
        graph[*vi].allele_info.bases == ai.bases)
      return vi;
  }
  return vi_end;
}

bool operator==(const Score& score1, const Score& score2)  {
      return score1.score == score2.score
              && std::equal(std::begin(score1.from), std::end(score1.from),
                            std::begin(score2.from))
              && std::equal(std::begin(score1.read_support),
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
                    {{"C", {"read1/0", "read2/0", "read4/0", "read5/0"}}}),
      MakeCandidate(110, 111,
                    {{"T", {"read1/0", "read2/0", "read3/0"}},  // SUB allele
                     {"G", {"read4/0", "read5/0"}}}             // SUB allele
                    )};

  // Create test reads.
  // NOTE: The read content "ACGTTGACTTGC" here isn't actually used in the logic
  // for this test, because the `candidates` are already created. You can ignore
  // "ACGTTGACTTGC" when reading this test.
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads({
          {"read1", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read2", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read3", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read4", "chr1", 109, "ACGTTGACTTGC", {"12M"}},
          {"read5", "chr1", 109, "ACGTTGACTTGC", {"12M"}},
      });

  direct_phasing.Build(candidates, reads);
  Vertex v_100_a = *FindVertex(direct_phasing.graph_,
                                 {AlleleType::SUBSTITUTION, 100, "A", {}});
  Vertex v_100_c = *FindVertex(direct_phasing.graph_,
                  {AlleleType::SUBSTITUTION, 100, "C", {}});
  Vertex v_105_c = *FindVertex(direct_phasing.graph_,
                  {AlleleType::SUBSTITUTION, 105, "C", {}});
  direct_phasing.UpdateStartingScore({v_100_a, v_100_c});
  Edge edge1, edge2;
  bool found = false;
  tie(edge1, found) = boost::edge(v_100_a, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);
  tie(edge2, found) = boost::edge(v_100_c, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);

  Score calculated_score = direct_phasing.CalculateScore(
      edge1, edge2);
  EXPECT_EQ(calculated_score,
             Score({
               .score = 4,
               .from = {v_100_a, v_100_c},
               .read_support = {{0, 1}, {3, 4}}
             }));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

TEST(DirectPhasingTest, CalculateScoreWirhPreviousScore) {
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

  // Create test reads.
  // NOTE: The read content "ACGTTGACTTGC" here isn't actually used in the logic
  // for this test, because the `candidates` are already created. You can ignore
  // "ACGTTGACTTGC" when reading this test.
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads =
      CreateTestReads({
          {"read1", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read2", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read3", "chr1", 99, "ACGTTGACTTGC", {"12M"}},
          {"read4", "chr1", 109, "ACGTTGACTTGC", {"12M"}},
          {"read5", "chr1", 109, "ACGTTGACTTGC", {"12M"}},
      });

  // Build graph
  direct_phasing.Build(candidates, reads);

  // Find all vertices.
  Vertex v_100_a = *FindVertex(direct_phasing.graph_,
                   {AlleleType::SUBSTITUTION, 100, "A", {}});
  Vertex v_100_c = *FindVertex(direct_phasing.graph_,
                  {AlleleType::SUBSTITUTION, 100, "C", {}});
  Vertex v_105_c = *FindVertex(direct_phasing.graph_,
                  {AlleleType::SUBSTITUTION, 105, "C", {}});
  Vertex v_110_t = *FindVertex(direct_phasing.graph_,
                  {AlleleType::SUBSTITUTION, 110, "T", {}});
  Vertex v_110_g = *FindVertex(direct_phasing.graph_,
                  {AlleleType::SUBSTITUTION, 110, "G", {}});

  // Update starting score.
  direct_phasing.UpdateStartingScore({v_100_a, v_100_c});
  Edge edge1, edge2;
  bool found = false;

  // Update the score for {edge1, edge2}
  tie(edge1, found) = boost::edge(v_100_a, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);
  tie(edge2, found) = boost::edge(v_100_c, v_105_c, direct_phasing.graph_);
  EXPECT_TRUE(found);
  direct_phasing.scores_[{v_105_c, v_105_c}] = direct_phasing.CalculateScore(
      edge1, edge2);

  // Verify scores for all combination of edge1, and edge2 (edge1, edge2
  // variables are reused).
  tie(edge1, found) = boost::edge(v_105_c, v_110_t, direct_phasing.graph_);
  EXPECT_TRUE(found);
  tie(edge2, found) = boost::edge(v_105_c, v_110_g, direct_phasing.graph_);
  EXPECT_TRUE(found);

  EXPECT_EQ(direct_phasing.CalculateScore(edge1, edge1),
            Score({
              .score = 4 + 2,
              .from = {v_105_c, v_105_c},
              .read_support = {{0, 1}, {}}
            }));
  EXPECT_EQ(direct_phasing.CalculateScore(edge2, edge2),
             Score({
              .score = 4 + 2,
              .from = {v_105_c, v_105_c},
              .read_support = {{}, {3, 4}}
            }));
  EXPECT_EQ(direct_phasing.CalculateScore(edge1, edge2),
           Score({
              .score = 4 + 4,
              .from = {v_105_c, v_105_c},
              .read_support = {{0, 1}, {3, 4}}
            }));
  EXPECT_EQ(direct_phasing.CalculateScore(edge2, edge1),
            Score({
              .score = 4 + 0,
              .from = {v_105_c, v_105_c},
              .read_support = {{}, {}}
            }));

  // Release memory.
  for (auto read : reads) {
    delete read.p_;
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
