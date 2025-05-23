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

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "boost/graph/graphviz.hpp"
#include "third_party/nucleus/core/statusor.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace learning {
namespace genomics {
namespace deepvariant {

const int kMinRefAlleleDepth = 3;
const int kMinAllelesToPhase = 2;
constexpr absl::string_view kRef = "REF";
const int kNumOfPhases = 2;
const float kMinMethylationThreshold = 0.4;
const float kMaxMethylationThreshold = 0.6;

std::string ReadKey(const nucleus::genomics::v1::Read& read) {
  return absl::StrCat(read.fragment_name(), "/", read.read_number());
}

nucleus::StatusOr<std::vector<int>> DirectPhasing::PhaseReads(
    absl::Span<const DeepVariantCall> candidates,
    absl::Span<const nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
        reads) {
  // Build graph from candidates.
  Build(candidates, reads);
  // Iterate positions in order. Calculate the score for each combination of
  // allele pairs.
  for (int i = 0; i < positions_.size(); i++) {
    // TODO Call UpdateStartingScore if score cannot be improved from
    // position to the next position. This happens with bad data where reads
    // are mismapped. Good example is chr1:143175001-143200000.
    // In this case we can break the graph and treat each piece as separate
    // graphs. This TODO has a lower impact because it happens with "bad" data
    // and DeepVariant will reject these candidates in most of the cases.
    // The work is tracked in internal
    if (i == 0) {
      UpdateStartingScore(vertices_by_position_[positions_[i]]);
      continue;
    }

    // If any of the vertices have no incoming edges we create zero-weighted
    // edges connecting to all vertices in the previous position. This is
    // needed so that we can consider a "broken" path.
    // Example:
    // ... ------- A -------- C   ------ G ------ ...
    // ... --------C        [ T ] ------ A ------ ...
    // This is a simplified example showing how a broken path may still
    // need to be considered. In this case we will create extra edges
    // connecting T with C.
    absl::btree_set<Edge> incoming_edges;
    if (i > 0 &&
        HasAtLeastOneIncomingEdge(vertices_by_position_[positions_[i]])) {
      for (const auto& v : vertices_by_position_[positions_[i]]) {
        auto [start, end] = boost::in_edges(v, graph_);

        // If there are no incoming edges for the vertex create zero weight
        // edges to all previous vertices to connect the graph.
        if (start == end && i > 0) {
          for (const auto& prev_v : vertices_by_position_[positions_[i - 1]]) {
            incoming_edges.insert(AddEdge(prev_v, v, 0));
          }
        }
        incoming_edges.insert(start, end);
      }
    } else {
      UpdateStartingScore(vertices_by_position_[positions_[i]]);
      continue;
    }

    absl::btree_map<std::pair<std::string, std::string>, Edge> keyed_edges;
    for (const auto& edge : incoming_edges) {
      std::string edge_source = graph_[edge.m_source].allele_info.bases;
      std::string edge_target = graph_[edge.m_target].allele_info.bases;
      keyed_edges[{edge_source, edge_target}] = edge;
    }

    // Enumerate all edge pairs
    for (const auto& edge_1 : keyed_edges) {
      for (const auto& edge_2 : keyed_edges) {
        const Vertex& to_1 = edge_1.second.m_target;
        const Vertex& to_2 = edge_2.second.m_target;
        Score score = CalculateScore(edge_1.second, edge_2.second);
        // If the score for the given vertices already exists then we update
        // it if the new score is higher.
        auto existing_score_it = scores_.find({to_1, to_2});
        // If the score for the given vertices does not exist then we create it.
        if (existing_score_it == scores_.end()) {
          scores_[{to_1, to_2}] = score;
        // If the score for the given vertices already exists then we update it.
        } else if (existing_score_it->second.score < score.score) {
          existing_score_it->second = score;
        // If the existing score is the same then we try to distinguish two
        // scores by comparing allele bases in "from" verteces.
        } else if (existing_score_it->second.score == score.score) {
          if (CompareVertexPairByBases(score.from[0], score.from[1],
                                      existing_score_it->second.from[0],
                                       existing_score_it->second.from[1])) {
            existing_score_it->second =  score;
          }
        }
      }  // for edge_2
    }    // for edge_1

    // If after assigning all score we find out that all scores are the same
    // we cannot phase this position. In that case the phasing is restarted from
    // the next position.
    if (AllScoresAreTheSame(keyed_edges)) {
      if (i < positions_.size() - 1) {
        i++;
        UpdateStartingScore(vertices_by_position_[positions_[i]]);
        continue;
      }
    }
  }  // for i in positions_
  // Backtrack from the last position. For each position where best partition is
  // not homozygous assign phases to vertices (alleles).
  AssignPhasesToVertices();

  // Phases are assigned to reads based on a set of alleles the read overlap.
  // If read overlaps more alleles of phase 1 then it is assigned a phase 1.
  // There 3 possible assignments: 0, 1, 2 where 0 is "phase unassigned".
  return AssignPhasesToReads(reads);
}

bool DirectPhasing::HasAtLeastOneIncomingEdge(
    const std::vector<Vertex>& vertecies) const {
  for (const auto& v : vertecies) {
    auto [start, end] = boost::in_edges(v, graph_);
    if (start != end) {
      return true;
      break;
    }
  }
  return false;
}

bool DirectPhasing::AllScoresAreTheSame(
    const absl::btree_map<std::pair<std::string, std::string>, Edge>&
        keyed_edges) const {
  int prev_score = -1;
  for (const auto& edge_1 : keyed_edges) {
    for (const auto& edge_2 : keyed_edges) {
      const Vertex& to_1 = edge_1.second.m_target;
      const Vertex& to_2 = edge_2.second.m_target;
      auto score_it = scores_.find({to_1, to_2});
      if (score_it != scores_.end()) {
        if (prev_score != -1 && score_it->second.score != prev_score) {
          return false;
        }
        prev_score = score_it->second.score;
      }
    }
  }
  return true;
}

bool DirectPhasing::CompareVertexPairByBases(
    const Vertex& v1_1, const Vertex& v1_2,
    const Vertex& v2_1, const Vertex& v2_2) const {
  if (v1_1 == nullptr || v1_2 == nullptr) {
    return false;
  }
  if (v2_1 == nullptr || v2_2 == nullptr) {
    return true;
  }
  return graph_[v1_1].allele_info.bases + graph_[v1_2].allele_info.bases >
      graph_[v2_1].allele_info.bases + graph_[v2_2].allele_info.bases;
}

absl::flat_hash_map<DirectPhasing::VertexPair,
                    DirectPhasing::Score>::const_iterator
DirectPhasing::MaxScore(int position) const {
  absl::flat_hash_map<VertexPair, Score>::const_iterator max_score_it =
      scores_.end();
  int max_score = 0;
  bool all_scores_equal = false;
  // Iterate all scores at positions_[position] and find the maximum.
  for (const Vertex& v1 : vertices_by_position_.at(positions_[position])) {
    for (const Vertex& v2 : vertices_by_position_.at(positions_[position])) {
      auto scores_it = scores_.find({v1, v2});
      if (scores_it == scores_.end()) {
        continue;
      }
      // TODO Add unit test for checking case where all scores are
      // equal for the candidate. This used to cause the non deterministic
      // behaviour and was fixed by adding allele bases comparison.
      if (scores_it->second.score > max_score) {
        max_score_it = scores_it;
        max_score = scores_it->second.score;
      } else if (scores_it->second.score == max_score) {
        // If scores are equal we will try to distinguish them by allele
        // bases
        if (max_score_it == scores_.end()
         || CompareVertexPairByBases(scores_it->first.phase_1_vertex,
                                     scores_it->first.phase_2_vertex,
                                     max_score_it->first.phase_1_vertex,
                                     max_score_it->first.phase_2_vertex)) {
          max_score_it = scores_it;
          max_score = scores_it->second.score;
        }
      }
    }
  }

  // If all the scores are the same at this position that means we couldn't
  // phase, move to the previous position.
  all_scores_equal = true;
  for (const Vertex& v1 : vertices_by_position_.at(positions_[position])) {
    for (const Vertex& v2 : vertices_by_position_.at(positions_[position])) {
      auto scores_it = scores_.find({v1, v2});
      if (scores_it != scores_.end()) {
        if (scores_it->second.score != max_score) {
          //  if (true) {
          all_scores_equal = false;
          break;
        }
      }
    }
    if (!all_scores_equal) break;
  }
  if (!all_scores_equal) {
    return max_score_it;
  } else {
    return scores_.end();
  }
}

void DirectPhasing::AssignPhasesToVertices() {
  // Assigning a random valid score. The max_score should be at least no less
  // then this.
  if (scores_.empty()) {
    return;
  }
  absl::flat_hash_map<VertexPair, Score>::const_iterator max_score_it =
      scores_.end();
  int i = positions_.size() - 1;
  // Going backward from the last positions find the first position that can
  // be phased.
  while (i >= 0) {
    // Skip position if max score cannot be found. In that case alleles at this
    // position can not be phased.
    while (i >= 0) {
      max_score_it = MaxScore(i);
      if (max_score_it == scores_.end()) {
        i--;
      } else {
        break;
      }
    }

    // Unwind from the first (from the last) position up and record phases.
    auto prev_score = max_score_it;
    while (max_score_it != scores_.end()) {
      if (max_score_it->first.phase_1_vertex !=
          max_score_it->first.phase_2_vertex) {
        graph_[max_score_it->first.phase_1_vertex].allele_info.phase = 1;
        graph_[max_score_it->first.phase_2_vertex].allele_info.phase = 2;
      }
      // Go to the next score.
      max_score_it = scores_.find(
          {max_score_it->second.from[0], max_score_it->second.from[1]});
      if (max_score_it == prev_score) {
        // Phasing cannot be continued. We have to break here.
        i--;
        break;
      }
      prev_score = max_score_it;
      i--;
    }
  }
}

std::vector<PhasedVariant> DirectPhasing::GetPhasedVariants() const {
  std::vector<PhasedVariant> phased_variants;
  for (int pos : positions_) {
    auto it = vertices_by_position_.find(pos);
    if (it == vertices_by_position_.end()) {
      continue;
    }

    std::array<std::string, 2> bases = {"", ""};
    for (const Vertex& v : it->second) {
      const auto& vertex = graph_[v];
      if (vertex.allele_info.phase == 1) {
        bases[0] = vertex.allele_info.bases;
      } else if (vertex.allele_info.phase == 2) {
        bases[1] = vertex.allele_info.bases;
      }
    }
    if (!bases[0].empty() && !bases[1].empty()) {
      phased_variants.push_back({.position = pos,
                                 .phase_1_bases = bases[0],
                                 .phase_2_bases = bases[1]});
    }
  }
  return phased_variants;
}

std::vector<int> DirectPhasing::AssignPhasesToReads(
    absl::Span<const nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
        reads) const {
  // Assign phase reads
  // 1. For each read find all allele the read overlaps.
  // 2. Assign the phase to the read based on the majority phase of all
  // overlapped alleles.
  // Each read is assigned a phase (1,2) or 0 if phase cannot be determined.
  std::vector<int> phases(reads.size(), 0);
  for (int i = 0; i < reads.size(); i++) {
    ReadIndex read_index = read_to_index_.at(ReadKey(*reads[i].p_));

    // Calculate the number of alleles of each phase the read overlaps.
    if (read_to_alleles_.contains(read_index)) {
      std::array<int, 3> read_phases = {0};

      for (auto allele_support : read_to_alleles_.at(read_index)) {
        const Vertex& v = allele_support.vertex;
        read_phases[graph_[v].allele_info.phase]++;
      }
      if (read_phases[1] > read_phases[2] &&
          read_phases[1] >= kMinAllelesToPhase) {
        phases[i] = 1;
      } else if (read_phases[2] > read_phases[1] &&
                 read_phases[2] >= kMinAllelesToPhase) {
        phases[i] = 2;
      } else {
        phases[i] = 0;
      }
    } else {
      phases[i] = 0;
    }
  }
  return phases;
}

void DirectPhasing::InitializeReadMaps(
    absl::Span<const nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
        reads) {
  size_t index = 0;
  for (const auto& read : reads) {
    read_to_index_[ReadKey(*read.p_)] = index;
    index++;
  }
}

// From <starting_score> we know the originating vertex. We need to find all the
// reads that support a connection between originating vertex in
// <starting_score> and a new <vertex>. In addition we count reads that start
// at <vertex>.
absl::flat_hash_set<ReadIndex> DirectPhasing::FindSupportingReads(
    const Vertex& vertex, const Score& starting_score, int phase) const {
  CHECK_GE(phase, 0);
  CHECK_LT(phase, kNumOfPhases);
  // Find all reads supporting <vertex> vertex
  absl::flat_hash_set<ReadIndex> reads;
  for (const ReadSupportInfo& rs : graph_[vertex].allele_info.read_support) {
    if (rs.is_first_allele ||
        starting_score.read_support[phase].contains(rs.read_index)) {
      reads.insert(rs.read_index);
    }
  }
  return reads;
}

DirectPhasing::Score DirectPhasing::CalculateScore(const Edge& edge1,
                                                   const Edge& edge2) const {
  Vertex from_vertices[2] = {edge1.m_source, edge2.m_source};
  Vertex to_vertices[2] = {edge1.m_target, edge2.m_target};

  // The function should not be called if preceding score does not exist.
  // TODO Replace with assert.
  if (!scores_.contains({from_vertices[0], from_vertices[1]})) {
    return Score();
  }

  // Getting a preceding score.
  const Score& prev_score = scores_.at({from_vertices[0], from_vertices[1]});

  // Get all reads that support a given path.
  absl::flat_hash_set<ReadIndex> supporting_reads_by_phase[kNumOfPhases];
  for (int phase = 0; phase < kNumOfPhases; phase++) {
    supporting_reads_by_phase[phase] =
        FindSupportingReads(to_vertices[phase], prev_score, phase);
  }

  absl::flat_hash_set<ReadIndex> all_reads;
  for (int phase = 0; phase < kNumOfPhases; phase++) {
    all_reads.insert(supporting_reads_by_phase[phase].begin(),
                     supporting_reads_by_phase[phase].end());
  }

  // New score is old score + number of all supporting reads.
  return Score{.score = static_cast<int>(prev_score.score + all_reads.size()),
               .from = {from_vertices[0], from_vertices[1]},
               .read_support = {supporting_reads_by_phase[0],
                                supporting_reads_by_phase[1]}};
}

void DirectPhasing::UpdateStartingScore(const std::vector<Vertex>& verts) {
  // Iterate all pairs of vertices.
  for (int i = 0; i < verts.size(); i++) {
    for (int j = i; j < verts.size(); j++) {
      const auto& v1 = verts[i];
      const auto& v2 = verts[j];
      absl::flat_hash_set<ReadIndex> cur1_support;
      for (auto rs : graph_[v1].allele_info.read_support) {
        cur1_support.insert(rs.read_index);
      }
      absl::flat_hash_set<ReadIndex> cur2_support;
      for (auto rs : graph_[v2].allele_info.read_support) {
        cur2_support.insert(rs.read_index);
      }
      // Score equals the total number of unique supporting reads. If candidate
      // is heterozygous then supporting reads are disjoint sets. If candidate
      // is homozygous then supporting reads are equal sets. With that in mind
      // we can optimzie the union of supporting reads with the following
      // expression.
      int score = (cur1_support == cur2_support)
                      ? cur1_support.size()
                      : cur1_support.size() + cur2_support.size();
      scores_[{v1, v2}] = Score{.score = score,
                                .from = {Vertex(), Vertex()},
                                .read_support = {cur1_support, cur2_support}};
    }
  }
}

std::vector<ReadSupportInfo> DirectPhasing::ReadSupportFromProto(
    const google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport>& read_support)
    const {
  std::vector<ReadSupportInfo> read_support_infos;
  read_support_infos.reserve(read_support.size());
  for (const auto& read_support_item : read_support) {
    auto it = read_to_index_.find(read_support_item.read_name());
    if (it != read_to_index_.end() && !read_support_item.is_low_quality()) {
      read_support_infos.push_back(ReadSupportInfo{
          .read_index = it->second,
          .is_low_quality = read_support_item.is_low_quality()});
    }
  }
  return read_support_infos;
}

DirectPhasing::Vertex DirectPhasing::AddVertex(
    int64_t position, AlleleType allele_type, absl::string_view bases,
    const google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport>& reads) {
  Vertex v = boost::add_vertex(
      VertexInfo{AlleleInfo{.type = allele_type,
                            .position = position,
                            .bases = std::string(bases),
                            .read_support = ReadSupportFromProto(reads)}},
      graph_);
  return v;
}

DirectPhasing::Edge DirectPhasing::AddEdge(const Vertex& in_vertex,
                                           const Vertex& out_vertex,
                                           float weight) {
  bool was_present;
  Edge edge;
  std::tie(edge, was_present) = boost::edge(in_vertex, out_vertex, graph_);
  if (!was_present) {
    std::tie(edge, std::ignore) =
        boost::add_edge(in_vertex, out_vertex, EdgeInfo{0}, graph_);
  }
  EdgeInfo& ei = graph_[edge];
  ei.weight += weight;
  return edge;
}

DirectPhasing::Edge DirectPhasing::AddEdge(const Vertex& in_vertex,
                                           bool is_low_quality_in,
                                           const Vertex& out_vertex,
                                           bool is_low_quality_out) {
  float edge_weight =
      (is_low_quality_in ? 0.25 : 0.5) + (is_low_quality_out ? 0.25 : 0.5);
  return AddEdge(in_vertex, out_vertex, edge_weight);
}

void DirectPhasing::UpdateReadToAllelesMap(const Vertex& v) {
  vertices_by_position_[graph_[v].allele_info.position].push_back(v);

  for (auto& read_support_info : graph_[v].allele_info.read_support) {
    bool is_first = (read_to_alleles_.find(read_support_info.read_index) ==
                     read_to_alleles_.end());
    read_support_info.is_first_allele = is_first;
    read_to_alleles_[read_support_info.read_index].push_back(
        AlleleSupport{.is_set = true,
                      .vertex = v,
                      .read_support = ReadSupportInfo{
                          .read_index = read_support_info.read_index,
                          .is_low_quality = read_support_info.is_low_quality,
                          .is_first_allele = is_first,
                      }});
  }
}

void DirectPhasing::AddMethylatedRefCandidate(
    const DeepVariantCall& candidate) {
  // Retrieve REF supporting reads.
  const google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport>& ref_reads =
      candidate.ref_support_ext().read_infos();

  google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport> methyl_reads;
  google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport> unmethyl_reads;

  // Extract MF from INFO fields.
  float ref_methylation = 0.0;
  if (candidate.variant().calls_size() > 0) {
    auto mf_it = candidate.variant().calls(0).info().find("MF");
    if (mf_it != candidate.variant().calls(0).info().end() &&
        !mf_it->second.values().empty()) {
      ref_methylation = mf_it->second.values(0).number_value();
    }
  }

  // Keep only potentially heterozygous methylation sites.
  if (ref_methylation <= kMinMethylationThreshold ||
      ref_methylation >= kMaxMethylationThreshold) {
    return;
  }

  // Separate methylated and unmethylated REF reads.
  for (const auto& read : ref_reads) {
    if (read.is_methylated()) {
      *methyl_reads.Add() = read;
    } else {
      *unmethyl_reads.Add() = read;
    }
  }

  std::vector<DeepVariantCall_ReadSupport> methyl_reads_vec(
      methyl_reads.begin(), methyl_reads.end());
  std::vector<DeepVariantCall_ReadSupport> unmethyl_reads_vec(
      unmethyl_reads.begin(), unmethyl_reads.end());

  std::sort(methyl_reads_vec.begin(), methyl_reads_vec.end(),
            [](const DeepVariantCall_ReadSupport& a,
               const DeepVariantCall_ReadSupport& b) {
              return a.read_name() < b.read_name();
            });

  std::sort(unmethyl_reads_vec.begin(), unmethyl_reads_vec.end(),
            [](const DeepVariantCall_ReadSupport& a,
               const DeepVariantCall_ReadSupport& b) {
              return a.read_name() < b.read_name();
            });

  // Add methylated REF reads.
  // Represent methylation state as "M" (methylated) and "U" (unmethylated),
  // analogous to how alleles are encoded.
  if (!methyl_reads.empty()) {
    UpdateReadToAllelesMap(AddVertex(candidate.variant().start(),
                                     AlleleType::REFERENCE, "M",
                                     methyl_reads));
  }

  // Add unmethylated REF reads.
  if (!unmethyl_reads.empty()) {
    UpdateReadToAllelesMap(AddVertex(candidate.variant().start(),
                                     AlleleType::REFERENCE, "U",
                                     unmethyl_reads));
  }
}

void DirectPhasing::AddCandidate(const DeepVariantCall& candidate) {
  // Rreference sites are included as candidates in
  // CallVariantPosition() only if they have methylated reads.
  // This check for reference site (i.e., the alternate
  // base is ".") helps identify positions where we need to perform phasing for
  // reference sites with methylation-specific signals separately.
  bool is_ref_site = candidate.variant().alternate_bases().size() == 1 &&
                     candidate.variant().alternate_bases(0) == ".";

  if (is_ref_site) {
    AddMethylatedRefCandidate(candidate);
    return;
  }

  // Add REF if it has read support.
  const google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport>& ref_reads =
      candidate.ref_support_ext().read_infos();
  // Add REF allele.
  if (ref_reads.size() >= kMinRefAlleleDepth) {
    UpdateReadToAllelesMap(AddVertex(candidate.variant().start(),
                                     AlleleType::REFERENCE, kRef, ref_reads));
  }

  // Add alt alleles.
  using AlleleSupportItem =
      std::pair<std::string, DeepVariantCall_SupportingReadsExt>;
  // We need alleles sorted in order to make the algorithm deterministic.
  // Without it alleles order (and therefore phase assignment)
  // is random, but phasing is still correct.
  std::vector<AlleleSupportItem> alleles(candidate.allele_support_ext().begin(),
                                         candidate.allele_support_ext().end());

  std::sort(
      alleles.begin(), alleles.end(),
      [](const AlleleSupportItem& allele1, const AlleleSupportItem& allele2) {
        return allele1.first < allele2.first;
      });
  for (const auto& [allele, read_support] : alleles) {
    UpdateReadToAllelesMap(AddVertex(candidate.variant().start(),
                                     AlleleTypeFromCandidate(allele, candidate),
                                     allele, read_support.read_infos()));
  }
}

// Filters out all homozygous candidates and candidates containing indels.
bool CandidateFilter(const DeepVariantCall& candidate, uint32_t* indel_end)  {
  // If there is only one allele and not enough support for the ref then
  // empirically we can consider this candidate homozygous.
  if (candidate.allele_support_ext().size() <= 1 &&
       candidate.ref_support_ext().read_infos_size() < kMinRefAlleleDepth) {
    return false;
  }
  // The test filters out all candidates containing indels.
  for (const auto& [allele, read_support] : candidate.allele_support_ext()) {
    // Allele must not be overlapped by an INDEL and allele has to be a SNP.
    if (candidate.variant().end() <= *indel_end || allele.size() !=
        candidate.variant().end() - candidate.variant().start()) {
      if (*indel_end < candidate.variant().end()) {
        *indel_end = candidate.variant().end();
      }
      return false;
    }
  }
  return true;
}

void DirectPhasing::Clear() {
  hom_positions_.clear();
  positions_.clear();
  vertices_by_position_.clear();
  vertex_index_map_.clear();
  scores_.clear();
  read_to_alleles_.clear();
  read_to_index_.clear();
  graph_.clear();
}

// Iterate through all candidates in the region. For each potentially
// heterozygous SNP candidate create a graph vertex corresponding to each
// allele. Candidate is heterozygous if there is a ref allele, or there are
// multiple distinctive alt alleles.
void DirectPhasing::Build(
    absl::Span<const DeepVariantCall> candidates,
    absl::Span<const nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
        reads) {
  Clear();
  InitializeReadMaps(reads);

  // Iterate all candidates and create graph nodes.
  // It is assumed that candidates are processed in the position order.
  uint32_t indel_end = 0;
  for (int i = 0; i < candidates.size(); i++) {
    const auto& candidate = candidates[i];
    if (i > 0) {
      CHECK_LT(candidates[i - 1].variant().start(),
               candidate.variant().start());
    }
    if (CandidateFilter(candidate, &indel_end)) {
      AddCandidate(candidate);
      // Keep an ordered vector of positions.
      positions_.push_back(candidate.variant().start());
    }
  }  // for candidates

  // Add edges. Edges are created only between consecutive positions.
  // read_to_vert contains a vector of alleles that the read supports. Alleles
  // are sorted by position.
  for (const auto& read_to_vert : read_to_alleles_) {
    bool is_first = true;
    AlleleSupport prev_allele_support;
    for (const auto& allele_support : read_to_vert.second) {
      if (is_first) {
        is_first = false;
        prev_allele_support = allele_support;
        continue;
      }
      CHECK(prev_allele_support.is_set);
      auto pos_it =
          std::find(positions_.begin(), positions_.end(),
                    graph_[allele_support.vertex].allele_info.position);
      int prev_pos = *(--pos_it);
      int prev_allele_pos =
          graph_[prev_allele_support.vertex].allele_info.position;
      if (pos_it == positions_.begin() || prev_pos == prev_allele_pos) {
        AddEdge(prev_allele_support.vertex,
                prev_allele_support.read_support.is_low_quality,
                allele_support.vertex,
                allele_support.read_support.is_low_quality);
      }
      prev_allele_support = allele_support;
    }
  }

  // TODO Control Pruning with parameter. It should be off for testing.
  // Also, investigate if it helps the algorithm.
  //  Prune();
  RebuildIndexMap();
}

void DirectPhasing::RebuildIndexMap() {
  RawVertexIndexMap table;
  VertexIterator vi, vend;
  std::tie(vi, vend) = boost::vertices(graph_);
  int index = 0;
  for (; vi != vend; ++vi) {
    table[*vi] = index;
    ++index;
  }
  vertex_index_map_ = table;
}

// Helper functions.
AlleleType AlleleTypeFromCandidate(std::string_view bases,
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
  return AlleleType::UNSPECIFIED;
}

int NumOfSubstitutionAlleles(const DeepVariantCall& candidate) {
  return std::count_if(
      candidate.allele_support_ext().begin(),
      candidate.allele_support_ext().end(),
      [candidate](
          std::pair<std::string, DeepVariantCall_SupportingReadsExt> it) {
        return (it.first != kUncalledAllele &&
                AlleleTypeFromCandidate(it.first, candidate) ==
                    AlleleType::SUBSTITUTION);
      });
}

int NumOfIndelAlleles(const DeepVariantCall& candidate) {
  return std::count_if(
      candidate.allele_support_ext().begin(),
      candidate.allele_support_ext().end(),
      [candidate](
          std::pair<std::string, DeepVariantCall_SupportingReadsExt> it) {
        return (it.first != kUncalledAllele &&
                (AlleleTypeFromCandidate(it.first, candidate) ==
                     AlleleType::DELETION ||
                 AlleleTypeFromCandidate(it.first, candidate) ==
                     AlleleType::INSERTION));
      });
}

int SubstitutionAllelesDepth(const DeepVariantCall& candidate) {
  int count = 0;
  for (const auto& allele_info_it : candidate.allele_support_ext()) {
    if (allele_info_it.first != kUncalledAllele &&
        AlleleTypeFromCandidate(allele_info_it.first, candidate) ==
            AlleleType::SUBSTITUTION) {
      // TODO Low quality reads are included here. To be verified.
      count += allele_info_it.second.read_infos_size();
    }
  }
  return count;
}

template <class BoostGraph>
class EdgeLabelWriter {
 public:
  explicit EdgeLabelWriter(const BoostGraph& g) : g_(g) {}

  void operator()(std::ostream& out, const DirectPhasing::Edge e) const {
    DirectPhasing::EdgeInfo ei = g_[e];
    out << "[label=" << ei.weight << "]";
  }

 private:
  const BoostGraph& g_;
};

template <class BoostGraph>
class VertexLabelWriter {
 public:
  explicit VertexLabelWriter(const BoostGraph& g) : g_(g) {}

  void operator()(std::ostream& out, const DirectPhasing::Vertex v) const {
    DirectPhasing::VertexInfo vi = g_[v];
    out << "[label=\"" << vi.allele_info.position
        << " " << vi.allele_info.bases << "\"]";
  }

 private:
  const BoostGraph& g_;
};

DirectPhasing::VertexIndexMap DirectPhasing::IndexMap() const {
  boost::const_associative_property_map<RawVertexIndexMap> vmap(
      vertex_index_map_);
  return vmap;
}

std::string DirectPhasing::GraphViz() const {
  std::stringstream graphviz;
  boost::write_graphviz(
      graphviz,
      graph_,
      VertexLabelWriter<BoostGraph>(graph_),
      EdgeLabelWriter<BoostGraph>(graph_),
      boost::default_writer(),
      IndexMap());
  return graphviz.str();
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
