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
#include <string>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/graph/graphviz.hpp"
#include "third_party/nucleus/protos/variants.pb.h"
#include "tensorflow/core/platform/logging.h"

namespace learning {
namespace genomics {
namespace deepvariant {

const int kMinRefAlleleDepth = 3;
const float kMinEdgeWeight = 2.0;
constexpr absl::string_view kRef = "REF";

absl::StatusOr<std::vector<int>> DirectPhasing::PhaseReads(
    const std::vector<DeepVariantCall>& candidates,
    const std::vector<
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>& reads) {
  // Not Implemented.
  return absl::Status(absl::StatusCode::kUnimplemented, "Not Implemented");
}

std::string ReadKey(const nucleus::genomics::v1::Read& read) {
  return absl::StrCat(read.fragment_name(), "/", read.read_number());
}

void DirectPhasing::InitializeReadMaps(
    const std::vector<
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>& reads) {
  size_t index = 0;
  for (const auto& read : reads) {
    read_to_index_[ReadKey(*read.p_)] = index;
    index++;
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

Vertex DirectPhasing::AddVertex(
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

Edge DirectPhasing::AddEdge(const std::pair<Vertex, bool>& in_vertex,
                            const std::pair<Vertex, bool>& out_vertex,
                            float weight) {
  bool was_present;
  Edge edge;
  std::tie(edge, was_present) =
      boost::edge(in_vertex.first, out_vertex.first, graph_);
  if (!was_present) {
    std::tie(edge, std::ignore) =
        boost::add_edge(in_vertex.first, out_vertex.first, EdgeInfo{0}, graph_);
  }
  EdgeInfo& ei = graph_[edge];
  ei.weight += weight;
  return edge;
}

Edge DirectPhasing::AddEdge(const std::pair<Vertex, bool>& in_vertex,
                            const std::pair<Vertex, bool>& out_vertex) {
  float edge_weight =
      (in_vertex.second ? 0.25 : 0.5) + (out_vertex.second ? 0.25 : 0.5);
  return AddEdge(in_vertex, out_vertex, edge_weight);
}

void DirectPhasing::UpdateReadToAllelesMap(const Vertex& v) {
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

void DirectPhasing::AddCandidate(const DeepVariantCall& candidate) {
  // Add REF if it has read support.
  const google::protobuf::RepeatedPtrField<DeepVariantCall_ReadSupport>& ref_reads =
      candidate.ref_support_ext().read_infos();
  // Add REF allele.
  if (ref_reads.size() >= kMinRefAlleleDepth) {
    UpdateReadToAllelesMap(AddVertex(candidate.variant().start(),
                                  AlleleType::REFERENCE, kRef, ref_reads));
  }

  // Add alt alleles.
  for (const auto& allele : candidate.allele_support_ext()) {
    UpdateReadToAllelesMap(
        AddVertex(candidate.variant().start(),
                  AlleleTypeFromCandidate(allele.first, candidate),
                  allele.first, allele.second.read_infos()));
  }
}

// Filter out all homozygious and low quality candidates. The graph should be
// built with only heteragious candidates.
// redacted
bool CandidateFilter(const DeepVariantCall& candidate) { return true; }

// Iterate through all candidates in the region. For each potentially
// heterozygious SNP candidate create a graph vertex coressponding to each
// allele. Candidate is heterozygous if there is a ref allele, or there are
// multiple distintive alt alleles.
void DirectPhasing::Build(
    const std::vector<DeepVariantCall>& candidates,
    const std::vector<
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>& reads) {
  InitializeReadMaps(reads);

  // Iterate all candidates and create graph nodes.
  for (const auto& candidate : candidates) {
    if (CandidateFilter(candidate)) {
      AddCandidate(candidate);
      // Keep an ordered vector of positions.
      if (positions_.empty() ||
          positions_.back() != candidate.variant().start()) {
        positions_.push_back(candidate.variant().start());
      }
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
        AddEdge(std::pair(prev_allele_support.vertex,
                          prev_allele_support.read_support.is_low_quality),
                std::pair(allele_support.vertex,
                          allele_support.read_support.is_low_quality));
      }
      prev_allele_support = allele_support;
    }
  }

  Prune();
  RebuildIndexMap();
}

void DirectPhasing::Prune() {
  // Remove low-weight edges.
  boost::remove_edge_if(
      [this](const Edge& e) { return graph_[e].weight < kMinEdgeWeight; },
      graph_);
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
