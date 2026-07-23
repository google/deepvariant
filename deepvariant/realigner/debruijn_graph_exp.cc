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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/node_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/ascii.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/reverse_graph.hpp"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

constexpr int kMaxNumPaths = 50000;
constexpr int kMaxNumPathsToKeep = 64;
constexpr int kDefaultKmerSize = 43;

using VertexIndexMap = DeBruijnGraphExp::VertexIndexMap;
using Read = nucleus::genomics::v1::Read;
using absl::string_view;

namespace {

// Visitor classes we will use to run boost algorithms.  N.B.: these classes
// operate by side effect, modifying pointers that are passed in.  This is not
// an optimal design, rather it is to work around the copying done by
// boost::visitor which makes it difficult for our code to retain a pointer to
// the visitor object that is actually used.

class CycleDetector : public boost::dfs_visitor<> {
 public:
  explicit CycleDetector(bool* has_cycle) : has_cycle(has_cycle) {}

  template <class EdgeExp, class Graph>
  void back_edge(EdgeExp, const Graph&) {
    *has_cycle = true;
  }

 private:
  bool* has_cycle;
};

template <class BoostGraph>
class EdgeLabelWriter {
 public:
  explicit EdgeLabelWriter(const BoostGraph& g) : g_(g) {}

  void operator()(std::ostream& out, const EdgeExp e) const {
    EdgeInfoExp ei = g_[e];
    out << "[label=" << std::to_string(ei.weight)
        << (ei.is_ref ? " color=red" : "") << "]";
  }

 private:
  const BoostGraph& g_;
};

class ReachableVertexVisitor : public boost::dfs_visitor<> {
 public:
  explicit ReachableVertexVisitor(std::set<VertexExp>* reachable_vertices)
      : reachable_vertices(reachable_vertices) {}

  template <class EdgeExp, class Graph>
  void tree_edge(EdgeExp e, const Graph& g) {
    VertexExp from = boost::source(e, g);
    if (reachable_vertices->find(from) != reachable_vertices->end()) {
      VertexExp to = boost::target(e, g);
      reachable_vertices->insert(to);
    }
  }

 private:
  std::set<VertexExp>* reachable_vertices;
};

}  // namespace

VertexExp DeBruijnGraphExp::EnsureVertex(string_view kmer) {
  VertexExp v;
  auto vertex_find = kmer_to_vertex_.find(kmer);
  if (vertex_find != kmer_to_vertex_.end()) {
    v = (*vertex_find).second;
  } else {
    string kmer_copy(kmer);
    v = boost::add_vertex(VertexInfoExp{kmer_copy, 0}, g_);
    // N.B.: must use the long-lived string in the map key as the referent of
    // the string_view key.
    kmer_to_vertex_[g_[v].kmer] = v;
  }
  g_[v].frequency++;

  return v;
}

VertexExp DeBruijnGraphExp::VertexForKmer(string_view kmer) const {
  return kmer_to_vertex_.at(kmer);
}

void DeBruijnGraphExp::RebuildIndexMap() {
  absl::flat_hash_map<VertexExp, int> table;
  VertexIterator vi, vend;
  std::tie(vi, vend) = boost::vertices(g_);
  int index = 0;
  for (; vi != vend; ++vi) {
    table[*vi] = index;
    ++index;
  }
  vertex_index_map_ = table;
}

VertexIndexMap DeBruijnGraphExp::IndexMap() const {
  boost::const_associative_property_map<RawVertexIndexMap> vmap(
      vertex_index_map_);
  return vmap;
}

bool DeBruijnGraphExp::HasCycle() const {
  bool has_cycle = false;
  CycleDetector cycle_detector(&has_cycle);
  boost::depth_first_search(
      g_, boost::visitor(cycle_detector).vertex_index_map(IndexMap()));
  return has_cycle;
}

DeBruijnGraphExp::DeBruijnGraphExp(
    absl::string_view ref,
    absl::Span<const nucleus::ConstProtoPtr<const Read>> reads,
    const Options& options, int k)
    : options_(options), k_(k) {
  CHECK_GT(k, 0);  // k should always be a positive integer.
  CHECK(static_cast<uint32_t>(k) < ref.size());
  // TODO: We may want to try adding the reference sequence to the graph.
  // Even though we don't use reference to assign source and sink nodes
  // reference k-mers can still help to bridge gaps between disjoint subtgraphs.
  for (const nucleus::ConstProtoPtr<const Read>& read_ptr : reads) {
    const Read& read = *read_ptr.p_;
    if (read.alignment().mapping_quality() >= options.min_mapq()) {
      AddEdgesForRead(read);
    }
  }
  RebuildIndexMap();
}

// Indicates that we couldn't find a minimum k that can be used.
constexpr int kBoundsNoWorkingK = -1;
struct KBounds {
  int min_k;  // Minimum k to consider (inclusive).
  int max_k;  // Maximum k to consider (inclusive).
};


KBounds KMinMaxFromReferenceExp(const string_view ref,
                             const DeBruijnGraphExp::Options& options) {
  KBounds bounds;
  bounds.min_k = options.min_k();
  bounds.max_k = std::min(options.max_k(), static_cast<int>(ref.size()) - 1);

  for (int k = options.min_k(); k <= bounds.max_k; k += options.step_k()) {
    bool has_cycle = false;
    absl::btree_set<string_view> kmers;

    for (int i = 0; i < ref.size() - k + 1; i++) {
      string_view kmer = ref.substr(i, k);
      if (kmers.insert(kmer).second == false) {
        // No insertion took place because the kmer already exists. This implies
        // that there's a cycle in the graph.
        has_cycle = true;
        break;
      }
    }

    if (!has_cycle) {
      bounds.min_k = k;
      break;
    }
  }

  return bounds;
}

std::unique_ptr<DeBruijnGraphExp> DeBruijnGraphExp::Build(
    absl::string_view ref,
    const std::vector<nucleus::ConstProtoPtr<const Read>>& reads,
    const DeBruijnGraphExp::Options& options) {
  KBounds bounds = KMinMaxFromReferenceExp(ref, options);

  std::unique_ptr<DeBruijnGraphExp> last_graph;
  for (int k = bounds.min_k; k <= bounds.max_k; k += options.step_k()) {
    last_graph = std::unique_ptr<DeBruijnGraphExp>(
        new DeBruijnGraphExp(ref, reads, options, k));
    if (!last_graph->HasCycle()) {
        last_graph->PruneLite();
      return last_graph;
    }
  }

  if (last_graph) {
    last_graph = std::unique_ptr<DeBruijnGraphExp>(
        new DeBruijnGraphExp(ref, reads, options, kDefaultKmerSize));

    last_graph->PruneLite();
    return last_graph;
  }

  return nullptr;
}

EdgeExp DeBruijnGraphExp::AddEdge(
    VertexExp from_vertex, VertexExp to_vertex, bool is_ref) {
  bool was_present;
  EdgeExp edge;
  std::tie(edge, was_present) = boost::edge(from_vertex, to_vertex, g_);
  if (!was_present) {
    std::tie(edge, std::ignore) = boost::add_edge(from_vertex, to_vertex,
                                                  EdgeInfoExp{0, false}, g_);
  }
  EdgeInfoExp& ei = g_[edge];
  ei.weight++;
  ei.is_ref |= is_ref;
  return edge;
}

void DeBruijnGraphExp::AddKmersAndEdges(string_view bases, int start, int end,
                                     bool is_ref) {
  CHECK_GE(start, 0);
  CHECK_LE(start + k_, bases.size());
  CHECK_LE(end + k_, bases.size());

  // End can be less than 0, in which case we return without doing any work.
  if (end > 0) {
    VertexExp vertex_prev = EnsureVertex(bases.substr(start, k_));
    for (int i = start + 1; i <= end; ++i) {
      VertexExp vertex_cur = EnsureVertex(bases.substr(i, k_));
      AddEdge(vertex_prev, vertex_cur, is_ref);
      vertex_prev = vertex_cur;
    }
  }
}

void DeBruijnGraphExp::AddEdgesForReference(string_view ref) {
  AddKmersAndEdges(ref, 0, ref.size() - k_, true /* is_ref */);
}


void DeBruijnGraphExp::AddEdgesForRead(
    const nucleus::genomics::v1::Read& read) {
  const string bases = absl::AsciiStrToUpper(read.aligned_sequence());

  // Lambda function to find the next bad position in the read, if one exists,
  // starting from offset `start` in the read. If all remains bases/quals are
  // good, returns bases.size().
  auto NextBadPosition = [&read, &bases, this](int start) -> int {
    for (int i = start; i < bases.size(); ++i) {
      if (!IsCanonicalBase(bases[i], nucleus::CanonicalBases::ACGT) ||
          (int)read.aligned_quality()[i] < options_.min_base_quality()) {
        return i;
      }
    }
    return bases.size();
  };

  // TODO: This comment is outdated.
  // This algorithm is simple and fast, but it isn't the most straightforward
  // implementation so it merits a few comments.
  //
  // Suppose I have the following data:
  //
  // offset: 01234567
  // bases:  ACGTAACC
  // bad? :  00010000
  // k_   :  2 <= using a kmer size of 2
  //
  // The algorithm below loops over positions (variable `i`), pulling kmers of
  // length k from positions `i` and `i + 1` to add as edges. The key
  // calculation is NextBadPosition that searches from the current `i` position
  // for the next position that is bad. In the above example, this would be the
  // 3 position. We then loop from i until `next_bad_position - k`, to create
  // our edges, since we know that everything from i to next_bad_position is
  // good but we cannot construct a valid kmer that overlaps next_bad_position
  // so it invalidates all kmer starts from `next_bad_position - k`. Finally, we
  // set i to `next_bad_position + 1`, which is the very next starting position
  // after the last bad position, and the algorithm repeats.
  //
  // This algorithm has many important properties for performance:
  //
  //   * It doesn't allocate any data structures to support the calculation.
  //   * It only examines whether a given position is good/bad once.
  //   * The loop to add edges is streamlined, without any unnecessary checks.
  //
  const string_view bases_view(bases);
  // Note that this SIGNED int type declaration is key to avoid
  // bases.size() - k_ underflowing.
  const int stop = bases.size() - k_;
  int i = 0;
  while (i < stop) {
    int next_bad_position = NextBadPosition(i);
    AddKmersAndEdges(bases_view, i, next_bad_position - k_, false);
    i = next_bad_position + 1;
  }
}

void DeBruijnGraphExp::CandidatePathsRankedHelper(
    VertexExp u,
    const absl::flat_hash_set<VertexExp>& sink_nodes,
    PathExp& current_path,
    std::priority_queue<PathExp>& pq, int& num_paths) const {
  int visit_count = 0;
  for (VertexExp v : current_path.path) {
    if (v == u) {
      visit_count++;
    }
  }

  bool is_cycle = (visit_count > 1);
  if (sink_nodes.contains(u) || boost::out_degree(u, g_) == 0
      || is_cycle) {
    num_paths++;
    pq.push(current_path);

    if (pq.size() > kMaxNumPathsToKeep) {
      pq.pop();
    }

    // If it was a cycle, stop recursing further down this path.
    if (is_cycle) {
      return;
    }
  }

  if (num_paths > kMaxNumPaths) return;

  double total_edge_frequency = 0;
  int out_degree = 0;
  AdjacencyIterator vi, vend;
  std::tie(vi, vend) = boost::adjacent_vertices(u, g_);
  std::vector<VertexExp> successors(vi, vend);
  if (successors.size() > 1) {
    // Successors share a (k-1)-prefix. Comparing only kmer[k_-1] suffices
    // to order them deterministically, even after Collapse() extends kmer
    // strings beyond length k_.
    std::sort(successors.begin(), successors.end(),
              [this](VertexExp v1, VertexExp v2) {
                return g_[v1].kmer[k_ - 1] < g_[v2].kmer[k_ - 1];
              });
  }

  for (VertexExp v : successors) {
    total_edge_frequency += g_[v].frequency;
    out_degree++;
  }

  double current_path_score = current_path.score;
  for (VertexExp v : successors) {
    if (out_degree > 1) {
      current_path.score = current_path_score + std::log10(g_[v].frequency) -
                           std::log10(total_edge_frequency);
    }
    current_path.path.push_back(v);
    CandidatePathsRankedHelper(v, sink_nodes, current_path, pq, num_paths);
    current_path.path.pop_back();
    if (num_paths > kMaxNumPaths) {
      // In a complex graph number of paths can explode exponentially. If we
      // went over the kMaxNumPaths we stop the search early.
      return;
    }
  }
  current_path.score = current_path_score;
}

std::vector<PathExp> DeBruijnGraphExp::CandidatePathsRanked() const {
  std::priority_queue<PathExp> pq;
  int num_paths = 0;
  std::vector<VertexExp> start_nodes = StartNodes();
  std::sort(
      start_nodes.begin(), start_nodes.end(),
      [this](VertexExp v1, VertexExp v2) { return g_[v1].kmer < g_[v2].kmer; });
  std::vector<VertexExp> sinks = SinkNodes();
  absl::flat_hash_set<VertexExp> sink_nodes(sinks.begin(), sinks.end());

  for (VertexExp source : start_nodes) {
    PathExp current_path = {{source}, 0.0};
    CandidatePathsRankedHelper(source, sink_nodes, current_path, pq, num_paths);
    if (num_paths > kMaxNumPaths) {
      LOG(INFO) << "Num paths > " << kMaxNumPaths << ", cutting off";
      break;
    }
  }

  // Deterministic traversal order (sorted start_nodes and successors) ensures
  // paths are inserted into the priority queue in a deterministic order.
  // std::stable_sort preserves that order among equal-scored paths.
  std::vector<PathExp> sorted_paths;
  sorted_paths.reserve(pq.size());
  while (!pq.empty()) {
    sorted_paths.push_back(pq.top());
    pq.pop();
  }
  std::stable_sort(
      sorted_paths.begin(), sorted_paths.end(),
      [](const PathExp& a, const PathExp& b) { return a.score > b.score; });
  return sorted_paths;
}

string DeBruijnGraphExp::HaplotypeForPath(const PathExp& path) const {
  std::stringstream haplotype;
  for (size_t i = 0; i < path.path.size(); ++i) {
    VertexExp v = path.path[i];
    if (i < path.path.size() - 1) {
      haplotype << g_[v].kmer.substr(0, g_[v].kmer.size() - (k_ - 1));
    } else {
      haplotype << g_[v].kmer;
    }
  }
  return haplotype.str();
}

std::vector<std::string> DeBruijnGraphExp::CandidateHaplotypesRanked(
    int min_haplotype_len) const {
  std::vector<std::string> haplotypes;
  for (const PathExp& path : CandidatePathsRanked()) {
    std::string haplotype = HaplotypeForPath(path);
    if (haplotype.size() >= min_haplotype_len) {
      haplotypes.push_back(haplotype);
    }
  }
  return haplotypes;
}

string DeBruijnGraphExp::GraphViz() const {
  std::stringstream graphviz;
  auto vertex_label_writer = boost::make_label_writer(
      boost::get(&VertexInfoExp::kmer, g_));
  boost::write_graphviz(
      graphviz,
      g_,
      vertex_label_writer,
      EdgeLabelWriter<BoostGraphExp>(g_),
      boost::default_writer(),
      IndexMap());
  return graphviz.str();
}

void DeBruijnGraphExp::Collapse() {
  bool changed = true;
  while (changed) {
    changed = false;
    VertexIterator vi, vend;
    std::tie(vi, vend) = boost::vertices(g_);
    for (; vi != vend; ++vi) {
      VertexExp v = *vi;
      if (boost::in_degree(v, g_) != 1) continue;

      EdgeExp in_edge = *boost::in_edges(v, g_).first;
      VertexExp u = boost::source(in_edge, g_);
      if (boost::out_degree(u, g_) != 1) continue;

      // Merge node v into u: node u's sequence will be extended to cover u->v,
      // node v will be removed, and out-edges of v will become out-edges of u.
      // Remove old sequence of u from map.
      kmer_to_vertex_.erase(g_[u].kmer);
      // Extend u's sequence by appending non-overlapping suffix of v's
      // sequence and update map.
      g_[u].kmer += g_[v].kmer.substr(k_ - 1);
      kmer_to_vertex_[g_[u].kmer] = u;

      // Collect out-edges of v.
      std::vector<std::pair<VertexExp, EdgeInfoExp>> out_edges;
      AdjacencyIterator ai, aend;
      std::tie(ai, aend) = boost::adjacent_vertices(v, g_);
      for (; ai != aend; ++ai) {
        EdgeExp e = boost::edge(v, *ai, g_).first;
        out_edges.push_back({*ai, g_[e]});
      }
      // Add edges from u to successors of v.
      for (const auto& oe : out_edges) {
        boost::add_edge(u, oe.first, oe.second, g_);
      }

      // Remove v from graph.
      kmer_to_vertex_.erase(g_[v].kmer);
      boost::clear_vertex(v, g_);
      boost::remove_vertex(v, g_);
      changed = true;
      break;
    }
  }
  RebuildIndexMap();
}

std::vector<VertexExp> DeBruijnGraphExp::StartNodes() const {
  std::vector<VertexExp> start_nodes;
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  for (; vbegin != vend; ++vbegin) {
    VertexExp v = *vbegin;
    if (boost::in_degree(v, g_) == 0) {
      start_nodes.push_back(v);
    }
  }
  return start_nodes;
}

std::vector<VertexExp> DeBruijnGraphExp::SinkNodes() const {
  std::vector<VertexExp> sink_nodes;
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  for (; vbegin != vend; ++vbegin) {
    VertexExp v = *vbegin;
    if (boost::in_degree(v, g_) > 0 && boost::out_degree(v, g_) == 0) {
      sink_nodes.push_back(v);
    }
  }
  return sink_nodes;
}

void DeBruijnGraphExp::PruneLite() {
  // Remove all edges with weight < 2.
  boost::remove_edge_if([this](const EdgeExp& e) {
     return g_[e].weight < 2; }, g_);

  // Remove vertices that have zero incoming and zero outgoing edges.
  std::vector<VertexExp> to_remove;
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  // We create a copy of the vertices because boost::remove_vertex invalidates
  // vertex descriptors and iterators.
  std::vector<VertexExp> vertices(vbegin, vend);
  for (VertexExp v : vertices) {
    if (boost::in_degree(v, g_) == 0 && boost::out_degree(v, g_) == 0) {
      to_remove.push_back(v);
    }
  }

  for (VertexExp v : to_remove) {
    kmer_to_vertex_.erase(g_[v].kmer);
    boost::clear_vertex(v, g_);
    boost::remove_vertex(v, g_);
  }
  RebuildIndexMap();
}


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
