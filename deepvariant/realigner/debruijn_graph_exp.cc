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
#include <map>
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

using Vertex = learning::genomics::deepvariant::Vertex;
using VertexIndexMap = DeBruijnGraphExp::VertexIndexMap;
using Edge = learning::genomics::deepvariant::Edge;
using Path = learning::genomics::deepvariant::Path;

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

  template <class Edge, class Graph>
  void back_edge(Edge, const Graph&) {
    *has_cycle = true;
  }

 private:
  bool* has_cycle;
};

template <class BoostGraph>
class EdgeLabelWriter {
 public:
  explicit EdgeLabelWriter(const BoostGraph& g) : g_(g) {}

  void operator()(std::ostream& out, const Edge e) const {
    EdgeInfo ei = g_[e];
    out << "[label=" << std::to_string(ei.weight)
        << (ei.is_ref ? " color=red" : "") << "]";
  }

 private:
  const BoostGraph& g_;
};

class ReachableVertexVisitor : public boost::dfs_visitor<> {
 public:
  explicit ReachableVertexVisitor(std::set<Vertex>* reachable_vertices)
      : reachable_vertices(reachable_vertices) {}

  template <class Edge, class Graph>
  void tree_edge(Edge e, const Graph& g) {
    Vertex from = boost::source(e, g);
    if (reachable_vertices->find(from) != reachable_vertices->end()) {
      Vertex to = boost::target(e, g);
      reachable_vertices->insert(to);
    }
  }

 private:
  std::set<Vertex>* reachable_vertices;
};

template <class BoostGraphT, class VertexIndexMapT>
std::set<Vertex> VerticesReachableFrom(
    Vertex v, const BoostGraphT& g, const VertexIndexMapT& vertex_index_map) {
  std::set<Vertex> reachable_vertices{v};
  ReachableVertexVisitor vis(&reachable_vertices);
  boost::depth_first_search(
      g, boost::visitor(vis).root_vertex(v).vertex_index_map(vertex_index_map));
  return reachable_vertices;
}

}  // namespace

Vertex DeBruijnGraphExp::EnsureVertex(string_view kmer) {
  Vertex v;
  auto vertex_find = kmer_to_vertex_.find(kmer);
  if (vertex_find != kmer_to_vertex_.end()) {
    v = (*vertex_find).second;
  } else {
    string kmer_copy(kmer);
    v = boost::add_vertex(VertexInfo{kmer_copy, 0}, g_);
    // N.B.: must use the long-lived string in the map key as the referent of
    // the string_view key.
    kmer_to_vertex_[g_[v].kmer] = v;
  }
  g_[v].frequency++;

  return v;
}

Vertex DeBruijnGraphExp::VertexForKmer(string_view kmer) const {
  return kmer_to_vertex_.at(kmer);
}

void DeBruijnGraphExp::RebuildIndexMap() {
  std::map<Vertex, int> table;
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


KBounds KMinMaxFromReference(const string_view ref,
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
  KBounds bounds = KMinMaxFromReference(ref, options);

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
        new DeBruijnGraphExp(ref, reads, options, 43));

    last_graph->PruneLite();
    return last_graph;
  }

  return nullptr;
}

Edge DeBruijnGraphExp::AddEdge(
    Vertex from_vertex, Vertex to_vertex, bool is_ref) {
  bool was_present;
  Edge edge;
  std::tie(edge, was_present) = boost::edge(from_vertex, to_vertex, g_);
  if (!was_present) {
    std::tie(edge, std::ignore) = boost::add_edge(from_vertex, to_vertex,
                                                  EdgeInfo{0, false}, g_);
  }
  EdgeInfo& ei = g_[edge];
  ei.weight++;
  ei.is_ref |= is_ref;
  return edge;
}

void DeBruijnGraphExp::AddKmersAndEdges(string_view bases, int start, int end,
                                     bool is_ref, bool is_debug) {
  CHECK_GE(start, 0);
  CHECK_LE(start + k_, bases.size());
  CHECK_LE(end + k_, bases.size());

  // End can be less than 0, in which case we return without doing any work.
  if (end > 0) {
    Vertex vertex_prev = EnsureVertex(bases.substr(start, k_));
    for (int i = start + 1; i <= end; ++i) {
      Vertex vertex_cur = EnsureVertex(bases.substr(i, k_));
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

  bool is_debug = false;
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
    AddKmersAndEdges(bases_view, i, next_bad_position - k_, false, is_debug);
    i = next_bad_position + 1;
  }
}

void DeBruijnGraphExp::CandidatePathsRankedHelper(
    Vertex u, const absl::flat_hash_set<Vertex>& sink_nodes, Path& current_path,
    std::priority_queue<Path>& pq, int& num_paths) const {
  int visit_count = 0;
  for (Vertex v : current_path.path) {
    if (v == u) {
      visit_count++;
    }
  }

  if (sink_nodes.contains(u) || boost::out_degree(u, g_) == 0
      || visit_count > 1) {
    num_paths++;
    pq.push(current_path);

    if (pq.size() > kMaxNumPathsToKeep) {
      pq.pop();
    }
  }

  if (num_paths > kMaxNumPaths) return;

  double total_edge_frequency = 0;
  int out_degree = 0;
  AdjacencyIterator vi, vend;
  std::tie(vi, vend) = boost::adjacent_vertices(u, g_);
  for (auto it = vi; it != vend; ++it) {
    total_edge_frequency += g_[*it].frequency;
    out_degree++;
  }

  double current_path_score = current_path.score;
  for (; vi != vend; ++vi) {
    Vertex v = *vi;
    if (out_degree > 1) {
      current_path.score = current_path_score + std::log10(g_[v].frequency) -
                           std::log10(total_edge_frequency);
    }
    current_path.path.push_back(v);
    CandidatePathsRankedHelper(v, sink_nodes, current_path, pq, num_paths);
    current_path.path.pop_back();
    if (num_paths > kMaxNumPaths) {
      LOG(INFO) << "Num paths > " << kMaxNumPaths << ", cutting off";
      return;
    }
  }
  current_path.score = current_path_score;
}

std::vector<Path> DeBruijnGraphExp::CandidatePathsRanked() const {
  std::priority_queue<Path> pq;
  int num_paths = 0;
  std::vector<Vertex> start_nodes = StartNodes();
  std::vector<Vertex> sinks = SinkNodes();
  absl::flat_hash_set<Vertex> sink_nodes(sinks.begin(), sinks.end());

  for (Vertex source : start_nodes) {
    Path current_path = {{source}, 0.0};
    CandidatePathsRankedHelper(source, sink_nodes, current_path, pq, num_paths);
    if (num_paths > kMaxNumPaths) break;
  }

  std::vector<Path> sorted_paths;
  while (!pq.empty()) {
    sorted_paths.push_back(pq.top());
    pq.pop();
  }
  std::sort(sorted_paths.begin(), sorted_paths.end(),
            [](const Path& a, const Path& b) { return a.score > b.score; });
  return sorted_paths;
}

string DeBruijnGraphExp::HaplotypeForPath(const Path& path) const {
  std::stringstream haplotype;
  for (size_t i = 0; i < path.path.size(); ++i) {
    Vertex v = path.path[i];
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
  for (const Path& path : CandidatePathsRanked()) {
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
      boost::get(&VertexInfo::kmer, g_));
  boost::write_graphviz(
      graphviz,
      g_,
      vertex_label_writer,
      EdgeLabelWriter<BoostGraph>(g_),
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
      Vertex v = *vi;
      if (boost::in_degree(v, g_) != 1) continue;

      Edge in_edge = *boost::in_edges(v, g_).first;
      Vertex u = boost::source(in_edge, g_);
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
      std::vector<std::pair<Vertex, EdgeInfo>> out_edges;
      AdjacencyIterator ai, aend;
      std::tie(ai, aend) = boost::adjacent_vertices(v, g_);
      for (; ai != aend; ++ai) {
        Edge e = boost::edge(v, *ai, g_).first;
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

std::vector<Vertex> DeBruijnGraphExp::StartNodes() const {
  std::vector<Vertex> start_nodes;
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  for (; vbegin != vend; ++vbegin) {
    Vertex v = *vbegin;
    if (boost::in_degree(v, g_) == 0) {
      start_nodes.push_back(v);
    }
  }
  return start_nodes;
}

std::vector<Vertex> DeBruijnGraphExp::SinkNodes() const {
  std::vector<Vertex> sink_nodes;
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  for (; vbegin != vend; ++vbegin) {
    Vertex v = *vbegin;
    if (boost::in_degree(v, g_) > 0 && boost::out_degree(v, g_) == 0) {
      sink_nodes.push_back(v);
    }
  }
  return sink_nodes;
}

void DeBruijnGraphExp::PruneLite() {
  // Remove all edges with weight < 2.
  boost::remove_edge_if([this](const Edge& e) { return g_[e].weight < 2; }, g_);

  // Remove vertices that have zero incoming and zero outgoing edges.
  std::vector<Vertex> to_remove;
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  // We create a copy of the vertices because boost::remove_vertex invalidates
  // vertex descriptors and iterators.
  std::vector<Vertex> vertices(vbegin, vend);
  for (Vertex v : vertices) {
    if (boost::in_degree(v, g_) == 0 && boost::out_degree(v, g_) == 0) {
      to_remove.push_back(v);
    }
  }

  for (Vertex v : to_remove) {
    kmer_to_vertex_.erase(g_[v].kmer);
    boost::clear_vertex(v, g_);
    boost::remove_vertex(v, g_);
  }
  RebuildIndexMap();
}


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
