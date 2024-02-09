/*
 * Copyright 2017 Google LLC.
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

#include "deepvariant/realigner/debruijn_graph.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/node_hash_set.h"
#include "absl/log/check.h"
#include "absl/strings/ascii.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/reverse_graph.hpp"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Vertex = DeBruijnGraph::Vertex;
using VertexIndexMap = DeBruijnGraph::VertexIndexMap;
using Edge = DeBruijnGraph::Edge;
using Path = DeBruijnGraph::Path;

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

Vertex DeBruijnGraph::EnsureVertex(string_view kmer) {
  Vertex v;
  auto vertex_find = kmer_to_vertex_.find(kmer);
  if (vertex_find != kmer_to_vertex_.end()) {
    v = (*vertex_find).second;
  } else {
    string kmer_copy(kmer);
    v = boost::add_vertex(VertexInfo{kmer_copy}, g_);
    // N.B.: must use the long-lived string in the map key as the referent of
    // the string_view key.
    kmer_to_vertex_[g_[v].kmer] = v;
  }

  return v;
}

Vertex DeBruijnGraph::VertexForKmer(string_view kmer) const {
  return kmer_to_vertex_.at(kmer);
}

void DeBruijnGraph::RebuildIndexMap() {
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

VertexIndexMap DeBruijnGraph::IndexMap() const {
  boost::const_associative_property_map<RawVertexIndexMap> vmap(
      vertex_index_map_);
  return vmap;
}

bool DeBruijnGraph::HasCycle() const {
  bool has_cycle = false;
  CycleDetector cycle_detector(&has_cycle);
  boost::depth_first_search(
      g_, boost::visitor(cycle_detector).vertex_index_map(IndexMap()));
  return has_cycle;
}

DeBruijnGraph::DeBruijnGraph(
    absl::string_view ref,
    absl::Span<const nucleus::ConstProtoPtr<const Read>> reads,
    const Options& options, int k)
    : options_(options), k_(k) {
  CHECK_GT(k, 0);  // k should always be a positive integer.
  CHECK(static_cast<uint32_t>(k) < ref.size());
  AddEdgesForReference(ref);
  source_ = VertexForKmer(ref.substr(0, k_));
  sink_ = VertexForKmer(ref.substr(ref.size() - k_, k_));
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
                             const DeBruijnGraph::Options& options) {
  KBounds bounds;
  bounds.min_k = kBoundsNoWorkingK;
  bounds.max_k  = std::min(options.max_k(), static_cast<int>(ref.size()) - 1);

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

std::unique_ptr<DeBruijnGraph> DeBruijnGraph::Build(
    const string& ref,
    const std::vector<nucleus::ConstProtoPtr<const Read>>& reads,
    const DeBruijnGraph::Options& options) {
  KBounds bounds = KMinMaxFromReference(ref, options);
  if (bounds.min_k == kBoundsNoWorkingK) return nullptr;

  for (int k = bounds.min_k; k <= bounds.max_k; k += options.step_k()) {
    std::unique_ptr<DeBruijnGraph> graph = std::unique_ptr<DeBruijnGraph>(
        new DeBruijnGraph(ref, reads, options, k));
    if (graph->HasCycle()) {
      continue;
    } else {
      graph->Prune();
      return graph;
    }
  }
  return nullptr;
}

Edge DeBruijnGraph::AddEdge(Vertex from_vertex, Vertex to_vertex, bool is_ref) {
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

void DeBruijnGraph::AddKmersAndEdges(string_view bases, int start, int end,
                                     bool is_ref) {
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

void DeBruijnGraph::AddEdgesForReference(string_view ref) {
  AddKmersAndEdges(ref, 0, ref.size() - k_, true /* is_ref */);
}


void DeBruijnGraph::AddEdgesForRead(const nucleus::genomics::v1::Read& read) {
  const string bases = absl::AsciiStrToUpper(read.aligned_sequence());

  // Lambda function to find the next bad position in the read, if one exists,
  // starting from offset `start` in the read. If all remains bases/quals are
  // good, returns bases.size().
  auto NextBadPosition = [&read, &bases, this](int start) -> int {
    for (int i = start; i < bases.size(); ++i) {
      if (!IsCanonicalBase(bases[i], nucleus::CanonicalBases::ACGT) ||
          read.aligned_quality()[i] < options_.min_base_quality()) {
        return i;
      }
    }
    return bases.size();
  };

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
    AddKmersAndEdges(bases_view, i, next_bad_position - k_, false /* is_ref */);
    i = next_bad_position + 1;
  }
}

std::vector<Path> DeBruijnGraph::CandidatePaths() const {
  std::vector<Path> terminated_paths;
  std::queue<Path> extendable_paths;

  CHECK_GT(boost::out_degree(source_, g_), 0);
  extendable_paths.push({source_});

  // Inefficient.
  while (!extendable_paths.empty()) {
    // Some windows can have an extremely branchy graph.  Ideally windows would
    // be chosen to avoid this.  We give up if we encounter too many paths.
    int n_total_paths = terminated_paths.size() + extendable_paths.size();
    if (n_total_paths > options_.max_num_paths()) {
      return {};
    }

    Path path = extendable_paths.front();
    extendable_paths.pop();
    Vertex last_v = path.back();
    // For each successor of last_v, add path::successor to the
    // appropriate queue.
    AdjacencyIterator vi, vend;
    std::tie(vi, vend) = boost::adjacent_vertices(last_v, g_);
    for (; vi != vend; ++vi) {
      Path extended_path(path);
      extended_path.push_back(*vi);
      if (*vi == sink_ || boost::out_degree(*vi, g_) == 0) {
        terminated_paths.push_back(extended_path);
      } else {
        extendable_paths.push(extended_path);
      }
    }
  }
  return terminated_paths;
}

string DeBruijnGraph::HaplotypeForPath(const Path& path) const {
  std::stringstream haplotype;
  for (Vertex v : path) {
    haplotype << g_[v].kmer[0];
  }
  if (!path.empty()) {
    haplotype << g_[path.back()].kmer.substr(1, k_ - 1);
  }
  return haplotype.str();
}

std::vector<std::string> DeBruijnGraph::CandidateHaplotypes() const {
  std::vector<std::string> haplotypes;
  for (const Path& path : CandidatePaths()) {
    haplotypes.push_back(HaplotypeForPath(path));
  }
  std::sort(haplotypes.begin(), haplotypes.end());
  return haplotypes;
}

string DeBruijnGraph::GraphViz() const {
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

void DeBruijnGraph::Prune() {
  // Remove low-weight edges not in the reference.
  boost::remove_edge_if(
      [this](const Edge& e) {
        return !g_[e].is_ref && g_[e].weight < options_.min_edge_weight();
      },
      g_);

  // Remove vertices not reachable forward from src or backward from sink.
  VertexIterator vbegin, vend;
  std::tie(vbegin, vend) = boost::vertices(g_);
  absl::flat_hash_set<Vertex> all_vertices(vbegin, vend);

  std::set<Vertex> fwd_reachable_vertices, rev_reachable_vertices;
  fwd_reachable_vertices = VerticesReachableFrom(
      source_, g_, IndexMap());
  rev_reachable_vertices = VerticesReachableFrom(
      sink_, boost::make_reverse_graph(g_), IndexMap());

  absl::flat_hash_set<Vertex> reachable_vertices;
  std::set_intersection(
      fwd_reachable_vertices.begin(), fwd_reachable_vertices.end(),
      rev_reachable_vertices.begin(), rev_reachable_vertices.end(),
      std::inserter(reachable_vertices, reachable_vertices.end()));
  for (Vertex v : all_vertices) {
    if (reachable_vertices.find(v) == reachable_vertices.end()) {
      kmer_to_vertex_.erase(g_[v].kmer);
      boost::clear_vertex(v, g_);
      boost::remove_vertex(v, g_);
    }
  }
  RebuildIndexMap();
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
