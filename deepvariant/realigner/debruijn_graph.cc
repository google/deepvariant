/*
 * Copyright 2017 Google Inc.
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
#include <memory>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/reverse_graph.hpp"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/strings/str_util.h"
#include "tensorflow/core/platform/logging.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using Vertex = DeBruijnGraph::Vertex;
using VertexIndexMap = DeBruijnGraph::VertexIndexMap;
using Edge = DeBruijnGraph::Edge;
using Path = DeBruijnGraph::Path;

using Read = nucleus::genomics::v1::Read;

using tensorflow::string;
using tensorflow::StringPiece;

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
    out << "[label=" << ei.weight << (ei.is_ref ? " color=red" : "") << "]";
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

Vertex DeBruijnGraph::EnsureVertex(StringPiece kmer) {
  Vertex v;
  auto vertex_find = kmer_to_vertex_.find(kmer);
  if (vertex_find != kmer_to_vertex_.end()) {
    v = (*vertex_find).second;
  } else {
    // N.B. TensorFlow StringPiece lacks explicit string conversion func.
    string kmer_copy = string(kmer.data(), kmer.size());
    v = boost::add_vertex(VertexInfo{kmer_copy}, g_);
    // N.B.: must use the long-lived string in the map key as the referent of
    // the StringPiece key.
    kmer_to_vertex_[StringPiece(g_[v].kmer)] = v;
  }
  return v;
}

Vertex DeBruijnGraph::VertexForKmer(StringPiece kmer) const {
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

DeBruijnGraph::DeBruijnGraph(const string& ref,
                             const std::vector<Read>& reads,
                             const Options& options,
                             int k)
    : options_(options), k_(k)
{
  CHECK_GT(k, 0);  // k should always be a positive integer.
  CHECK(static_cast<uint32_t>(k) < ref.size());
  AddEdgesForReference(ref);
  source_ = VertexForKmer(ref.substr(0, k_));
  sink_ = VertexForKmer(ref.substr(ref.size() - k_, k_));
  for (const Read& read : reads) {
    if (read.alignment().mapping_quality() >= options.min_mapq()) {
      AddEdgesForRead(read);
    }
  }
  RebuildIndexMap();
}

std::unique_ptr<DeBruijnGraph> DeBruijnGraph::Build(
    const string& ref, const std::vector<Read>& reads,
    const DeBruijnGraph::Options& options) {

  std::unique_ptr<DeBruijnGraph> graph, trial_graph;
  int max_k  = std::min(options.max_k(), static_cast<int>(ref.size()) - 1);

  for (int k = options.min_k(); k <= max_k; k += options.step_k()) {
    // If we can't get an acyclic graph from just the reference, we should go on
    // to the next k -- an optimization.
    // N.B.: MakeUnique doesn't work with private constructors.
    trial_graph = std::unique_ptr<DeBruijnGraph>(
        new DeBruijnGraph(ref, {}, options, k));
    if (trial_graph->HasCycle()) {
      continue;
    }
    graph = std::unique_ptr<DeBruijnGraph>(
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

Edge DeBruijnGraph::AddEdge(StringPiece from, StringPiece to, bool is_ref) {
  Vertex from_vertex = EnsureVertex(from);
  Vertex to_vertex = EnsureVertex(to);
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

void DeBruijnGraph::AddEdgesForReference(StringPiece ref) {
  StringPiece kmer_prev, kmer_cur;
  const signed int ref_length = ref.size();
  for (int i = 0; i < ref_length - k_ + 1; i++) {
    kmer_prev = kmer_cur;
    kmer_cur = ref.substr(i, k_);
    if (i > 0) {
      AddEdge(kmer_prev, kmer_cur, true);
    }
  }
}

void DeBruijnGraph::AddEdgesForRead(const nucleus::genomics::v1::Read& read) {
  string bases = tensorflow::str_util::Uppercase(read.aligned_sequence());
  StringPiece bases_view(bases);
  std::vector<int> qual(read.aligned_quality().begin(),
                        read.aligned_quality().end());
  CHECK(qual.size() == bases.size());

  const signed int read_length = bases.size();

  // This set maintains the QC-failing positions among [i..i+k].
  std::set<int> recent_qc_fail_positions;

  for (int i = 0; i < read_length - k_; ++i) {
    // Update QC fail set: remove (i-1), add (i+k) if it fails QC.
    recent_qc_fail_positions.erase(i - 1);
    if (!IsCanonicalBase(bases[i + k_], nucleus::CanonicalBases::ACGT) ||
        qual[i + k_] < options_.min_base_quality()) {
      recent_qc_fail_positions.insert(i + k_);
    }
    if (recent_qc_fail_positions.empty()) {
      AddEdge(bases_view.substr(i, k_), bases_view.substr(i + 1, k_), false);
    }
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

std::vector<string> DeBruijnGraph::CandidateHaplotypes() const {
  std::vector<string> haplotypes;
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
  std::unordered_set<Vertex> all_vertices(vbegin, vend);

  std::set<Vertex> fwd_reachable_vertices, rev_reachable_vertices;
  fwd_reachable_vertices = VerticesReachableFrom(
      source_, g_, IndexMap());
  rev_reachable_vertices = VerticesReachableFrom(
      sink_, boost::make_reverse_graph(g_), IndexMap());

  std::unordered_set<Vertex> reachable_vertices;
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
