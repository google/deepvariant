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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_DEBRUIJN_GRAPH_H_
#define LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_DEBRUIJN_GRAPH_H_

#include <map>
#include <memory>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "third_party/nucleus/protos/reads.pb.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/lib/hash/hash.h"
#include "tensorflow/core/platform/types.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using tensorflow::string;
using tensorflow::StringPiece;
using tensorflow::StringPieceHasher;

struct VertexInfo {
  string kmer;
};

struct EdgeInfo {
  int weight;   // The # of multiedges this edge represents.
  bool is_ref;  // True iff this edge is reflected by the reference sequence.
};


class DeBruijnGraph {
 private:
  using BoostGraph = boost::adjacency_list<
    boost::setS,            // Out edge list type.
    boost::listS,           // Vertex list type.
    boost::bidirectionalS,  // Directed graph.
    VertexInfo,             // Vertex label.
    EdgeInfo>;              // Edge label.

  using VertexIterator = boost::graph_traits<BoostGraph>::vertex_iterator;
  using EdgeIterator = boost::graph_traits<BoostGraph>::edge_iterator;
  using AdjacencyIterator = boost::graph_traits<BoostGraph>::adjacency_iterator;

 public:
  using Vertex = boost::graph_traits<BoostGraph>::vertex_descriptor;
  using Edge = boost::graph_traits<BoostGraph>::edge_descriptor;
  using Path = std::vector<Vertex>;

  using RawVertexIndexMap = std::map<Vertex, int>;
  using VertexIndexMap =
      boost::const_associative_property_map<RawVertexIndexMap>;

  using Options = RealignerOptions::DeBruijnGraphOptions;

 private:
  // Convenience method for rebuilding a table usable as the vertex_index_t
  // property that algorithms require.
  void RebuildIndexMap();

  // Accessor for the vertex index table.
  VertexIndexMap IndexMap() const;

  // Ensure a vertex with label kmer is present--adding if necessary.
  Vertex EnsureVertex(StringPiece kmer);

  // Look up the vertex with this kmer label.
  Vertex VertexForKmer(StringPiece kmer) const;

  // Is this graph cyclic?
  bool HasCycle() const;

  // Private constructor.  Public interface via factory only allows access to
  // acyclic DeBruijn graphs.  Argument `k` is used to construct the graph;
  // filtering settings are taken from options.
  DeBruijnGraph(const string& ref,
                const std::vector<nucleus::genomics::v1::Read>& reads,
                const Options& options,
                int k);

  // Add edge, implicitly adding the vertices if needed.  If such an edge is
  // already present, we merely increment its weight to reflect its "multiedge"
  // degree.
  Edge AddEdge(StringPiece from, StringPiece to, bool is_ref);

  // Add all the edges implied by the given reference string.
  void AddEdgesForReference(StringPiece ref);

  // Add all the edges implied by the given read (and according to our edge
  // filtering criteria).
  void AddEdgesForRead(const nucleus::genomics::v1::Read& read);

  // Returns candidate haplotype paths through the graph.  If more that
  // options.max_num_paths paths are found, this will return an empty vector.
  std::vector<Path> CandidatePaths() const;

  // Returns the string traced by a path through the graph.
  string HaplotypeForPath(const Path& path) const;

  // Removes low weight non-ref edges from the graph.
  void Prune();

 public:
  // We attempt to build acyclic graphs with increasing kmer size until we
  // achieve an acyclic graph---kmer size starts with options.min_k and goes
  // linearly up to options.max_k, stepping by options.step_k.  If we are able
  // to construct an acyclic DeBruijn graph in this manner, it is returned;
  // otherwise we return nullptr.
  static std::unique_ptr<DeBruijnGraph> Build(
      const string& ref,
      const std::vector<nucleus::genomics::v1::Read>& reads,
      const Options& options);

  // Gets all the candidate haplotypes defined by paths through the graph.  If
  // more than options.max_num_paths() haplotypes are identified, returns an
  // empty vector, to preempt excessive computation.
  std::vector<string> CandidateHaplotypes() const;

  // Gets a GraphViz representation of the graph.
  string GraphViz() const;

  // Gets the kmer size used in this graph.
  int KmerSize() const { return k_; }

 private:
  BoostGraph g_;
  Options options_;
  int k_;
  Vertex source_;
  Vertex sink_;
  // N.B.: kmer strings are owned by VertexInfo objects;
  // map keys are merely pointers.
  std::unordered_map<StringPiece, Vertex, StringPieceHasher>
      kmer_to_vertex_;
  RawVertexIndexMap vertex_index_map_;
};


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_REALIGNER_DEBRUIJN_GRAPH_H_
