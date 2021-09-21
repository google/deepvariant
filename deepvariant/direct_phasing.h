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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_DIRECT_PHASING_H_
#define LEARNING_GENOMICS_DEEPVARIANT_DIRECT_PHASING_H_

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/flat_hash_map.h"
#include "absl/strings/string_view.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "absl/status/statusor.h"

namespace learning {
namespace genomics {
namespace deepvariant {


inline constexpr absl::string_view kUncalledAllele = "UNCALLED_ALLELE";

// redacted
using ReadIndex = uint16_t;

// Data type storing read id and mapping quality. It is used in AlleleInfo.
struct ReadSupportInfo {
  ReadIndex read_index;
  bool is_low_quality;
};

// Data type associated with graph nodes. It uniquely defines an allele by its
// type and bases along with the vector of supporting read ids.
struct AlleleInfo {
  int position = 0;
  std::string bases = "";
  std::vector<ReadSupportInfo> read_support;
};

struct VertexInfo {
  AlleleInfo allele_info;
};

struct EdgeInfo {
  float weight;
};

using BoostGraph =
    boost::adjacency_list<boost::setS,            // Out edge list type.
                          boost::listS,           // Vertex list type.
                          boost::bidirectionalS,  // Directed graph.
                          VertexInfo,             // Vertex label.
                          EdgeInfo>;              // Edge label.

using Vertex = boost::graph_traits<BoostGraph>::vertex_descriptor;
using Edge = boost::graph_traits<BoostGraph>::edge_descriptor;
using RawVertexIndexMap = absl::flat_hash_map<Vertex, int>;
using VertexIndexMap =
    boost::const_associative_property_map<RawVertexIndexMap>;

using VertexIterator = boost::graph_traits<BoostGraph>::vertex_iterator;
using EdgeIterator = boost::graph_traits<BoostGraph>::edge_iterator;
using AdjacencyIterator = boost::graph_traits<BoostGraph>::adjacency_iterator;

struct AlleleSupport {
  Vertex vertex;
  bool is_low_quality;
};

// Class that implements Direct Phasing algorithm. This class is only used by
// make_examples.py. There are two exported methods:
// * PhaseReads - called for each region and returns read phases calculated from
//                candidates.
// * GraphViz - auxiliary methods to create graphviz output for debugging
//              purposes.
class DirectPhasing {
 public:
  // Function returns read phases for each read in the input reads preserving
  // the order. Python wrapper will be used to add phases to read protos in
  // order to avoid copying gigabytes of memory.
  absl::StatusOr<std::vector<int>> PhaseReads(
      const std::vector<DeepVariantCall>& candidates,
      const std::vector<
          nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>& reads);

  // Helper function to output graph into graphviz for debugging. This function
  // is exported to Python.
  std::string GraphViz() const;

 private:
  // Build graph from candidates.
  void Build(
      const std::vector<DeepVariantCall>& candidates,
      const std::vector<
          nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>& reads);

  BoostGraph graph_;
  Vertex source_;
  Vertex sink_;
  RawVertexIndexMap vertex_index_map_;  // This is needed for GraphViz.

  // Allele support for each read. Map is keyed by read id. Alleles are sorted
  // by position. This map allows to quickly query all alleles that a read
  // supports. Boolean variable designates if read to allele support is
  // low_quality. If true then read supports the allele with low quality.
  absl::flat_hash_map<std::string, std::vector<AlleleSupport>>
      read_to_alleles_;
  // Map read name to read id.
  absl::flat_hash_map<std::string, ReadIndex> read_to_index_;
};

// Helper functions.

// Calculate AlleleType by comparing alt allele size and candidate interval.
AlleleType AlleleTypeFromCandidate(
    std::string_view bases,
    const DeepVariantCall& candidate);

// Calculate number of alt alleles that are SUBs.
int NumOfSubstitutionAlleles(const DeepVariantCall& candidate);

// Calculate number of alt alleles that are INDELs.
int NumOfIndelAlleles(const DeepVariantCall& candidate);

// Calculate the depth of all SUB alt alleles. This is done by enumerating all
// supporting reads for all SUB alleles.
int SubstitutionAllelesDepth(const DeepVariantCall& candidate);

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_DIRECT_PHASING_H_
