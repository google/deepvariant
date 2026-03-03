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

#include "deepvariant/realigner/debruijn_graph.h"

#include <memory>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using ::testing::ElementsAre;

class DeBruijnGraphTest : public ::testing::Test {
 protected:
  DeBruijnGraphOptions options() {
    DeBruijnGraphOptions options;
    options.set_min_k(10);
    options.set_max_k(100);
    options.set_step_k(1);
    options.set_min_mapq(20);
    options.set_min_base_quality(20);
    options.set_min_edge_weight(1);
    options.set_max_num_paths(10);
    return options;
  }
};

TEST_F(DeBruijnGraphTest, TestCollapse) {
  std::string ref = "GATTACA";
  // Use k=3.
  // GAT -> ATT -> TTA -> TAC -> ACA
  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads;
  auto dbg = DeBruijnGraph::Build(ref, reads, opts);
  ASSERT_NE(dbg, nullptr);

  EXPECT_THAT(dbg->CandidateHaplotypes(), ElementsAre("GATTACA"));

  // Before collapse, it should have 5 vertices.
  // Actually, I can't easily check vertex count from public API, but I can
  // check haplotypes.

  dbg->Collapse();

  // After collapse, it should still have the same haplotype.
  EXPECT_THAT(dbg->CandidateHaplotypes(), ElementsAre("GATTACA"));
}

TEST_F(DeBruijnGraphTest, TestCollapseWithBranch) {
  std::string ref = "GATTACA";
  // GAT -> ATT -> TTA -> TAC -> ACA

  // Create a read that branches and then rejoins.
  // GAT -> ATG -> TGA -> GAC -> ACA
  nucleus::genomics::v1::Read read;
  read.set_aligned_sequence("GATGACA");
  read.set_aligned_quality(std::string(7, 30));
  read.mutable_alignment()->set_mapping_quality(30);

  DeBruijnGraphOptions opts = options();
  opts.set_min_k(3);
  opts.set_max_k(3);

  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>> reads;
  reads.push_back(nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(
    &read));

  auto dbg = DeBruijnGraph::Build(ref, reads, opts);
  ASSERT_NE(dbg, nullptr);

  EXPECT_THAT(dbg->CandidateHaplotypes(), ElementsAre("GATGACA", "GATTACA"));

  dbg->Collapse();

  EXPECT_THAT(dbg->CandidateHaplotypes(), ElementsAre("GATGACA", "GATTACA"));
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
