/*
 * Copyright 2023 Google LLC.
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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_MERGE_PHASED_READS_H_
#define LEARNING_GENOMICS_DEEPVARIANT_MERGE_PHASED_READS_H_

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Structure to hold merged reads with phasing.
struct MergedPhaseRead {
  std::string fragment_name;                 // Uniquely identifies a read.
  int phase;                                 // Phasing {0, 1, 2}.
  absl::flat_hash_map<int, int> phase_dist;  // Different phases the read was
                                             // assigned after merging. This is
                                             // needed to count number of reads
                                             // with inconsistent phasing.
};

// Input reads are grouped by shard and region order. At each step the merging
// is done between two groups with adjacent shards and the same region order.
// For example group (shard_2, region_1) is merged with (shard_1, region_1).
// Each group contains the map of unique read id to the read's phase.
// When groups are merged we need to find reads with the same
// fragment_name. To make it faster numeric IDs are used instead of string ids.
struct Group {
  absl::flat_hash_map<int, int> read_id_to_phase;
};

struct ShardRegion {
  int shard;
  int region;

  template <typename H>
  friend H AbslHashValue(H h, const ShardRegion& sr) {
    return H::combine(std::move(h), sr.shard, sr.region);
  }
};

inline bool operator==(const ShardRegion& l, const ShardRegion& r) {
  return l.shard == r.shard && l.region == r.region;
}

// Implementation of phased reads merging algorithm.
class Merger {
 public:
  // Loads input files.
  void LoadFromFiles(absl::string_view input_path);

  // Main API entry. Call it to merge reads.
  void MergeReads();

  // Scans reads for inconsistent phasing, correct where possible and print out
  // the results.
  void CorrectAndPrintout(const std::string_view& output_path);

 private:
  // Groups reads.
  void GroupReads();

  // Helper function to compare two reads.
  bool ComparePhasing(const ShardRegion& group_1,
                      const ShardRegion& group_2) const;
  void ReversePhasing(const ShardRegion& group);
  void MergeGroup(const ShardRegion& group);
  int UpdateReadsMap(const std::string& fragment_name);

  // Merged reads with phasing data.
  std::vector<MergedPhaseRead> phased_reads_;

  // Map from read_id to phased_reads index.
  absl::flat_hash_map<std::string, int> phased_reads_map_;

  // Map from (shard, region_num) to Group.
  absl::flat_hash_map<ShardRegion, Group> groups_;
  int num_shards_;
  int num_groups_;
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_MERGE_PHASED_READS_H_
