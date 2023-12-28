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
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Structure to hold input reads.
struct UnmergedRead {
  std::string fragment_name;  // Uniquely identifies a read.
  int phase = 0;              // Phasing {0, 1, 2}.
  int region_order = 0;       // Each make_examples shard
                              // runs over a set of regions. This
                              // field contains an order of a
                              // region.
  int shard = 0;              // Shard.
  int id = 0;                 // id in merged_reads_
  std::vector<int> consensus_phases;
};

// Structure to hold merged reads with phasing.
struct MergedPhaseRead {
  std::string fragment_name;                 // Uniquely identifies a read.
  int phase = 0;                             // Phasing {0, 1, 2}.
  absl::flat_hash_map<int, int> phase_dist;  // Different phases the read was
                                             // assigned after merging. This is
                                             // needed to count number of reads
                                             // with inconsistent phasing.
};

// Group of related reads.
struct Group {
  // Key is a merged read id, value is unmerged read id.
  absl::flat_hash_map<int, int> merged_id_to_unmerged_id;
};

struct ShardRegion {
  int shard = 0;
  int region = 0;

  template <typename H>
  friend H AbslHashValue(H h, const ShardRegion& sr) {
    return H::combine(std::move(h), sr.shard, sr.region);
  }
};

inline bool operator==(const ShardRegion& l, const ShardRegion& r) {
  return l.shard == r.shard && l.region == r.region;
}

// The components of a sharded file spec.
struct ShardedFileSpec {
  std::string basename;
  int nshards = 0;
  std::string suffix;
};

// Parses a simple sharded file spec into its components.
absl::StatusOr<ShardedFileSpec> parse_sharded_file_spec(
    absl::string_view file_spec);

// Implementation of phased reads merging algorithm.
class Merger {
 public:
  // Loads input files.
  void LoadFromFiles(absl::string_view input_path);

  // Main API entry. Call it to merge reads.
  void MergeReads();

  // Corrects phasing of reads that have an inconsistent phasing.
  // Returns the number of corrected reads. As a result of correction the phase
  // can be reversed or reset to zero.
  int CorrectPhasing();

  // Scans reads for inconsistent phasing, correct where possible and print out
  // the results.
  void CorrectAndPrintReadStats(const std::string& output_path);

 private:
  friend class MergerPeer;

  // Groups reads.
  void GroupReads();

  // Helper function to compare two reads.
  bool CompareGroups(const ShardRegion& group_1,
                     const ShardRegion& group_2) const;
  void ReversePhasing(const ShardRegion& group);
  void MergeGroup(const ShardRegion& group);
  int UpdateReadsMap(absl::string_view fragment_name);

  std::vector<UnmergedRead> unmerged_reads_;
  // Merged reads with phasing data.
  std::vector<MergedPhaseRead> merged_reads_;

  // Map from read name to merged_reads_ index.
  absl::flat_hash_map<std::string, int> merged_reads_map_;

  // Input reads are grouped by shard and region order. At each step the merging
  // is done between two groups with adjacent shards and the same region order.
  // For example group (shard_2, region_1) is merged with (shard_1, region_1).
  // Each group contains the map of unique read id to the read's phase.
  // When groups are merged we need to find reads with the same
  // fragment_name. To make it faster numeric IDs are used instead of string
  // ids.
  absl::flat_hash_map<ShardRegion, Group> groups_;
  int num_shards_;
  int num_groups_;
};

// Peer class for unit testing.
class MergerPeer {
 public:
  // Populates unmerged_reads_. This is needed for unit testing.
  static void SetUnmergedReads(
      Merger& merger, const std::vector<UnmergedRead>& unmerged_reads);

  // Returns merged_reads_. This is needed for unit testing.
  static const std::vector<MergedPhaseRead>& merged_reads(
      const Merger& merger) {
    return merger.merged_reads_;
  }
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_MERGE_PHASED_READS_H_
