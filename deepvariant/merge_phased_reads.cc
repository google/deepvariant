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

// This utility is used to merge phased reads from different shards.
// We can find a consistent phasing if there are reads that overlap multiple
// shards. Please note, that input file must be local, this utility does not
// support Google paths.

#include "deepvariant/merge_phased_reads.h"

#include <sys/types.h>
#include <unistd.h>

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/ascii.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "re2/re2.h"

namespace learning {
namespace genomics {
namespace deepvariant {

absl::StatusOr<ShardedFileSpec>
parse_sharded_file_spec(absl::string_view file_spec) {
  ShardedFileSpec output;
  std::string nshards_str;
  static const absl::string_view kShardSpecPattern =
      R"((.*)\@(\d*[1-9]\d*)(?:\.(.+))?)";
  if (!RE2::FullMatch(file_spec, kShardSpecPattern, &output.basename,
                      &nshards_str, &output.suffix)) {
    // No @<N> found.
    return absl::InvalidArgumentError(
        absl::StrCat("'", file_spec, "' is not a valid sharded file spec."));
  }

  int nshards_int;
  if (!absl::SimpleAtoi(nshards_str, &nshards_int)) {
    return absl::OutOfRangeError(absl::StrCat(
        "'", file_spec, "' is not a valid sharded file spec. '", nshards_str,
        "' is too large to represent as an integer."));
  }
  output.nshards = nshards_int;

  return output;
}

int shard_with(int num_shards) {
  if (num_shards < 100000) return 5;
  if (num_shards < 1000000) return 6;
  if (num_shards < 10000000) return 7;
  if (num_shards < 100000000) return 8;
  if (num_shards < 1000000000) return 9;
  LOG(FATAL) << "num_shards == " << num_shards << ": Unsupported";
}

// Generates a sharded file name from ShardedFileSpec and shard number.
// Format: <basename>-<shard>-of-<num_shards>[.<extension>]
std::string generate_sharded_filename(const ShardedFileSpec& spec, int shard) {
  const int num_shards = spec.nshards;
  DCHECK_LE(0, shard);
  DCHECK_LE(0, num_shards);
  return absl::StrCat(spec.basename, "-",
                      absl::StrFormat("%0*d", shard_with(num_shards), shard),
                      "-of-", absl::StrFormat("%05d", num_shards), ".",
                      spec.suffix);
}

// Loads input files from a sharded path.
void Merger::LoadFromFiles(absl::string_view input_path) {
  absl::StatusOr<ShardedFileSpec> sharded_input =
      parse_sharded_file_spec(input_path);
  if (!sharded_input.ok()) {
    LOG(FATAL) << "Could not read " << input_path;
  }
  num_shards_ = sharded_input->nshards;
  LOG(INFO) << "basename=" << sharded_input->basename << ", " << num_shards_
            << " shards";

  for (int shard = 0; shard < num_shards_; ++shard) {
    const std::string filename =
        generate_sharded_filename(sharded_input.value(), shard);
    LOG(INFO) << "Loading " << filename;

    std::ifstream csv_file;
    csv_file.open(filename);
    if (!csv_file.is_open()) {
      LOG(FATAL) << "Could not open " << filename;
    }
    std::string line;
    bool is_first_line = true;
    while (std::getline(csv_file, line)) {
      absl::string_view stripped = absl::StripAsciiWhitespace(line);
      if (stripped.empty()) {
        continue;
      }
      if (is_first_line) {
        is_first_line = false;
        continue;
      }

      std::vector<absl::string_view> tokens = absl::StrSplit(stripped, '\t');
      if (tokens.size() < 3) {
        LOG(WARNING) << "Skipping malformed line with " << tokens.size()
                     << " columns (expected >= 3): " << line;
        continue;
      }
      int region;
      int phase;
      if (!absl::SimpleAtoi(tokens[2], &region) ||
          !absl::SimpleAtoi(tokens[1], &phase)) {
        LOG(WARNING) << "Skipping line with invalid integer field: " << line;
        continue;
      }
      if (region <= 0) {
        LOG(WARNING) << "Skipping line with invalid region: " << line;
        continue;
      }
      int id = UpdateReadsMap(tokens[0]);
      unmerged_reads_.push_back({
          .fragment_name = std::string(tokens[0]),
          .phase = phase,
          .region_order = region,
          .shard = shard,
          .id = id,
      });
    }
  }
  LOG(INFO) << "Total records loaded: " << unmerged_reads_.size();
}

int Merger::UpdateReadsMap(absl::string_view fragment_name) {
  auto it = merged_reads_map_.find(fragment_name);
  if (it == merged_reads_map_.end()) {
    merged_reads_.push_back({.fragment_name = std::string(fragment_name),
                             .phase = 0,
                             .phase_dist = {}});
    merged_reads_map_[fragment_name] = merged_reads_.size() - 1;
  }
  return merged_reads_map_[fragment_name];
}

void Merger::GroupReads() {
  for (auto index = 0; index < unmerged_reads_.size(); index++) {
    Group& read_group =
        groups_[{.shard = unmerged_reads_[index].shard,
                 .region = unmerged_reads_[index].region_order}];
    size_t merged_index =
        merged_reads_map_[unmerged_reads_[index].fragment_name];
    read_group.merged_id_to_unmerged_id[merged_index] = index;
  }
  num_groups_ = groups_.size();
}

// Returns true if number of reads with mismatched phases are greater than
// number of reads with matching phases.
Merger::ComparisonResult Merger::CompareGroups(
    const ShardRegion& group_1, const ShardRegion& group_2) const {
  int num_reads_not_matching_phase = 0;
  int num_reads_matching_phase = 0;
  auto group_1_it = groups_.find(group_1);
  auto group_2_it = groups_.find(group_2);
  if (group_1_it == groups_.end() || group_2_it == groups_.end()) {
    return ComparisonResult::NOT_ENOUGH_OVERLAP;
  }
  // Iterate read ids in group_2.
  for (auto [merged_reads_idx_2, unmerged_reads_idx2] :
       group_2_it->second.merged_id_to_unmerged_id) {
    // Find a matching read id in group_1.
    auto group1_index_map_it =
        group_1_it->second.merged_id_to_unmerged_id.find(merged_reads_idx_2);
    // If read is not found in group_1 then do nothing.
    if (group1_index_map_it ==
        group_1_it->second.merged_id_to_unmerged_id.end()) {
      continue;
    }
    // Only consider pairs that have different phases. If one of the reads have
    // phase zero it means it is unphased and we cannot compare it to another
    // one.
    int unmerged_reads_idx1 = group1_index_map_it->second;
    if (unmerged_reads_[unmerged_reads_idx2].phase == 0 ||
        unmerged_reads_[unmerged_reads_idx1].phase == 0) {
      continue;
    }
    // Count number of reads that have matching and unmatching phases.
    if (unmerged_reads_[unmerged_reads_idx2].phase !=
        unmerged_reads_[unmerged_reads_idx1].phase) {
      num_reads_not_matching_phase++;
    } else {
      num_reads_matching_phase++;
    }
  }
  if (std::abs(num_reads_not_matching_phase - num_reads_matching_phase) < 2) {
    return ComparisonResult::NOT_ENOUGH_OVERLAP;
  }
  return num_reads_not_matching_phase > num_reads_matching_phase
             ? ComparisonResult::SWITCH
             : ComparisonResult::MATCH;
}

// Reverses phase for the group. Phases are reversed as follow:
// Phase 1 -> Phase 2, Phase 2 -> Phase 1, Phase 0 -> Phase 0.
void Merger::ReversePhasing(const ShardRegion& group) {
  for (auto [phased_reads_idx, unphased_reads_idx] :
       groups_[group].merged_id_to_unmerged_id) {
    if (unmerged_reads_[unphased_reads_idx].phase > 0) {
      unmerged_reads_[unphased_reads_idx].phase =
          3 - unmerged_reads_[unphased_reads_idx].phase;
    }
  }
}

// Merge reads from the group into merged_reads_ vector. If read already exist
// in the merged_reads_ vector it's phase is not changed unless it is 0.
void Merger::MergeGroup(const ShardRegion& group) {
  for (auto [phased_read_index, unphased_reads_index] :
       groups_[group].merged_id_to_unmerged_id) {
    // If merged_reads_ already contains the read we keep its phase and update
    // phase distribution for the read.
    std::string read_id = unmerged_reads_[unphased_reads_index].fragment_name;
    auto& merged_read = merged_reads_[phased_read_index];
    if (merged_read.phase == 0) {
      merged_read.phase = unmerged_reads_[unphased_reads_index].phase;
    }
    merged_read.phase_dist[unmerged_reads_[unphased_reads_index].phase]++;
  }
}

// Main entry point function. Reads are merged one group at a time iterating
// groups in the same order they were processed by make_examples.
// 1. reads are grouped by shard and region.
// 2. Read phases are compared between last merged group and the group being
//    merged. If most phases are not matched then phase is reversed for the
//    group.
// 3. Group is merged into merged_reads_.
void Merger::MergeReads(absl::string_view switches_output_path) {
  // Write to a temporary file first, then atomically rename to prevent
  // readers from seeing partially-written data.
  std::string tmp_path;
  std::ofstream csv_file_switches;
  if (!switches_output_path.empty()) {
    tmp_path = absl::StrCat(switches_output_path, ".tmp.", getpid());
    csv_file_switches.open(tmp_path);
    if (!csv_file_switches.is_open()) {
      LOG(FATAL) << "Failed to open temporary switches file: " << tmp_path;
    }
  }
  GroupReads();
  int cur_region = 1;
  int processed_groups = 0;
  ShardRegion prev_group;
  while (processed_groups < num_groups_) {
    for (int shard = 0; shard < num_shards_; shard++) {
      const auto& cur_group_it = groups_.find({shard, cur_region});
      if (cur_group_it == groups_.end()) {
        continue;
      }
      const ComparisonResult comparison_result =
          CompareGroups(prev_group, {shard, cur_region});
      if (comparison_result == ComparisonResult::SWITCH) {
        ReversePhasing({shard, cur_region});
      }
      if (csv_file_switches.is_open()) {
        csv_file_switches << shard << "\t" << cur_region << "\t"
                          << static_cast<int>(comparison_result) << "\n";
      }
      MergeGroup({shard, cur_region});
      processed_groups++;
      LOG_EVERY_N(INFO, 1000) << "Processed " << processed_groups << " groups";
      prev_group = cur_group_it->first;
    }
    cur_region++;
  }
  if (csv_file_switches.is_open()) {
    csv_file_switches.close();
    if (csv_file_switches.fail()) {
      LOG(FATAL) << "Failed to write out phase switches.";
    }
  }
  if (!tmp_path.empty()) {
    if (std::rename(tmp_path.c_str(),
                    std::string(switches_output_path).c_str()) != 0) {
      LOG(FATAL) << "Failed to rename " << tmp_path << " to "
                 << switches_output_path;
    }
  }
}

void MergerPeer::SetUnmergedReads(
    Merger& merger, absl::Span<const UnmergedRead> unmerged_reads) {
  // Cannot use absl::btree_set as the rbegin() iterator is not provided.
  std::set<int> shards;
  for (const auto& unmerged_read : unmerged_reads) {
    shards.insert(unmerged_read.shard);
    merger.unmerged_reads_.push_back(unmerged_read);
    merger.UpdateReadsMap(unmerged_read.fragment_name);
  }
  if (!shards.empty()) {
    merger.num_shards_ = *shards.rbegin() + 1;
  } else {
    merger.num_shards_ = 0;
  }
}

int Merger::CorrectPhasing() {
  int count_reads_corrected = 0;
  for (auto& read_info : merged_reads_) {
    std::array<int, 3> phase_counts{};
    for (int phase : {1, 2}) {
      auto it = read_info.phase_dist.find(phase);
      phase_counts[phase] = it == read_info.phase_dist.end()
                        ? 0
                        : it->second;
    }
    // Correct phasing taking the majority phase. If there are equal number of
    // phase 1 and phase 2 hits then reads become unphased.
    int old_phase = read_info.phase;
    if (phase_counts[1] == phase_counts[2]) {
      read_info.phase = 0;
    } else {
      read_info.phase = phase_counts[1] > phase_counts[2] ? 1 : 2;
    }
    if (old_phase != read_info.phase) {
      count_reads_corrected++;
    }
  }
  return count_reads_corrected;
}

void Merger::CorrectAndPrintReadStats(absl::string_view output_path) {
  std::string tmp_path;
  std::ofstream csv_file;
  if (!output_path.empty()) {
    tmp_path = absl::StrCat(output_path, ".tmp.", getpid());
    csv_file.open(tmp_path);
    if (!csv_file.is_open()) {
      LOG(FATAL) << "Failed to open temporary corrected stats file: "
                 << tmp_path;
    }
  }
  int count_reads_corrected = CorrectPhasing();
  std::array<int, 3> phase_counts = {0, 0, 0};
  int n_reads = 0;
  for (auto& read_info : merged_reads_) {
    phase_counts[read_info.phase]++;
    if (csv_file.is_open()) {
      LOG_EVERY_N(INFO, 20000) << "Written " << n_reads << " reads";
      csv_file << read_info.fragment_name << "\t" << read_info.phase << "\n";
    }
    n_reads++;
  }
  if (csv_file.is_open()) {
    csv_file.close();
    if (csv_file.fail()) {
      LOG(FATAL) << "Failed to write out corrected read stats.";
    }
  }
  if (!tmp_path.empty()) {
    if (std::rename(tmp_path.c_str(), std::string(output_path).c_str()) != 0) {
      LOG(FATAL) << "Failed to rename " << tmp_path << " to " << output_path;
    }
  }
  LOG(INFO) << "Count of reads not phased " << phase_counts[0];
  LOG(INFO) << "Count of reads with phase 1 " << phase_counts[1];
  LOG(INFO) << "Count of reads with phase 2 " << phase_counts[2];
  LOG(INFO) << "Count of corrected reads " << count_reads_corrected;
  LOG(INFO) << "Total reads processed: " << n_reads;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
