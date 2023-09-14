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
//
// Usage:
// blaze-bin/learning/genomics/deepvariant/merge_phased_reads_cpp \
// --input_path <Path to sharded tsv file> \
// --output_path <Path to output file> \
// --logtostderr

#include "deepvariant/merge_phased_reads.h"

#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "re2/re2.h"

ABSL_FLAG(std::string, input_path, "", "Sharded input.");
ABSL_FLAG(std::string, output_path, "", "Output path.");

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
                      "-of-", absl::StrFormat("%05d", num_shards), spec.suffix,
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
    if (!sharded_input.status().ok()) {
      LOG(FATAL) << sharded_input.status();
    }
    const std::string filename =
        generate_sharded_filename(sharded_input.value(), shard);
    LOG(INFO) << "Loading " << filename;

    std::ifstream csv_file;
    csv_file.open(filename);
    std::string line;
    bool is_first_line = true;
    while (std::getline(csv_file, line)) {
      if (is_first_line) {
        is_first_line = false;
        continue;
      }

      std::istringstream iss(line);
      std::vector<std::string> tokens(5);
      MergedPhaseRead phased_read;
      int i = 0;
      while (std::getline(iss, tokens[i], '\t')) {
        i++;
      }
      int id = UpdateReadsMap(tokens[0]);
      unmerged_reads_.push_back({
          .fragment_name = tokens[0],
          .phase = std::stoi(tokens[1]),
          .region_order = std::stoi(tokens[2]),
          .shard = shard,
          .id = id,
      });
    }
  }
  LOG(INFO) << "Total records loaded: " << merged_reads_.size();
}

int Merger::UpdateReadsMap(const std::string& fragment_name) {
  auto it = merged_reads_map_.find(fragment_name);
  if (it == merged_reads_map_.end()) {
    merged_reads_.push_back(
        {.fragment_name = fragment_name, .phase = 0, .phase_dist = {}});
    merged_reads_map_[fragment_name] = merged_reads_.size() - 1;
  }
  return merged_reads_map_[fragment_name];
}

void Merger::GroupReads() {
  for (auto index = 0; index < unmerged_reads_.size(); index++) {
    Group& read_group = groups_[{.shard = merged_reads_[index].shard,
                                 .region = merged_reads_[index].region_order}];
    size_t merged_index =
        merged_reads_map_[unmerged_reads_[index].fragment_name];
    read_group.merged_id_to_unmerged_id[merged_index] = index;
  }
  num_groups_ = groups_.size();
}

bool Merger::CompareGroups(const ShardRegion& group_1,
                           const ShardRegion& group_2) const {
  int num_reads_not_matching_phase = 0;
  int num_reads_matching_phase = 0;
  auto group_1_it = groups_.find(group_1);
  auto group_2_it = groups_.find(group_2);
  if (group_1_it == groups_.end() || group_2_it == groups_.end()) {
    LOG(FATAL) << "CompareGroups() called with invalid arguments";
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
  return num_reads_not_matching_phase > num_reads_matching_phase;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main(int argc, char** argv) {
  QCHECK(FLAGS_input_path.CurrentValue().empty() ||
      FLAGS_output_path.CurrentValue().empty()) <<
    "ERROR: --input_path and --output_path flags must be set.";

  return 0;
}
