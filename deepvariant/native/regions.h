#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"

namespace deepvariant {

// Parse a region string like "chr1", "chr1:100-200" (1-based, inclusive).
// Populates *out with a 0-based half-open [start, end) Range proto.
// Returns false on parse failure.
bool ParseRegionString(
    const std::string& s,
    const std::unordered_map<std::string, int64_t>& contig_lengths,
    nucleus::genomics::v1::Range* out);

// Build calling regions from reference contigs intersected with
// user-specified region strings. An empty regions_to_include means
// "all contigs". Regions in regions_to_exclude are subtracted.
std::vector<nucleus::genomics::v1::Range> BuildCallingRegions(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::vector<std::string>& regions_to_include,
    const std::vector<std::string>& regions_to_exclude);

// Return the subset of calling_regions assigned to shard task_id
// (0-based) out of num_shards total. Regions are partitioned by
// cumulative base-pair count to produce balanced shards.
std::vector<nucleus::genomics::v1::Range> ShardRegions(
    const std::vector<nucleus::genomics::v1::Range>& calling_regions,
    int task_id, int num_shards);

// Split each calling region into chunks of at most `partition_size`
// basepairs. Mirrors upstream's `regions.partition()` — required for
// realigner parity, since the WindowSelector + DBG run independently
// on each chunk and adjacent chunks emit overlapping windows at the
// chunk boundary.
std::vector<nucleus::genomics::v1::Range> PartitionRegions(
    const std::vector<nucleus::genomics::v1::Range>& calling_regions,
    int64_t partition_size);

}  // namespace deepvariant
