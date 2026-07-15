#include "deepvariant/native/regions.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "absl/log/log.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_replace.h"
#include "absl/strings/str_split.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"

namespace deepvariant {

bool ParseRegionString(
    const std::string& s,
    const std::unordered_map<std::string, int64_t>& contig_lengths,
    nucleus::genomics::v1::Range* out) {
  // Find last ':' to split "chrom:start-end".
  auto colon = s.rfind(':');
  if (colon == std::string::npos) {
    // Whole contig.
    auto it = contig_lengths.find(s);
    if (it == contig_lengths.end()) {
      LOG(ERROR) << "Unknown contig: " << s;
      return false;
    }
    out->set_reference_name(s);
    out->set_start(0);
    out->set_end(it->second);
    return true;
  }

  std::string chrom = s.substr(0, colon);
  std::string range_part = s.substr(colon + 1);

  auto dash = range_part.find('-');
  int64_t start_1based, end_1based;
  // Strip comma grouping from numeric substrings (e.g. "1,000,000") before
  // parsing, matching upstream ranges.parse_literal which accepts [0-9,]+.
  if (dash == std::string::npos) {
    // Single position: "chr1:1000" → half-open [999, 1000)
    const std::string pos_str = absl::StrReplaceAll(range_part, {{",", ""}});
    if (!absl::SimpleAtoi(pos_str, &start_1based)) {
      LOG(ERROR) << "Cannot parse position in region: " << s;
      return false;
    }
    end_1based = start_1based;
  } else {
    const std::string start_str =
        absl::StrReplaceAll(range_part.substr(0, dash), {{",", ""}});
    const std::string end_str =
        absl::StrReplaceAll(range_part.substr(dash + 1), {{",", ""}});
    if (!absl::SimpleAtoi(start_str, &start_1based) ||
        !absl::SimpleAtoi(end_str, &end_1based)) {
      LOG(ERROR) << "Cannot parse range in region: " << s;
      return false;
    }
  }

  auto it = contig_lengths.find(chrom);
  if (it == contig_lengths.end()) {
    LOG(ERROR) << "Unknown contig: " << chrom;
    return false;
  }

  out->set_reference_name(chrom);
  out->set_start(start_1based - 1);  // 1-based → 0-based
  out->set_end(std::min(end_1based, it->second));
  return true;
}

std::vector<nucleus::genomics::v1::Range> BuildCallingRegions(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::vector<std::string>& regions_to_include,
    const std::vector<std::string>& regions_to_exclude) {
  // Build a length map.
  std::unordered_map<std::string, int64_t> lengths;
  for (const auto& c : contigs) {
    lengths[c.name()] = c.n_bases();
  }

  std::vector<nucleus::genomics::v1::Range> regions;

  if (regions_to_include.empty()) {
    // All contigs in reference order.
    for (const auto& c : contigs) {
      nucleus::genomics::v1::Range r;
      r.set_reference_name(c.name());
      r.set_start(0);
      r.set_end(c.n_bases());
      regions.push_back(std::move(r));
    }
  } else {
    for (const auto& spec : regions_to_include) {
      nucleus::genomics::v1::Range r;
      if (ParseRegionString(spec, lengths, &r)) {
        regions.push_back(std::move(r));
      }
    }
  }

  if (regions_to_exclude.empty()) return regions;

  // Build exclude set as sorted intervals per contig.
  std::vector<nucleus::genomics::v1::Range> exclusions;
  for (const auto& spec : regions_to_exclude) {
    nucleus::genomics::v1::Range r;
    if (ParseRegionString(spec, lengths, &r)) {
      exclusions.push_back(std::move(r));
    }
  }

  // For each inclusion region, subtract all overlapping exclusion regions.
  std::vector<nucleus::genomics::v1::Range> result;
  for (auto& inc : regions) {
    std::vector<std::pair<int64_t, int64_t>> fragments = {
        {inc.start(), inc.end()}};
    for (const auto& exc : exclusions) {
      if (exc.reference_name() != inc.reference_name()) continue;
      std::vector<std::pair<int64_t, int64_t>> next;
      for (auto& [s, e] : fragments) {
        if (exc.end() <= s || exc.start() >= e) {
          next.push_back({s, e});
        } else {
          if (s < exc.start()) next.push_back({s, exc.start()});
          if (exc.end() < e) next.push_back({exc.end(), e});
        }
      }
      fragments = std::move(next);
    }
    for (auto& [s, e] : fragments) {
      if (s < e) {
        nucleus::genomics::v1::Range r;
        r.set_reference_name(inc.reference_name());
        r.set_start(s);
        r.set_end(e);
        result.push_back(std::move(r));
      }
    }
  }
  return result;
}

std::vector<nucleus::genomics::v1::Range> ShardRegions(
    const std::vector<nucleus::genomics::v1::Range>& calling_regions,
    int task_id, int num_shards) {
  if (num_shards <= 1) return calling_regions;

  // Compute total bp and target bp per shard.
  int64_t total_bp = 0;
  for (const auto& r : calling_regions) total_bp += r.end() - r.start();
  int64_t target = (total_bp + num_shards - 1) / num_shards;

  int64_t accumulated = 0;
  int current_shard = 0;
  std::vector<nucleus::genomics::v1::Range> shard_regions;

  for (const auto& r : calling_regions) {
    int64_t r_start = r.start();
    int64_t r_end = r.end();

    while (r_start < r_end) {
      int64_t shard_end_bp = (current_shard + 1) * target;
      int64_t this_end =
          std::min(r_end, r_start + (shard_end_bp - accumulated));
      if (this_end <= r_start) this_end = r_end;  // safety

      if (current_shard == task_id) {
        nucleus::genomics::v1::Range chunk;
        chunk.set_reference_name(r.reference_name());
        chunk.set_start(r_start);
        chunk.set_end(this_end);
        shard_regions.push_back(std::move(chunk));
      }

      accumulated += this_end - r_start;
      r_start = this_end;

      if (accumulated >= (current_shard + 1) * target) {
        ++current_shard;
        if (current_shard > task_id) return shard_regions;
      }
    }
  }
  return shard_regions;
}

std::vector<nucleus::genomics::v1::Range> PartitionRegions(
    const std::vector<nucleus::genomics::v1::Range>& calling_regions,
    int64_t partition_size) {
  std::vector<nucleus::genomics::v1::Range> out;
  if (partition_size <= 0) return calling_regions;
  for (const auto& r : calling_regions) {
    int64_t s = r.start();
    while (s < r.end()) {
      const int64_t e = std::min<int64_t>(s + partition_size, r.end());
      nucleus::genomics::v1::Range chunk;
      chunk.set_reference_name(r.reference_name());
      chunk.set_start(s);
      chunk.set_end(e);
      out.push_back(std::move(chunk));
      s = e;
    }
  }
  return out;
}

}  // namespace deepvariant
