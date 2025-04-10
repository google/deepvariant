/*
 * Copyright 2025 Google LLC.
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

#include "deepvariant/methylation_aware_phasing.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "boost/math/distributions/normal.hpp"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

const float kPThreshold = 0.05;
const float kRankSumVarianceDenominator = 12.0;

// This function checks whether a methylated site is differentially methylated
// between haplotypes.
//
// It checks if the site's methylation values from haplotype 1 and
// haplotype 2 originate from significantly different distributions. It is used
// to identify sites where methylation is haplotype-specific, which can be
// leveraged to phase previously unphased reads based on their methylation
// patterns. It returns true if the p-value of the Wilcoxon Rank-Sum Test is
// less than the p-value threshold.
bool IsDifferentiallyMethylated(const std::vector<double>& hap1_methyl,
                                const std::vector<double>& hap2_methyl) {
  double p_val = WilcoxonRankSumTest(hap1_methyl, hap2_methyl);
  return p_val >= 0 && p_val < kPThreshold;
}

// Performs two-sided Wilcoxon Rank-Sum Test (Mann-Whitney U Test) between two
// haplotypes’ methylation values. The calculated p-value is returned. Returns
// -1.0 if either haplotype input vector is empty, in which case the test cannot
// be performed. This may occur if no reads were phased to a given haplotype at
// this site.
double WilcoxonRankSumTest(const std::vector<double>& hap1_methyl,
                           const std::vector<double>& hap2_methyl) {
  size_t n1 = hap1_methyl.size(), n2 = hap2_methyl.size();
  if (n1 == 0 || n2 == 0) {
    return -1.0;
  }

  struct RankedValue {
    double value;
    int group;  // 0 for hap1, 1 for hap2
  };

  std::vector<RankedValue> combined;
  combined.reserve(n1 + n2);

  for (double val : hap1_methyl) combined.push_back({val, 0});
  for (double val : hap2_methyl) combined.push_back({val, 1});

  std::sort(combined.begin(), combined.end(),
            [](const RankedValue& a, const RankedValue& b) {
              return a.value < b.value;
            });

  std::vector<double> ranks(combined.size());
  for (size_t i = 0; i < combined.size(); ++i) {
    size_t j = i;
    while (j + 1 < combined.size() &&
           combined[j + 1].value == combined[i].value) {
      ++j;
    }
    double avg_rank = (i + j + 2) / 2.0;  // ranks are 1-based
    for (size_t k = i; k <= j; ++k) {
      ranks[k] = avg_rank;
    }
    i = j;
  }

  double rank_sum_1 = 0.0;
  for (size_t i = 0; i < combined.size(); ++i) {
    if (combined[i].group == 0) {
      rank_sum_1 += ranks[i];
    }
  }

  double U1 = rank_sum_1 - (n1 * (n1 + 1)) / 2.0;
  double U2 = n1 * n2 - U1;
  double U = std::min(U1, U2);

  double mean_U = (n1 * n2) / 2.0;
  double std_U =
      std::sqrt((n1 * n2 * (n1 + n2 + 1)) / kRankSumVarianceDenominator);
  double z = (U - mean_U) / std_U;

  boost::math::normal_distribution<> normal_dist(0.0, 1.0);
  double p_value = 2 * (1 - boost::math::cdf(normal_dist, std::fabs(z)));

  return p_value;
}

// Assigns an unphased read to a haplotype based on methylation similarity.
//
// For each informative methylated reference site, this function compares the
// unphased read's methylation level to the average methylation levels of
// haplotype 1 and haplotype 2 reads at that site. It tallies a "vote" for the
// haplotype whose average methylation is closer to that of the unphased read.
// After all sites are evaluated, the haplotype with the majority vote (with a
// minimum of 3 votes) is returned as the predicted haplotype for the read.
//
// Args:
//   unphased_read: A DeepVariantCall_ReadSupport proto representing the
//                  unphased read.
//   informative_calls: Vector of DeepVariantCall objects containing methylated
//                      reference sites that show differential methylation.
//   hap1_reads: Pointers to ReadSupport entries currently assigned to
//                haplotype 1.
//   hap2_reads: Pointers to ReadSupport entries currently assigned to
//                haplotype 2.
//
// Returns:
//   1 if the read is assigned to haplotype 1,
//   2 if assigned to haplotype 2,
//   or 0 if the read cannot be confidently assigned to either haplotype.
int HaplotypeVoteWithMethylation(
  const DeepVariantCall_ReadSupport& unphased_read,
  const std::vector<DeepVariantCall>& informative_calls,
  const std::vector<const DeepVariantCall_ReadSupport*>& hap1_reads,
  const std::vector<const DeepVariantCall_ReadSupport*>& hap2_reads) {
  const std::string unphased_name = unphased_read.read_name();
  int hap1_votes = 0, hap2_votes = 0;

  for (const auto& call : informative_calls) {
    // Get methylation level for the current site's unphased read
    double read_methyl = -1.0;
    for (const auto& support : call.ref_support_ext().read_infos()) {
      if (support.read_name() == unphased_name) {
        read_methyl = GetMethylationLevelAtSite(support);
        break;
      }
    }
    if (read_methyl < 0) continue;  // no valid methylation info for this read

    // Gather hap1/hap2 methylation values at this site
    absl::flat_hash_map<std::string, double> hap1_map, hap2_map;
    for (const auto& support : call.ref_support_ext().read_infos()) {
      double m = GetMethylationLevelAtSite(support);
      if (m < 0) continue;
      std::string name = support.read_name();
      if (std::find_if(hap1_reads.begin(), hap1_reads.end(),
                       [&](const auto* r) { return r->read_name() == name; }) !=
          hap1_reads.end()) {
        hap1_map[name] = m;
      } else if (std::find_if(hap2_reads.begin(), hap2_reads.end(),
                              [&](const auto* r) {
                                return r->read_name() == name;
                              }) != hap2_reads.end()) {
        hap2_map[name] = m;
      }
    }

    std::vector<double> hap1_methyl, hap2_methyl;
    for (const auto& p : hap1_map) hap1_methyl.push_back(p.second);
    for (const auto& p : hap2_map) hap2_methyl.push_back(p.second);
    if (hap1_methyl.empty() || hap2_methyl.empty()) continue;

    // Compare the unphased read's methylation level to the closest haplotype
    double hap1_mean =
        std::accumulate(hap1_methyl.begin(), hap1_methyl.end(), 0.0) /
        hap1_methyl.size();
    double hap2_mean =
        std::accumulate(hap2_methyl.begin(), hap2_methyl.end(), 0.0) /
        hap2_methyl.size();

    if (std::abs(read_methyl - hap1_mean) < std::abs(read_methyl - hap2_mean)) {
      ++hap1_votes;
    } else {
      ++hap2_votes;
    }
  }

  if (hap1_votes >=3 && hap1_votes > hap2_votes) return 1;
  if (hap2_votes >= 3 && hap2_votes > hap1_votes) return 2;
  return 0;
}

// Returns the normalized methylation level of a read.
//
// The methylation level is stored as an integer from 0 to 255 in the
// `methylation_level` field of the DeepVariantCall_ReadSupport proto.
// This function normalizes it to a probability in the range [0.0, 1.0].
//
// A value of 0 indicates either no methylation or that the field is unset.
// To avoid treating unset values as valid data, this function returns -1.0
// when the methylation level is 0.
//
// Args:
//   read: A DeepVariantCall_ReadSupport proto.
//
// Returns:
//   A normalized methylation level in [0.0, 1.0] or -1.0 if unavailable.
double GetMethylationLevelAtSite(const DeepVariantCall_ReadSupport& read) {
  if (read.methylation_level() == 0) {
    return -1.0;
  }
  return static_cast<double>(read.methylation_level()) / 255.0;
}

// Identifies methylated reference sites that show statistically significant
// differences in methylation levels between reads assigned to haplotype 1 and
// haplotype 2. For each site, the function collects methylation levels from
// hap1/hap2 reads (based on read name matches in the call’s ReadSupport list),
// and runs a Wilcoxon rank-sum test to determine if the distributions differ.
//
// Args:
//   methylated_calls: Vector of DeepVariantCall objects for methylated ref
//                     sites.
//   hap1_reads: Pointers to ReadSupport entries assigned to haplotype 1.
//   hap2_reads: Pointers to ReadSupport entries assigned to haplotype 2.
//
// Returns:
//   A vector of methylated sites that are considered informative
//   for distinguishing haplotypes based on methylation signal.
std::vector<DeepVariantCall> IdentifyInformativeSites(
  const std::vector<DeepVariantCall>& methylated_calls,
  const std::vector<const DeepVariantCall_ReadSupport*>& hap1_reads,
  const std::vector<const DeepVariantCall_ReadSupport*>& hap2_reads) {
  std::vector<DeepVariantCall> informative_sites;
  absl::flat_hash_set<std::string> hap1_names, hap2_names;
  for (const auto* r : hap1_reads) hap1_names.insert(r->read_name());
  for (const auto* r : hap2_reads) hap2_names.insert(r->read_name());

  for (const auto& call : methylated_calls) {
    std::vector<double> hap1_methyl;
    std::vector<double> hap2_methyl;

    for (const auto& support : call.ref_support_ext().read_infos()) {
      std::string name = support.read_name();
      double m = GetMethylationLevelAtSite(support);
      if (m < 0) continue;
      if (hap1_names.contains(name)) {
        hap1_methyl.push_back(m);
      } else if (hap2_names.contains(name)) {
        hap2_methyl.push_back(m);
      }
    }

    // ------------------------- Filtering Block -------------------------
    // Adjust accordingly
    const int total_reads = hap1_methyl.size() + hap2_methyl.size();

    // Skip sites with low coverage on either haplotype
    if (hap1_methyl.size() < 2 || hap2_methyl.size() < 2) {
      continue;
    }

    // Skip sites with low total read support (avoid sparse data)
    if (total_reads < 6) {
      continue;
    }

    // Compute per-haplotype methylation means
    double hap1_mean =
        std::accumulate(hap1_methyl.begin(), hap1_methyl.end(), 0.0) /
        hap1_methyl.size();
    double hap2_mean =
        std::accumulate(hap2_methyl.begin(), hap2_methyl.end(), 0.0) /
        hap2_methyl.size();

    // Skip sites with low methylation mean difference between haplotypes
    if (std::abs(hap1_mean - hap2_mean) < 0.25) {
      continue;
    }

    // Skip sites where methylation is highly variable within a haplotype
    auto stddev = [](const std::vector<double>& values, double mean) {
      double sum_sq = 0.0;
      for (double v : values) sum_sq += (v - mean) * (v - mean);
      return std::sqrt(sum_sq / values.size());
    };
    if (stddev(hap1_methyl, hap1_mean) > 0.2 ||
        stddev(hap2_methyl, hap2_mean) > 0.2) {
      continue;
    }
    // --------------------------------------------------------------------

    if (IsDifferentiallyMethylated(hap1_methyl, hap2_methyl)) {
      informative_sites.push_back(call);
    }
  }
  return informative_sites;
}

// Returns a unique key identifying a read based on its fragment name and
// read number.
//
// This key is used to match reads across different representations (e.g.
// nucleus Read and DeepVariantCall_ReadSupport), and to ensure correct mapping
// of phasing information.
//
// Args:
//   read: A nucleus::genomics::v1::Read object.
//
// Returns:
//   A string in the format "<fragment_name>/<read_number>".
std::string ReadKeyForMethylationAwarePhasing(
    const nucleus::genomics::v1::Read& read) {
  return absl::StrCat(read.fragment_name(), "/", read.read_number());
}

// Extracts reads assigned to a specific haplotype phase.
//
// This function filters the read support objects based on the current phase
// assignments, and returns pointers to the reads with the specified haplotype
// phase.
//
// Args:
//   reads: Vector of DeepVariantCall_ReadSupport objects representing all
//          reads.
//   phases: Vector of phase assignments corresponding to each read
//           (0 = unphased, 1 or 2).
//   target_phase: The haplotype phase to extract (typically 1 or 2).
//
// Returns:
//   A vector of pointers to the reads assigned to the specified phase.
std::vector<const DeepVariantCall_ReadSupport*> ExtractReadsByPhase(
  const std::vector<DeepVariantCall_ReadSupport>& reads,
  const std::vector<int>& phases,
  int target_phase) {
  std::vector<const DeepVariantCall_ReadSupport*> filtered;
  for (size_t i = 0; i < reads.size(); ++i) {
    if (phases[i] == target_phase) {
      filtered.push_back(&reads[i]);
    }
  }
  return filtered;
}

// Performs iterative methylation-aware phasing on a set of reads.
//
// This function assigns haplotype phases to unphased reads based on their
// methylation profiles at informative CpG sites. It uses already-phased reads
// from direct phasing to define the haplotype-specific methylation signatures,
// identifies differentially methylated sites, and votes for each unphased
// read's haplotype assignment by comparing its methylation levels to those
// signatures.
//
// The phasing is done iteratively. In each iteration:
//   1. Reads phased to haplotype 1 and 2 are used to identify informative
//      CpG sites.
//   2. Unphased reads are assigned to the haplotype with more similar
//      methylation patterns.
//   3. The loop terminates when either convergence is reached (no new reads
//      were phased), or a maximum number of iterations is exceeded.
//
// Args:
//   reads_to_phase: Vector of all reads in the region, including phased and
//                   unphased.
//   initial_read_phases: Vector of phase assignments for each read
//                        (0 = unphased) from direct phasing.
//   methylated_ref_sites: List of candidate CpG sites used for
//                        methylation-aware phasing.
//   max_iter: Maximum number of phasing iterations to perform.
//
// Returns:
//   A vector of integer phase calls (0, 1, or 2) for each read, same size and
//   order as reads_to_phase.
std::vector<int> PerformMethylationAwarePhasing(
  const std::vector<nucleus::genomics::v1::Read>& reads_to_phase,
  const std::vector<int>& initial_read_phases,
  const std::vector<DeepVariantCall>& methylated_ref_sites,
  int max_iter) {
  // Build a map of read support for each read.
  std::unordered_map<std::string, DeepVariantCall_ReadSupport> read_support_map;
  for (const auto& call : methylated_ref_sites) {
    for (const auto& support : call.ref_support_ext().read_infos()) {
      read_support_map[support.read_name()] = support;
    }
  }

  std::vector<DeepVariantCall_ReadSupport> full_read_supports;
  full_read_supports.reserve(reads_to_phase.size());
  for (const auto& r : reads_to_phase) {
    std::string key = ReadKeyForMethylationAwarePhasing(r);
    if (read_support_map.count(key)) {
      full_read_supports.push_back(read_support_map[key]);
    } else {
      DeepVariantCall_ReadSupport blank;
      blank.set_read_name(key);
      full_read_supports.push_back(blank);
    }
  }

  std::vector<int> current_phases = initial_read_phases;

  LOG(INFO) << "Starting methylation-aware phasing";
  LOG(INFO) << "Total methylated reference sites: "
            << methylated_ref_sites.size();

  for (int iter = 0; iter < max_iter; ++iter) {
    LOG(INFO) << "--- Iteration " << iter + 1 << " ---";
    auto hap1_reads =
        ExtractReadsByPhase(full_read_supports, current_phases, 1);
    auto hap2_reads =
        ExtractReadsByPhase(full_read_supports, current_phases, 2);

    const int num_phased = hap1_reads.size() + hap2_reads.size();
    const int num_unphased = full_read_supports.size() - num_phased;

    LOG(INFO) << "Phased reads: " << num_phased
              << " (Hap1: " << hap1_reads.size()
              << ", Hap2: " << hap2_reads.size() << ")"
              << ", Unphased: " << num_unphased;

    // If there are no unphased reads, we are done.
    if (num_unphased == 0) {
      LOG(INFO) << "No unphased reads found, exiting methylation-aware phasing";
      break;
    }

    auto informative_sites =
        IdentifyInformativeSites(methylated_ref_sites, hap1_reads, hap2_reads);

    LOG(INFO) << "Informative methylated sites: " << informative_sites.size();

    int newly_phased = 0;
    for (size_t i = 0; i < full_read_supports.size(); ++i) {
      if (current_phases[i] == 0) {
        int vote = HaplotypeVoteWithMethylation(
            full_read_supports[i], informative_sites, hap1_reads, hap2_reads);
        if (vote > 0) {
          current_phases[i] = vote;
          ++newly_phased;
        }
      }
    }

    LOG(INFO) << "Newly phased reads: " << newly_phased;
    if (newly_phased == 0) {
      LOG(INFO) << "No new reads phased, exiting methylation-aware phasing";
      break;
    }
  }

  return current_phases;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
