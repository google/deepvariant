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
#include <vector>

#include "boost/math/distributions/normal.hpp"

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

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
