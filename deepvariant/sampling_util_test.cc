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
#include "deepvariant/sampling_util.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <utility>
#include <vector>

#include "deepvariant/distribution_functor.h"
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/container/btree_set.h"

namespace learning::genomics::deepvariant::sampling::internal {

typedef absl::btree_set<int> int_set;

std::vector<int_set> all_subsets(int_set univ, size_t size) {
  std::vector<int> univ_vec(univ.begin(), univ.end());
  std::vector<bool> keep(univ.size(), false);
  for (int i = 0; i != size; i++) {
    keep[i] = true;
  }

  std::vector<int_set> subsets;
  do {
    int_set subset;
    for (int i = 0; i != univ.size(); i++) {
      if (keep[i]) {
        subset.insert(univ_vec[i]);
      }
    }
    subsets.push_back(subset);
  } while (std::prev_permutation(keep.begin(), keep.end()));
  return subsets;
}

TEST(SamplingUtilTest, ReservoirSampleIsUniform) {
  // The population is the set of all integers from 0 to 6.
  int_set population = {0, 1, 2, 3, 4, 5, 6};

  // A function that provides a family of uniform distributions over integers.
  distribution_functor::DistributionGenerator<size_t, size_t>
      uniform_integer_distribution_provider([](size_t max) {
        std::vector<size_t> domain(max + 1);
        std::iota(domain.begin(), domain.end(), 0);
        return distribution_functor::uniform<size_t>(domain);
      });

  // Function under test.
  auto fut = [&population](std::function<size_t(size_t)> index_provider) {
    return ReservoirSampleImpl(3, index_provider, population);
  };

  auto final_distribution = distribution_functor::dist_map(
      uniform_integer_distribution_provider, fut);

  // Check the final distribution is the uniform distribution over all subsets
  // of size 3.
  EXPECT_EQ(final_distribution,
            distribution_functor::uniform<int_set>(all_subsets(population, 3)));
}

int_set complement(const int_set& inner, const int_set& outer) {
  int_set outer_copy(outer.begin(), outer.end());
  absl::erase_if(outer_copy, [&inner](int i) { return inner.contains(i); });
  return outer_copy;
}

typedef std::pair<std::vector<int>, std::vector<int>> vector_pair;

vector_pair get_vector_pair(const int_set& left, const int_set& right) {
  return std::make_pair(std::vector<int>(left.begin(), left.end()),
                        std::vector<int>(right.begin(), right.end()));
}

TEST(SamplingUtilTest, CheckSampleWithPartitionMinsDistribution) {
  // Consider the case of two partitions, each with 3 elements, and we want to
  // sample 4 elements from their union, with at least 1 element from each
  // partition.
  absl::btree_set<int_set> partition = {{0, 1, 2}, {3, 4, 5}};
  size_t sample_size = 4;
  int min_per_partition = 1;
  int_set first = {0, 1, 2};
  int_set second = {3, 4, 5};

  // A function that provides a family of uniform distributions over subsets.
  distribution_functor::DistributionGenerator<int_set, int_set, size_t>
      uniform_subset_distribution_provider([](int_set set, size_t size) {
        return distribution_functor::uniform<int_set>(all_subsets(set, size));
      });

  // The function under test.
  auto fut = [&partition, &sample_size, &min_per_partition](
                 std::function<int_set(int_set, size_t)> subset_provider) {
    return *SampleWithPartitionMinsImpl<int>(
        partition, sample_size, min_per_partition, subset_provider);
  };

  auto final_distribution =
      distribution_functor::dist_map(uniform_subset_distribution_provider, fut);

  // Project the distribution into a bernoulli representing the chance the set
  // is balanced across the partition.
  auto is_balanced_distribution =
      distribution_functor::dist_map(final_distribution, [](auto v) {
        int c = 0;
        for (int i = 0; i < 3; ++i) {
          if (v.contains(i)) {
            c++;
          }
        }
        return (c % 2) == 0;
      });

  // The distribution should be balanced 2/3 of the time.
  EXPECT_EQ(is_balanced_distribution,
            distribution_functor::Distribution<bool>::FromWeightMap(
                {{true, 2}, {false, 1}}));
}

}  // namespace learning::genomics::deepvariant::sampling::internal
