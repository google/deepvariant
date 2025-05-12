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

#include <functional>
#include <numeric>
#include <vector>

#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

namespace learning::genomics::deepvariant::sampling::internal {
using ::testing::_;
using ::testing::Each;
using ::testing::Pair;

bool next_returns(const std::vector<int>& limits, std::vector<int>& current) {
  int index = current.size() - 1;
  while (index >= 0 && current[index] == limits[index]) {
    current[index] = 0;
    index--;
  }
  if (index < 0) {
    return false;
  }
  current[index]++;
  return true;
}

TEST(SamplingUtilTest, InPlaceReservoirSampleIsUniform) {
  // Test picking 3 elements from {0, 1, 2, 3, 4, 5, 6}
  std::vector<int> population(7);
  std::vector<int> returns{0, 0, 0, 0};
  const std::vector<int> limits = {3, 4, 5, 6};

  // We iterate over all possible functions that given n in [3, 6], return an
  // element from [0, n]. And we count how many times each subset occurs as a
  // sample.
  absl::flat_hash_map<absl::flat_hash_set<int>, int> counter;
  do {
    std::iota(population.begin(), population.end(), 0);
    std::function f = [&returns](int i) {
      if (i >= 3 && i < 7) {
        return returns[i - 3];
      } else {
        return 0;
      }
    };
    int sample_size = InPlaceReservoirSampleImpl(3, f, population);
    absl::flat_hash_set<int> sample(population.begin(),
                                    population.begin() + sample_size);
    counter[sample]++;
  } while (next_returns(limits, returns));

  // Choose 3 elements from 7 elements: (7 choose 3) = 35 combinations.
  EXPECT_EQ(counter.size(), 35);
  // Since there are 4*5*6*7 = 35*24 possible functions, each sample should
  // occur exactly 24 times.
  EXPECT_THAT(counter, Each(Pair(_, 24)));
}

}  // namespace learning::genomics::deepvariant::sampling::internal
