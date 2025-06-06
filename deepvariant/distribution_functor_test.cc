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
#include "deepvariant/distribution_functor.h"

#include <functional>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace learning::genomics::deepvariant::distribution_functor {

using ::testing::Pair;
using ::testing::UnorderedElementsAre;

TEST(DistributionFunctorTest, DistributionConstructionAndAccessors) {
  Distribution<int>::WeightMap weight_map = {{1, 2}, {2, 3}, {3, 1}};
  Distribution<int> dist = Distribution<int>::FromWeightMap(weight_map);

  EXPECT_EQ(dist.get_weight_map(), weight_map);
  EXPECT_EQ(dist.get_total_weight(), 6);
}

TEST(DistributionFunctorTest, UnitFactory) {
  auto dist = unit(5);
  EXPECT_THAT(dist.get_weight_map(), UnorderedElementsAre(Pair(5, 1)));
  EXPECT_EQ(dist.get_total_weight(), 1);
}

TEST(DistributionFunctorTest, UniformFactory) {
  auto dist = uniform<int>({1, 2, 3});
  EXPECT_THAT(dist.get_weight_map(),
              UnorderedElementsAre(Pair(1, 1), Pair(2, 1), Pair(3, 1)));
  EXPECT_EQ(dist.get_total_weight(), 3);
}

TEST(DistributionFunctorTest, DistributionEquality) {
  auto dist1 = Distribution<int>::FromWeightMap({{1, 2}, {2, 3}});
  auto dist2 = Distribution<int>::FromWeightMap({{1, 2}, {2, 3}});
  auto dist3 = Distribution<int>::FromWeightMap({{1, 2}, {3, 3}});

  EXPECT_EQ(dist1, dist2);
  EXPECT_NE(dist1, dist3);
}

TEST(DistributionFunctorTest, DistMapSimple) {
  auto dist = uniform<int>({1, 2, 3});
  auto mapped_dist = dist_map(dist, [](int x) { return x * 2; });
  EXPECT_THAT(mapped_dist.get_weight_map(),
              UnorderedElementsAre(Pair(2, 1), Pair(4, 1), Pair(6, 1)));
  EXPECT_EQ(mapped_dist.get_total_weight(), 3);
}

TEST(DistributionFunctorTest, DistMapIdentity) {
  auto dist = uniform<int>({1, 2, 3});
  auto mapped_dist = dist_map(dist, [](int x) { return x; });
  EXPECT_EQ(mapped_dist, dist);
}

TEST(DistributionFunctorTest, DistMapDifferentCardinality) {
  auto dist = uniform<int>({1, 2});
  auto mapped_dist = dist_map(dist, [](int x) {
    if (x == 1) {
      return "one";
    } else {
      return "two_a_two_b";
    }
  });
  EXPECT_THAT(mapped_dist.get_weight_map(),
              UnorderedElementsAre(Pair("one", 1), Pair("two_a_two_b", 1)));
  EXPECT_EQ(mapped_dist.get_total_weight(), 2);
}

TEST(DistributionFunctorTest, DistBindSimple) {
  auto dist = uniform<int>({1, 2});
  auto bound_dist =
      dist_bind(dist, [](int x) { return uniform<int>({x, x + 1}); });
  EXPECT_THAT(bound_dist.get_weight_map(),
              UnorderedElementsAre(Pair(1, 1), Pair(2, 2), Pair(3, 1)));
  EXPECT_EQ(bound_dist.get_total_weight(), 4);
}

TEST(DistributionFunctorTest, DistBindDifferentCardinality) {
  auto dist = uniform<int>({1, 2});
  auto bound_dist = dist_bind(dist, [](int x) {
    if (x == 1) {
      return uniform<std::string>({"a", "b", "c"});
    } else {
      return uniform<std::string>({"d"});
    }
  });
  EXPECT_THAT(bound_dist.get_weight_map(),
              UnorderedElementsAre(Pair("a", 1), Pair("b", 1), Pair("c", 1),
                                   Pair("d", 3)));
  EXPECT_EQ(bound_dist.get_total_weight(), 6);
}

TEST(DistributionFunctorTest, DistributionGeneratorSingleParam) {
  DistributionGenerator<int, int> gen(
      [](int x) { return uniform<int>({x, x + 1}); });
  auto traced_dist = dist_map(gen, [](std::function<int(int)> provider) {
    return provider(1) + provider(2);
  });
  EXPECT_THAT(traced_dist.get_weight_map(),
              UnorderedElementsAre(Pair(3, 1), Pair(4, 2), Pair(5, 1)));
  EXPECT_EQ(traced_dist.get_total_weight(), 4);
}

TEST(DistributionFunctorTest, DistributionGeneratorMultiParam) {
  DistributionGenerator<int, int, int> gen(
      [](int x, int y) { return uniform<int>({x + y}); });
  auto traced_dist = dist_map(gen, [](std::function<int(int, int)> provider) {
    return provider(1, 2) + provider(2, 1);
  });
  EXPECT_THAT(traced_dist.get_weight_map(), UnorderedElementsAre(Pair(6, 1)));
  EXPECT_EQ(traced_dist.get_total_weight(), 1);
}

TEST(DistributionFunctorTest, DistributionGeneratorDifferentReturnType) {
  DistributionGenerator<std::string, int> gen(
      [](int x) { return uniform<std::string>({"a", "b"}); });
  auto traced_dist =
      dist_map(gen, [](std::function<std::string(int)> provider) {
        return provider(1) + provider(2);
      });
  EXPECT_THAT(traced_dist.get_weight_map(),
              UnorderedElementsAre(Pair("aa", 1), Pair("ab", 1), Pair("ba", 1),
                                   Pair("bb", 1)));
  EXPECT_EQ(traced_dist.get_total_weight(), 4);
}

TEST(DistributionFunctorTest, DistributionGeneratorNonUniformSingleParam) {
  // Generator that returns a non-uniform distribution based on input.
  DistributionGenerator<int, int> gen([](int x) {
    return Distribution<int>::FromWeightMap({{x, 2}, {x + 1, 1}});
  });
  auto traced_dist = dist_map(gen, [](std::function<int(int)> provider) {
    return provider(1) + provider(2);
  });
  // When x=1, we have {1:2, 2:1}, and when x=2 we have {2:2, 3:1}.
  // The possible outcomes are:
  // 1+2 = 3 (2*2 = 4 times)
  // 1+3 = 4 (2*1 = 2 times)
  // 2+2 = 4 (1*2 = 2 times)
  // 2+3 = 5 (1*1 = 1 time)
  EXPECT_THAT(traced_dist.get_weight_map(),
              UnorderedElementsAre(Pair(3, 4), Pair(4, 4), Pair(5, 1)));
  EXPECT_EQ(traced_dist.get_total_weight(), 9);
}

TEST(DistributionFunctorTest, DistributionGeneratorNonUniformMultiParam) {
  // Generator that returns a non-uniform distribution based on two inputs.
  DistributionGenerator<int, int, int> gen([](int x, int y) {
    return Distribution<int>::FromWeightMap(
        {{x + y, 3}, {x + y + 1, 2}, {x + y + 2, 1}});
  });
  auto traced_dist = dist_map(gen, [](std::function<int(int, int)> provider) {
    return provider(1, 1) + provider(1, 2);
  });
  // When x=1, y=1 we have {2:3, 3:2, 4:1}
  // When x=1, y=2 we have {3:3, 4:2, 5:1}
  // The possible outcomes are:
  // 2+3 = 5 (3*3 = 9 times)
  // 2+4 = 6 (3*2 = 6 times)
  // 2+5 = 7 (3*1 = 3 times)
  // 3+3 = 6 (2*3 = 6 times)
  // 3+4 = 7 (2*2 = 4 times)
  // 3+5 = 8 (2*1 = 2 times)
  // 4+3 = 7 (1*3 = 3 times)
  // 4+4 = 8 (1*2 = 2 times)
  // 4+5 = 9 (1*1 = 1 time)
  EXPECT_THAT(traced_dist.get_weight_map(),
              UnorderedElementsAre(Pair(5, 9), Pair(6, 12), Pair(7, 10),
                                   Pair(8, 4), Pair(9, 1)));
  EXPECT_EQ(traced_dist.get_total_weight(), 36);
}

}  // namespace learning::genomics::deepvariant::distribution_functor
