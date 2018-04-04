/*
 * Copyright 2018 Google Inc.
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
 *
 */

#include "third_party/nucleus/util/samplers.h"

#include "third_party/nucleus/testing/test_utils.h"
#include "tensorflow/core/lib/random/simple_philox.h"

#include "tensorflow/core/platform/test.h"

namespace nucleus {

using ::testing::DoubleNear;

template<typename RNG>
void VerifySampler(const FractionalSampler<RNG>& sampler, double fraction) {
  int n_kept = 0;
  int n_trials = 1000000;
  for (int i = 0; i < n_trials; ++i) {
    if (sampler.Keep()) {
      n_kept++;
    }
  }
  const double actual_fraction = n_kept / (1.0 * n_trials);
  EXPECT_THAT(actual_fraction, DoubleNear(fraction, 0.001));
}

class FractionalSamplerTest : public ::testing::TestWithParam<double> {};

TEST_P(FractionalSamplerTest, TestFractionalSampler) {
  // Test that the fractional sampler produces approximately fraction * n_trials
  // Keep()=true values over a large number of trials.
  const double fraction = GetParam();
  tensorflow::random::PhiloxRandom gen(123456 /* random seed */);
  tensorflow::random::SimplePhilox rng(&gen);
  FractionalSampler<tensorflow::random::SimplePhilox> sampler(fraction, rng);
  VerifySampler(sampler, fraction);
}

TEST_P(FractionalSamplerTest, TestPhiloxFractionalSampler) {
  // Test that the fractional sampler produces approximately fraction * n_trials
  // Keep()=true values over a large number of trials.
  const double fraction = GetParam();
  PhiloxFractionalSampler sampler(fraction, 123456 /* random seed */);
  VerifySampler(sampler, fraction);
}


INSTANTIATE_TEST_CASE_P(FractionalSamplerTest1, FractionalSamplerTest,
                        ::testing::Values(0.9, 0.1, 0.01, 0.05));

}  // namespace nucleus
