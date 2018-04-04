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

#ifndef THIRD_PARTY_NUCLEUS_UTIL_SAMPLERS_H_
#define THIRD_PARTY_NUCLEUS_UTIL_SAMPLERS_H_

#include "tensorflow/core/lib/random/simple_philox.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/types.h"

using tensorflow::uint64;

namespace nucleus {

// Helper class for randomly sampling a fraction of values.
//
// API is simple: only a fraction_to_keep calls to Keep() will return true.
//
// So keeping a 10% fraction of the values in a vector<int> x is:
//
// FractionalSampler sampler(0.10, rng);
// for( int v : x ) {
//  if (sampler.Keep()) {
//    ...
//  }
//
// This is a templated class so it can use any Random Number Generator (RNG)
// that (1) can be copied and (2) provides a RandDouble() function that emits
// a value between 0.0 and 1.0.
template<typename RNG>
class FractionalSampler {
 public:
  // Creates a new FractionalSampler that keeps fraction_to_keep elements on
  // average among N calls to Keep(). RNG is be a random number generator
  // object that can be copied into a member variable of this class.
  explicit FractionalSampler(double fraction_to_keep, const RNG& rng)
      : fraction_to_keep_(fraction_to_keep), rng_(rng) {
    CHECK_GE(fraction_to_keep, 0.0) << "Must be between 0.0 and 1.0";
    CHECK_LE(fraction_to_keep, 1.0) << "Must be between 0.0 and 1.0";
  }

  // Randomly return true approximately fraction_to_keep of the time.
  bool Keep() const { return rng_.RandDouble() <= fraction_to_keep_; }

  // Gets the fraction of elements that will be kept.
  double FractionKept() const { return fraction_to_keep_; }

 private:
  const double fraction_to_keep_;
  mutable RNG rng_;
};

class PhiloxFractionalSampler
    : public FractionalSampler<tensorflow::random::SimplePhilox> {
 public:
  // Creates a new PhiloxFractionalSampler that keeps fraction_to_keep elements
  // on average among N calls to Keep() using a PhiloxRandom initialized with
  // with random_seed.
  explicit PhiloxFractionalSampler(double fraction_to_keep, uint64 random_seed)
      : FractionalSampler(fraction_to_keep,
                          tensorflow::random::SimplePhilox(&rng_gen_)),
        rng_gen_(random_seed) {}

  // Fraction of non-variant sites to emit as DeepVariantCalls.
  mutable tensorflow::random::PhiloxRandom rng_gen_;
};

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_UTIL_SAMPLERS_H_
