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

#include "third_party/nucleus/util/math.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>


#include "tensorflow/core/platform/test.h"

namespace nucleus {

using ::testing::Eq;
using ::testing::DoubleEq;
using ::testing::DoubleNear;
using ::testing::ElementsAreArray;

constexpr double TOL = 1e-4;

TEST(PhredToPError, HandlesValidInputs) {
  EXPECT_THAT(PhredToPError(0), Eq(1));
  EXPECT_THAT(PhredToPError(10), Eq(0.1));
  EXPECT_THAT(PhredToPError(20), Eq(0.01));
  EXPECT_THAT(PhredToPError(30), Eq(0.001));
  EXPECT_THAT(PhredToPError(40), Eq(0.0001));
}

TEST(PhredToLog10PError, HandlesValidInputs) {
  EXPECT_THAT(PhredToLog10PError(0), DoubleEq(0));
  EXPECT_THAT(PhredToLog10PError(10), DoubleEq(-1.0));
  EXPECT_THAT(PhredToLog10PError(20), DoubleEq(-2.0));
  EXPECT_THAT(PhredToLog10PError(30), DoubleEq(-3.0));
  EXPECT_THAT(PhredToLog10PError(40), DoubleEq(-4.0));
}

TEST(PErrorToPhred, HandlesValidInputs) {
  EXPECT_THAT(PErrorToPhred(0.1), DoubleNear(10.0, TOL));
  EXPECT_THAT(PErrorToPhred(0.01), DoubleNear(20.0, TOL));
  EXPECT_THAT(PErrorToPhred(0.001), DoubleNear(30.0, TOL));
  EXPECT_THAT(PErrorToPhred(0.0001), DoubleNear(40.0, TOL));
  EXPECT_THAT(PErrorToPhred(0.00015), DoubleNear(38.2391874, TOL));
}

TEST(PErrorToRoundedPhred, HandlesValidInputs) {
  EXPECT_THAT(PErrorToRoundedPhred(0.1), Eq(10));
  EXPECT_THAT(PErrorToRoundedPhred(0.01), Eq(20));
  EXPECT_THAT(PErrorToRoundedPhred(0.001), Eq(30));
  EXPECT_THAT(PErrorToRoundedPhred(0.0001), Eq(40));
  EXPECT_THAT(PErrorToRoundedPhred(0.00015), Eq(38));
}

TEST(PErrorToLog10PError, HandlesValidInputs) {
  EXPECT_THAT(PErrorToLog10PError(0.1), DoubleNear(-1.0, TOL));
  EXPECT_THAT(PErrorToLog10PError(0.01), DoubleNear(-2.0, TOL));
  EXPECT_THAT(PErrorToLog10PError(0.001), DoubleNear(-3.0, TOL));
  EXPECT_THAT(PErrorToLog10PError(0.0001), DoubleNear(-4.0, TOL));
}

TEST(Log10PErrorToPhred, HandlesValidInputs) {
  EXPECT_THAT(Log10PErrorToPhred(0.0), DoubleEq(0.0));
  EXPECT_THAT(Log10PErrorToPhred(-1.0), DoubleEq(10.0));
  EXPECT_THAT(Log10PErrorToPhred(-1.52), DoubleEq(15.2));
  EXPECT_THAT(Log10PErrorToPhred(-2.0), DoubleEq(20.0));
  EXPECT_THAT(Log10PErrorToPhred(-3.0), DoubleEq(30.0));
  EXPECT_THAT(Log10PErrorToPhred(-4.0), DoubleEq(40.0));
}

TEST(Log10PErrorToRoundedPhred, HandlesValidInputs) {
  EXPECT_THAT(Log10PErrorToRoundedPhred(0.0), Eq(0));
  EXPECT_THAT(Log10PErrorToRoundedPhred(-1.0), Eq(10));
  EXPECT_THAT(Log10PErrorToRoundedPhred(-1.52), Eq(15));
  EXPECT_THAT(Log10PErrorToRoundedPhred(-2.0), Eq(20));
  EXPECT_THAT(Log10PErrorToRoundedPhred(-3.0), Eq(30));
  EXPECT_THAT(Log10PErrorToRoundedPhred(-4.0), Eq(40));
}

TEST(Log10PErrorToRoundedPhred, RoundsInsteadOfFloors) {
  EXPECT_THAT(Log10PErrorToRoundedPhred(-0.000585209927521646), Eq(0));
}

/* R expectations
expected <- function(log10p) {
  paste(-10 * log10(1 - 10^log10p))
}
expected(-0.01)
expected(-0.25)
expected(-0.5)
expected(-1)
expected(-2)
expected(-3)

expected(-1e-9)
expected(-1e-15)
expected(-1e-20)
*/
TEST(Log10PTrueToPhred, HandlesInputsBelowMax) {
  EXPECT_THAT(Log10PTrueToPhred(-0.01, 1000),
              DoubleNear(16.4277471723837, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-0.25, 1000),
              DoubleNear(3.58864458983019, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-0.5, 1000),
              DoubleNear(1.6508853862677, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-1, 1000),
              DoubleNear(0.457574905606751, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-2, 1000),
              DoubleNear(0.0436480540245009, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-3, 1000),
              DoubleNear(0.00434511774017692, 1e-5));
}

TEST(Log10PTrueToPhred, HandlesInfiniteValues) {
  EXPECT_THAT(Log10PTrueToPhred(-1e-9, 1000.0),
              DoubleNear(86.3778430572196, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-1e-15, 1000.0),
              DoubleNear(146.323704754571, 1e-5));
  EXPECT_THAT(Log10PTrueToPhred(-1e-20, 1000.0),
              DoubleEq(1000.0));
}

TEST(Log10ToReal, HandlesValidInputs) {
  EXPECT_THAT(Log10ToReal(0.0), Eq(1.0));
  EXPECT_THAT(Log10ToReal(-1.0), Eq(0.1));
  EXPECT_THAT(Log10ToReal(-2.0), Eq(0.01));
  EXPECT_THAT(Log10ToReal(-3.0), Eq(0.001));
  EXPECT_THAT(Log10ToReal(-4.0), Eq(0.0001));
}

TEST(ZeroShiftLikelihoods, HandlesValidInputs) {
  std::vector<double> test_data {-3.0, -100.7, -88.0};
  EXPECT_THAT(ZeroShiftLikelihoods(test_data),
              ElementsAreArray({0.0, -97.7, -85.0}));
}

}  // namespace nucleus
