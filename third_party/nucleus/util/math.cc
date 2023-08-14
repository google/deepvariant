/*
 * Copyright 2018 Google LLC.
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

#include <algorithm>
#include <cmath>

#include "absl/log/check.h"
#include "third_party/nucleus/util/math.h"

namespace nucleus {

double PhredToPError(const int phred) {
  CHECK_GE(phred, 0);
  return std::pow(10.0, -static_cast<double>(phred) / 10.0);
}

double PhredToLog10PError(const int phred) {
  CHECK_GE(phred, 0);
  return -static_cast<double>(phred) / 10;
}

double PErrorToPhred(const double perror) {
  CHECK_GT(perror, 0);
  CHECK_LE(perror, 1);
  return Log10PErrorToPhred(PErrorToLog10PError(perror));
}

int PErrorToRoundedPhred(const double perror) {
  CHECK_GT(perror, 0);
  CHECK_LE(perror, 1);
  return Log10PErrorToRoundedPhred(PErrorToLog10PError(perror));
}

double PErrorToLog10PError(const double perror) {
  CHECK_GT(perror, 0);
  CHECK_LE(perror, 1);
  return std::log10(perror);
}

double Log10PErrorToPhred(const double log10_perror) {
  CHECK_LE(log10_perror, 0);
  return -10 * log10_perror;
}

int Log10PErrorToRoundedPhred(const double log10_perror) {
  return std::abs(std::round(Log10PErrorToPhred(log10_perror)));
}

double Log10PTrueToPhred(const double log10_ptrue,
                         const double value_if_not_finite) {
  const double ptrue = Log10ToReal(log10_ptrue);
  const double perror = std::log10(1 - ptrue);
  return std::isfinite(perror) ? -10 * perror : value_if_not_finite;
}

double Log10ToReal(const double log10_probability) {
  CHECK_LE(log10_probability, 0.0);
  return std::pow(10, log10_probability);
}

std::vector<double> ZeroShiftLikelihoods(
    const std::vector<double>& likelihoods) {
  std::vector<double> normalized(likelihoods.size());
  double max = *std::max_element(likelihoods.cbegin(), likelihoods.cend());
  std::transform(likelihoods.cbegin(), likelihoods.cend(),
                 normalized.begin(),
                 [max](double x) { return x - max; });
  return normalized;
}

}  // namespace nucleus
