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

// Core mathematical routines.
//
// A quick note on terminology here.
//
// There are a bunch kinds of probabilities used commonly in genomics:
//
// -- pError: the probability of being wrong.
// -- pTrue: the probability of being correct.
//
// Normalized probabilities vs. unnormalized likelihoods:
//
// -- Normalized probabilities: p_1, ..., p_n such that sum(p_i) == 1 are said
//    said to be normalized because they represent a valid probability
//    distribution over the states 1 ... n.
// -- Unnormalized likelihoods: p_1, ..., p_n where sum(p_i) != 1. These are not
//    normalized and so aren't a valid probabilities distribution.
//
// To add even more complexity, probabilities are often represented in three
// semi-equivalent spaces:
//
// -- Real-space: the classic space, with values ranging from [0.0, 1.0]
//    inclusive.
// -- log10-space: If p is the real-space value, in log10-space this would be
//    represented as log10(p). How the p == 0 case is handled is often function
//    dependent, which may accept/return -Inf or not handle the case entirely.
// -- Phred-scaled: See https://en.wikipedia.org/wiki/Phred_quality_score for
//    more information. But briefly, the Phred-scale maintains resolution in the
//    lower parts of the probability space using integer quality scores (though
//    using ints is optional, really). The phred-scale is defined as
//
//      `phred(p) = -10 * log10(p)`
//
//    where p is a real-space probability.
//
// The code in math.h dealing with probabilities is very explicit about what
// kind probability and representation is expects and returns, as unfortunately
// these are all represented commonly as doubles in C++. Though tempting to
// address this issue with classic software engineering practices like creating
// a Probability class, in practice this is extremely difficult to do as this
// code is often performance critical and the low-level mathematical operations
// used in this code (e.g., log10) don't distiguish themselves among the types
// of probabilities.
#ifndef THIRD_PARTY_NUCLEUS_UTIL_MATH_H_
#define THIRD_PARTY_NUCLEUS_UTIL_MATH_H_

#include <vector>

namespace nucleus {

// Converts Phred scale to probability scale. Phred value must be >= 0.
double PhredToPError(int phred);

// Converts Phred scale to log10 scale. Phred value must be >= 0.
double PhredToLog10PError(int phred);

// Converts Phred scale to pError probability scale.
// Note: There is no Phred Scale equivalent for PError = 0 (would be
// infinity), so this function does not accept PError == 0.
double PErrorToPhred(double perror);
int PErrorToRoundedPhred(double perror);

// Converts probability space to Log10 space.
// Note: There is no Phred Scale equivalent for probability = 0 (would be
// infinity), so this function does not accept probability == 0.
double PErrorToLog10PError(double perror);

// Converts Log10 scale to Phred scale.
double Log10PErrorToPhred(double log10_perror);
int Log10PErrorToRoundedPhred(double log10_perror);

// Converts a Log10(ptrue) value into a phred-scaled value of 1 - 10^log10p.
//
// This operation is common when you've got a probability of an event occurring,
// p, and you want to emit the Phred-equivalent of it being wrong, which is
// -10 * log10(1 - p). The operation 1 - p can easily underflow, causing the us
// to evaluate log10(0), leading to an infinite value. In that case, the
// function returns value_if_not_finite.
double Log10PTrueToPhred(double log10_ptrue, double value_if_not_finite);

// Converts Log10 scale to real scale.
double Log10ToReal(double log10_probability);

// Takes the maximum value (remember, likelihoods are in log10 space and are all
// negative values) and subtract it from all genotype likelihoods so that the
// most likely likelihood is 0. This gives a bit more resolution in the
// conversion.
std::vector<double> ZeroShiftLikelihoods(
    const std::vector<double>& likelihoods);

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_UTIL_MATH_H_
