# Copyright 2017 Google Inc.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""Mathematics functions.

A quick note on terminology here.

There are a bunch kinds of probabilities used commonly in genomics:

-- pError: the probability of being wrong.
-- pTrue: the probability of being correct.

Normalized probabilities vs. unnormalized likelihoods:

-- Normalized probabilities: p_1, ..., p_n such that sum(p_i) == 1 are said
   said to be normalized because they represent a valid probability
   distribution over the states 1 ... n.
-- Unnormalized likelihoods: p_1, ..., p_n where sum(p_i) != 1. These are not
   normalized and so aren't a valid probabilities distribution.

To add even more complexity, probabilities are often represented in three
semi-equivalent spaces:

-- Real-space: the classic space, with values ranging from [0.0, 1.0]
   inclusive.
-- log10-space: If p is the real-space value, in log10-space this would be
   represented as log10(p). How the p == 0 case is handled is often function
   dependent, which may accept/return -Inf or not handle the case entirely.
-- Phred-scaled: See https://en.wikipedia.org/wiki/Phred_quality_score for
   more information. But briefly, the Phred-scale maintains resolution in the
   lower parts of the probability space using integer quality scores (though
   using ints is optional, really). The phred-scale is defined as

     `phred(p) = -10 * log10(p)`

   where p is a real-space probability.

The code in math.h dealing with probabilities is very explicit about what
kind probability and representation is expects and returns, as unfortunately
these are all represented commonly as doubles in C++. Though tempting to
address this issue with classic software engineering practices like creating
a Probability class, in practice this is extremely difficult to do as this
code is often performance critical and the low-level mathematical operations
used in this code (e.g., log10) don't distiguish themselves among the types
of probabilities.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



import numpy as np
import scipy.misc
import scipy.stats

from deepvariant.core.python import math as math_

# Scipy dependency:
#
# Note that we depend on scipy here, which is a very heavy-weight dependency
# if it needs to be built from source. scipy is liberally licensed:
# https://github.com/scipy/scipy/blob/v0.19.0/LICENSE.txt
# So we might consider vendoring in a few functions from that package to avoid
# the dependency. Here are links to source for a few that we use:
#
# logsumexp [relatively simple]:
# https://github.com/scipy/scipy/blob/v0.19.0/scipy/special/_logsumexp.py
#
# logpmf [much more complex]:
# https://github.com/scipy/scipy/blob/v0.19.0/scipy/stats/_discrete_distns.py

# C++ CLIF functions:
#
# We are enumerating the C++ functions exported by python/math.clif here, so
# it's clear to people what functions are available in python without digging
# into the raw python/C++ CLIF code.
log10_ptrue_to_phred = math_.log10_ptrue_to_phred
phred_to_perror = math_.phred_to_perror
phred_to_log10_perror = math_.phred_to_log10_perror
perror_to_log10_perror = math_.perror_to_log10_perror
perror_to_phred = math_.perror_to_phred
log10_perror_to_phred = math_.log10_perror_to_phred
perror_to_rounded_phred = math_.perror_to_rounded_phred
log10_perror_to_rounded_phred = math_.log10_perror_to_rounded_phred
log10_perror_to_perror = math_.log10_perror_to_perror
zero_shift_log10_probs = math_.zero_shift_log10_probs

# Maximum confidence in a variant call. Used to prevent overflow with log10.
_MAX_CONFIDENCE = 1.0 - 1e-15

LOG_10_OF_E = np.log10(np.e)
LOG_E_OF_10 = np.log(10.0)


def perror_to_bounded_log10_perror(perror, min_prob=1.0 - _MAX_CONFIDENCE):
  """Computes log10(p) for the given probability.

  The log10 probability is capped by -_MAX_CONFIDENCE.

  Args:
    perror: the probability to log10
    min_prob: minimum allowed probability

  Returns:
    log10(p).

  Raises:
    ValueError: If probability is outside of [0.0, 1.0].
  """
  if not 0 <= perror <= 1:
    raise ValueError('perror must be between zero and one:', perror)
  return perror_to_log10_perror(max(perror, min_prob))


def ptrue_to_bounded_phred(ptrue, max_prob=_MAX_CONFIDENCE):
  """Computes the phred scaled confidence from the given ptrue probability.

  See https://en.wikipedia.org/wiki/Phred_quality_score for more information.
  The quality score is capped by _MAX_CONFIDENCE.

  Args:
    ptrue: the ptrue probability to phred scale.
    max_prob: maximum allowed probability

  Returns:
    Phred scaled version of the 1 - ptrue.

  Raises:
    ValueError: If ptrue is outside of [0.0, 1.0].
  """
  if not 0 <= ptrue <= 1:
    raise ValueError('ptrue must be between zero and one:', ptrue)
  return perror_to_phred(1.0 - min(ptrue, max_prob))


def log10_binomial(k, n, p):
  """Calculates numerically-stable value of log10(binomial(k, n, p)).

  Returns the log10 of the binomial density for k successes in n trials where
  each success has a probability of occurring of p.

  In real-space, we would calculate:

     result = (n choose k) * (1-p)^(n-k) * p^k

  This function computes the log10 of result, which is:

     log10(result) = log10(n choose k) + (n-k) * log10(1-p) + k * log10(p)

  This is equivalent to invoking the R function:
    dbinom(x=k, size=n, prob=p, log=TRUE)

  See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Binomial.html
  for more details on the binomial.

  Args:
    k: int >= 0. Number of successes.
    n: int >= 0. Number of trials.
    p: 0.0 <= float <= 1.0. Probability of success.

  Returns:
    log10 probability of seeing k successes in n trials with p.
  """
  return scipy.stats.binom.logpmf(k, n, p) / LOG_E_OF_10


def log10sumexp(log10_probs):
  """Returns log10(sum(10^log10_probs)) computed in a numerically-stable way.

  Args:
    log10_probs: array-like of floats. An array of log10 probabilties.

  Returns:
    Float.
  """
  return LOG_10_OF_E * scipy.misc.logsumexp(np.array(log10_probs) * LOG_E_OF_10)


def normalize_log10_probs(log10_probs):
  """Approximately normalizes log10 probabilities.

  This function normalizes log10 probabilities. What this means is that we
  return an equivalent array of probabilities but whereas sum(10^log10_probs) is
  not necessarily 1.0, the resulting array is sum(10^result) ~= 1.0. The ~=
  indicates that the result is not necessarily == 1.0 but very close.

  This function is a fast and robust approximation of the true normalization of
  a log10 transformed probability vector. To understand the approximation,
  let's start with the exact calculation. Suppose I have three models, each
  emitting a probability that some data was generated by that model:

    data = {0.1, 0.01, 0.001} => probabilities from models A, B, and C

  These probabilities are unnormalized, in the sense that the total probability
  over the vector doesn't sum to 1 (sum(data) = 0.111). In many applications we
  want to normalize this vector so that sum(normalized(data)) = 1 and the
  relative magnitudes of the original probabilities are preserved (i.e,:

    data[i] / data[j] = normalized(data)[i] / normalized(data)[j]

  for all pairs of values indexed by i and j. For much of the work we do in
  genomics, we have so much data that representing these raw probability
  vectors in real-space risks numeric underflow/overflow, so we instead
  represent our probability vectors in log10 space:

    log10_data = log10(data) = {-1, -2, -3}

  Given that we expect numeric problems in real-space, normalizing this log10
  vector is hard, because the standard way you'd do the normalization is via:

    data[i] = data[i] / sum(data)
    log10_data[i] = log10_data[i] - log10(sum(10^data))

  But computing the sum of log10 values this way is dangerous because the naive
  implementation converts back to real-space to do the sum, the very operation
  we're trying to avoid due to numeric instability.

  This function implements an approximate normalization, which relaxes the need
  for an exact calculation of the sum. This function ensures that the
  normalization is numerically safe at the expense of the sum not being exactly
  equal to 1 but rather just close.

  Args:
    log10_probs: array-like of floats. An array of log10 probabilties.

  Returns:
    np.array with the same shape as log10_probs but where sum(10^result) ~= 1.0.

  Raises:
    ValueError: if any log10_probs > 0.0
  """
  log10_probs = np.array(log10_probs)
  if np.max(log10_probs) > 0.0:
    raise ValueError('log10_probs all must be <= 0', log10_probs)
  lse = log10sumexp(log10_probs)
  # np.minimum protects us from producing values slightly > 0.0 (e.g., 1e-16).
  return np.minimum(log10_probs - lse, 0.0)
