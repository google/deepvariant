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
"""Tests for deepvariant .core.math."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import numpy.testing as npt

from deepvariant.core import math


class MathTests(parameterized.TestCase):

  @parameterized.parameters(
      (0.9000000, None, 10.0),
      (0.9900000, None, 20.0),
      (0.9990000, None, 30.0),
      (0.9999000, None, 40.0),
      (0.9999900, None, 50.0),
      (0.9999990, None, 60.0),
      (0.9999999, None, 70.0),
      # Check that bounding works.
      (0.9999999, 1 - 1e-1, 10.0),
      (0.9999999, 1 - 1e-2, 20.0),
      (0.9999999, 1 - 1e-3, 30.0),
      (0.9999999, 1 - 1e-9, 70.0),
  )
  def test_phred_scale(self, prob, bound, expected):
    if bound:
      actual = math.ptrue_to_bounded_phred(prob, bound)
    else:
      actual = math.ptrue_to_bounded_phred(prob)
    self.assertAlmostEqual(actual, expected, places=6)

  @parameterized.parameters(
      (1.000000, None, 0.0),
      (0.100000, None, -1.0),
      (0.010000, None, -2.0),
      (0.001000, None, -3.0),
      (0.000100, None, -4.0),
      (0.000010, None, -5.0),
      (0.000001, None, -6.0),
      # Check that bounding works.
      (0.000100, 1e-1, -1.0),
      (0.000100, 1e-2, -2.0),
      (0.000100, 1e-3, -3.0),
      (0.000100, 1e-4, -4.0),
      (0.000100, 1e-5, -4.0),
      (0.000100, 1e-6, -4.0),
  )
  def test_log10_prob(self, prob, bound, expected):
    if bound:
      actual = math.perror_to_bounded_log10_perror(prob, bound)
    else:
      actual = math.perror_to_bounded_log10_perror(prob)
    self.assertAlmostEqual(actual, expected, places=6)

  @parameterized.parameters(
      (np.log10(0.900000), -1.0, 10.0),
      (np.log10(0.990000), -1.0, 20.0),
      (np.log10(0.999000), -1.0, 30.0),
      # A huge negative value is handled correctly.
      (-10000000.0, -1.0, 0.0),
      # Check that we can pass in a 0.0 probability and get a good value.
      (0.0, -1.0, -1.0),
      # This probability is still calculated correctly, included for safety.
      (0 - 1e-16, -1.0, 156.53559774527022),
      # Pass in a prob close to one, making sure we get bound value back.
      (0 - 1e-32, -1.0, -1.0),
  )
  def test_log10_ptrue_to_phred(self, prob, value_if_not_finite, expected):
    actual = math.log10_ptrue_to_phred(prob, value_if_not_finite)
    self.assertAlmostEqual(actual, expected, places=6)

  # R code to produce the expectation table.
  # expected <- function(k, n, p) {
  #   pbin <- dbinom(k, n, p, log=T) * log10(exp(1))
  #   likelihoods = paste(sprintf("%.6f", pbin), collapse=", ")
  #   result = paste(k, n, p, pbin, sep=", ")
  #   cat(paste("(", result, "),\n", sep=""))
  # }
  #
  # for (n in c(0, 5, 10)) {
  #   for (k in seq(0, n)) {
  #     for (p in c(0.01, 0.5)) {
  #       expected(k, n, p)
  #     }
  #   }
  # }
  # expected(0, 1000, 0.5)
  # expected(0, 10000, 0.5)
  # expected(100, 10000, 0.5)
  @parameterized.parameters(
      (0, 0, 0.01, 0),
      (0, 0, 0.5, 0),
      (0, 5, 0.01, -0.0218240270122504),
      (0, 5, 0.5, -1.50514997831991),
      (1, 5, 0.01, -1.31848921727378),
      (1, 5, 0.5, -0.806179973983887),
      (2, 5, 0.01, -3.01309441620735),
      (2, 5, 0.5, -0.505149978319906),
      (3, 5, 0.01, -5.0087296108049),
      (3, 5, 0.5, -0.505149978319906),
      (4, 5, 0.01, -7.30539480106643),
      (4, 5, 0.5, -0.806179973983887),
      (5, 5, 0.01, -10),
      (5, 5, 0.5, -1.50514997831991),
      (0, 10, 0.01, -0.0436480540245008),
      (0, 10, 0.5, -3.01029995663981),
      (1, 10, 0.01, -1.03928324862205),
      (1, 10, 0.5, -2.01029995663981),
      (2, 10, 0.01, -2.38170592944426),
      (2, 10, 0.5, -1.35708744286447),
      (3, 10, 0.01, -3.95137239176953),
      (3, 10, 0.5, -0.931118710592187),
      (4, 10, 0.01, -5.70396953768078),
      (4, 10, 0.5, -0.688080661905893),
      (5, 10, 0.01, -7.62042348623071),
      (5, 10, 0.5, -0.608899415858268),
      (6, 10, 0.01, -9.69523992687588),
      (6, 10, 0.5, -0.688080661905893),
      (7, 10, 0.01, -11.9339131701597),
      (7, 10, 0.5, -0.931118710592187),
      (8, 10, 0.01, -14.3555170970296),
      (8, 10, 0.5, -1.35708744286447),
      (9, 10, 0.01, -17.0043648054024),
      (9, 10, 0.5, -2.01029995663981),
      (10, 10, 0.01, -20),
      (10, 10, 0.5, -3.01029995663981),
      (0, 1000, 0.5, -301.029995663981),
      (0, 10000, 0.5, -3010.29995663981),
      (100, 10000, 0.5, -2768.48565263445),
  )
  def test_log10_binomial(self, k, n, p, expected):
    self.assertAlmostEqual(math.log10_binomial(k, n, p), expected)

  @parameterized.parameters(
      ([0], 0.0),
      ([0.0], 0.0),
      ([0.0, -10000.0], 0.0),
      ([-1000.0, -10000.0], -1000.0),
      ([-1, -10, -100], -1.0),
      ([-1, -10, -1], -0.69897),
      ([-1, -1, -1], -0.5228787),
      ([-1, -1, -1, -100], -0.5228787),
      ([-1, -1, -1, -100, -1000], -0.5228787),
      ([-1, -1, -1, -100, -1000, -10000], -0.5228787),
      ([-1, -1, -1, -100, -1000, -10000, -100000], -0.5228787),
  )
  def test_log10sumexp(self, log10_probs, expected):
    self.assertAlmostEqual(math.log10sumexp(log10_probs), expected)

  # R code to compute expected results.
  # expected <- function(lprobs) {
  #   result = lprobs - log10(sum(10^lprobs))
  #   lprob_str = paste("[", paste(sprintf("%.6f", lprobs), collapse=", "),
  #                     "]", sep="")
  #   result_str = paste("[", paste(sprintf("%.6f", result), collapse=", "),
  #                     "]", sep="")
  #   cat(paste("(", lprob_str, ", ", result_str, "),\n", sep=""))
  # }
  #
  # expected(c(0))
  # expected(c(-1, -10))
  # expected(c(-1, -100))
  # expected(c(-1, -1000))
  # expected(c(-1, -2))
  # expected(c(-1, -2, -3))
  # expected(c(-1, -2, -3, -100))
  # expected(c(-1, -2, -100))
  # expected(c(-1, -2, -100, -100))
  @parameterized.parameters(
      ([0.000000], [0.000000]),
      ([-1.000000, -10.000000], [-0.000000, -9.000000]),
      ([-1.000000, -100.000000], [0.000000, -99.000000]),
      ([-1.000000, -1000.000000], [0.000000, -999.000000]),
      ([-1.000000, -2.000000], [-0.041393, -1.041393]),
      ([-1.000000, -2.000000, -3.000000], [-0.045323, -1.045323, -2.045323]),
      ([-1.000000, -2.000000, -3.000000, -100.000000],
       [-0.045323, -1.045323, -2.045323, -99.045323]),
      ([-1.000000, -2.000000, -100.000000], [-0.041393, -1.041393, -99.041393]),
      ([-1.000000, -2.000000, -100.000000, -100.000000],
       [-0.041393, -1.041393, -99.041393, -99.041393]),
  )
  def test_normalize_log10_probs(self, log10_probs, expected):
    npt.assert_allclose(
        math.normalize_log10_probs(log10_probs), expected, atol=1e-6)


if __name__ == '__main__':
  absltest.main()
