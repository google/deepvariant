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
"""Tests for Math CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest

from deepvariant.core.python import math


class MathWrapTest(absltest.TestCase):

  def test_one_minus_log10_prob_to_phred(self):
    self.assertAlmostEqual(16.4277471723837,
                           math.log10_ptrue_to_phred(-0.01, 1000))

  def test_phred_to_prob(self):
    self.assertEqual(0.1, math.phred_to_perror(10))

  def test_phred_to_log10_prob(self):
    self.assertEqual(-1, math.phred_to_log10_perror(10))

  def test_prob_to_phred(self):
    self.assertEqual(10.0, math.perror_to_phred(0.1))

  def test_prob_to_rounded_phred(self):
    self.assertEqual(10, math.perror_to_rounded_phred(0.1))

  def test_prob_to_log10_prob(self):
    self.assertEqual(-1, math.perror_to_log10_perror(0.1))

  def test_log10_prob_to_phred(self):
    self.assertEqual(10, math.log10_perror_to_phred(-1))

  def test_log10_prob_to_rounded_phred(self):
    self.assertEqual(10, math.log10_perror_to_rounded_phred(-1))

  def test_log10_prob_to_prob(self):
    self.assertEqual(0.1, math.log10_perror_to_perror(-1))

  def test_zero_shift_log10_probs(self):
    self.assertSequenceEqual([0, -1, -2],
                             math.zero_shift_log10_probs([-1, -2, -3]))


if __name__ == '__main__':
  absltest.main()
