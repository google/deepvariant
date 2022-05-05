# Copyright 2018 Google LLC.
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

"""Tests for third_party.nucleus.util.utils."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import collections

from absl.testing import absltest

from absl.testing import parameterized
import numpy as np
import numpy.testing as npt

from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils


class UtilsTest(parameterized.TestCase):

  def test_read_range(self):
    """Tests reads have their ranges calculated correctly."""
    start = 10000001
    read = test_utils.make_read(
        'AAACAG',
        chrom='chrX',
        start=start,
        cigar='2M1I3M',
        quals=range(10, 16),
        name='read1')
    self.assertEqual(
        ranges.make_range('chrX', start, start + 5), utils.read_range(read))
    read = test_utils.make_read(
        'AAACAG',
        chrom='chrX',
        start=start,
        cigar='2M16D3M',
        quals=range(10, 16),
        name='read1')
    self.assertEqual(
        ranges.make_range('chrX', start, start + 5 + 16),
        utils.read_range(read))

  def test_read_end(self):
    """Tests reads have their ends calculated correctly."""
    start = 10000001
    read = test_utils.make_read(
        'AAACAG',
        chrom='chrX',
        start=start,
        cigar='2M1I3M',
        quals=range(10, 16),
        name='read1')
    self.assertEqual(
        start + 5, utils.read_end(read))
    read = test_utils.make_read(
        'AAACAG',
        chrom='chrX',
        start=start,
        cigar='2M16D3M',
        quals=range(10, 16),
        name='read1')
    self.assertEqual(
        start + 5 + 16,
        utils.read_end(read))

  def test_reservoir_sample_length(self):
    """Tests samples have expected length."""
    first_ten_ints = range(10)
    # Test sampling with k > len(iterable).
    self.assertEqual(len(utils.reservoir_sample(first_ten_ints, 11)), 10)
    # Test sampling with k == len(iterable).
    self.assertEqual(len(utils.reservoir_sample(first_ten_ints, 10)), 10)
    # Test sampling with k < len(iterable).
    self.assertEqual(len(utils.reservoir_sample(first_ten_ints, 9)), 9)
    # Test sampling with k == 0.
    self.assertEqual(len(utils.reservoir_sample(first_ten_ints, 0)), 0)
    # Test sampling with k < 0 (bad args).
    with self.assertRaises(ValueError):
      utils.reservoir_sample(first_ten_ints, -1)

  @parameterized.parameters(
      (10, 0),
      (1, 1),
      (10, 1),
      (1, 3),
      (3, 3),
      (6, 3),
      (10, 3),
  )
  def test_reservoir_sample_frequency(self, iterable_size, k):
    """Tests observed frequency is close to expected frequency."""
    # Use a fixed random number so our test is deterministic.
    random = np.random.RandomState(123456789)
    n_replicates = 100000
    counts = collections.Counter(
        item
        for _ in range(n_replicates)
        for item in utils.reservoir_sample(range(iterable_size), k, random))
    expected_frequency = min(k / float(iterable_size), 1.0)
    for c in counts.values():
      observed_frequency = c / float(n_replicates)
      npt.assert_allclose(observed_frequency, expected_frequency, atol=0.01)

  @parameterized.parameters(
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=4, e2=10, expected=False),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=3, e2=10, expected=False),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=2, e2=10, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=1, e2=10, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=0, e2=10, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=0, e2=1, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=0, e2=2, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=0, e2=3, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=1, e2=2, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=1, e2=3, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=2, e2=3, expected=True),
      # dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=3, e2=3, expected=False),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=0, e2=4, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr1', s2=1, e2=4, expected=True),
      dict(ref1='chr1', s1=0, e1=3, ref2='chr2', s2=1, e2=4, expected=False),
  )
  def test_read_overlaps_region(self, ref1, s1, e1, ref2, s2, e2, expected):

    def check_overlaps(chr1, start1, end1, chr2, start2, end2, expected):
      nbp = end1 - start1
      read = test_utils.make_read(
          'A' * nbp, chrom=chr1, start=start1, cigar='{}M'.format(nbp))
      region = ranges.make_range(chr2, start2, end2)
      self.assertEqual(utils.read_overlaps_region(read, region), expected)
      # This check ensures we get the same result calling ranges.ranges_overlap.
      self.assertEqual(
          ranges.ranges_overlap(region, utils.read_range(read)), expected)

    check_overlaps(ref1, s1, e1, ref2, s2, e2, expected)
    check_overlaps(ref2, s2, e2, ref1, s1, e1, expected)


if __name__ == '__main__':
  absltest.main()
