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

"""Tests for third_party.nucleus.util.sequence_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest

from absl.testing import parameterized

from third_party.nucleus.util import sequence_utils


class SequenceUtilsTests(parameterized.TestCase):

  @parameterized.parameters(
      dict(seq='', expected=''),
      dict(seq='A', expected='T'),
      dict(seq='T', expected='A'),
      dict(seq='C', expected='G'),
      dict(seq='G', expected='C'),
      dict(seq='GGGCAGATT', expected='AATCTGCCC'),
      dict(
          seq='GGGCAGANN',
          expected='NNTCTGCCC',
          complement_dict=sequence_utils.DNA_COMPLEMENT_UPPER),
      dict(
          seq='accgt',
          expected='acggt',
          complement_dict=sequence_utils.DNA_COMPLEMENT),
      dict(
          seq='ATCGRYSWKMBVDHN',
          expected='NDHBVKMWSRYCGAT',
          complement_dict=sequence_utils.IUPAC_DNA_COMPLEMENT_UPPER),
      dict(
          seq='ATCGRYSWKMBVDHNatcgryswkmbvdhn',
          expected='ndhbvkmwsrycgatNDHBVKMWSRYCGAT',
          complement_dict=sequence_utils.IUPAC_DNA_COMPLEMENT),
  )
  def test_reverse_complement(self, seq, expected, complement_dict=None):
    """Tests canonical DNA sequences are reverse complemented correctly."""
    self.assertEqual(
        sequence_utils.reverse_complement(seq, complement_dict), expected)

  @parameterized.parameters(
      dict(seq='GGGCAGANN'),
      dict(seq='accgt'),
      dict(
          seq='ATCGRYSWKMBVDHNatcgryswkmbvdhn',
          complement_dict=sequence_utils.IUPAC_DNA_COMPLEMENT_UPPER),
      dict(seq='X', complement_dict=sequence_utils.IUPAC_DNA_COMPLEMENT),
  )
  def test_bad_reverse_complement(self, seq, complement_dict=None):
    """Tests error is raised when complement_dict does not cover given seq."""
    with self.assertRaisesRegexp(sequence_utils.Error, 'Unknown base in'):
      sequence_utils.reverse_complement(seq, complement_dict)

  @parameterized.parameters(
      dict(
          bases_set=sequence_utils.STRICT_DNA_BASES_UPPER,
          complement_dict=sequence_utils.STRICT_DNA_COMPLEMENT_UPPER),
      dict(
          bases_set=sequence_utils.STRICT_DNA_BASES,
          complement_dict=sequence_utils.STRICT_DNA_COMPLEMENT),
      dict(
          bases_set=sequence_utils.DNA_BASES_UPPER,
          complement_dict=sequence_utils.DNA_COMPLEMENT_UPPER),
      dict(
          bases_set=sequence_utils.DNA_BASES,
          complement_dict=sequence_utils.DNA_COMPLEMENT),
  )
  def test_base_set_definitions(self, bases_set, complement_dict):
    """Tests that base set and complement dict definitions are consistent."""
    self.assertEqual(bases_set, frozenset(complement_dict.keys()))


if __name__ == '__main__':
  absltest.main()
