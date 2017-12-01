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
"""Tests for cigar."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools



from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.core import cigar
from deepvariant.core.genomics import cigar_pb2

_CIGAR_TUPLES_AND_CIGAR_UNITS = [
    ((1, 'M'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.ALIGNMENT_MATCH, operation_length=1)),
    ((2, 'I'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.INSERT, operation_length=2)),
    ((3, 'D'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.DELETE, operation_length=3)),
    ((4, 'N'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.SKIP, operation_length=4)),
    ((5, 'S'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.CLIP_SOFT, operation_length=5)),
    ((6, 'H'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.CLIP_HARD, operation_length=6)),
    ((7, 'P'),
     cigar_pb2.CigarUnit(operation=cigar_pb2.CigarUnit.PAD,
                         operation_length=7)),
    ((8, '='),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.SEQUENCE_MATCH, operation_length=8)),
    ((9, 'X'),
     cigar_pb2.CigarUnit(
         operation=cigar_pb2.CigarUnit.SEQUENCE_MISMATCH, operation_length=9)),
]


def _example_cigar_string_and_units(repeat=3):
  examples = {}
  for x in itertools.product(_CIGAR_TUPLES_AND_CIGAR_UNITS, repeat=repeat):
    lengths_and_opstrs, cigar_units = zip(*x)
    cigar_str = ''.join(str(l) + opstr for l, opstr in lengths_and_opstrs)
    examples[cigar_str] = list(cigar_units)
  return examples


class CigarTests(parameterized.TestCase):

  @parameterized.parameters((cigar_units, cigar_str)
                            for cigar_str, cigar_units in
                            _example_cigar_string_and_units(3).iteritems())
  def test_format_cigar_units(self, cigar_units, expected):
    self.assertEqual(cigar.format_cigar_units(cigar_units), expected)

  @parameterized.parameters(
      ('10M', 10),
      ('10=', 10),
      ('10X', 10),
      ('10M2I3M', 13),
      ('10M2D3M', 15),
      ('10M2N3M', 15),
      ('1S10M2D3M', 15),
      ('1S10M2D3M1S', 15),
      ('1S10M2D3M1S5H', 15),
      ('8H1S10M2D3M1S5H', 15),
      ('8H1S10M2N3M1S5H', 15),
  )
  def test_alignment_length(self, cigar_str, expected):
    cigar_units = cigar.parse_cigar_string(cigar_str)
    self.assertEqual(cigar.alignment_length(cigar_units), expected)

  @parameterized.parameters(
      (length, opstr, expected)
      for (length, opstr), expected in _CIGAR_TUPLES_AND_CIGAR_UNITS)
  def test_to_cigar_unit(self, length, opstr, expected):
    # Check we can convert a tuple and list of length and opstr.
    self.assertEqual(cigar.to_cigar_unit((length, opstr)), expected)
    self.assertEqual(cigar.to_cigar_unit([length, opstr]), expected)

    # Check that we can convert a string version len+opstr.
    self.assertEqual(cigar.to_cigar_unit(str(length) + opstr), expected)

    # Check that we can pass a CigarUnit through without modification.
    self.assertEqual(cigar.to_cigar_unit(expected), expected)

    # Check that we can pass length as a string.
    self.assertEqual(cigar.to_cigar_unit((str(length), opstr)), expected)

  @parameterized.parameters(
      '-1M',
      '0M',
      '',
      'M',
      'M12',
      '4m',
      # Have to be wrapped in a list to stop parameterized from treating the
      # tuple as the positional arguments to the test function.
      [()],
      [(4)],
      [(4, 'M', 'X')],
      [(4, 'M', 10)],
      [{4, 'M'}],
      # This integer is too large to fit in an int64 cigar, make sure an
      # exception is thrown. Max int64 is 9,223,372,036,854,775,807, so we try
      # one more.
      '9223372036854775808M',
  )
  def test_to_cigar_unit_detects_bad_args(self, bad):
    with self.assertRaises(ValueError):
      cigar.to_cigar_unit(bad)

  @parameterized.parameters(
      zip(*to_convert)
      for to_convert in itertools.product(
          _CIGAR_TUPLES_AND_CIGAR_UNITS, repeat=3))
  def test_to_cigar_units(self, to_convert, expected):
    # We can convert the raw form.
    to_convert = list(to_convert)
    expected = list(expected)

    actual = cigar.to_cigar_units(to_convert)
    self.assertEqual(actual, expected)

    # We can also convert the string form by formatting actual.
    self.assertEqual(
        cigar.to_cigar_units(cigar.format_cigar_units(actual)), expected)

  @parameterized.parameters(
      (str(length) + opstr, [expected])
      for (length, opstr), expected in _CIGAR_TUPLES_AND_CIGAR_UNITS)
  def test_parse_cigar_string_single(self, cigar_str, expected):
    self.assertEqual(cigar.parse_cigar_string(cigar_str), expected)

  @parameterized.parameters(_example_cigar_string_and_units(3).iteritems())
  def test_parse_cigar_string_three_pieces(self, cigar_str, expected):
    self.assertEqual(cigar.parse_cigar_string(cigar_str), expected)

  @parameterized.parameters(
      '',
      '12',
      '12m',
      '12?',
      'M12',
      '12M1',
      '12MI',
      '12M-1I',
      '12.0M',
  )
  def test_parse_cigar_string_detects_bad_inputs(self, bad_cigar_str):
    with self.assertRaises(ValueError):
      cigar.parse_cigar_string(bad_cigar_str)


if __name__ == '__main__':
  absltest.main()
