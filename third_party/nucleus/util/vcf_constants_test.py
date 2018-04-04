# Copyright 2018 Google Inc.
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

"""Tests for third_party.nucleus.util.vcf_constants."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from absl.testing import absltest
from absl.testing import parameterized
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import struct_utils
from third_party.nucleus.util import vcf_constants


class VcfConstantsTest(parameterized.TestCase):

  def test_unique_reserved_filter(self):
    num_reserved_filter = len(vcf_constants.RESERVED_FILTER_FIELDS)
    unique_filt_ids = {filt.id for filt in vcf_constants.RESERVED_FILTER_FIELDS}
    self.assertLen(unique_filt_ids, num_reserved_filter)

  def test_unique_reserved_info(self):
    num_reserved_info = len(vcf_constants.RESERVED_INFO_FIELDS)
    unique_info_ids = {info.id for info in vcf_constants.RESERVED_INFO_FIELDS}
    self.assertLen(unique_info_ids, num_reserved_info)

  def test_unique_reserved_format(self):
    num_reserved_format = len(vcf_constants.RESERVED_FORMAT_FIELDS)
    unique_format_ids = {f.id for f in vcf_constants.RESERVED_FORMAT_FIELDS}
    self.assertLen(unique_format_ids, num_reserved_format)

  def test_get_reserved_filter(self):
    filt = vcf_constants.reserved_filter_field('PASS')
    self.assertIsInstance(filt, variants_pb2.VcfFilterInfo)
    self.assertEqual(filt.id, 'PASS')
    self.assertEqual(filt.description, 'All filters passed')

  @parameterized.parameters(
      'RefCall',
      'LowQual',
      'AD',
      'DP',
      'GT',
      'GQ',
  )
  def test_invalid_get_reserved_filter(self, field_id):
    with self.assertRaisesRegexp(ValueError, 'No reserved field with id'):
      vcf_constants.reserved_filter_field(field_id)

  @parameterized.parameters(
      'AA',
      'AC',
      'AD',
      'ADF',
      'END',
      'H2',
  )
  def test_get_reserved_info(self, field_id):
    info = vcf_constants.reserved_info_field(field_id)
    self.assertIsInstance(info, variants_pb2.VcfInfo)
    self.assertEqual(info.id, field_id)

  @parameterized.parameters(
      'PASS',
      'GT',
      'GQ',
      'GL',
      'FT',
  )
  def test_invalid_get_reserved_info(self, field_id):
    with self.assertRaisesRegexp(ValueError, 'No reserved field with id'):
      vcf_constants.reserved_info_field(field_id)

  @parameterized.parameters(
      'AD',
      'ADF',
      'DP',
      'GT',
      'GQ',
      'GL',
      'FT',
      'PL',
  )
  def test_get_reserved_format(self, field_id):
    fmt = vcf_constants.reserved_format_field(field_id)
    self.assertIsInstance(fmt, variants_pb2.VcfFormatInfo)
    self.assertEqual(fmt.id, field_id)

  @parameterized.parameters(
      'PASS',
      'AN',
      '1000G',
      'END',
      'H2',
  )
  def test_invalid_get_reserved_format(self, field_id):
    with self.assertRaisesRegexp(ValueError, 'No reserved field with id'):
      vcf_constants.reserved_format_field(field_id)

  @parameterized.parameters(
      dict(
          value_type=vcf_constants.CHARACTER_TYPE,
          values=['a'],
          number='1',
          expected='a'),
      dict(
          value_type=vcf_constants.CHARACTER_TYPE,
          values=['b'],
          number='.',
          expected=['b']),
      dict(
          value_type=vcf_constants.CHARACTER_TYPE,
          values=['c', 'd'],
          number='R',
          expected=['c', 'd']),
      dict(
          value_type=vcf_constants.FLAG_TYPE,
          values=[True],
          number='0',
          expected=True),
      dict(
          value_type=vcf_constants.FLOAT_TYPE,
          values=[2.5],
          number='1',
          expected=2.5),
      dict(
          value_type=vcf_constants.FLOAT_TYPE,
          values=[2.5],
          number='.',
          expected=[2.5]),
      dict(
          value_type=vcf_constants.FLOAT_TYPE,
          values=[2.5, 3.5],
          number='A',
          expected=[2.5, 3.5]),
      dict(
          value_type=vcf_constants.INTEGER_TYPE,
          values=[2, 3, 4],
          number='G',
          expected=[2, 3, 4]),
      dict(
          value_type=vcf_constants.STRING_TYPE,
          values=['a', 'bc'],
          number='.',
          expected=['a', 'bc']),
  )
  def test_create_get_fn(self, value_type, values, number, expected):
    info = variants_pb2.Variant().info
    set_fn = vcf_constants.SET_FN_LOOKUP[value_type]
    set_fn(info, 'field', values)
    get_fn = vcf_constants.create_get_fn(value_type, number)
    actual = get_fn(info, 'field')
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(field='CIGAR', expected=struct_utils.set_string_field),
      dict(field='DP', expected=struct_utils.set_int_field),
      dict(field='MQ', expected=struct_utils.set_number_field),
      dict(field='SOMATIC', expected=struct_utils.set_bool_field),
  )
  def test_reserved_info_field_set_fn(self, field, expected):
    actual = vcf_constants.reserved_info_field_set_fn(field)
    self.assertIs(actual, expected)

  @parameterized.parameters(
      dict(field='INVALID'),
      dict(field='EC'),
      dict(field='HQ'),
  )
  def test_invalid_reserved_info_field_set_fn(self, field):
    with self.assertRaisesRegexp(ValueError, 'Unknown reserved INFO field:'):
      vcf_constants.reserved_info_field_set_fn(field)

  def test_reserved_info_field_get_fn(self):
    info = variants_pb2.Variant().info
    values = ['C']
    struct_utils.set_string_field(info, 'AA', values)
    get_fn = vcf_constants.reserved_info_field_get_fn('AA')
    actual = get_fn(info, 'AA')
    self.assertEqual(actual, values[0])

  @parameterized.parameters(
      dict(field='INVALID'),
      dict(field='EC'),
      dict(field='HQ'),
  )
  def test_invalid_reserved_info_field_get_fn(self, field):
    with self.assertRaisesRegexp(ValueError,
                                 'Unknown reserved INFO field to get:'):
      vcf_constants.reserved_info_field_get_fn(field)

  @parameterized.parameters(
      dict(field='AD', expected=struct_utils.set_int_field),
      dict(field='GL', expected=struct_utils.set_number_field),
      dict(field='FT', expected=struct_utils.set_string_field),
  )
  def test_reserved_format_field_set_fn(self, field, expected):
    actual = vcf_constants.reserved_format_field_set_fn(field)
    self.assertIs(actual, expected)

  @parameterized.parameters(
      dict(field='INVALID'),
      dict(field='CIGAR'),
      dict(field='H2'),
  )
  def test_invalid_reserved_format_field_set_fn(self, field):
    with self.assertRaisesRegexp(ValueError, 'Unknown reserved FORMAT field:'):
      vcf_constants.reserved_format_field_set_fn(field)

  def test_reserved_format_field_get_fn(self):
    info = variants_pb2.VariantCall().info
    expected = [0.2, 0.5, 0.3]
    struct_utils.set_number_field(info, 'GP', expected[:])
    get_fn = vcf_constants.reserved_format_field_get_fn('GP')
    actual = get_fn(info, 'GP')
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(field='INVALID'),
      dict(field='CIGAR'),
      dict(field='H2'),
  )
  def test_invalid_reserved_format_field_get_fn(self, field):
    with self.assertRaisesRegexp(ValueError,
                                 'Unknown reserved FORMAT field to get:'):
      vcf_constants.reserved_format_field_get_fn(field)


if __name__ == '__main__':
  absltest.main()
