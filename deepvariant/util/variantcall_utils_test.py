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
"""Tests for variantcall_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import functools

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.util.genomics import struct_pb2
from deepvariant.util.genomics import variants_pb2
from deepvariant.util import struct_utils
from deepvariant.util import variantcall_utils


class DummyVcfHeader(object):
  """Dummy VcfHeader class that implements format_field_set_fn."""

  def format_field_set_fn(self, field_name):
    del field_name  # Unused.
    return struct_utils.set_string_field

  def format_field_get_fn(self, field_name):
    return functools.partial(
        struct_utils.get_string_field, is_single_field=True)


class VariantcallUtilsTests(parameterized.TestCase):

  def _assert_struct_lists_equal(self, actual, expected):
    self.assertEqual(len(actual), len(expected))
    for actual_elem, expected_elem in zip(actual, expected):
      self.assertEqual(actual_elem, expected_elem)

  @parameterized.parameters(
      dict(
          field_name='GP',
          value=[.1, .2, .7],
          header=None,
          expected=[struct_pb2.Value(number_value=v) for v in [.1, .2, .7]]),
      dict(
          field_name='AD',
          value=[23],
          header=None,
          expected=[struct_pb2.Value(int_value=23)]),
      dict(
          field_name='FT',
          value=['PASS'],
          header=None,
          expected=[struct_pb2.Value(string_value='PASS')]),
      dict(
          field_name='FT',
          value=['PASS'],
          header=DummyVcfHeader(),
          expected=[struct_pb2.Value(string_value='PASS')]),
  )
  def test_set_format(self, field_name, value, header, expected):
    call = variants_pb2.VariantCall()
    variantcall_utils.set_format(call, field_name, value, header)
    actual = call.info[field_name].values
    self._assert_struct_lists_equal(actual, expected)

  @parameterized.parameters(
      dict(field_name='GP', header=None, expected=[.1, .2, .7]),
      dict(field_name='AD', header=None, expected=[55, 3]),
      dict(field_name='DP', header=None, expected=58),
      dict(field_name='GL', header=None, expected=[-1, -3, -5.5]),
      dict(field_name='GT', header=None, expected=[0, 1]),
      dict(field_name='FT', header=None, expected='LowQual'),
      dict(field_name='FT', header=DummyVcfHeader(), expected='LowQual'),
  )
  def test_get_format(self, field_name, header, expected):
    call = variants_pb2.VariantCall()
    variantcall_utils.set_format(call, 'GP', [.1, .2, .7])
    variantcall_utils.set_format(call, 'AD', [55, 3])
    variantcall_utils.set_format(call, 'DP', 58)
    variantcall_utils.set_format(call, 'GL', [-1, -3, -5.5])
    variantcall_utils.set_format(call, 'GT', [0, 1])
    variantcall_utils.set_format(call, 'FT', ['LowQual'])
    actual = variantcall_utils.get_format(call, field_name, vcf_header=header)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(
          field_name='AD',
          setter=variantcall_utils.set_ad,
          getter=variantcall_utils.get_ad,
          values=[[1, 5], [30, 29]]),
      dict(
          field_name='GL',
          setter=variantcall_utils.set_gl,
          getter=variantcall_utils.get_gl,
          values=[[-1, -2, -3.3], [-0.001, -3, -10]]),
      dict(
          field_name='GQ',
          setter=variantcall_utils.set_gq,
          getter=variantcall_utils.get_gq,
          values=range(10)),
      dict(
          field_name='GT',
          setter=variantcall_utils.set_gt,
          getter=variantcall_utils.get_gt,
          values=[[0, 1], [1, 1], [1, 2]]),
      dict(
          field_name='MIN_DP',
          setter=variantcall_utils.set_min_dp,
          getter=variantcall_utils.get_min_dp,
          values=range(10)),
  )
  def test_variantcall_format_roundtrip(self, field_name, setter, getter,
                                        values):
    vc = variants_pb2.VariantCall()
    self.assertNotIn(field_name, vc.info)
    for value in values:
      setter(vc, value)
      if field_name not in ['GT', 'GL']:
        self.assertIn(field_name, vc.info)
      actual = getter(vc)
      self.assertEqual(actual, value)


if __name__ == '__main__':
  absltest.main()
