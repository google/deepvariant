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

"""Tests for variantcall_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import functools

from absl.testing import absltest
from absl.testing import parameterized
import mock

from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import struct_utils
from third_party.nucleus.util import variantcall_utils


class VariantcallUtilsTests(parameterized.TestCase):

  def _assert_struct_lists_equal(self, actual, expected):
    self.assertEqual(len(actual), len(expected))
    for actual_elem, expected_elem in zip(actual, expected):
      self.assertEqual(actual_elem, expected_elem)

  @parameterized.parameters(
      dict(
          field_name='GP',
          value=[.1, .2, .7],
          reader=None,
          expected=[struct_pb2.Value(number_value=v) for v in [.1, .2, .7]]),
      dict(
          field_name='AD',
          value=[23],
          reader=None,
          expected=[struct_pb2.Value(int_value=23)]),
      dict(
          field_name='FT',
          value=['PASS'],
          reader=None,
          expected=[struct_pb2.Value(string_value='PASS')]),
      dict(
          field_name='FT',
          value=['PASS'],
          reader=True,
          expected=[struct_pb2.Value(string_value='PASS')]),
  )
  def test_set_format(self, field_name, value, reader, expected):
    if reader is not None:
      reader = mock.Mock()
      reader.field_access_cache.format_field_set_fn.return_value = (
          struct_utils.set_string_field)
    call = variants_pb2.VariantCall()
    variantcall_utils.set_format(call, field_name, value, reader)
    actual = call.info[field_name].values
    self._assert_struct_lists_equal(actual, expected)

  @parameterized.parameters(
      dict(field_name='GP', reader=None, expected=[.1, .2, .7]),
      dict(field_name='AD', reader=None, expected=[55, 3]),
      dict(field_name='DP', reader=None, expected=58),
      dict(field_name='GL', reader=None, expected=[-1, -3, -5.5]),
      dict(field_name='GT', reader=None, expected=[0, 1]),
      dict(field_name='FT', reader=None, expected='LowQual'),
      dict(field_name='FT', reader=True, expected='LowQual'),
      dict(field_name='MT', reader=None, expected='0/0'),
  )
  def test_get_format(self, field_name, reader, expected):
    if reader is not None:
      reader = mock.Mock()
      reader.field_access_cache.format_field_get_fn.return_value = (
          functools.partial(
              struct_utils.get_string_field, is_single_field=True))

    call = variants_pb2.VariantCall()
    variantcall_utils.set_format(call, 'GP', [.1, .2, .7])
    variantcall_utils.set_format(call, 'AD', [55, 3])
    variantcall_utils.set_format(call, 'DP', 58)
    variantcall_utils.set_format(call, 'GL', [-1, -3, -5.5])
    variantcall_utils.set_format(call, 'GT', [0, 1])
    variantcall_utils.set_format(call, 'FT', ['LowQual'])
    variantcall_utils.set_format(call, 'MT', ['0/0'])
    actual = variantcall_utils.get_format(call, field_name, vcf_object=reader)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(
          field_name='AD',
          setter=variantcall_utils.set_ad,
          getter=variantcall_utils.get_ad,
          values=[[1, 5], [30, 29]],
      ),
      dict(
          field_name='GL',
          setter=variantcall_utils.set_gl,
          getter=variantcall_utils.get_gl,
          values=[[-1, -2, -3.3], [-0.001, -3, -10]],
      ),
      dict(
          field_name='GQ',
          setter=variantcall_utils.set_gq,
          getter=variantcall_utils.get_gq,
          values=range(10),
      ),
      dict(
          field_name='GT',
          setter=variantcall_utils.set_gt,
          getter=variantcall_utils.get_gt,
          values=[[0, 1], [1, 1], [1, 2]],
      ),
      dict(
          field_name='MIN_DP',
          setter=variantcall_utils.set_min_dp,
          getter=variantcall_utils.get_min_dp,
          values=range(10),
      ),
      dict(
          field_name='MED_DP',
          setter=variantcall_utils.set_med_dp,
          getter=variantcall_utils.get_med_dp,
          values=range(10),
      ),
      dict(
          field_name='MT',
          setter=variantcall_utils.set_mt,
          getter=variantcall_utils.get_mt,
          values=['0/0', '0/1'],
      ),
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

  def test_set_ps(self):
    vc = variants_pb2.VariantCall()
    variantcall_utils.set_ps(vc, '1-81')
    self.assertEqual(vc.info['PS'].values[0].int_value, 100_081)

  @parameterized.parameters(
      dict(genotype=[], expected=False),
      dict(genotype=[-1], expected=False),
      dict(genotype=[-1, -1], expected=False),
      dict(genotype=[-1, -1073741825], expected=False),
      dict(genotype=[-1, 0], expected=True),
      dict(genotype=[0, 0], expected=True),
      dict(genotype=[0, 1], expected=True),
      dict(genotype=[0, 1, -1], expected=True),
      dict(genotype=[-1, 0, -1073741825], expected=True),
  )
  def test_has_genotypes(self, genotype, expected):
    call = variants_pb2.VariantCall(genotype=genotype)
    actual = variantcall_utils.has_genotypes(call)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(genotype=[], expected=True),
      dict(genotype=[-1], expected=False),
      dict(genotype=[-1, -1073741825], expected=False),
      dict(genotype=[-1, 0], expected=False),
      dict(genotype=[0, 0], expected=True),
      dict(genotype=[0, 1], expected=True),
      dict(genotype=[0, 1, -1], expected=False),
      dict(genotype=[1, 0, -1073741825], expected=False),
  )
  def test_has_full_genotypes(self, genotype, expected):
    call = variants_pb2.VariantCall(genotype=genotype)
    actual = variantcall_utils.has_full_genotypes(call)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(genotype=[], expected=0),
      dict(genotype=[-1], expected=1),
      dict(genotype=[-1, -1], expected=2),
      dict(genotype=[-1, -1073741825], expected=1),
      dict(genotype=[-1, 0], expected=2),
      dict(genotype=[0, 0], expected=2),
      dict(genotype=[0, 1], expected=2),
      dict(genotype=[0, 1, -1], expected=3),
      dict(genotype=[-1, 0, -1073741825], expected=2),
  )
  def test_ploidy(self, genotype, expected):
    call = variants_pb2.VariantCall(genotype=genotype)
    actual = variantcall_utils.ploidy(call)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(genotype=[], expected=False),
      dict(genotype=[-1], expected=False),
      dict(genotype=[-1, -1], expected=False),
      dict(genotype=[-1, -1073741825], expected=False),
      dict(genotype=[-1, 0], expected=False),
      dict(genotype=[0, 0], expected=False),
      dict(genotype=[0, 1], expected=True),
      dict(genotype=[0, 1, -1], expected=True),
      dict(genotype=[-1, 0, -1073741825], expected=False),
  )
  def test_has_variation(self, genotype, expected):
    call = variants_pb2.VariantCall(genotype=genotype)
    actual = variantcall_utils.has_variation(call)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(genotype=[], expected=False),
      dict(genotype=[-1], expected=False),
      dict(genotype=[-1, -1], expected=False),
      dict(genotype=[-1, -1073741825], expected=False),
      dict(genotype=[-1, 0], expected=False),
      dict(genotype=[0, 0], expected=False),
      dict(genotype=[0, 1], expected=True),
      dict(genotype=[0, 2], expected=True),
      dict(genotype=[1, 1], expected=False),
      dict(genotype=[2, 2], expected=False),
      dict(genotype=[1, 2], expected=True),
      dict(genotype=[0, 1, -1], expected=True),
      dict(genotype=[0, 1, 2], expected=True),
      dict(genotype=[-1, 0, -1073741825], expected=False),
  )
  def test_is_heterozygous(self, genotype, expected):
    call = variants_pb2.VariantCall(genotype=genotype)
    actual = variantcall_utils.is_heterozygous(call)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(mf=[0.0, 0.0], expected='0/0'),
      dict(mf=[0.2, 0.2], expected='0/0'),
      dict(mf=[0.1, 0.9], expected='0/1'),
      dict(mf=[0.5, 0.5], expected='0/0'),
      dict(mf=[0.9, 0.1], expected='0/1'),
      dict(mf=[0.8, 0.8], expected='1/1'),
      dict(mf=[1.0, 1.0], expected='1/1'),
      dict(mf=[], expected=''),
      dict(mf=None, expected=''),
  )
  def test_determine_methylation_type(self, mf, expected):
    actual = variantcall_utils.determine_methylation_type(mf)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(
          variant_call=variants_pb2.VariantCall(
              info={
                  'MT': struct_pb2.ListValue(
                      values=[struct_pb2.Value(string_value='0/0')]
                  )
              }
          ),
          expected='0/0',
      ),
      dict(
          variant_call=variants_pb2.VariantCall(
              info={
                  'MT': struct_pb2.ListValue(
                      values=[struct_pb2.Value(string_value='0/1')]
                  )
              }
          ),
          expected='0/1',
      ),
      dict(
          variant_call=variants_pb2.VariantCall(
              info={
                  'MT': struct_pb2.ListValue(
                      values=[struct_pb2.Value(string_value='1/1')]
                  )
              }
          ),
          expected='1/1',
      ),
      dict(
          variant_call=variants_pb2.VariantCall(info={}),
          expected='',
      ),
      dict(
          variant_call=variants_pb2.VariantCall(),
          expected='',
      ),
  )
  def test_get_mt(self, variant_call, expected):
    got = variantcall_utils.get_mt(variant_call=variant_call)
    self.assertEqual(expected, got)

  @parameterized.parameters(
      dict(
          variant_call=variants_pb2.VariantCall(
              info={
                  'MF': struct_pb2.ListValue(
                      values=[
                          struct_pb2.Value(number_value=0.1),
                          struct_pb2.Value(number_value=0.2),
                      ]
                  )
              }
          ),
          expected=[0.1, 0.2],
      ),
      dict(
          variant_call=variants_pb2.VariantCall(
              info={
                  'MF': struct_pb2.ListValue(
                      values=[
                          struct_pb2.Value(number_value=0.5),
                          struct_pb2.Value(number_value=0.8),
                      ]
                  )
              }
          ),
          expected=[0.5, 0.8],
      ),
      dict(
          variant_call=variants_pb2.VariantCall(
              info={
                  'MF': struct_pb2.ListValue(
                      values=[struct_pb2.Value(number_value=1.0)]
                  )
              }
          ),
          expected=[1.0],
      ),
      dict(
          variant_call=variants_pb2.VariantCall(info={}),
          expected=[],
      ),
      dict(
          variant_call=variants_pb2.VariantCall(),
          expected=[],
      ),
  )
  def test_get_mf(self, variant_call, expected):
    got = variantcall_utils.get_mf(variant_call=variant_call)
    self.assertEqual(expected, got)


if __name__ == '__main__':
  absltest.main()
