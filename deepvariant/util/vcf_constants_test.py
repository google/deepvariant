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
"""Tests for deepvariant.util.vcf_constants."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from absl.testing import absltest
from absl.testing import parameterized
from deepvariant.util import struct_utils
from deepvariant.util import vcf_constants


class VcfConstantsTest(parameterized.TestCase):

  def test_unique_reserved_info(self):
    num_reserved_info = len(vcf_constants.RESERVED_INFO_FIELDS)
    unique_info_ids = {info.id for info in vcf_constants.RESERVED_INFO_FIELDS}
    self.assertLen(unique_info_ids, num_reserved_info)

  def test_unique_reserved_format(self):
    num_reserved_format = len(vcf_constants.RESERVED_FORMAT_FIELDS)
    unique_format_ids = {f.id for f in vcf_constants.RESERVED_FORMAT_FIELDS}
    self.assertLen(unique_format_ids, num_reserved_format)

  @parameterized.parameters(
      dict(field='CIGAR', expected=struct_utils.set_string_field),
      dict(field='DP', expected=struct_utils.set_int_field),
      dict(field='MQ', expected=struct_utils.set_number_field),
      dict(field='SOMATIC', expected=struct_utils.set_bool_field),
  )
  def test_get_reserved_info_field(self, field, expected):
    actual = vcf_constants.reserved_info_field_set_fn(field)
    self.assertIs(actual, expected)

  @parameterized.parameters(
      dict(field='INVALID'),
      dict(field='EC'),
      dict(field='HQ'),
  )
  def test_invalid_reserved_info_field(self, field):
    with self.assertRaisesRegexp(ValueError, 'Unknown reserved INFO field:'):
      vcf_constants.reserved_info_field_set_fn(field)

  @parameterized.parameters(
      dict(field='AD', expected=struct_utils.set_int_field),
      dict(field='GL', expected=struct_utils.set_number_field),
      dict(field='FT', expected=struct_utils.set_string_field),
  )
  def test_get_reserved_format_field(self, field, expected):
    actual = vcf_constants.reserved_format_field_set_fn(field)
    self.assertIs(actual, expected)

  @parameterized.parameters(
      dict(field='INVALID'),
      dict(field='CIGAR'),
      dict(field='H2'),
  )
  def test_invalid_reserved_format_field(self, field):
    with self.assertRaisesRegexp(ValueError, 'Unknown reserved FORMAT field:'):
      vcf_constants.reserved_format_field_set_fn(field)


if __name__ == '__main__':
  absltest.main()
