# Copyright 2019 Google LLC.
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
"""Tests for third_party.nucleus.io.tabix."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import os
import shutil

from absl.testing import absltest

from tensorflow.python.platform import gfile
from third_party.nucleus.io import tabix
from third_party.nucleus.io import vcf
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges


class TabixTest(absltest.TestCase):
  """Test the functionality of tabix.build_index."""

  def setUp(self):
    super(TabixTest, self).setUp()
    self.input_file = test_utils.genomics_core_testdata('test_samples.vcf.gz')
    self.output_file = test_utils.test_tmpfile('test_samples.vcf.gz')
    shutil.copyfile(self.input_file, self.output_file)
    self.tbx_index_file = self.output_file + '.tbi'
    self.csi_index_file = self.output_file + '.csi'

  def tearDown(self):
    super(TabixTest, self).tearDown()
    os.remove(self.output_file)
    try:
      os.remove(self.tbx_index_file)
    except OSError:
      pass
    try:
      os.remove(self.csi_index_file)
    except OSError:
      pass

  def test_build_index_tbx(self):
    self.assertFalse(gfile.Exists(self.tbx_index_file))
    tabix.build_index(self.output_file)
    self.assertTrue(gfile.Exists(self.tbx_index_file))

  def test_build_index_csi(self):
    min_shift = 14
    self.assertFalse(gfile.Exists(self.csi_index_file))
    tabix.build_csi_index(self.output_file, min_shift)
    self.assertTrue(gfile.Exists(self.csi_index_file))

  def test_vcf_query_tbx(self):
    tabix.build_index(self.output_file)
    self.input_reader = vcf.VcfReader(self.input_file)
    self.output_reader = vcf.VcfReader(self.output_file)

    range1 = ranges.parse_literal('chr3:100,000-500,000')
    self.assertEqual(
        list(self.input_reader.query(range1)),
        list(self.output_reader.query(range1)))

  def test_vcf_query_csi(self):
    min_shift = 14
    tabix.build_csi_index(self.output_file, min_shift)
    self.input_reader = vcf.VcfReader(self.input_file)
    self.output_reader = vcf.VcfReader(self.output_file)

    range1 = ranges.parse_literal('chr3:100,000-500,000')
    self.assertEqual(
        list(self.input_reader.query(range1)),
        list(self.output_reader.query(range1)))

if __name__ == '__main__':
  absltest.main()
