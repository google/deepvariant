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

"""Tests for third_party.nucleus.io.gfile."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest

from third_party.nucleus.io import gfile
from third_party.nucleus.testing import test_utils


class GfileTest(absltest.TestCase):

  def test_exists(self):
    self.assertTrue(gfile.Exists(
        test_utils.genomics_core_testdata('test_regions.bedgraph')))
    self.assertFalse(gfile.Exists(
        test_utils.genomics_core_testdata('does_not_exist')))

  def test_glob(self):
    g1 = gfile.Glob(test_utils.genomics_core_testdata('test*'))
    self.assertGreater(len(g1), 1)
    self.assertIn(
        test_utils.genomics_core_testdata('test.bam'), g1)
    g2 = gfile.Glob(test_utils.genomics_core_testdata('does_not_exist*'))
    self.assertEqual([], g2)

  def test_reading(self):
    with gfile.Open(test_utils.genomics_core_testdata('headerless.sam')) as f:
      for line in f:
        self.assertTrue(line.startswith('SRR3656745'))

  def test_writing(self):
    path = test_utils.test_tmpfile('test_gfile')
    with gfile.Open(path, 'w') as f:
      f.write('test\n')
      f.write('end\n')

    with gfile.Open(path, 'r') as f2:
      lines = f2.readlines()

    self.assertEqual(['test\n', 'end\n'], lines)


if __name__ == '__main__':
  absltest.main()
