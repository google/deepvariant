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
"""Tests for gff_reader CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import clif_postproc
from third_party.nucleus.io.python import gff_reader
from third_party.nucleus.protos import gff_pb2
from third_party.nucleus.testing import test_utils


class GffReaderTest(parameterized.TestCase):

  def setUp(self):
    self.options = gff_pb2.GffReaderOptions()
    self.first = gff_pb2.GffRecord()
    self.first.range.reference_name = 'ctg123'
    self.first.range.start = 999
    self.first.range.end = 9000
    self.first.source = 'GenBank'
    self.first.type = 'gene'
    self.first.score = 2.5
    self.first.strand = gff_pb2.GffRecord.FORWARD_STRAND
    self.first.phase = 0
    self.first.attributes['ID'] = 'gene00001'
    self.first.attributes['Name'] = 'EDEN'

    self.second = gff_pb2.GffRecord()
    self.second.range.reference_name = 'ctg123'
    self.second.range.start = 999
    self.second.range.end = 1012
    self.second.phase = -1
    self.second.score = -float('inf')

  @parameterized.parameters('test_features.gff', 'test_features.gff.gz')
  def test_gff_iterate(self, test_features_gff_filename):
    file_path = test_utils.genomics_core_testdata(test_features_gff_filename)
    with gff_reader.GffReader.from_file(file_path, self.options) as reader:
      iterable = reader.iterate()
      self.assertIsInstance(iterable, clif_postproc.WrappedCppIterable)
      actual = list(iterable)
      self.assertLen(actual, 2)
      self.assertEqual(actual[0], self.first)
      self.assertEqual(actual[1], self.second)

  def test_from_file_raises_with_missing_gff(self):
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(ValueError, 'Could not open missing.gff'):
      gff_reader.GffReader.from_file('missing.gff', self.options)

  def test_ops_on_closed_reader_raise(self):
    file_path = test_utils.genomics_core_testdata('test_features.gff')
    reader = gff_reader.GffReader.from_file(file_path, self.options)
    with reader:
      pass
    # At this point the reader is closed.
    with self.assertRaisesRegexp(ValueError, 'Cannot Iterate a closed'):
      reader.iterate()


if __name__ == '__main__':
  absltest.main()
