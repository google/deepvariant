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
"""Tests for third_party.nucleus.io.vcf."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
import six

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.protos import gff_pb2
from third_party.nucleus.testing import test_utils


class DummyReader(genomics_reader.GenomicsReader):
  """A GenomicsReader that produces consecutive integers."""

  def __init__(self, input_path):
    self.limit = int(input_path)
    super(DummyReader, self).__init__()

  def iterate(self):
    for i in range(self.limit):
      yield i

  def query(self, region):
    raise NotImplementedError('Can not query DummyReader')

  def __exit__(self, exit_type, exit_value, exit_traceback):
    pass


class GenomicsReaderTests(absltest.TestCase):
  """Tests for GenomicsReader."""

  def testIteration(self):
    dreader = DummyReader('10')
    self.assertEqual(list(range(10)), list(dreader))

  def testTwoIteratorsAtTheSameTime(self):
    dreader = DummyReader('100')
    iter2 = iter(dreader)
    for i in range(100):
      self.assertEqual(i, six.next(dreader))
      self.assertEqual(i, six.next(iter2))


class TFRecordReaderTests(absltest.TestCase):
  """Tests for TFRecordReader."""

  def testUncompressed(self):
    reader = genomics_reader.TFRecordReader(
        test_utils.genomics_core_testdata('test_features.gff.tfrecord'),
        gff_pb2.GffRecord(),
    )
    records = list(reader.iterate())
    self.assertEqual('GenBank', records[0].source)
    self.assertEqual('ctg123', records[1].range.reference_name)
    self.assertNotEqual(reader.c_reader, 0)

  def testUncompressedExplicit(self):
    reader = genomics_reader.TFRecordReader(
        test_utils.genomics_core_testdata('test_features.gff.tfrecord'),
        gff_pb2.GffRecord(),
        compression_type=''
    )
    records = list(reader.iterate())
    self.assertEqual('GenBank', records[0].source)
    self.assertEqual('ctg123', records[1].range.reference_name)

  def testCompressed(self):
    reader = genomics_reader.TFRecordReader(
        test_utils.genomics_core_testdata('test_features.gff.tfrecord.gz'),
        gff_pb2.GffRecord(),
    )
    records = list(reader.iterate())
    self.assertEqual('GenBank', records[0].source)
    self.assertEqual('ctg123', records[1].range.reference_name)

  def testCompressedExplicit(self):
    reader = genomics_reader.TFRecordReader(
        test_utils.genomics_core_testdata('test_features.gff.tfrecord.gz'),
        gff_pb2.GffRecord(),
        compression_type='GZIP'
    )
    records = list(reader.iterate())
    self.assertEqual('GenBank', records[0].source)
    self.assertEqual('ctg123', records[1].range.reference_name)


if __name__ == '__main__':
  absltest.main()
