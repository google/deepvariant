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
"""Tests for third_party.nucleus.io.vcf."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest

from third_party.nucleus.io import genomics_reader
from tensorflow.python.lib.io import python_io


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
    values = list(dreader)
    self.assertEqual(range(10), values)

  def testTwoIteratorsAtTheSameTime(self):
    dreader = DummyReader('100')
    iter2 = iter(dreader)
    for i in range(100):
      self.assertEqual(i, dreader.next())
      self.assertEqual(i, iter2.next())


class DummyProto(object):
  """A pretend protocol buffer class that only provides FromString."""

  def FromString(self, buf):
    return buf


class TFRecordReaderTests(absltest.TestCase):
  """Tests for TFRecordReader."""

  def setUp(self):
    python_io.tf_record_iterator = lambda x, y: x.split(',')

  def testMock(self):
    reader = genomics_reader.TFRecordReader('a,b,c,d,e', DummyProto())
    self.assertEqual(['a','b','c','d','e'], list(reader))

  def testTwoIteratorsAtTheSameTime(self):
    dreader = genomics_reader.TFRecordReader('0,1,2,3,4,5', DummyProto())
    iter2 = iter(dreader)
    for i in range(6):
      self.assertEqual(str(i), dreader.next())
      self.assertEqual(str(i), iter2.next())


if __name__ == '__main__':
  absltest.main()
