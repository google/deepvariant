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

"""Tests for third_party.nucleus.io.tfrecord."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import types

from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import tfrecord

from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.testing import test_utils


class IOTest(parameterized.TestCase):

  def write_test_protos(self, filename):
    protos = [reference_pb2.ContigInfo(name=str(i)) for i in range(10)]
    path = test_utils.test_tmpfile(filename)
    tfrecord.write_tfrecords(protos, path)
    return protos, path

  @parameterized.parameters('foo.tfrecord', 'foo@2.tfrecord', 'foo@3.tfrecord')
  def test_read_write_tfrecords(self, filename):
    protos, path = self.write_test_protos(filename)

    # Create our generator of records from read_tfrecords.
    reader = tfrecord.read_tfrecords(path, reference_pb2.ContigInfo)

    # Make sure it's actually a generator.
    self.assertEqual(type(reader), types.GeneratorType)

    # Check the round-trip contents.
    if '@' in filename:
      # Sharded outputs are striped across shards, so order isn't preserved.
      self.assertCountEqual(protos, reader)
    else:
      self.assertEqual(protos, list(reader))

  @parameterized.parameters((filename, max_records)
                            for max_records in [None, 0, 1, 3, 100]
                            for filename in ['foo.tfrecord', 'foo@2.tfrecord'])
  def test_read_tfrecords_max_records(self, filename, max_records):
    protos, path = self.write_test_protos(filename)

    # Create our generator of records from read_tfrecords.
    if max_records is None:
      expected_n = len(protos)
    else:
      expected_n = min(max_records, len(protos))
    actual = tfrecord.read_tfrecords(
        path, reference_pb2.ContigInfo, max_records=max_records)
    self.assertLen(list(actual), expected_n)

  @parameterized.parameters('foo.tfrecord', 'foo@2.tfrecord', 'foo@3.tfrecord')
  def test_shard_sorted_tfrecords(self, filename):
    protos, path = self.write_test_protos(filename)

    # Create our generator of records.
    key = lambda x: int(x.name)
    reader = tfrecord.read_shard_sorted_tfrecords(
        path, key=key, proto=reference_pb2.ContigInfo)

    # Make sure it's actually a generator.
    self.assertEqual(type(reader), types.GeneratorType)

    # Check the round-trip contents.
    contents = list(reader)
    self.assertEqual(protos, contents)
    self.assertEqual(contents, sorted(contents, key=key))

  @parameterized.parameters((filename, max_records)
                            for max_records in [None, 0, 1, 3, 100]
                            for filename in ['foo.tfrecord', 'foo@2.tfrecord'])
  def test_shard_sorted_tfrecords_max_records(self, filename, max_records):
    protos, path = self.write_test_protos(filename)

    if max_records is None:
      expected_n = len(protos)
    else:
      expected_n = min(max_records, len(protos))
    # Create our generator of records from read_tfrecords.
    actual = tfrecord.read_shard_sorted_tfrecords(
        path,
        key=lambda x: int(x.name),
        proto=reference_pb2.ContigInfo,
        max_records=max_records)
    self.assertLen(list(actual), expected_n)


if __name__ == '__main__':
  absltest.main()
