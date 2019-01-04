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

"""Tests for third_party.nucleus.io.tfrecord."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import types

from absl.testing import absltest
from absl.testing import parameterized
import mock

from third_party.nucleus.io import tfrecord

from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.testing import test_utils
from tensorflow.python.lib.io import python_io


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

  @parameterized.parameters(
      ('foo.tfrecord', ''),
      ('foo.tfrecord.gz', 'GZIP'),
      (['foo.tfrecord', 'bar.tfrecord'], ''),
      (['foo.tfrecord.gz', 'bar.tfrecord.gz'], 'GZIP'),
  )
  def test_make_tfrecord_options(self, filenames, expected_compression_type):
    compression_type = python_io.TFRecordOptions.get_compression_type_string(
        tfrecord.make_tfrecord_options(filenames))
    self.assertEqual(compression_type, expected_compression_type)

  @parameterized.parameters(
      (['foo.tfrecord', 'bar.tfrecord.gz'],),
      (['foo.tfrecord', 'bar.tfrecord', 'baz.tfrecord.gz'],),
  )
  def test_make_tfrecord_options_with_bad_inputs(self, filenames):
    with self.assertRaisesRegexp(
        ValueError,
        'Incorrect value: {}. Filenames need to be all of the same type: '
        'either all with .gz or all without .gz'.format(','.join(filenames))):
      tfrecord.make_tfrecord_options(filenames)

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


class RawProtoWriterAdaptorTests(parameterized.TestCase):

  def setUp(self):
    self.proto1 = reference_pb2.ContigInfo(
        name='p1', n_bases=10, pos_in_fasta=0)
    self.proto2 = reference_pb2.ContigInfo(
        name='p2', n_bases=20, pos_in_fasta=1)
    self.protos = [self.proto1, self.proto2]

  @parameterized.parameters(
      dict(take_ownership=True),
      dict(take_ownership=False),
  )
  def test_adaptor_with_ownership(self, take_ownership):
    mock_writer = mock.MagicMock()
    adaptor = tfrecord.RawProtoWriterAdaptor(
        mock_writer, take_ownership=take_ownership)

    # Write out protos to our adaptor.
    with adaptor as enter_return_value:
      # Make sure that __enter__ returns the adaptor itself.
      self.assertIs(adaptor, enter_return_value)
      adaptor.write(self.proto1)
      adaptor.write(self.proto2)

    if take_ownership:
      # If we took ownership, mock_writer __enter__ and __exit__ should have
      # been called.
      mock_writer.__enter__.assert_called_once_with()
      test_utils.assert_called_once_workaround(mock_writer.__exit__)
    else:
      # If not, they shouldn't have been called.
      test_utils.assert_not_called_workaround(mock_writer.__enter__)
      test_utils.assert_not_called_workaround(mock_writer.__exit__)

    self.assertEqual(mock_writer.write.call_args_list,
                     [mock.call(r.SerializeToString()) for r in self.protos])


if __name__ == '__main__':
  absltest.main()
