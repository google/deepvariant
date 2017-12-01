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
"""Tests for deepvariant .core.io_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import Queue
import types


from absl.testing import absltest
from absl.testing import parameterized
import mock

import tensorflow as tf
from deepvariant.core import io_utils as io
from deepvariant.core import test_utils
from deepvariant.core.protos import core_pb2


class IOTest(parameterized.TestCase):

  @parameterized.parameters(
      # Unsharded outputs pass through as expected.
      (0, ['foo.txt'], [None, 'foo.txt']),
      (0, ['foo.txt', 'bar.txt'], [None, 'foo.txt', 'bar.txt']),
      (0, ['bar.txt', 'foo.txt'], [None, 'bar.txt', 'foo.txt']),
      # It's ok to have False values for other bindings.
      (0, ['foo.txt', None], [None, 'foo.txt', None]),
      (0, ['foo.txt', ''], [None, 'foo.txt', '']),
      (0, ['foo@10.txt', None], [10, 'foo-00000-of-00010.txt', None]),
      (0, ['foo@10.txt', ''], [10, 'foo-00000-of-00010.txt', '']),
      # Simple check that master behaves as expected.
      (0, ['foo@10.txt', None], [10, 'foo-00000-of-00010.txt', None]),
      (0, ['foo@10', None], [10, 'foo-00000-of-00010', None]),
      (1, ['foo@10', None], [10, 'foo-00001-of-00010', None]),
      (9, ['foo@10', None], [10, 'foo-00009-of-00010', None]),
      # Make sure we handle sharding of multiple outputs.
      (0, ['foo@10', 'bar@10', 'baz@10'],
       [10, 'foo-00000-of-00010', 'bar-00000-of-00010', 'baz-00000-of-00010']),
      (9, ['foo@10', 'bar@10', 'baz@10'],
       [10, 'foo-00009-of-00010', 'bar-00009-of-00010', 'baz-00009-of-00010']),
  )
  def test_resolve_filespecs(self, task_id, outputs, expected):
    self.assertEqual(io.resolve_filespecs(task_id, *outputs), expected)

  @parameterized.parameters(
      # shard >= num_shards.
      (10, ['foo@10']),
      # shard > 0 but master isn't sharded.
      (1, ['foo']),
      # Inconsistent sharding.
      (0, ['foo@10', 'bad@11']),
      # master isn't sharded but bad is.
      (0, ['foo', 'bad@11']),
  )
  def test_resolve_filespecs_raises_with_bad_inputs(self, task_id, outputs):
    with self.assertRaises(ValueError):
      io.resolve_filespecs(task_id, *outputs)

  @parameterized.parameters(
      # Unsharded files work.
      ('foo.txt', ['foo.txt']),
      ('foo-00000-of-00010.txt', ['foo-00000-of-00010.txt']),
      # Sharded file patterns work.
      ('foo@3.txt', [
          'foo-00000-of-00003.txt', 'foo-00001-of-00003.txt',
          'foo-00002-of-00003.txt'
      ]),
      ('foo@3',
       ['foo-00000-of-00003', 'foo-00001-of-00003', 'foo-00002-of-00003']),
  )
  def test_maybe_generate_sharded_filenames(self, filespec, expected):
    self.assertEqual(io.maybe_generate_sharded_filenames(filespec), expected)

  def write_test_protos(self, filename):
    protos = [core_pb2.ContigInfo(name=str(i)) for i in range(10)]
    path = test_utils.test_tmpfile(filename)
    io.write_tfrecords(protos, path)
    return protos, path

  @parameterized.parameters('foo.tfrecord', 'foo@2.tfrecord', 'foo@3.tfrecord')
  def test_read_write_tfrecords(self, filename):
    protos, path = self.write_test_protos(filename)

    # Create our generator of records from read_tfrecords.
    reader = io.read_tfrecords(path, core_pb2.ContigInfo)

    # Make sure it's actually a generator.
    self.assertEqual(type(reader), types.GeneratorType)

    # Check the round-trip contents.
    if '@' in filename:
      # Sharded outputs are striped across shards, so order isn't perserved.
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
    compression_type = tf.python_io.TFRecordOptions.get_compression_type_string(
        io.make_tfrecord_options(filenames))
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
      io.make_tfrecord_options(filenames)

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
    actual = io.read_tfrecords(
        path, core_pb2.ContigInfo, max_records=max_records)
    self.assertLen(list(actual), expected_n)


class RawProtoWriterAdaptorTests(parameterized.TestCase):

  def setUp(self):
    self.proto1 = core_pb2.ContigInfo(name='p1', n_bases=10, pos_in_fasta=0)
    self.proto2 = core_pb2.ContigInfo(name='p2', n_bases=20, pos_in_fasta=1)
    self.protos = [self.proto1, self.proto2]

  @parameterized.parameters(
      dict(take_ownership=True),
      dict(take_ownership=False),
  )
  def test_adaptor_with_ownership(self, take_ownership):
    mock_writer = mock.MagicMock()
    adaptor = io.RawProtoWriterAdaptor(
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


class AsyncWriterTest(parameterized.TestCase):

  def test_thread_state(self):
    writer = io.AsyncWriter(lambda _: 0, maxsize=100)

    # Thread isn't started at creation.
    self.assertIsNone(writer._thread)

    with writer:
      # Once we enter the with block, our thread should exist and be flagged as
      # a daemon thread.
      self.assertIsNotNone(writer._thread)
      self.assertTrue(writer._thread.daemon)

    # Thread is None.
    self.assertIsNone(writer._thread)

  @parameterized.parameters(
      # Everything works even without any objects written to the queue.
      (
          [],),
      # Three simple checks for writing 1, 2, or 10 ints in order.
      (
          [1],),
      ([1, 2],),
      (list(range(10)),),
      # Checks that we can pass in multiple arguments to write_fn.
      (
          [(1, 2), (3, 4), (5, 6)],),
  )
  def test_writing(self, objs):
    # Make sure everything is a tuple for the *args of write.
    objs = [(obj,) if not isinstance(obj, tuple) else obj for obj in objs]

    # Create a writer and write each of objs into it.
    write_fn = mock.MagicMock()
    with io.AsyncWriter(write_fn, maxsize=100) as writer:
      for obj in objs:
        writer.write(*obj)

    # Tests that the calls to our write_fn is exactly those we expect in their
    # exact order.
    self.assertEqual([mock.call(*obj) for obj in objs], write_fn.call_args_list)

  def test_write_blocking_and_timeout(self):
    write_fn = mock.MagicMock()
    writer = io.AsyncWriter(write_fn, maxsize=1)

    # Queue is empty, so we succeed immediately.
    writer.write(1)

    # With a non-blocking write we'll immediately fail with a Full exception.
    with self.assertRaises(Queue.Full):
      writer.write(2, block=False)

    # We block and still get the Full exception. This test will fail if the
    # timeout isn't respected as we'll time out in the unittest itself.
    with self.assertRaises(Queue.Full):
      writer.write(3, block=True, timeout=1)

    # This construct starts the writer thread, so after we exit the with block
    # we'll have written all of the Queue elements.
    with writer:
      pass
    write_fn.assert_called_once_with(1)

  def test_flush_is_called_when_provided(self):
    write_fn = mock.MagicMock()
    flush_fn = mock.MagicMock()

    with io.AsyncWriter(write_fn, flush_fn=flush_fn, maxsize=1) as writer:
      writer.write(1)
      self.assertFalse(flush_fn.called)

    write_fn.assert_called_once_with(1)
    flush_fn.assert_called_once_with()

  @parameterized.parameters(
      (KeyError, 'write_fn'),
      (NameError, 'os is not defined'),
  )
  def test_exceptions_are_propogated(self, exception, text):
    write_fn = mock.MagicMock(side_effect=exception(text))

    # Just creating the writer doesn't actually trigger the exception, since
    # that only occurs in write_fn, which isn't run until we are in the
    # queue.
    writer = io.AsyncWriter(write_fn, maxsize=1)

    # We can write the writer, without being in the with block, and we don't
    # except out.
    writer.write(1)

    # We've already got 1 in the queue, so as soon as the thread starts it'll
    # except out. This code, with just a pass in the body, tests that we raise
    # an exception properly from the __exit__ function.
    with self.assertRaisesRegexp(exception, text):
      with writer:
        pass

    # The thread is no longer running, but the exception was already triggered.
    # So we want to make sure that a new writer.write call triggers the
    # exception again. This isn't a natural test, in the sense that we really
    # want a producer thread to get the exception from a write call while the
    # thread is actually running, but this is close enough without doing locking
    # etc gymnastics just for testing.
    with self.assertRaisesRegexp(exception, text):
      writer.write(2)

    # We already raised the exception to the producer via writer.write(2) call
    # above. This means that our next use of writer as a context manager should
    # not signal another exception from the __exit__ function.
    with writer:
      pass

  def test_exception_with_non_empty_queue_and_exit(self):
    # This test is for a very specific bug I observed using this code. If an
    # exception occurs during write, and we have elements in the queue, the
    # previous implement would hang when trying to exit the with block. This
    # checks for that exact situation.
    write_fn = mock.MagicMock(side_effect=KeyError('write_fn'))
    writer = io.AsyncWriter(write_fn)
    writer.write(1)
    writer.write(2)

    # There are now two element in the queue. We'll remove one and throw an
    # exception, leaving one element in the queue and the exception. We should
    # still except out immediately in the with block, rather than hang.
    with self.assertRaisesRegexp(KeyError, 'write_fn'):
      with writer:
        pass

  @parameterized.parameters(
      (KeyError, 'flush_fn'),
      (ImportError, 'bad import'),
  )
  def test_exception_from_flush_is_thrown_to_producer(self, exception, text):
    """Tests that we get an exception thrown by flush_fn."""
    write_fn = mock.MagicMock()
    flush_fn = mock.MagicMock(side_effect=exception(text))

    with self.assertRaisesRegexp(exception, text):
      with io.AsyncWriter(write_fn, flush_fn=flush_fn):
        pass


class ShardsTest(parameterized.TestCase):

  @parameterized.named_parameters(
      ('no_suffix', '/dir/foo/bar@3', '/dir/foo/bar', 3, None),
      ('suffix', '/dir/foo/bar@3.txt', '/dir/foo/bar', 3, 'txt'),
  )
  def testParseShardedFileSpec(self, spec, expected_basename,
                               expected_num_shards, expected_suffix):

    basename, num_shards, suffix = io.ParseShardedFileSpec(spec)
    self.assertEquals(basename, expected_basename)
    self.assertEquals(num_shards, expected_num_shards)
    self.assertEquals(suffix, expected_suffix)

  def testParseShardedFileSpecInvalid(self):
    self.assertRaises(io.ShardError, io.ParseShardedFileSpec, '/dir/foo/bar@0')

  @parameterized.named_parameters(
      ('no_suffix', '/dir/foo/bar@3', [
          '/dir/foo/bar-00000-of-00003', '/dir/foo/bar-00001-of-00003',
          '/dir/foo/bar-00002-of-00003'
      ]),
      ('suffix', '/dir/foo/bar@3.txt', [
          '/dir/foo/bar-00000-of-00003.txt', '/dir/foo/bar-00001-of-00003.txt',
          '/dir/foo/bar-00002-of-00003.txt'
      ]),
  )
  def testGenerateShardedFilenames(self, spec, expected):
    names = io.GenerateShardedFilenames(spec)
    self.assertEquals(names, expected)

  def testGenerateShardedFilenamesManyShards(self):
    names = io.GenerateShardedFilenames('/dir/foo/bar@100000')
    self.assertEquals(len(names), 100000)
    self.assertEquals(names[99999], '/dir/foo/bar-099999-of-100000')

  @parameterized.named_parameters(
      ('no_spec', '/dir/foo/bar'),
      ('zero_shards', '/dir/foo/bar@0'),
  )
  def testGenerateShardedFilenamesError(self, spec):
    self.assertRaises(io.ShardError, io.GenerateShardedFilenames, spec)

  @parameterized.named_parameters(
      ('basic', '/dir/foo/bar@3', True),
      ('suffix', '/dir/foo/bar@3,txt', True),
      ('many_shards', '/dir/foo/bar@123456', True),
      ('invalid_spec', '/dir/foo/bar@0', False),
      ('not_spec', '/dir/foo/bar', False),
  )
  def testIsShardedFileSpec(self, spec, expected):
    actual = io.IsShardedFileSpec(spec)
    self.assertEquals(actual, expected,
                      'io.IshShardedFileSpec({0}) is {1} expected {2}'.format(
                          spec, actual, expected))

  @parameterized.named_parameters(
      ('no_suffix', '/dir/foo/bar', 3, None, '/dir/foo/bar-?????-of-00003'),
      ('suffix', '/dir/foo/bar', 3, 'txt', '/dir/foo/bar-?????-of-00003.txt'),
      ('many', '/dir/foo/bar', 1234567, 'txt',
       '/dir/foo/bar-???????-of-1234567.txt'),
  )
  def testGenerateShardedFilePattern(self, basename, num_shards, suffix,
                                     expected):

    self.assertEquals(
        io.GenerateShardedFilePattern(basename, num_shards, suffix), expected)

  @parameterized.named_parameters(
      ('no_spec', '/dir/foo/bar', '/dir/foo/bar'),
      ('suffix', '/dir/foo/bar@3.txt', '/dir/foo/bar-?????-of-00003.txt'),
      ('no_suffix', '/dir/foo/bar@3', '/dir/foo/bar-?????-of-00003'),
      ('1000', '/dir/foo/bar@1000', '/dir/foo/bar-?????-of-01000'),
      ('many', '/dir/foo/bar@12345678', '/dir/foo/bar-????????-of-12345678'),
  )
  def testNormalizeToShardedFilePattern(self, spec, expected):
    self.assertEquals(expected, io.NormalizeToShardedFilePattern(spec))


if __name__ == '__main__':
  absltest.main()
