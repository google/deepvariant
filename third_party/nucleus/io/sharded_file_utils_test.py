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

"""Tests for third_party.nucleus.io.sharded_file_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import sharded_file_utils as io
from third_party.nucleus.testing import test_utils


class IOTest(parameterized.TestCase):

  @parameterized.parameters(
      # Unsharded outputs pass through as expected.
      dict(task_id=0, filespecs=['foo.txt'], expected=[0, 'foo.txt']),
      dict(
          task_id=0,
          filespecs=['foo.txt', 'bar.txt'],
          expected=[0, 'foo.txt', 'bar.txt']),
      dict(
          task_id=0,
          filespecs=['bar.txt', 'foo.txt'],
          expected=[0, 'bar.txt', 'foo.txt']),
      # It's ok to have False values for other bindings.
      dict(
          task_id=0, filespecs=['foo.txt', None], expected=[0, 'foo.txt',
                                                            None]),
      dict(task_id=0, filespecs=['foo.txt', ''], expected=[0, 'foo.txt', '']),
      dict(
          task_id=0,
          filespecs=['foo@10.txt', None],
          expected=[10, 'foo-00000-of-00010.txt', None]),
      dict(
          task_id=0,
          filespecs=['foo@10.txt', ''],
          expected=[10, 'foo-00000-of-00010.txt', '']),
      # Simple check that master behaves as expected.
      dict(
          task_id=0,
          filespecs=['foo@10.txt', None],
          expected=[10, 'foo-00000-of-00010.txt', None]),
      dict(
          task_id=0,
          filespecs=['foo@10', None],
          expected=[10, 'foo-00000-of-00010', None]),
      dict(
          task_id=1,
          filespecs=['foo@10', None],
          expected=[10, 'foo-00001-of-00010', None]),
      dict(
          task_id=9,
          filespecs=['foo@10', None],
          expected=[10, 'foo-00009-of-00010', None]),
      # Make sure we handle sharding of multiple filespecs.
      dict(
          task_id=0,
          filespecs=['foo@10', 'bar@10', 'baz@10'],
          expected=[
              10, 'foo-00000-of-00010', 'bar-00000-of-00010',
              'baz-00000-of-00010'
          ]),
      dict(
          task_id=9,
          filespecs=['foo@10', 'bar@10', 'baz@10'],
          expected=[
              10, 'foo-00009-of-00010', 'bar-00009-of-00010',
              'baz-00009-of-00010'
          ]),
  )
  def test_resolve_filespecs(self, task_id, filespecs, expected):
    self.assertEqual(io.resolve_filespecs(task_id, *filespecs), expected)

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


class ShardsTest(parameterized.TestCase):

  @parameterized.named_parameters(
      ('no_suffix', '/dir/foo/bar@3', '/dir/foo/bar', 3, ''),
      ('suffix-dot', '/dir/foo/bar@3.txt', '/dir/foo/bar', 3, '.txt'),
  )
  def testParseShardedFileSpec(self, spec, expected_basename,
                               expected_num_shards, expected_suffix):

    basename, num_shards, suffix = io.parse_sharded_file_spec(spec)
    self.assertEqual(basename, expected_basename)
    self.assertEqual(num_shards, expected_num_shards)
    self.assertEqual(suffix, expected_suffix)

  def testParseShardedFileSpecInvalid(self):
    self.assertRaises(io.ShardError,
                      io.parse_sharded_file_spec, '/dir/foo/bar@0')

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
    names = io.generate_sharded_filenames(spec)
    self.assertEqual(names, expected)

  def testGenerateShardedFilenamesManyShards(self):
    names = io.generate_sharded_filenames('/dir/foo/bar@100000')
    self.assertEqual(len(names), 100000)
    self.assertEqual(names[99999], '/dir/foo/bar-099999-of-100000')

  @parameterized.named_parameters(
      ('no_spec', '/dir/foo/bar'),
      ('zero_shards', '/dir/foo/bar@0'),
  )
  def testGenerateShardedFilenamesError(self, spec):
    self.assertRaises(io.ShardError, io.generate_sharded_filenames, spec)

  @parameterized.named_parameters(
      ('basic', '/dir/foo/bar@3', True),
      ('suffix', '/dir/foo/bar@3,txt', True),
      ('many_shards', '/dir/foo/bar@123456', True),
      ('invalid_spec', '/dir/foo/bar@0', False),
      ('not_spec', '/dir/foo/bar', False),
  )
  def testIsShardedFileSpec(self, spec, expected):
    actual = io.is_sharded_file_spec(spec)
    self.assertEqual(actual, expected,
                     'io.IshShardedFileSpec({0}) is {1} expected {2}'.format(
                         spec, actual, expected))

  @parameterized.named_parameters(
      ('basic', '/dir/foo/bar-00001-of-00003', True),
      ('suffix', '/dir/foo/bar-00001-of-00003,txt', True),
      ('many_shards', '/dir/foo/bar-00000-of-12345', True),
      ('invalid_spec', '/dir/foo/bar-00000-of-00000', False),
      ('not_spec', '/dir/foo/bar', False),
      ('directory_not_match', '/dir/foo/bar-00001-of-10000/baz.txt', False),
      ('file_match', '/dir/foo/bar-00001-of-10000.baz.txt', True),
  )
  def testIsShardedFilename(self, filename, expected):
    actual = io.is_sharded_filename(filename)
    self.assertEqual(actual, expected,
                     'io.is_sharded_filename({0}) is {1} expected {2}'.format(
                         filename, actual, expected))

  @parameterized.named_parameters(
      ('no_suffix', '/dir/foo/bar', 3, '', '/dir/foo/bar-?????-of-00003'),
      ('suffix', '/dir/foo/bar', 3, '.txt', '/dir/foo/bar-?????-of-00003.txt'),
      ('many', '/dir/foo/bar', 1234567, '.txt',
       '/dir/foo/bar-???????-of-1234567.txt'),
  )
  def testGenerateShardedFilePattern(self, basename, num_shards, suffix,
                                     expected):

    self.assertEqual(io.generate_sharded_file_pattern(
        basename, num_shards, suffix), expected)

  @parameterized.named_parameters(
      ('no_spec', '/dir/foo/bar', '/dir/foo/bar'),
      ('suffix', '/dir/foo/bar@3.txt', '/dir/foo/bar-?????-of-00003.txt'),
      ('no_suffix', '/dir/foo/bar@3', '/dir/foo/bar-?????-of-00003'),
      ('1000', '/dir/foo/bar@1000', '/dir/foo/bar-?????-of-01000'),
      ('many', '/dir/foo/bar@12345678', '/dir/foo/bar-????????-of-12345678'),
  )
  def testNormalizeToShardedFilePattern(self, spec, expected):
    self.assertEqual(expected, io.normalize_to_sharded_file_pattern(spec))

  @parameterized.named_parameters(
      ('no_spec', 'no_spec', ['no_spec']),
      ('sharded', 'sharded@3', ['sharded-00000-of-00003',
                                'sharded-00001-of-00003',
                                'sharded-00002-of-00003']),
      ('wildcard1', '*.ext', ['cat.ext', 'dog.ext']),
      ('wildcard2', 'fo?bar', ['foobar']),
      ('comma_list', 'file1,file2,file3', ['file1', 'file2', 'file3']),
      ('mixed_list', 'mixed.*txt,mixed@1,mixed_file',
       ['mixed.1txt', 'mixed.2txt', 'mixed-00000-of-00001', 'mixed_file']),
      ('with_dups', 'with_dups*',
       ['with_dups.1txt', 'with_dups.2txt', 'with_dups-00000-of-00001',
        'with_dups']),
  )
  def testGlobListShardedFilePatterns(self, specs, expected_files):
    # First, create all expected_files so Glob will work later.
    expected_full_files = [test_utils.test_tmpfile(f, '')
                           for f in expected_files]
    # Create the full spec names. This one doesn't create the files.
    full_specs = ','.join(
        [test_utils.test_tmpfile(spec) for spec in specs.split(',')])
    self.assertEqual(sorted(set(expected_full_files)),
                     io.glob_list_sharded_file_patterns(full_specs))

  @parameterized.parameters(
      ('name3-00000-of-00001', 'name3', 0, 1, ''),
      ('name4-12-of-20.foo.bar', 'name4', 12, 20, '.foo.bar'),
      ('name5-123456-of-999999.baz', 'name5', 123456, 999999, '.baz'),
      ('dir/name6.xxx-01111-of-02222.yyy', 'dir/name6.xxx', 1111, 2222, '.yyy'),
      ('/dir/foo/bar-00001-of-10000.baz.txt', '/dir/foo/bar', 1, 10000,
       '.baz.txt'),
  )
  def testParseShardedFilename(self, filename, expected_basename,
                               expected_shard, expected_num_shards,
                               expected_suffix):
    (basename, shard, num_shards,
     suffix) = io.parse_sharded_filename(filename)
    self.assertEqual(expected_basename, basename)
    self.assertEqual(expected_shard, int(shard))
    self.assertEqual(expected_num_shards, int(num_shards))
    self.assertEqual(expected_suffix, suffix)

if __name__ == '__main__':
  absltest.main()
