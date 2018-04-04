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

"""Generic utilities for IO.

Important: Please keep this module free of TensorFlow c++ extensions.
This makes it easy to build pure python packages for training that work with
CMLE.

redacted
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import heapq
import math
import os
import re



import contextlib2

from tensorflow.core.example import example_pb2
from tensorflow.python.lib.io import python_io

SHARD_SPEC_PATTERN = re.compile(R'((.*)\@(\d*[1-9]\d*)(?:\.(.+))?)')


class ShardError(Exception):
  """An IO error."""
  pass


def ParseShardedFileSpec(spec):  # pylint:disable=invalid-name
  """Parse a sharded file specification.

  Args:
    spec: The sharded file specification. A sharded file spec is one like
       gs://some/file@200.txt. In this case @200 specifies the number of shards.

  Returns:
    basename: The basename for the files.
    num_shards: The number of shards.
    suffix: The suffix if there is one, or '' if not.
  Raises:
    ShardError: If the spec is not a valid sharded specification.
  """
  m = SHARD_SPEC_PATTERN.match(spec)
  if not m:
    raise ShardError(('The file specification {0} is not a sharded file '
                      'specification because it did not match the regex '
                      '{1}').format(spec, SHARD_SPEC_PATTERN.pattern))

  # If there's a non-empty suffix, we need to prepend '.' so we get files like
  # foo@20.ext instead of foo@ext. The original C++ parser version has:
  # string ext = StrCat(suff.empty() ? "" : ".", suff);
  suffix = '.' + m.group(4) if m.group(4) else ''

  return m.group(2), int(m.group(3)), suffix


def _ShardWidth(num_shards):  # pylint:disable=invalid-name
  """Return the width of the shard matcher based on the number of shards."""
  return max(5, int(math.floor(math.log10(num_shards)) + 1))


def GenerateShardedFilenames(spec):  # pylint:disable=invalid-name
  """Generate the list of filenames corresponding to the sharding path.

  Args:
    spec: Sharding specification.

  Returns:
    List of filenames.

  Raises:
    ShardError: If spec is not a valid sharded file specification.
  """
  basename, num_shards, suffix = ParseShardedFileSpec(spec)
  files = []
  width = _ShardWidth(num_shards)
  format_str = '{{0}}-{{1:0{0}}}-of-{{2:0{0}}}{{3}}'.format(width)
  for i in range(num_shards):
    files.append(format_str.format(basename, i, num_shards, suffix))

  return files


def GenerateShardedFilePattern(basename, num_shards, suffix):  # pylint:disable=invalid-name
  """Generate a sharded file pattern.

  Args:
    basename: str; The basename for the files.
    num_shards: int; The number of shards.
    suffix: str; The suffix if there is one or ''.
  Returns:
    pattern:
  """
  width = _ShardWidth(num_shards)
  specifier = '?' * width
  format_str = '{{0}}-{{1}}-of-{{2:0{0}}}{{3}}'.format(width)
  return format_str.format(basename, specifier, num_shards, suffix)


def NormalizeToShardedFilePattern(spec_or_pattern):  # pylint:disable=invalid-name
  """Take a sharding spec or sharding file pattern and return a sharded pattern.

  The input can be a sharding spec(e.g '/some/file@10') or a sharded file
  pattern (e.g. '/some/file-?????-of-00010)

  Args:
    spec_or_pattern: A sharded file specification or sharded file pattern.

  Returns:
    sharded file pattern.
  """
  try:
    basename, num_shards, suffix = ParseShardedFileSpec(spec_or_pattern)
  except ShardError:
    return spec_or_pattern
  return GenerateShardedFilePattern(basename, num_shards, suffix)


def IsShardedFileSpec(spec):  # pylint:disable=invalid-name
  """Returns true if spec is a sharded file specification."""
  m = SHARD_SPEC_PATTERN.match(spec)
  return m is not None


VCF_EXTENSIONS = frozenset(['.vcf', '.vcf.gz'])


# redacted
def sharded_filename(spec, i):
  """Gets a path appropriate for writing the ith file of a sharded spec."""
  return GenerateShardedFilenames(spec)[i]


# redacted
# readability when there are multiple input filespecs.
def resolve_filespecs(shard, *filespecs):
  """Transforms potentially sharded filespecs into their paths for single shard.

  This function takes a shard number and a varargs potentially sharded
  filespecs, and returns a list where the filespecs have been resolved into
  concrete file paths for a single shard.

  This function has a concept of a master filespec, which is used to constrain
  and check the validity of other filespecs in filespecs. The first filespec is
  considered the master, and it cannot be None. For example, if master is not
  sharded, none of the other specs can be sharded, and vice versa. They
  must all also have a consistent sharding (e.g., master is @10, then all others
  must be @10).

  Note that filespecs (except the master) may be None or any other False value,
  which are returned as is in the output list.

  Args:
    shard: int >= 0. Our shard number.
    *filespecs: list[str]. Contains all of the filespecs we want to resolve into
      shard-specific file paths.

  Returns:
    A list. The first element is the number of shards, followed by the
    shard-specific paths for each filespec, in order.

  Raises:
    ValueError: if any filespecs are consistent.
  """
  if not filespecs:
    raise ValueError('filespecs must have at least one element.')

  master = filespecs[0]
  master_is_sharded = IsShardedFileSpec(master)

  master_num_shards = None
  if master_is_sharded:
    _, master_num_shards, _ = ParseShardedFileSpec(master)
    if shard >= master_num_shards or shard < 0:
      raise ValueError('Invalid shard={} value with master={} sharding', shard,
                       master)
  elif shard > 0:
    raise ValueError('Output is not sharded but shard > 0', shard)

  def resolve_one(filespec):
    """Resolves a single filespec into a concrete filepath."""
    if not filespec:
      return filespec

    is_sharded = IsShardedFileSpec(filespec)
    if master_is_sharded != is_sharded:
      raise ValueError('Master={} and {} have inconsistent sharding'.format(
          master, filespec))

    if not is_sharded:  # Not sharded => filespec is the actual filename.
      return filespec

    _, filespec_num_shards, _ = ParseShardedFileSpec(filespec)
    if filespec_num_shards != master_num_shards:
      raise ValueError('Master={} and {} have inconsistent sharding'.format(
          master, filespec))
    # return GenerateShardedFilename(basename, shard, filespec_shards)
    return sharded_filename(filespec, shard)

  return [master_num_shards] + [resolve_one(spec) for spec in filespecs]


def maybe_generate_sharded_filenames(filespec):
  """Potentially expands sharded filespec into a list of paths.

  This function takes in a potentially sharded filespec and expands it into a
  list containing the full set of corresponding concrete sharded file paths. If
  the input filespec is not sharded then a list containing just that file path
  is returned. This function is useful, for example, when the input to a binary
  can either be sharded or not.

  Args:
    filespec: String. A potentially sharded filespec to expand.

  Returns:
    A list of file paths.

  Raises:
    TypeError: if filespec is not valid basestring type.
  """
  if not isinstance(filespec, basestring):
    raise TypeError('Invalid filespec: %s' % filespec)
  if IsShardedFileSpec(filespec):
    return GenerateShardedFilenames(filespec)
  else:
    return [filespec]


def read_tfrecords(path, proto=None, max_records=None, options=None):
  """Yields the parsed records in tfrecord formatted file path.

  Note that path can be sharded filespec (path@N) in which case this function
  will read each shard in order.

  Args:
    path: String. A path to a tfrecord formatted file containing protos.
    proto: A proto class. proto.FromString() will be called on each serialized
      record in path to parse it.
    max_records: int >= 0 or None. Maximum number of records to read from path.
      If None, the default, all records will be read.
    options: A python_io.TFRecordOptions object for the reader.

  Yields:
    proto.FromString() values on each record in path in order.
  """
  if not proto:
    proto = example_pb2.Example

  if not options:
    options = make_tfrecord_options(path)

  if IsShardedFileSpec(path):
    paths = GenerateShardedFilenames(path)
  else:
    paths = [path]

  i = 0
  for path in paths:
    for buf in python_io.tf_record_iterator(path, options):
      i += 1
      if max_records is not None and i > max_records:
        return
      yield proto.FromString(buf)


def read_shard_sorted_tfrecords(path,
                                key,
                                proto=None,
                                max_records=None,
                                options=None):
  """Yields the parsed records in TFRecord-formatted file path in sorted order.

  The input TFRecord file must have each shard already in sorted order when
  using the key function for comparison (but elements can be interleaved across
  shards). Under those constraints, the elements will be yielded in a global
  sorted order.

  Args:
    path: String. A path to a TFRecord-formatted file containing protos.
    key: Callable. A function that takes as input a single instance of the proto
      class and returns a value on which the comparison for sorted ordering is
      performed.
    proto: A proto class. proto.FromString() will be called on each serialized
      record in path to parse it.
    max_records: int >= 0 or None. Maximum number of records to read from path.
      If None, the default, all records will be read.
    options: A python_io.TFRecordOptions object for the reader.

  Yields:
    proto.FromString() values on each record in path in sorted order.
  """
  if proto is None:
    proto = example_pb2.Example

  if options is None:
    options = make_tfrecord_options(path)

  if IsShardedFileSpec(path):
    paths = GenerateShardedFilenames(path)
  else:
    paths = [path]

  keyed_iterables = []
  for path in paths:
    protos = (
        proto.FromString(buf)
        for buf in python_io.tf_record_iterator(path, options))
    keyed_iterables.append(((key(elem), elem) for elem in protos))

  for i, (_, value) in enumerate(heapq.merge(*keyed_iterables)):
    if max_records is not None and i >= max_records:
      return
    yield value


def write_tfrecords(protos, output_path, options=None):
  """Writes protos to output_path.

  This function writes serialized strings of each proto in protos to output_path
  in their original order. If output_path is a sharded file (e.g., foo@2), this
  function will write the protos spread out as evenly as possible among the
  individual componented of the sharded spec (e.g., foo-00000-of-00002 and
  foo-00001-of-00002). Note that the order of records in the sharded files may
  differ from the order in protos due to the striping.

  Args:
    protos: An iterable of protobufs. The objects we want to write out.
    output_path: str. The filepath where we want to write protos.
    options: A python_io.TFRecordOptions object for the writer.
  """
  if not options:
    options = make_tfrecord_options(output_path)

  if IsShardedFileSpec(output_path):
    with contextlib2.ExitStack() as stack:
      _, n_shards, _ = ParseShardedFileSpec(output_path)
      writers = [
          stack.enter_context(
              make_tfrecord_writer(sharded_filename(output_path, i), options))
          for i in range(n_shards)
      ]
      for i, proto in enumerate(protos):
        writers[i % n_shards].write(proto.SerializeToString())
  else:
    with make_tfrecord_writer(output_path, options) as writer:
      for proto in protos:
        writer.write(proto.SerializeToString())


def make_tfrecord_options(filenames):
  """Creates a python_io.TFRecordOptions for the specified filename.

  Args:
    filenames: str or list[str]. A path or a list of paths where we'll
      read/write our TFRecord.
  Returns:
    A python_io.TFRecordOptions object.
  Raises:
    ValueError: If the filenames contain inconsistent file types.
  """
  # Backward compatibility: if input is one string, make it a list first.
  if not isinstance(filenames, list):
    filenames = [filenames]

  # If there are multiple file patterns, they have to be all the same type.
  extensions = set(os.path.splitext(filename)[1] for filename in filenames)
  if len(extensions) != 1:
    raise ValueError(
        'Incorrect value: {}. Filenames need to be all of the same type: '
        'either all with .gz or all without .gz'.format(','.join(filenames)))
  if extensions == {'.gz'}:
    compression_type = python_io.TFRecordCompressionType.GZIP
  else:
    compression_type = python_io.TFRecordCompressionType.NONE
  return python_io.TFRecordOptions(compression_type)


def make_tfrecord_writer(outfile, options=None):
  """Creates a python_io.TFRecordWriter for the specified outfile.

  Args:
    outfile: str. A path where we'll write our TFRecords.
    options: python_io.TFRecordOptions or None. If None, one
      will be inferred from the filename.
  Returns:
    A python_io.TFRecordWriter object.
  """
  if not options:
    options = make_tfrecord_options(outfile)
  return python_io.TFRecordWriter(outfile, options)


def make_proto_writer(outfile):
  """Creates a write to outfile writing general Protos.

  Args:
    outfile: A path to a file where we want to write protos.

  Returns:
    An writer object and a write_fn accepting a proto that writes to writer.
  """
  writer = make_tfrecord_writer(outfile)
  write_fn = lambda proto: writer.write(proto.SerializeToString())

  return writer, write_fn


# redacted
# this class wrapper.
class RawProtoWriterAdaptor(object):
  """Adaptor class wrapping a low-level bytes writer with a write(proto) method.

  This class provides a simple wrapper around low-level writers that accept
  serialized protobufs (via the SerializeToString()) for their write() methods.
  After wrapping this low-level writer will have a write(proto) method that
  accepts a protocol message `proto` and calls the low-level writer with
  `proto.SerializeToString()`. Given that many C++ writers require the proto
  to write properly (e.g., VCF writer) this allows us to provide a uniform API
  to clients who call write(proto) and either have that write call go directly
  to a type-specific writer or to a low-level writer via this
  RawProtoWriterAdaptor.
  """

  def __init__(self, raw_writer, take_ownership=True):
    """Creates a new RawProtoWriterAdaptor.

    Arguments:
      raw_writer: A low-level writer with a write() method accepting a
        serialized protobuf. Must also support __enter__ and __exit__ if
        take_ownership is True.
      take_ownership: bool. If True, we will all __enter__ and __exit__ on the
        raw_writer if/when this object's __enter__ and __exit__ are called. If
        False, no calls to these methods will be invoked on raw_writer.
    """
    self.raw_writer = raw_writer
    self.take_ownership = take_ownership

  def __enter__(self):
    if self.take_ownership:
      self.raw_writer.__enter__()
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    if self.take_ownership:
      self.raw_writer.__exit__(exc_type, exc_value, traceback)

  def write(self, proto):
    """Writes `proto.SerializeToString` to raw_writer."""
    self.raw_writer.write(proto.SerializeToString())
