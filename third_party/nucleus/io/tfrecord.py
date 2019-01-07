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
"""I/O for TFRecord files.

Utilities for reading and writing TFRecord files, especially those containing
serialized TensorFlow Example protocol buffers.
"""

# Important: Please keep this module free of TensorFlow c++ extensions.
# This makes it easy to build pure python packages for training that work with
# CMLE.

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import heapq
import os

import contextlib2

from third_party.nucleus.io import sharded_file_utils
from tensorflow.core.example import example_pb2
from tensorflow.python.lib.io import python_io


def read_tfrecords(path, proto=None, max_records=None, options=None):
  """Yields the parsed records in a TFRecord file path.

  Note that path can be sharded filespec (path@N) in which case this function
  will read each shard in order; i.e. shard 0 will read each entry in order,
  then shard 1, ...

  Args:
    path: String. A path to a TFRecord file containing protos.
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

  if sharded_file_utils.is_sharded_file_spec(path):
    paths = sharded_file_utils.generate_sharded_filenames(path)
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
  """Yields the parsed records in a TFRecord file path in sorted order.

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

  if sharded_file_utils.is_sharded_file_spec(path):
    paths = sharded_file_utils.generate_sharded_filenames(path)
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
  individual components of the sharded spec (e.g., foo-00000-of-00002 and
  foo-00001-of-00002). Note that the order of records in the sharded files may
  differ from the order in protos due to the striping.

  Args:
    protos: An iterable of protobufs. The objects we want to write out.
    output_path: str. The filepath where we want to write protos.
    options: A python_io.TFRecordOptions object for the writer.
  """
  if not options:
    options = make_tfrecord_options(output_path)

  if sharded_file_utils.is_sharded_file_spec(output_path):
    with contextlib2.ExitStack() as stack:
      _, n_shards, _ = sharded_file_utils.parse_sharded_file_spec(output_path)
      writers = [
          stack.enter_context(
              make_tfrecord_writer(sharded_file_utils.sharded_filename(
                  output_path, i), options))
          for i in range(n_shards)
      ]
      for i, proto in enumerate(protos):
        writers[i % n_shards].write(proto.SerializeToString())
  else:
    with make_tfrecord_writer(output_path, options) as writer:
      for proto in protos:
        writer.write(proto.SerializeToString())


def make_tfrecord_options(filenames):
  """Returns a python_io.TFRecordOptions for the specified filename.

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

  # If there are multiple file patterns, they all have to be the same type.
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
  """Returns a python_io.TFRecordWriter for the specified outfile.

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
  """Returns a writer capable of writing general Protos to outfile.

  Args:
    outfile: str. A path to a file where we want to write protos.

  Returns:
    A writer object and a write_fn accepting a proto that writes to writer.
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
  to write properly (e.g., VCF writer), this allows us to provide a uniform API
  to clients who call write(proto) and either have that write call go directly
  to a type-specific writer or to a low-level writer via this
  RawProtoWriterAdaptor.
  """

  def __init__(self, raw_writer, take_ownership=True):
    """Creates a new RawProtoWriterAdaptor.

    Arguments:
      raw_writer: A low-level writer with a write() method that accepts a
        serialized protobuf. Must also support __enter__ and __exit__ if
        take_ownership is True.
      take_ownership: bool. If True, we will call __enter__ and __exit__ on the
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
