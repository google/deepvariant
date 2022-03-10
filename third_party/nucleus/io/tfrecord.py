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

# Important: Please keep this module free of TensorFlow C++ extensions.
# This makes it easy to build pure python packages for training that work with
# CMLE.

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import heapq

import contextlib2

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io import sharded_file_utils
from tensorflow.core.example import example_pb2


# pylint: disable=invalid-name
def Reader(path, proto=None, compression_type=None):
  """A TFRecordReader that defaults to tf.Example protos."""
  if not proto:
    proto = example_pb2.Example

  return genomics_reader.TFRecordReader(
      path, proto, compression_type=compression_type)


def Writer(path, compression_type=None):
  """A convenience wrapper around genomics_writer.TFRecordWriter."""
  return genomics_writer.TFRecordWriter(path, compression_type=compression_type)


# pylint: enable=invalid-name


# TODO: Refactor all of the following (internal).
def read_tfrecords(path, proto=None, max_records=None, compression_type=None):
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
    compression_type: 'GZIP', 'ZLIB', '' (uncompressed), or None to autodetect
      based on file extension.

  Yields:
    proto.FromString() values on each record in path in order.
  """
  if sharded_file_utils.is_sharded_file_spec(path):
    paths = sharded_file_utils.generate_sharded_filenames(path)
  else:
    paths = [path]

  i = 0
  for path in paths:
    with Reader(path, proto, compression_type=compression_type) as reader:
      for record in reader.iterate():
        i += 1
        if max_records is not None and i > max_records:
          return
        yield record


def expanded_paths_if_sharded(path):
  """Returns all file paths for the given (optionally sharded) tfrecord path.

  Args:
      path: String. Path to the (optionally sharded) TFRecords.  Returns a lits
        with the original path if the tfrecord file is not sharded, or a list of
        paths to all shards if it is.
  """
  if sharded_file_utils.is_sharded_file_spec(path):
    paths = sharded_file_utils.generate_sharded_filenames(path)
  else:
    paths = [path]
  return paths


def read_shard_sorted_tfrecords(path,
                                key,
                                proto=None,
                                max_records=None,
                                compression_type=None):
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
    compression_type: 'GZIP', 'ZLIB', '' (uncompressed), or None to autodetect
      based on file extension.

  Yields:
    proto.FromString() values on each record in path in sorted order.
  """
  paths = expanded_paths_if_sharded(path)

  keyed_iterables = []
  for path in paths:
    protos = Reader(path, proto, compression_type=compression_type).iterate()
    keyed_iterables.append(((key(elem), elem) for elem in protos))

  for i, (_, value) in enumerate(heapq.merge(*keyed_iterables)):
    if max_records is not None and i >= max_records:
      return
    yield value


def write_tfrecords(protos, output_path, compression_type=None):
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
    compression_type: 'GZIP', 'ZLIB', '' (uncompressed), or None to autodetect
      based on file extension.
  """
  if sharded_file_utils.is_sharded_file_spec(output_path):
    with contextlib2.ExitStack() as stack:
      _, n_shards, _ = sharded_file_utils.parse_sharded_file_spec(output_path)
      writers = [
          stack.enter_context(
              Writer(
                  sharded_file_utils.sharded_filename(output_path, i),
                  compression_type=compression_type)) for i in range(n_shards)
      ]
      for i, proto in enumerate(protos):
        writers[i % n_shards].write(proto)
  else:
    with Writer(output_path, compression_type=compression_type) as writer:
      for proto in protos:
        writer.write(proto)
