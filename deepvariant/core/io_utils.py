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
"""Generic utilities for IO.

Important: Please keep this module free of TensorFlow c++ extensions.
This makes it easy to build pure python packages for training that work with
CMLE.

redacted
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math
import os
import Queue
import re
import threading



import contextlib2
import tensorflow as tf

from absl import logging

from tensorflow.core.example import example_pb2

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
    suffix: The suffix if there is one or None.
  Raises:
    ShardError: If the spec is not a valid sharded specification.
  """
  m = SHARD_SPEC_PATTERN.match(spec)
  if not m:
    raise ShardError(('The file specification {0} is not a sharded file '
                      'specification because it did not match the regex '
                      '{1}').format(spec, SHARD_SPEC_PATTERN.pattern))

  return m.group(2), int(m.group(3)), m.group(4)


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

  if not suffix:
    suffix = ''
  else:
    suffix = '.' + suffix
  files = []
  width = _ShardWidth(num_shards)
  format_str = '{{0}}-{{1:0{0}}}-of-{{2:0{0}}}{{3}}'.format(width)
  for i in range(num_shards):
    files.append(format_str.format(basename, i, num_shards, suffix))

  return files


def GenerateShardedFilePattern(basename, num_shards, suffix):  # pylint:disable=invalid-name
  """Generate a sharded file pattern.

  Args:
    basename: The basename for the files.
    num_shards: The number of shards.
    suffix: The suffix if there is one or None. Should not include '.'.
  Returns:
    pattern:
  """
  if suffix:
    suffix = '.' + suffix
  else:
    suffix = ''
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


# redacted
class OutputsWriter(object):
  """Manages all of the outputs of make_examples in a single place."""

  def __init__(self, options):
    self._writers = {}
    if options.candidates_filename:
      logging.info('Writing candidates to %s', options.candidates_filename)
      self._add_writer('candidates',
                       make_tfrecord_writer(options.candidates_filename))
    if options.examples_filename:
      logging.info('Writing examples to %s', options.examples_filename)
      self._add_writer('examples',
                       make_tfrecord_writer(options.examples_filename))
    if options.gvcf_filename:
      logging.info('Writing a gvcf to %s', options.gvcf_filename)
      # redacted
      # we can write out a VCF or a TFRecords as requested.
      self._add_writer('gvcfs', make_tfrecord_writer(options.gvcf_filename))

  def _add_writer(self, name, writer):
    self._writers[name] = writer

  def __enter__(self):
    """API function to support with syntax."""
    for writer in self._writers.itervalues():
      writer.__enter__()
    return self

  def __exit__(self, exception_type, exception_value, traceback):
    for writer in self._writers.itervalues():
      writer.__exit__(exception_type, exception_value, traceback)

  def write(self, writer_name, *protos):
    writer = self._writers.get(writer_name, None)
    if writer:
      for proto in protos:
        writer.write(proto.SerializeToString())


class AsyncWriter(object):
  """Asynchronous writer.

  An asynchronous writer is a wrapper around a standard writer that supports a
  simple producer + single output consumer thread for writing. The writer
  takes a write_fn argument, which will be called over and over with the
  elements written to the writer by the producer(s). When the writer is put into
  a with context, it starts up a daemon consumer thread that calls write_fn
  on each of the elements "enqueued" by calls to write(), in the same order that
  they were written by the producer.

  For example, we'd expect to use this code like:

  # Use an async writer to writer out serialized protobufs.
  with tf.python_io.TFRecordWriter(path) as sync_writer:
    with AsyncWriter(sync_writer.write) as writer:
      for pb in protos:
        writer.write(pb.SerializeToString())

  In another example, the writer below is just like the one above, but does the
  serialization of the protobuf in the writer thread, saving more main thread
  cycles. It illustrates why we take a function (write_fn) instead of a concrete
  writer class, as it allows us to conveniently move pre-processing code needed
  into the consumer thread:

  with tf.python_io.TFRecordWriter(path) as sync_writer:
    def write_candidate(candidate):
      sync_writer.write(candidate.SerializeToString())

    with AsyncWriter(write_candidate) as writer:
      for candidate in candidates:
        writer.write(candidate)

  An AsyncWriter also supports a flush_fn() that is called when the AsyncWriter
  exits from its with context block. For example, we can use the flush_fn to
  write out buffered elements:

  with tf.python_io.TFRecordWriter(path) as sync_writer:
    pending = []

    def write_pending()
      for p in pending:
        sync_writer.write(p)
      pending = []

    def write_candidate(item):
      if pending and pending[-1] != item:
        write_pending()
      pending.append(item)

    with AsyncWriter(write_candidate, flush_fn=write_pending) as write:
      for candidate in candidates:
        writer.write(candidate)

  Which 'uniq' the items being written and uses the flush_fn to write out the
  final pending items from the running pending list.

  Return values and exceptions:

  The return values of write_fn and flush_fn are ignored by AsyncWriter and so
  are not accessible to producer queues.

  The AsyncWriter does its best to communicate exceptions from write_fn and
  flush_fn as soon as possible. If write_fn raises an exception, this exception
  will be raised immediately on any subsequent calls to write(). In addition,
  if an exception in write_fn occurs with no calls to write() after, this
  exception will be raised by __exit__() when the context manager exits.
  Similiarly, an exception raised by flush_fn will be raised when it is called
  in __exit__() to the producer.

  Once an exception occurs in write_fn, the AsyncWriter is effectively dead.
  Subsequent calls to write(*args) will raise the write_fn() exception without
  ever adding *args to the queue. The write thread will continue to remove
  elements from the queue until empty, allowing the writer to close the with
  block.
  """

  def __init__(self, write_fn, flush_fn=None, maxsize=0):
    """Creates a new AsyncWriter.

    Args:
      write_fn: A function to called on each element(s) added with write().
      flush_fn: optional. Function taking no arguments that will be called once
        when after all elements have been added and written with write_fn
        after the __exit__ call to AsyncWriter.
      maxsize: int >= 0. Maximum number of elements to store in the Queue. If
        0, then the queue will have infinite capacity.

    Returns:
      A new AsyncWriter.
    """
    self._write_fn = write_fn
    self._flush_fn = flush_fn
    self._queue = Queue.Queue(maxsize=maxsize)  # Blocking FifoQueue.
    self._thread = None
    self._exception = None
    self._exception_raised_to_producer = False

  def write(self, *args, **kwargs):
    """Writes *args to this writer.

    This function enqueues the objects *args into this writer's queue for future
    asynchronous writing of *args via a write_fn(*args) call.

    Args:
      *args: A tuple of objects to enqueue for writing.
      **kwargs: keyword arguments. Due to constraints with keyword args and
        positional varargs, this argument has to be a generic dictionary. The
        only keyword args allows here are: 'block', which should be set to True,
        the default, to block on the enqueue until space if free in our queue,
        or False which won't block but can drop *args if the queue is full; and
        'timeout', the length of time in seconds we will block before failing to
        enqueue *args before we give up and raise a Queue.Full exception. If
        timeout is None, the default, we will block forever until the queue is
        full.

    Raises:
      Queue.Full: If either block is False and the queue is full, or block is
        True, timeout is >= 0, and timeout seconds have passed without a slot
        freeing up in the Queue.
      Exception: Can raise an exception of arbitrary type that occurred in a
        previous invocation of the write_fn. This exception may not apply to the
        most recently added *args, but to some previously written *args.
    """
    if self._exception:
      # If an exception occurred, raise it in write so the producer sees it as
      # soon as possible.
      self._exception_raised_to_producer = True
      raise self._exception  # pylint: disable=raising-bad-type
    # Python 2 doesn't support keyword args with defaults after *args.
    block = kwargs.pop('block', True)
    timeout = kwargs.pop('timeout', None)
    self._queue.put(args, block=block, timeout=timeout)

  def __enter__(self):
    """API function to support with syntax."""
    self._thread = threading.Thread(target=self._loop)
    self._thread.daemon = True
    self._thread.start()
    return self

  def __exit__(self, exception_type, exception_value, traceback):
    """API function to support with syntax.

    API docs are:
    https://docs.python.org/2/reference/datamodel.html#with-statement-context-managers

    Args:
      exception_type: see API docs.
      exception_value: see API docs.
      traceback: see API docs.

    Exception: Can raise an exception of arbitrary type that occurred in a
      previous invocation of the write_fn.
    """
    self._queue.join()
    self._thread = None  # Remove our thread.

    if self._exception and not self._exception_raised_to_producer:
      # We only raise the exception if we haven't already communicated it to the
      # producer.
      raise self._exception  # pylint: disable=raising-bad-type

    if self._flush_fn:
      self._flush_fn()

  def _loop(self):
    """Runs an endless loop calling write_fn on elements from our queue."""
    while True:
      # Keep grabbing elements from the queue, even if we have an exception.
      objs = self._queue.get()
      if not self._exception:
        # Only process the element if we don't have a pending exception.
        try:
          self._write_fn(*objs)
        except Exception as e:  # pylint: disable=broad-except
          self._exception = e
      # Must mark the task as done even if we have an exception so that future
      # calls to queue.join() will not deadlock.
      self._queue.task_done()


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
    options: A tf.python_io.TFRecordOptions object for the reader.

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
    for buf in tf.python_io.tf_record_iterator(path, options):
      i += 1
      if max_records is not None and i > max_records:
        return
      yield proto.FromString(buf)


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
    options: A tf.python_io.TFRecordOptions object for the writer.
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
  """Creates a tf.python_io.TFRecordOptions for the specified filename.

  Args:
    filenames: str or list[str]. A path or a list of paths where we'll
      read/write our TFRecord.
  Returns:
    A tf.python_io.TFRecordOptions object.
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
    compression_type = tf.python_io.TFRecordCompressionType.GZIP
  else:
    compression_type = tf.python_io.TFRecordCompressionType.NONE
  return tf.python_io.TFRecordOptions(compression_type)


def make_tfrecord_writer(outfile, options=None):
  """Creates a tf.python_io.TFRecordWriter for the specified outfile.

  Args:
    outfile: str. A path where we'll write our TFRecords.
    options: tf.python_io.TFRecordOptions or None. If None, one
      will be inferred from the filename.
  Returns:
    A tf.python_io.TFRecordWriter object.
  """
  if not options:
    options = make_tfrecord_options(outfile)
  return tf.python_io.TFRecordWriter(outfile, options)


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
