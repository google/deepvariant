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

"""Abstract base class for objects reading genomics data.

Most users will want to use a subclass of GenomicsReader, with code like
the following:

  with GenomicsReaderSubClass(output_path, options) as reader:
    for proto in reader:
      do_something(reader.header, proto)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import abc


from absl import logging

from tensorflow.python.lib.io import python_io


class GenomicsReader(object):
  """Abstract base class for reading genomics data.

  In addition to the abstractmethods defined below, sub-classes should
  also set a .header member variable in their objects.
  """

  __metaclass__ = abc.ABCMeta

  @abc.abstractmethod
  def iterate(self):
    """Returns an iterator for going through the file's records."""

  @abc.abstractmethod
  def query(self, region):
    """Returns an iterator for going through the records in the region.

    Args:
      region:  A nucleus.genomics.v1.Range.

    Returns:
      An iterator.
    """

  def __enter__(self):
    """Enter a `with` block."""
    return self

  def __exit__(self, unused_type, unused_value, unused_traceback):
    """Exit a `with` block.  Typically, this will close the file."""
    pass

  def __init__(self):
    """Allows users to use the object as an iterator."""
    # Some readers can only support one iterator at a time, so don't
    # create one now.  Rather, create it when needed in next().
    self.iterator = None

  def __iter__(self):
    """Allows users to use the object as an iterator."""
    return self.iterate()

  def next(self):
    """Allows users to use the object as an iterator."""
    if self.iterator is None:
      self.iterator = self.iterate()
    return self.iterator.next()


class TFRecordReader(GenomicsReader):
  """A GenomicsReader that reads from a TFRecord file.

  Example usage:
    reader = TFRecordReader('/tmp/my_file.tfrecords.gz',
                            proto=tensorflow.Example)
    for example in reader:
      process(example)

  Note that TFRecord files do not have headers, and do not need
  to be wrapped in a "with" block.
  """

  def __init__(self, input_path, proto, tf_options=None):
    """Initializes the TFRecordReader.

    Args:
      input_path:  The filename of the file to read.
      proto:  The protocol buffer type the TFRecord file is expected to
        contain.  For example, variants_pb2.Variant or reads_pb2.Read.
      tf_options:  A python_io.TFRecordOptions object.  If not set,
        __init__ will create one with the compression type based on
        whether input_path ends in '.gz' or not.
    """
    super(TFRecordReader, self).__init__()

    self.input_path = input_path
    self.proto = proto
    self.header = None

    if not tf_options:
      compressed = input_path.endswith('.gz')
      tf_options = python_io.TFRecordOptions(
          python_io.TFRecordCompressionType.GZIP if compressed else
          python_io.TFRecordCompressionType.NONE)
    self.tf_options = tf_options

  def iterate(self):
    # redacted
    for buf in python_io.tf_record_iterator(self.input_path, self.tf_options):
      yield self.proto.FromString(buf)

  def query(self, region):
    raise NotImplementedError('Can not query TFRecord file')

  def __exit__(self, exit_type, exit_value, exit_traceback):
    # tf_record_iterator closes the file when out of records.
    pass


class DispatchingGenomicsReader(GenomicsReader):
  """A GenomicsReader that dispatches based on the file extension.

  If '.tfrecord' is present in the filename, a TFRecordReader is used,
  otherwise a native reader is.

  Sub-classes of DispatchingGenomicsReader must define the following methods:
    * _native_reader()
    * _record_proto()
  """

  def __init__(self, input_path, **kwargs):
    super(DispatchingGenomicsReader, self).__init__()

    if '.tfrecord' in input_path:
      self._reader = TFRecordReader(input_path, proto=self._record_proto(),
                                    tf_options=kwargs.get('tf_options', None))
    else:
      # Remove tf_options, if present, from the arguments we pass to the
      # native reader.
      kwargs.pop('tf_options', None)
      self._reader = self._native_reader(input_path, **kwargs)
    logging.info('Reading %s with %s',
                 input_path, self._reader.__class__.__name__)
    self.header = getattr(self._reader, 'header', None)
    self._post_init_hook()

  @abc.abstractmethod
  def _native_reader(self, output_path, **kwargs):
    """Returns a GenomicsReader for writing the records `natively`.

    Args:
      output_path: The path to write the records to.
      **kwargs:  Zero or more keyword arguments.

    Returns:
      A GenomicsReader.
    """

  @abc.abstractmethod
  def _record_proto(self):
    """Returns the protocol buffer type used by this reader."""

  def iterate(self):
    return self._reader.iterate()

  def query(self, region):
    return self._reader.query(region)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)

  def _post_init_hook(self):
    """Hook for sub-classes to run code at the end of __init__."""
    pass
