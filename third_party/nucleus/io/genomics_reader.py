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
"""Classes that provide the interface for reading genomics data.

`GenomicsReader` defines the core API supported by readers, and is subclassed
directly or indirectly (via `DispatchingGenomicsReader`) for all concrete
implementations.

`TFRecordReader` is an implementation of the `GenomicsReader` API for reading
`TFRecord` files. This is usable for all data types when encoding data in
protocol buffers.

`DispatchingGenomicsReader` is an abstract class defined for convenience on top
of `GenomicsReader` that supports reading from either the native file format or
from `TFRecord` files of the corresponding protocol buffer used to encode data
of that file type. The input format assumed is dependent upon the filename of
the input data.

Concrete implementations for individual file types (e.g. BED, SAM, VCF, etc.)
reside in type-specific modules in this package. The instantiation of readers
may have reader-specific requirements documented there. General examples of the
`iterate()` and `query()` functionality are shown below.

```python
# Equivalent ways to iterate through all elements in a reader.
# 1. Using the reader itself as an iterable object.
kwargs = ...  # Reader-specific keyword arguments.
with GenomicsReaderSubClass(output_path, **kwargs) as reader:
  for proto in reader:
    do_something(reader.header, proto)

# 2. Calling the iterate() method of the reader explicitly.
with GenomicsReaderSubClass(output_path, **kwargs) as reader:
  for proto in reader.iterate():
    do_something(reader.header, proto)

# Querying for all elements within a specific region of the genome.
from third_party.nucleus.protos import range_pb2
region = range_pb2.Range(reference_name='chr1', start=10, end=20)

with GenomicsReaderSubClass(output_path, **kwargs) as reader:
  for proto in reader.query(region):
    do_something(reader.header, proto)
```
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import abc
import errno

from absl import logging
import six

from third_party.nucleus.io.python import tfrecord_reader


class GenomicsReader(six.Iterator):
  """Abstract base class for reading genomics data.

  In addition to the abstractmethods defined below, subclasses should
  also set a `header` member variable in their objects.
  """

  __metaclass__ = abc.ABCMeta

  @abc.abstractmethod
  def iterate(self):
    """Returns an iterator for going through all the file's records."""

  @abc.abstractmethod
  def query(self, region):
    """Returns an iterator for going through the records in the region.

    Args:
      region:  A nucleus.genomics.v1.Range.

    Returns:
      An iterator containing all and only records within the specified region.
    """

  def __enter__(self):
    """Enter a `with` block."""
    return self

  def __exit__(self, unused_type, unused_value, unused_traceback):
    """Exit a `with` block.  Typically, this will close the file."""

  def __init__(self):
    """Initializer."""
    # Some readers can only support one iterator at a time, so don't
    # create one now.  Rather, create it when needed in next().
    self.iterator = None

  def __iter__(self):
    """Allows users to use the object itself as an iterator."""
    return self.iterate()

  def __next__(self):
    """Allows users to use the object itself as an iterator."""
    if self.iterator is None:
      self.iterator = self.iterate()
    return six.next(self.iterator)


class TFRecordReader(GenomicsReader):
  """A GenomicsReader that reads protocol buffers from a TFRecord file.

  Example usage:
    reader = TFRecordReader('/tmp/my_file.tfrecords.gz',
                            proto=tensorflow.Example)
    for example in reader:
      process(example)

  Note that TFRecord files do not have headers, and do not need
  to be wrapped in a "with" block.
  """

  def __init__(self, input_path, proto, compression_type=None):
    """Initializes the TFRecordReader.

    Args:
      input_path:  The filename of the file to read.
      proto:  The protocol buffer type the TFRecord file is expected to
        contain.  For example, variants_pb2.Variant or reads_pb2.Read.
      compression_type:  Either 'ZLIB', 'GZIP', '' (uncompressed), or
        None.  If None, __init__ will guess the compression type based on
        the input_path's suffix.

    Raises:
      IOError: if there was any problem opening input_path for reading.
    """
    super(TFRecordReader, self).__init__()

    self.input_path = input_path
    self.proto = proto
    self.header = None

    if compression_type is None:
      compression_type = 'GZIP' if input_path.endswith('.gz') else ''

    self.reader = tfrecord_reader.TFRecordReader.from_file(
        input_path, compression_type)
    if self.reader is None:
      raise IOError(errno.EIO,
                    'Error trying to open %s for reading' % input_path)

  def iterate(self):
    """Returns an iterator for going through all the file's records."""
    while self.reader.get_next():
      yield self.proto.FromString(self.reader.get_record())

  def query(self, region):
    """Returns an iterator for going through the records in the region.

    NOTE: This function is not currently implemented by TFRecordReader as the
    TFRecord format does not provide a general mechanism for fast random access
    to elements in genome order.
    """
    raise NotImplementedError('Can not query TFRecord file')

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self.reader.close()

  @property
  def c_reader(self):
    """Returns the underlying C++ reader."""
    return self.reader


class DispatchingGenomicsReader(GenomicsReader):
  """A GenomicsReader that dispatches based on the file extension.

  If '.tfrecord' is present in the filename, a TFRecordReader is used.
  Otherwise, a native reader is.

  Subclasses of DispatchingGenomicsReader must define the following methods:
    * _native_reader()
    * _record_proto()
  """

  def __init__(self, input_path, **kwargs):
    super(DispatchingGenomicsReader, self).__init__()

    if '.tfrecord' in input_path:
      self._reader = TFRecordReader(
          input_path, proto=self._record_proto(),
          compression_type=kwargs.get('compression_type', None))
    else:
      # Remove compression_type, if present, from the arguments we pass to the
      # native reader.
      kwargs.pop('compression_type', None)
      self._reader = self._native_reader(input_path, **kwargs)
    logging.info('Reading %s with %s',
                 input_path, self._reader.__class__.__name__)
    self.header = getattr(self._reader, 'header', None)
    self._post_init_hook()

  @abc.abstractmethod
  def _native_reader(self, input_path, **kwargs):
    """Returns a GenomicsReader for reading the records `natively`.

    Args:
      input_path: The path to the native file to read.
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
    """Hook for subclasses to run code at the end of __init__."""
