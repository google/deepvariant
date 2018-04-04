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
"""Abstract base class for objects writing genomics data.

Most users will want to use a subclass of GenomicsWriter, with code like
the following:

  with GenomicsWriterSubClass(output_path, options) as writer:
    for proto in records:
      writer.write(proto)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import abc


from absl import logging

from tensorflow.python.lib.io import python_io


class GenomicsWriter(object):
  """Abstract base class for writing genomics data.

  A GenomicsWriter only has one method, write, which writes a single
  protocol buffer to a file.
  """

  __metaclass__ = abc.ABCMeta

  @abc.abstractmethod
  def write(self, proto):
    """Writes proto to the file.

    Args:
      proto:  A protocol buffer.
    """
    pass

  def __enter__(self):
    """Enter a `with` block."""
    return self

  @abc.abstractmethod
  def __exit__(self, unused_type, unused_value, unused_traceback):
    """Exit a `with` block.  Typically, this will close the file."""
    pass


class TFRecordWriter(GenomicsWriter):
  """A GenomicsWriter that writes to a TFRecord file."""

  def __init__(self, output_path, header=None):
    super(TFRecordWriter, self).__init__()

    compressed = output_path.endswith('.gz')
    options = python_io.TFRecordOptions(
        python_io.TFRecordCompressionType.GZIP if compressed else
        python_io.TFRecordCompressionType.NONE)
    self._writer = python_io.TFRecordWriter(output_path, options=options)
    self.header = header

  def write(self, proto):
    self._writer.write(proto.SerializeToString())

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class DispatchingGenomicsWriter(GenomicsWriter):
  """A GenomicsWriter that dispatches based on the file extension.

  If '.tfrecord' is present in the filename, a TFRecordWriter is used.
  Otherwise, a native writer is.

  Sub-classes of DispatchingGenomicsWriter must define a _native_writer()
  method.
  """

  def __init__(self, output_path, **kwargs):
    super(DispatchingGenomicsWriter, self).__init__()
    self.header = kwargs.get('header', None)

    if '.tfrecord' in output_path:
      self._writer = TFRecordWriter(output_path, header=self.header)
    else:
      self._writer = self._native_writer(output_path, **kwargs)
    logging.info('Writing %s with %s',
                 output_path, self._writer.__class__.__name__)
    self._post_init_hook()

  @abc.abstractmethod
  def _native_writer(self, output_path, **kwargs):
    """Returns a GenomicsWriter for writing the records `natively`.

    Args:
      output_path: The path to write the records to.
      **kwargs:  Zero or more keyword arguments.

    Returns:
      A GenomicsWriter.
    """

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)

  def _post_init_hook(self):
    """Hook for sub-classes to run code at the end of __init__."""
    pass
