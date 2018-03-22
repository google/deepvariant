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

  def __init__(self, output_path):
    super(TFRecordWriter, self).__init__()

    compressed = output_path.endswith('.gz')
    options = python_io.TFRecordOptions(
        python_io.TFRecordCompressionType.GZIP if compressed else
        python_io.TFRecordCompressionType.NONE)
    self._writer = python_io.TFRecordWriter(output_path, options=options)

  def write(self, proto):
    self._writer.write(proto.SerializeToString())

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class DispatchingGenomicsWriter(GenomicsWriter):
  """A GenomicsWriter that dispatches based on the file extension.

  Sub-classes of DispatchingGenomicsWriter must define _get_extensions()
  and _native_writer() methods.  Using those methods, the object will then
  use the native writer if one of the extensions is present in the
  output filename, and use a TFRecordWriter otherwise.
  """

  def __init__(self, output_path, **kwargs):
    super(DispatchingGenomicsWriter, self).__init__()

    if '.tfrecord' not in output_path and any(
        ext in output_path for ext in self._get_extensions()):
      self._writer = self._native_writer(output_path, **kwargs)
    else:
      self._writer = TFRecordWriter(output_path)
    logging.info('Writing %s with %s',
                 output_path, self._writer.__class__.__name__)

  @abc.abstractmethod
  def _get_extensions(self):
    """Returns a list of file extensions to use the native writer for.

    Returns:
      A list of strings.
    """

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
