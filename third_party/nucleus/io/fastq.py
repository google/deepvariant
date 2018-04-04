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
"""Class for reading FASTQ files.

API for reading:
  with FastqReader(input_path) as reader:
    for record in reader:
      print(record)

API for writing:
  with FastqWriter(output_path) as writer:
    for record in records:
      writer.write(record)

where `record` is a nucleus.genomics.v1.FastqRecord protocol buffer.

If the path contains '.tfrecord' as an extension, a TFRecord file is
assumed.  Otherwise, it is treated as a true FASTQ file.  In either case,
an extension of '.gz' will cause the file to be treated as compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import fastq_reader
from third_party.nucleus.io.python import fastq_writer
from third_party.nucleus.protos import fastq_pb2


class NativeFastqReader(genomics_reader.GenomicsReader):
  """Class for reading from native FASTQ files.

  Most users will want to use FastqReader instead, because it dynamically
  dispatches between reading native FASTQ files and TFRecord files based on the
  filename's extension.
  """

  def __init__(self, input_path):
    """Initializes a NativeFastqReader.

    Args:
      input_path: string. A path to a resource containing FASTQ records.
    """
    super(NativeFastqReader, self).__init__()

    fastq_path = input_path.encode('utf8')
    if fastq_path.endswith('.gz'):
      options = fastq_pb2.FastqReaderOptions(
          compression_type=fastq_pb2.FastqReaderOptions.GZIP)
    else:
      options = fastq_pb2.FastqReaderOptions()
    self._reader = fastq_reader.FastqReader.from_file(fastq_path, options)
    self.header = None

  def query(self):
    raise NotImplementedError('Can not query a FASTQ file')

  def iterate(self):
    """Returns an iterable of FastqRecord protos in the file."""
    return self._reader.iterate()

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class FastqReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading FastqRecord protos from FASTQ or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeFastqReader(input_path, **kwargs)

  def _record_proto(self):
    return fastq_pb2.FastqRecord


class NativeFastqWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native FASTQ files.

  Most users will want FastqWriter, which will write to either native FASTQ
  files or TFRecord files, based on the output filename's extension.
  """

  def __init__(self, output_path):
    """Initializer for NativeFastqWriter.

    Args:
      output_path: str. The path to which to write the FASTQ file.
    """
    super(NativeFastqWriter, self).__init__()

    writer_options = fastq_pb2.FastqWriterOptions()
    self._writer = fastq_writer.FastqWriter.to_file(output_path, writer_options)

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class FastqWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing FastqRecord protos to FASTQ or TFRecord files."""

  def _native_writer(self, output_path):
    return NativeFastqWriter(output_path)
