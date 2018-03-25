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

where `record` is a nucleus.genomics.v1.FastqRecord protocol buffer.

If input_path ends with '.gz', it is assumed to be compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections


from deepvariant.util.io import genomics_reader
from deepvariant.util.genomics import fastq_pb2
from deepvariant.util.python import fastq_reader

_FASTQ_EXTENSIONS = frozenset(['.fq', '.fastq'])


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

  def _get_extensions(self):
    return _FASTQ_EXTENSIONS

  def _native_reader(self, input_path, **kwargs):
    return NativeFastqReader(input_path, **kwargs)

  def _record_proto(self):
    return fastq_pb2.FastqRecord
