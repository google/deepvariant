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
"""Classes for reading and writing GFF files.

The GFF format is described at
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md.

API for reading:

```python
from third_party.nucleus.io import gff

# Iterate through all records.
with gff.GffReader(input_path) as reader:
  for record in reader:
    print(record)
```

where `record` is a `nucleus.genomics.v1.GffRecord` protocol buffer.

API for writing:

```python
from third_party.nucleus.io import gff
from third_party.nucleus.protos import gff_pb2

# records is an iterable of nucleus.genomics.v1.GffRecord protocol buffers.
records = ...

header = gff_pb2.GffHeader()

# Write all records to the desired output path.
with gff.GffWriter(output_path, header) as writer:
  for record in records:
    writer.write(record)
```

For both reading and writing, if the path provided to the constructor contains
'.tfrecord' as an extension, a `TFRecord` file is assumed and attempted to be
read or written. Otherwise, the filename is treated as a true GFF file.

Files that end in a '.gz' suffix cause the file to be treated as compressed
(with BGZF if it is a true GFF file, and with gzip if it is a TFRecord file).
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import gff_reader
from third_party.nucleus.io.python import gff_writer
from third_party.nucleus.protos import gff_pb2


class NativeGffReader(genomics_reader.GenomicsReader):
  """Class for reading from native GFF files.

  Most users will want to use GffReader instead, because it dynamically
  dispatches between reading native GFF files and TFRecord files based on the
  filename's extension.
  """

  def __init__(self, input_path):
    """Initializes a NativeGffReader.

    Args:
      input_path: string. A path to a resource containing GFF records.
    """
    super(NativeGffReader, self).__init__()
    gff_path = input_path.encode('utf8')
    reader_options = gff_pb2.GffReaderOptions()
    self._reader = gff_reader.GffReader.from_file(gff_path, reader_options)
    self.header = self._reader.header

  def query(self):
    """Returns an iterator for going through the records in the region.

    NOTE: This function is not currently implemented by NativeGffReader though
    it could be implemented for sorted, tabix-indexed GFF files.
    """
    raise NotImplementedError('Can not currently query a GFF file')

  def iterate(self):
    """Returns an iterable of GffRecord protos in the file."""
    return self._reader.iterate()

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class GffReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading GffRecord protos from GFF or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeGffReader(input_path, **kwargs)

  def _record_proto(self):
    return gff_pb2.GffRecord


class NativeGffWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native GFF files.

  Most users will want GffWriter, which will write to either native GFF
  files or TFRecord files, based on the output filename's extension.
  """

  def __init__(self, output_path, header):
    """Initializer for NativeGffWriter.

    Args:
      output_path: str. The path to which to write the GFF file.
      header: nucleus.genomics.v1.GffHeader. The header that defines all
        information germane to the constituent GFF records.
    """
    super(NativeGffWriter, self).__init__()
    writer_options = gff_pb2.GffWriterOptions()
    self._writer = gff_writer.GffWriter.to_file(output_path, header,
                                                writer_options)

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class GffWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing GffRecord protos to GFF or TFRecord files."""

  def _native_writer(self, output_path, header):
    return NativeGffWriter(output_path, header=header)
