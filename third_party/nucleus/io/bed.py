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
"""Classes for reading and writing BED files.

The BED format is described at
https://genome.ucsc.edu/FAQ/FAQformat.html#format1

API for reading:

```python
from third_party.nucleus.io import bed

# Iterate through all records.
with bed.BedReader(input_path) as reader:
  for record in reader:
    print(record)
```

where `record` is a `nucleus.genomics.v1.BedRecord` protocol buffer.

API for writing:

```python
from third_party.nucleus.io import bed
from third_party.nucleus.protos import bed_pb2

# records is an iterable of nucleus.genomics.v1.BedRecord protocol buffers.
records = ...

# header defines how many fields to write out.
header = bed_pb2.BedHeader(num_fields=5)

# Write all records to the desired output path.
with bed.BedWriter(output_path, header) as writer:
  for record in records:
    writer.write(record)
```

For both reading and writing, if the path provided to the constructor contains
'.tfrecord' as an extension, a `TFRecord` file is assumed and attempted to be
read or written. Otherwise, the filename is treated as a true BED file.

Files that end in a '.gz' suffix cause the file to be treated as compressed
(with BGZF if it is a true BED file, and with gzip if it is a TFRecord file).
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import bed_reader
from third_party.nucleus.io.python import bed_writer
from third_party.nucleus.protos import bed_pb2


class NativeBedReader(genomics_reader.GenomicsReader):
  """Class for reading from native BED files.

  Most users will want to use BedReader instead, because it dynamically
  dispatches between reading native BED files and TFRecord files based on the
  filename's extension.
  """

  def __init__(self, input_path, num_fields=0):
    """Initializes a NativeBedReader.

    Args:
      input_path: string. A path to a resource containing BED records.
      num_fields: int. The number of fields to read in the BED. If unset or set
        to zero, all fields in the input are read.
    """
    super(NativeBedReader, self).__init__()

    bed_path = input_path.encode('utf8')
    options = bed_pb2.BedReaderOptions(num_fields=num_fields)
    self._reader = bed_reader.BedReader.from_file(bed_path, options)
    self.header = self._reader.header

  def query(self):
    """Returns an iterator for going through the records in the region.

    NOTE: This function is not currently implemented by NativeBedReader though
    it could be implemented for sorted, tabix-indexed BED files.
    """
    raise NotImplementedError('Can not currently query a BED file')

  def iterate(self):
    """Returns an iterable of BedRecord protos in the file."""
    return self._reader.iterate()

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class BedReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading BedRecord protos from BED or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeBedReader(input_path, **kwargs)

  def _record_proto(self):
    return bed_pb2.BedRecord


class NativeBedWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native BED files.

  Most users will want BedWriter, which will write to either native BED
  files or TFRecord files, based on the output filename's extension.
  """

  def __init__(self, output_path, header=None):
    """Initializer for NativeBedWriter.

    Args:
      output_path: str. The path to which to write the BED file.
      header: nucleus.genomics.v1.BedHeader. The header that defines all
        information germane to the constituent BED records.
    """
    super(NativeBedWriter, self).__init__()
    if header is None:
      header = bed_pb2.BedHeader(num_fields=3)
    writer_options = bed_pb2.BedWriterOptions()
    self._writer = bed_writer.BedWriter.to_file(output_path, header,
                                                writer_options)

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class BedWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing BedRecord protos to BED or TFRecord files."""

  def _native_writer(self, output_path, header):
    return NativeBedWriter(output_path, header=header)
