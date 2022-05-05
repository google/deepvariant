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
# pylint: disable=line-too-long
"""Classes for reading and writing SAM and BAM files.

The SAM/BAM/CRAM formats are described at
https://samtools.github.io/hts-specs/SAMv1.pdf
https://samtools.github.io/hts-specs/CRAMv3.pdf

API for reading:

```python
from third_party.nucleus.io import sam

with sam.SamReader(input_path) as reader:
  for read in reader:
    print(read)
```

where `read` is a `nucleus.genomics.v1.Read` protocol buffer. input_path will
dynamically decode the underlying records depending the file extension, with
`.sam` for SAM files, `.bam` for BAM files, and `.cram` for CRAM files. It will
also search for an appropriate index file to use to enable calls to the
`query()` method.

API for writing SAM/BAM:

```python
from third_party.nucleus.io import sam

# reads is an iterable of nucleus.genomics.v1.Read protocol buffers.
reads = ...

with sam.SamWriter(output_path, header=header) as writer:
  for read in reads:
    writer.write(read)
```

API for writing CRAM:

```python
# ref_path is required for writing CRAM files. If embed_ref, the output CRAM
# file will embed reference sequences.
with sam.SamWriter(output_path, header=header, ref_path=ref_path,
                   embed_ref=embed_ref) as writer:
  for read in reads:
    writer.write(read)
```

For both reading and writing, if the path provided to the constructor contains
'.tfrecord' as an extension, a `TFRecord` file is assumed and attempted to be
read or written. Otherwise, the filename is treated as a true SAM/BAM/CRAM file.

For `TFRecord` files, ending in a '.gz' suffix causes the file to be treated as
compressed with gzip.

Notes on using CRAM with SamReader
--------------------------------

Nucleus supports reading from CRAM files using the same API as for SAM/BAM:

```python
from third_party.nucleus.io import sam

with sam.SamReader("/path/to/sample.cram") as reader:
  for read in reader:
    print(read)
```

There is one type of CRAM file, though, that has a slightly more complicated
API. If the CRAM file uses read sequence compression with an external reference
file, and this reference file is no longer accessible in the location specified
by the CRAM file's "UR" tag and cannot be found in the local genome cache, its
location must be passed to SamReader via the ref_path parameter:

```python
from third_party.nucleus.io import sam

cram_path = "/path/to/sample.cram"
ref_path = "/path/to/genome.fasta"
with sam.SamReader(cram_path, ref_path=ref_path) as reader:
  for read in reader:
    print(read)
```

Unfortunately, htslib is unable to load the ref_path from anything other than a
POSIX filesystem. (htslib plugin filesystems like S3 or GCS buckets won't work).
For that reason, we don't recommend the use of CRAM files with external
reference files, but instead suggest using read sequence compression with
embedded reference data. (This has a minor impact on file size, but
significantly improves file access simplicity and safety.)

For more information about CRAM, see:
* The `samtools` documentation at http://www.htslib.org/doc/samtools.html
* The "Global Options" section of the samtools docs at http://www.htslib.org/doc/samtools.html#GLOBAL_OPTIONS
* How reference sequences are encoded in CRAM at http://www.htslib.org/doc/samtools.html#REFERENCE_SEQUENCES
* Finally, benchmarking of different CRAM options http://www.htslib.org/benchmarks/CRAM.html
"""
# pylint: enable=line-too-long

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import sam_reader
from third_party.nucleus.io.python import sam_writer
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils


class NativeSamReader(genomics_reader.GenomicsReader):
  """Class for reading from native SAM/BAM/CRAM files.

  Most users will want to use SamReader instead, because it dynamically
  dispatches between reading native SAM/BAM/CRAM files and TFRecord files based
  on the filename's extensions.
  """

  def __init__(self,
               input_path,
               ref_path=None,
               read_requirements=None,
               parse_aux_fields=False,
               hts_block_size=None,
               downsample_fraction=None,
               random_seed=None,
               use_original_base_quality_scores=False,
               aux_fields_to_keep=None):
    """Initializes a NativeSamReader.

    Args:
      input_path: str. A path to a resource containing SAM/BAM/CRAM records.
        Currently supports SAM text format, BAM binary format, and CRAM.
      ref_path: optional str or None. Only used for CRAM decoding, and only
        necessary if the UR encoded path in the CRAM itself needs to be
        overridden. If provided, we will tell the CRAM decoder to use this FASTA
        for the reference sequence.
      read_requirements: optional ReadRequirement proto. If not None, this proto
        is used to control which reads are filtered out by the reader before
        they are passed to the client.
      parse_aux_fields: optional bool, defaulting to False. If False, we do not
        parse the auxiliary fields of the SAM/BAM/CRAM records (see SAM spec for
        details). Parsing the aux fields is unnecessary for many applications,
        and adds a significant parsing cost to access. If you need these aux
        fields, set parse_aux_fields to True and these fields will be parsed and
        populate the appropriate Read proto fields (e.g., read.info).
      hts_block_size: int or None. If specified, this configures the block size
        of the underlying htslib file object. Larger values (e.g. 1M) may be
        beneficial for reading remote files. If None, the reader uses the
        default htslib block size.
      downsample_fraction: float in the interval [0.0, 1.0] or None. If
        specified as a positive float, the reader will only keep each read with
        probability downsample_fraction, randomly. If None or zero, all reads
        are kept.
      random_seed: None or int. The random seed to use with this sam reader, if
        needed. If None, a fixed random value will be assigned.
      use_original_base_quality_scores: optional bool, defaulting to False. If
        True, quality scores are read from OQ tag.
      aux_fields_to_keep: None or list[str]. If None, we keep all aux fields if
        they are parsed. If set, we only keep the aux fields with the names in
        this list.

    Raises:
      ValueError: If downsample_fraction is not None and not in the interval
        (0.0, 1.0].
      ImportError: If someone tries to load a tfbam file.
    """
    if input_path.endswith('.tfbam'):
      # Delayed loading of tfbam_lib.
      try:
        from tfbam_lib import tfbam_reader  # pylint: disable=g-import-not-at-top
        self._reader = tfbam_reader.make_sam_reader(
            input_path,
            read_requirements=read_requirements,
            unused_block_size=hts_block_size,
            downsample_fraction=downsample_fraction,
            random_seed=random_seed)
      except ImportError:
        raise ImportError(
            'tfbam_lib module not found, cannot read .tfbam files.')
    else:
      aux_field_handling = reads_pb2.SamReaderOptions.SKIP_AUX_FIELDS
      if parse_aux_fields:
        aux_field_handling = reads_pb2.SamReaderOptions.PARSE_ALL_AUX_FIELDS

      # We make 0 be a valid value that means "keep all reads" so that proto
      # defaults (=0) do not omit all reads.
      if downsample_fraction is not None and downsample_fraction != 0:
        if not 0.0 < downsample_fraction <= 1.0:
          raise ValueError(
              'downsample_fraction must be in the interval (0.0, 1.0]',
              downsample_fraction)

      if random_seed is None:
        # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
        random_seed = 2928130004

      self._reader = sam_reader.SamReader.from_file(
          input_path.encode('utf8'),
          ref_path.encode('utf8') if ref_path is not None else '',
          reads_pb2.SamReaderOptions(
              read_requirements=read_requirements,
              aux_field_handling=aux_field_handling,
              aux_fields_to_keep=aux_fields_to_keep,
              hts_block_size=(hts_block_size or 0),
              downsample_fraction=downsample_fraction,
              random_seed=random_seed,
              use_original_base_quality_scores=use_original_base_quality_scores)
      )

      self.header = self._reader.header

    super(NativeSamReader, self).__init__()

  def iterate(self):
    """Returns an iterable of Read protos in the file."""
    return self._reader.iterate()

  def query(self, region):
    """Returns an iterator for going through the reads in the region."""
    return self._reader.query(region)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class SamReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading Read protos from SAM/BAM/CRAM or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeSamReader(input_path, **kwargs)

  def _record_proto(self):
    return reads_pb2.Read


class NativeSamWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native SAM/BAM/CRAM files.

  Most users will want SamWriter, which will write to either native SAM/BAM/CRAM
  files or TFRecords files, based on the output filename's extensions.
  """

  def __init__(self, output_path, header, ref_path=None, embed_ref=False):
    """Initializer for NativeSamWriter.

    Args:
      output_path: str. A path where we'll write our SAM/BAM/CRAM file.
      ref_path: str. Path to the reference file. Required for CRAM file.
      embed_ref: bool. Whether to embed the reference sequences in CRAM file.
        Default is False.
      header: A nucleus.SamHeader proto.  The header is used both for writing
        the header, and to control the sorting applied to the rest of the file.
    """
    super(NativeSamWriter, self).__init__()
    self._writer = sam_writer.SamWriter.to_file(
        output_path,
        ref_path.encode('utf8') if ref_path is not None else '', embed_ref,
        header)

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class SamWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing Read protos to SAM or TFRecord files."""

  def _native_writer(self, output_path, **kwargs):
    return NativeSamWriter(output_path, **kwargs)


class InMemorySamReader(object):
  """Python interface class for in-memory SAM/BAM/CRAM reader.

  Attributes:
    reads: list[nucleus.genomics.v1.Read]. The list of in-memory reads.
    is_sorted: bool, True if reads are sorted.
  """

  def __init__(self, reads, is_sorted=False):
    self.replace_reads(reads, is_sorted=is_sorted)

  def replace_reads(self, reads, is_sorted=False):
    """Replace the reads stored by this reader."""
    self.reads = reads
    self.is_sorted = is_sorted

  def iterate(self):
    """Iterate over all records in the reads.

    Returns:
      An iterator over nucleus.genomics.v1.Read's.
    """
    return self.reads

  def query(self, region):
    """Returns an iterator for going through the reads in the region.

    Args:
      region: nucleus.genomics.v1.Range. The query region.

    Returns:
      An iterator over nucleus.genomics.v1.Read protos.
    """
    # TODO: Add a faster query version for sorted reads.
    return (
        read for read in self.reads if utils.read_overlaps_region(read, region))
