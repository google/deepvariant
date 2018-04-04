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
"""Classes for reading and writing SAM and BAM files.

API for reading:
  with SamReader(input_path) as reader:
    for read in reader:
      process(reader.header, read)

API for writing:

  with SamWriter(output_path) as writer:
    for read in reads:
      writer.write(read)

where read is a nucleus.genomics.v1.Read protocol buffer.

If the path contains '.tfrecord', a TFRecord file is assumed; otherwise
it is treated as a true SAM file.  Also, file names ending with '.gz'
are assumed to be compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import sam_reader
from third_party.nucleus.protos import index_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils


class NativeSamReader(genomics_reader.GenomicsReader):
  """Class for reading from native SAM files.

  Most users will want to use SamReader instead, because it dynamically
  dispatches between reading native SAM files and TFRecord files based
  on the filename's extensions.
  """

  def __init__(self, input_path,
               use_index=True,
               read_requirements=None,
               parse_aux_fields=False,
               hts_block_size=None,
               downsample_fraction=None,
               random_seed=None):
    """Initializes a NativeSamReader.

    Args:
      input_path: string. A path to a resource containing SAM/BAM records.
        Currently supports SAM text format and BAM binary format.
      use_index: optional bool, defaulting to True. If True, we will attempt to
        load an index file for reads_source to enable the query() API call. If
        True an index file must exist. If False, we will not attempt to load an
        index for reads_source, disabling the query() call.
      read_requirements: optional ReadRequirement proto. If not None, this proto
        is used to control which reads are filtered out by the reader before
        they are passed to the client.
      parse_aux_fields: optional bool. If False, the default, we will not parse
        the auxillary fields of the SAM/BAM records (see SAM spec for details).
        Parsing the aux fields is often unnecessary for many applications, and
        adds a significant parsing cost to access. If you need these aux fields,
        set parse_aux_fields to True and these fields will be parsed and
        populate the appropriate Read proto fields (e.g., read.info).
      hts_block_size: integer or None.  If None, will use the default htslib
        block size.  Otherwise, will configure the underlying block size of the
        underlying htslib file object.  Larger values (e.g. 1M) may be
        beneficial for reading remote files.
      downsample_fraction: None or float in the interval [0.0, 1.0]. If not
        None or 0.0, the reader will only keep each read with probability
        downsample_fraction, randomly.
      random_seed: None or int. The random seed to use with this sam reader, if
        needed. If None, a fixed random value will be assigned.

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
            use_index=use_index,
            unused_block_size=hts_block_size,
            downsample_fraction=downsample_fraction,
            random_seed=random_seed)
      except ImportError:
        raise ImportError(
            'tfbam_lib module not found, cannot read .tfbam files.')
    else:
      index_mode = index_pb2.INDEX_BASED_ON_FILENAME
      if not use_index:
        index_mode = index_pb2.DONT_USE_INDEX

      aux_field_handling = reads_pb2.SamReaderOptions.SKIP_AUX_FIELDS
      if parse_aux_fields:
        aux_field_handling = reads_pb2.SamReaderOptions.PARSE_ALL_AUX_FIELDS

      if downsample_fraction:
        if not 0.0 < downsample_fraction <= 1.0:
          raise ValueError(
              'downsample_fraction must be in the interval (0.0, 1.0]',
              downsample_fraction)

      if random_seed is None:
        # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
        random_seed = 2928130004

      self._reader = sam_reader.SamReader.from_file(
          input_path.encode('utf8'),
          reads_pb2.SamReaderOptions(
              read_requirements=read_requirements,
              index_mode=index_mode,
              aux_field_handling=aux_field_handling,
              hts_block_size=(hts_block_size or 0),
              downsample_fraction=downsample_fraction,
              random_seed=random_seed))

      self.header = self._reader.header

    super(NativeSamReader, self).__init__()

  def iterate(self):
    return self._reader.iterate()

  def query(self, region):
    return self._reader.query(region)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class SamReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading Read protos from SAM or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeSamReader(input_path, **kwargs)

  def _record_proto(self):
    return reads_pb2.Read


class NativeSamWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native SAM files.

  Most users will want SamWriter, which will write to either native SAM
  files or TFRecords files, based on the output filename's extensions.
  """

  def __init__(self, output_path, header):
    """Initializer for NativeSamWriter.

    Args:
      output_path: str. A path where we'll write our SAM/BAM file.
      header: A nucleus.SamHeader protobuf.  The header is used both
        for writing the header, and to control the sorting applied to
        the rest of the file.
    """
    raise NotImplementedError

  def write(self, proto):
    raise NotImplementedError

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class SamWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing Variant protos to SAM or TFRecord files."""

  def _native_writer(self, output_path, header):
    return NativeSamWriter(output_path, header)


class InMemorySamReader(object):
  """Python interface class for in-memory sam reader.

  Attributes:
    reads: [third_party.nucleus.protos.Read], the list of in-memory reads.
    is_sorted: bool, if reads are sorted.
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
      An iterator over third_party.nucleus.protos.Read.
    """
    return self.reads

  def query(self, region):
    """Iterate over records overlapping a query region.

    Args:
      region: third_party.nucleus.protos.Range, query region.

    Returns:
      An iterator over third_party.nucleus.protos.Read
    """
    # redacted
    return (read for read in self.reads
            if ranges.ranges_overlap(region, utils.read_range(read)))
