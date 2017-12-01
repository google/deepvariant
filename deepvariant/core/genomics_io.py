# Copyright 2017 Google Inc.
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
"""Higher-level APIs for reading and writing genomics data."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re



from deepvariant.core import io_utils
from deepvariant.core.protos import core_pb2
from deepvariant.core.python import reference_fai
from deepvariant.core.python import sam_reader as sam_reader_
from deepvariant.core.python import vcf_reader as vcf_reader_
from deepvariant.core.python import vcf_writer

SHARD_SPEC_PATTERN = re.compile(R'((.*)\@(\d*[1-9]\d*)(?:\.(.+))?)')

VCF_EXTENSIONS = frozenset(['.vcf', '.vcf.gz'])
SAM_EXTENSIONS = frozenset(['.sam', '.bam'])


def make_ref_reader(reference_filename):
  """Creates an indexed GenomeReference for reference_filename."""
  return reference_fai.GenomeReferenceFai.from_file(
      reference_filename.encode('utf8'),
      reference_filename.encode('utf8') + '.fai')


def make_sam_reader(reads_source,
                    read_requirements=None,
                    use_index=True,
                    parse_aux_fields=False,
                    hts_block_size=None,
                    downsample_fraction=None,
                    random_seed=None):
  """Creates a SamReader for reads_source.

  This function creates a SAM/BAM reader from reads_source, configured by the
  additional arguments to this function. This returned reader object supports
  the iterate() and query(range) APIs.

  Args:
    reads_source: string. A path to a resource containing SAM/BAM records.
      Currently supports SAM text format and BAM binary format.
    read_requirements: optional ReadRequirement proto. If not None, this proto
      is used to control which reads are filtered out by the reader before they
      are passed to the client.
    use_index: optional bool, defaulting to True. If True, we will attempt to
      load an index file for reads_source to enable the query() API call. If
      True an index file must exist. If False, we will not attempt to load an
      index for reads_source, disabling the query() call.
    parse_aux_fields: optional bool. If False, the default, we will not parse
      the auxillary fields of the SAM/BAM records (see SAM spec for details).
      Parsing the aux fields is often unnecessary for many applications, and
      adds a significant parsing cost to access. If you need these aux fields,
      set parse_aux_fields to True and these fields will be parsed and populate
      the appropriate Read proto fields (e.g., read.info).
    hts_block_size: integer or None.  If None, will use the default htslib
      block size.  Otherwise, will configure the underlying block size of the
      underlying htslib file object.  Larger values (e.g. 1M) may be beneficial
      for reading remote files.
    downsample_fraction: None or float in the interval (0.0, 1.0]. If not None,
      this sam_reader will only keep each read with probability
      downsample_fraction, randomly.
    random_seed: None or int. The random seed to use with this sam reader, if
      needed. If None, a fixed random value will be assigned.

  Returns:
    A sam_reader object. The exact class implementing this API is not specified.

  Raises:
    ValueError: If downsample_fraction is not None and not in the interval
      (0.0, 1.0].
  """
  if reads_source.endswith('.tfbam'):
    # delay load tfbam_lib. This is a simple plugin mechanism.
    try:
      from tfbam_lib import tfbam_reader  # pylint: disable=g-import-not-at-top
      return tfbam_reader.make_sam_reader(
          reads_source,
          read_requirements=read_requirements,
          use_index=use_index,
          unused_block_size=hts_block_size,
          downsample_fraction=downsample_fraction,
          random_seed=random_seed)
    except ImportError:
      raise ImportError('tfbam_lib module not found, cannot read .tfbam files.')
  aux_field_handling = core_pb2.SamReaderOptions.SKIP_AUX_FIELDS
  if parse_aux_fields:
    aux_field_handling = core_pb2.SamReaderOptions.PARSE_ALL_AUX_FIELDS
  index_mode = (
      core_pb2.INDEX_BASED_ON_FILENAME
      if use_index else core_pb2.DONT_USE_INDEX)

  if downsample_fraction:
    if not 0.0 < downsample_fraction <= 1.0:
      raise ValueError('downsample_fraction must be in the interval (0.0, 1.0]',
                       downsample_fraction)

  if not random_seed:
    # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
    random_seed = 2928130004

  return sam_reader_.SamReader.from_file(
      reads_source.encode('utf8'),
      core_pb2.SamReaderOptions(
          read_requirements=read_requirements,
          index_mode=index_mode,
          aux_field_handling=aux_field_handling,
          hts_block_size=(hts_block_size or 0),
          downsample_fraction=downsample_fraction,
          random_seed=random_seed))


def make_vcf_reader(variants_source, use_index=True, include_likelihoods=False):
  """Creates an indexed VcfReader for variants_source."""
  if use_index:
    index_mode = core_pb2.INDEX_BASED_ON_FILENAME
  else:
    index_mode = core_pb2.DONT_USE_INDEX

  if include_likelihoods:
    desired_vcf_fields = core_pb2.OptionalVariantFieldsToParse()
  else:
    desired_vcf_fields = core_pb2.OptionalVariantFieldsToParse(
        exclude_genotype_quality=True, exclude_genotype_likelihood=True)

  return vcf_reader_.VcfReader.from_file(
      variants_source.encode('utf8'),
      core_pb2.VcfReaderOptions(
          index_mode=index_mode, desired_format_entries=desired_vcf_fields))


def make_vcf_writer(outfile, contigs, samples, filters):
  """Creates a VcfWriter.

  Args:
    outfile: str. A path where we'll write our VCF file.
    contigs: Iterable of learning.genomics.deepvariant.core.ContigInfo protobufs
      used to populate the contigs info in the VCF header.
    samples: Iterable of str. The name of the samples we will be writing to this
      VCF file. The order of samples provided here must be the same as the order
      of VariantCall elements in any Variant written to this writer.
    filters: Iterable of learning.genomics.deepvariant.core.VcfFilterInfo
      protos. Description of the filter fields that may occur in Variant protos
      written to this writer. Filters can include filter descriptions that never
      occur in any Variant proto, but all filter field values among all of the
      written Variant protos must be provided here.

  Returns:
    vcf_writer.VcfWriter.
  """
  writer_options = core_pb2.VcfWriterOptions(
      contigs=contigs, sample_names=samples, filters=filters)
  return vcf_writer.VcfWriter.to_file(outfile, writer_options)


def make_variant_writer(outfile, contigs, samples=None, filters=None):
  """Creates a writer to outfile writing variant Protos.

  This function creates an writer that accepts Variant protos and writes
  them to outfile. The type of the writer is determined by the extension of
  outfile. If it's one of VCF_EXTENSIONS, we will write out VCF
  records via make_vcf_writer. Otherwise we will write out TFRecord file of
  serialized Variant protos.

  Args:
    outfile: A path to a file where we want to write our variant calls.
    contigs: A list of the reference genome contigs for writers that need contig
      information.
    samples: A list of sample names we will be writing out. Can be None if the
      list of samples isn't required for the intended output type. Will raise an
      exception is None and is required.
    filters: A list of VcfFilterInfo protos defining the filter fields present
      in the to-be-written records. Can be None if the intended writer doesn't
      require this information, but will raise an exception if None and required
      by the destination writer.

  Returns:
    An writer object and a write_fn accepting a Variant proto that writes to
    writer.

  Raises:
    ValueError: If any of the optional arguments needed for the specific output
      type of outfile are missing.
  """
  if any(outfile.endswith(ext) for ext in VCF_EXTENSIONS):
    if samples is None:
      raise ValueError('samples must be provided for vcf output')
    if filters is None:
      raise ValueError('filters must be provided for vcf output')
    writer = make_vcf_writer(
        outfile, contigs=contigs, samples=samples, filters=filters)
    write_fn = writer.write
  else:
    writer = io_utils.make_tfrecord_writer(outfile)
    write_fn = lambda variant: writer.write(variant.SerializeToString())

  return writer, write_fn


def make_read_writer(outfile, contigs=None):
  """Creates a writer to outfile writing Read protos.

  This function creates an writer that accepts Read protos and writes them to
  outfile. The type of the writer is determined by the extension of outfile. If
  it's one of SAM_EXTENSIONS, we will write out SAM records via make_sam_writer.
  Otherwise we will write out TFRecord file of serialized Read protos.

  Args:
    outfile: A path to a file where we want to write our variant calls.
    contigs: A list of the reference genome contigs for writers that need contig
      information.

  Returns:
    An writer object and a write_fn accepting a Read proto that writes to
    writer.

  Raises:
    ValueError: If any of the optional arguments needed for the specific output
      type of outfile are missing.
  """
  if any(outfile.endswith(ext) for ext in SAM_EXTENSIONS):
    if contigs is None:
      raise ValueError('contigs must be provided for SAM/BAM output')
    raise NotImplementedError
  else:
    return io_utils.RawProtoWriterAdaptor(
        io_utils.make_tfrecord_writer(outfile))
