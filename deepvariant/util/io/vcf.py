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
"""Classes for reading and writing VCF files.

API for reading:
  with VcfReader(output_path) as reader:
    for variant in reader:
      process(reader.header, variant)

API for writing:

  with VcfWriter(output_path, contigs, samples, filters) as writer:
    for variant in variants:
      writer.write(variant)

where variant is a nucleus.genomics.v1.Variant protocol buffer.

If output_path contains '.vcf' as an extension, then a true VCF file
will be output.  Otherwise, a TFRecord file will be output.  In either
case, an extension of '.gz' will cause the output to be compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections



from deepvariant.util.io import genomics_reader
from deepvariant.util.io import genomics_writer
from deepvariant.util.genomics import variants_pb2
from deepvariant.util.protos import core_pb2
from deepvariant.util.python import vcf_reader
from deepvariant.util.python import vcf_writer

_VCF_EXTENSIONS = frozenset(['.vcf'])

# redacted
VcfHeader = collections.namedtuple(
    'VcfHeader', ['contigs', 'filters', 'samples'])


class NativeVcfReader(genomics_reader.GenomicsReader):
  """Class for reading from native VCF files.

  Most users will want to use VcfReader instead, because it dynamically
  dispatches between reading native VCF files and TFRecord files based
  on the filename's extensions.
  """

  def __init__(self, output_path, use_index=True, include_likelihoods=False):
    index_mode = core_pb2.INDEX_BASED_ON_FILENAME
    if not use_index:
      index_mode = core_pb2.DONT_USE_INDEX

    # redacted
    # list of strings.
    desired_vcf_fields = core_pb2.OptionalVariantFieldsToParse()
    if not include_likelihoods:
      desired_vcf_fields.exclude_genotype_quality = True
      desired_vcf_fields.exclude_genotype_likelihood = True

    self._reader = vcf_reader.VcfReader.from_file(
        output_path.encode('utf8'),
        core_pb2.VcfReaderOptions(
            index_mode=index_mode,
            desired_format_entries=desired_vcf_fields))

    # redacted
    self.header = VcfHeader(contigs=self._reader.Contigs(),
                            filters=self._reader.Filters(),
                            samples=self._reader.Samples())

    genomics_reader.GenomicsReader.__init__(self)

  def iterate(self):
    return self._reader.iterate()

  def query(self, region):
    return self._reader.query(region)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class VcfReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading Variant protos from VCF or TFRecord files."""

  def _get_extensions(self):
    return _VCF_EXTENSIONS

  def _native_reader(self, output_path, **kwargs):
    return NativeVcfReader(output_path, **kwargs)

  def _record_proto(self):
    return variants_pb2.Variant


class NativeVcfWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native VCF files.

  Most users will want VcfWriter, which will write to either native VCF
  files or TFRecords files, based on the output filename's extensions.
  """

  def __init__(self, output_path,
               contigs, samples, filters, round_qualities=False):
    """Initializer for NativeVcfWriter.

    Args:
      output_path: str. A path where we'll write our VCF file.
      contigs: Iterable of nucleus.genomics.v1.ContigInfo protobufs
        used to populate the contigs info in the VCF header.
      samples: Iterable of str. The name of the samples we will be writing to
        this VCF file. The order of samples provided here must be the same as
        the order of VariantCall elements in any Variant written to this writer.
      filters: Iterable of third_party.nucleus.util.VcfFilterInfo
        protos. Description of the filter fields that may occur in Variant
        protos written to this writer. Filters can include filter descriptions
        that never occur in any Variant proto, but all filter field values
        among all of the written Variant protos must be provided here.
      round_qualities: bool. If True, the QUAL field is rounded to one point
        past the decimal.
    """
    writer_options = core_pb2.VcfWriterOptions(
        contigs=contigs,
        sample_names=samples,
        filters=filters,
        round_qual_values=round_qualities)
    self._writer = vcf_writer.VcfWriter.to_file(output_path, writer_options)

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class VcfWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing Variant protos to VCF or TFRecord files."""

  def _get_extensions(self):
    return _VCF_EXTENSIONS

  def _native_writer(self, output_path,
                     contigs, samples, filters, round_qualities=False):
    return NativeVcfWriter(
        output_path, contigs, samples, filters, round_qualities)
