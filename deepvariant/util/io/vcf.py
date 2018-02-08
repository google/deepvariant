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

API for writing:

  with VcfWriter(output_path, contigs, samples, filters) as writer:
    for proto in records:
      writer.write(proto)

where proto is a nucleus.genomics.v1.Variant protocol buffer.

If output_path contains '.vcf' as an extension, then a true VCF file
will be output.  Otherwise, a TFRecord file will be output.  In either
case, an extension of '.gz' will cause the output to be compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from deepvariant.util.io import genomics_writer
from deepvariant.util.protos import core_pb2
from deepvariant.util.python import vcf_writer


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
      contigs: Iterable of third_party.nucleus.util.ContigInfo protobufs
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
    return frozenset(['.vcf'])

  def _native_writer(self, output_path,
                     contigs, samples, filters, round_qualities=False):
    return NativeVcfWriter(
        output_path, contigs, samples, filters, round_qualities)
