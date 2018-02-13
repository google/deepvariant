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
"""Classes for reading and writing SAM and BAM files.

API for writing:

  with SamWriter(output_path) as writer:
    for read in reads:
      writer.write(read)

where read is a nucleus.genomics.v1.Read protocol buffer.

If output_path contains '.sam' as an extension, then a true SAM file
will be output.  Otherwise, a TFRecord file will be output.  In either
case, an extension of '.gz' will cause the output to be compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from deepvariant.util.io import genomics_writer


class NativeSamWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native SAM files.

  Most users will want SamWriter, which will write to either native SAM
  files or TFRecords files, based on the output filename's extensions.
  """

  def __init__(self, output_path, contigs):
    """Initializer for NativeSamWriter.

    Args:
      output_path: str. A path where we'll write our SAM/BAM file.
      contigs: Iterable of third_party.nucleus.util.ContigInfo protobufs
        used to populate the contigs info in the SAM header.
    """
    raise NotImplementedError

  def write(self, proto):
    raise NotImplementedError

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class SamWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing Variant protos to SAM or TFRecord files."""

  def _get_extensions(self):
    return frozenset(['.bam', '.sam'])

  def _native_writer(self, output_path, contigs):
    return NativeSamWriter(output_path, contigs)
