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
"""Classes for reading and writing VCF files.

API for reading:
  with VcfReader(input_path) as reader:
    for variant in reader:
      process(reader.header, variant)

API for writing:
  with VcfWriter(output_path, header) as writer:
    for variant in variants:
      writer.write(variant)

where variant is a nucleus.genomics.v1.Variant protocol buffer.

If the path contains '.tfrecord', then a TFRecord file is assumed.
Otherwise, it is treated as a true VCF file.  In either case, an
extension of '.gz' will cause the file to be treated as compressed.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import vcf_reader
from third_party.nucleus.io.python import vcf_writer
from third_party.nucleus.protos import index_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import vcf_constants


def _create_get_fn_cache(fields):
  """Returns a dictionary from field to a callable that extracts its value."""
  return {
      field.id: vcf_constants.create_get_fn(field.type, field.number)
      for field in fields
  }


def _create_set_fn_cache(fields):
  """Returns a dictionary from field to a callable that sets its value."""
  return {field.id: vcf_constants.SET_FN_LOOKUP[field.type] for field in fields}


class VcfHeaderCache(object):
  """This class creates a cache of accessors to structured fields in Variants.

  The INFO and FORMAT fields within Variant protos are structured and typed,
  with types defined by the corresponding VCF header. This cache object provides
  provides {info,format}_field_{get,set}_fn functions that can be used to
  extract information from the structured Variant protos based on the types
  defined therein.

  Note: Users should not need to interact with this class at all. It is used
  by the variant_utils.{get,set}_info and variantcall_utils.{get,set}_format
  functions for interacting with the INFO and FORMAT fields in a Variant proto.
  """

  def __init__(self, header):
    """Constructor.

    Args:
      header: nucleus.genomics.v1.VcfHeader proto. Used to define the accessor
        functions needed.
    """
    if header is None:
      header = variants_pb2.VcfHeader()
    self._info_get_cache = _create_get_fn_cache(header.infos)
    self._info_set_cache = _create_set_fn_cache(header.infos)
    self._format_get_cache = _create_get_fn_cache(header.formats)
    self._format_set_cache = _create_set_fn_cache(header.formats)

  def info_field_get_fn(self, field_name):
    """Returns a callable that extracts the given INFO field based on its type.

    Args:
      field_name: str. The INFO field name of interest, e.g. 'AA', 'DB', 'AF'.

    Returns:
      A callable used to extract the given INFO field from a Variant proto.
    """
    return self._info_get_cache[field_name]

  def info_field_set_fn(self, field_name):
    """Returns a callable that sets the given INFO field based on its type."""
    return self._info_set_cache[field_name]

  def format_field_get_fn(self, field_name):
    """Returns a callable that gets the given FORMAT field based on its type."""
    return self._format_get_cache[field_name]

  def format_field_set_fn(self, field_name):
    """Returns a callable that sets the given FORMAT field based on its type."""
    return self._format_set_cache[field_name]


class NativeVcfReader(genomics_reader.GenomicsReader):
  """Class for reading from native VCF files.

  Most users will want to use VcfReader instead, because it dynamically
  dispatches between reading native VCF files and TFRecord files based
  on the filename's extensions.
  """

  def __init__(self,
               input_path,
               use_index=True,
               excluded_info_fields=None,
               excluded_format_fields=None):
    """Initializer for NativeVcfReader.

    Args:
      input_path: str. The path to the VCF file to read.
      use_index: bool. If True, the input is assumed to be bgzipped and tabix
        indexed, and the `query` functionality is supported.
      excluded_info_fields: list(str). A list of INFO field IDs that should not
        be parsed into the Variants. If None, all INFO fields are included.
      excluded_format_fields: list(str). A list of FORMAT field IDs that should
        not be parsed into the Variants. If None, all FORMAT fields are
        included.
    """
    super(NativeVcfReader, self).__init__()

    index_mode = index_pb2.INDEX_BASED_ON_FILENAME
    if not use_index:
      index_mode = index_pb2.DONT_USE_INDEX

    self._reader = vcf_reader.VcfReader.from_file(
        input_path.encode('utf8'),
        variants_pb2.VcfReaderOptions(
            index_mode=index_mode,
            excluded_info_fields=excluded_info_fields,
            excluded_format_fields=excluded_format_fields))

    self.header = self._reader.header
    self.field_access_cache = VcfHeaderCache(self.header)

  def iterate(self):
    return self._reader.iterate()

  def query(self, region):
    return self._reader.query(region)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class VcfReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading Variant protos from VCF or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeVcfReader(input_path, **kwargs)

  def _record_proto(self):
    return variants_pb2.Variant

  def _post_init_hook(self):
    # Initialize field_access_cache.  If we are dispatching to a
    # NativeVcfReader, we use its field_access_cache.  Otherwise, we
    # need to create a new one.
    self.field_access_cache = getattr(
        self._reader, 'field_access_cache', VcfHeaderCache(self.header))


class NativeVcfWriter(genomics_writer.GenomicsWriter):
  """Class for writing to native VCF files.

  Most users will want VcfWriter, which will write to either native VCF
  files or TFRecords files, based on the output filename's extensions.
  """

  def __init__(self,
               output_path,
               header=None,
               round_qualities=False,
               excluded_info_fields=None,
               excluded_format_fields=None):
    """Initializer for NativeVcfWriter.

    Args:
      output_path: str. A path where we'll write our VCF file.
      header: nucleus.genomics.v1.VcfHeader. The header that defines all
        information germane to the constituent variants. This includes contigs,
        FILTER fields, INFO fields, FORMAT fields, samples, and all other
        structured and unstructured header lines.
      round_qualities: bool. If True, the QUAL field is rounded to one point
        past the decimal.
      excluded_info_fields: list(str). A list of INFO field IDs that should not
        be written to the output. If None, all INFO fields are included.
      excluded_format_fields: list(str). A list of FORMAT field IDs that should
        not be written to the output. If None, all FORMAT fields are included.
    """
    super(NativeVcfWriter, self).__init__()

    if header is None:
      header = variants_pb2.VcfHeader()
    writer_options = variants_pb2.VcfWriterOptions(
        round_qual_values=round_qualities,
        excluded_info_fields=excluded_info_fields,
        excluded_format_fields=excluded_format_fields,
    )
    self._writer = vcf_writer.VcfWriter.to_file(output_path, header,
                                                writer_options)
    self.field_access_cache = VcfHeaderCache(header)

  def write(self, proto):
    self._writer.write(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class VcfWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing Variant protos to VCF or TFRecord files."""

  def _native_writer(self,
                     output_path,
                     header,
                     round_qualities=False,
                     excluded_info_fields=None,
                     excluded_format_fields=None):
    return NativeVcfWriter(
        output_path,
        header=header,
        round_qualities=round_qualities,
        excluded_info_fields=excluded_info_fields,
        excluded_format_fields=excluded_format_fields)

  def _post_init_hook(self):
    # Initialize field_access_cache.  If we are dispatching to a
    # NativeVcfWriter, we use its field_access_cache.  Otherwise, we
    # need to create a new one.
    self.field_access_cache = getattr(
        self._writer, 'field_access_cache', VcfHeaderCache(self.header))


class InMemoryVcfReader(genomics_reader.GenomicsReader):
  """Class for "reading" Variant protos from an in-memory cache of variants.

  API:
    variants = [... Variant protos ...]
    header = variants_pb2.VcfHeader()
    with InMemoryVcfReader(variants, header) as reader:
      for variant in reader:
        process(reader.header, variant)

  This class accepts a collection of variants and optionally a header and
  provides all of the standard API functions of VcfReader but instead of
  fetching variants from a file the variants are queried from an in-memory cache
  of variant protos.

  Note that the input variants provided to this class aren't checked in any way,
  and their ordering determines the order of variants emitted by this class for
  iterate and query operation. This is intentional, to make this class easy to
  use for testing where you often want to use less-than-perfectly formed inputs.
  In order to fully meet the contract of a standard VcfReader, variants should
  be sorted by their contig ordering and then by their start and finally by
  their ends.

  Implementation note:
    The current implementation will be very slow for query() if the provided
    cache of variants is large, as we do a O(n) search to collect all of the
    overlapping variants for each query. There are several straightforward
    optimizations to do if we need/want to scale this up. (a) sort the variants
    and use a binary search to find overlapping variants (b) partition the
    variants by contig, so we have dict[contig] => [variants on contig], which
    allows us to completely avoid considering any variants on any other contigs.
    Neither of these optimizations are worth it if len(variants) is small, but
    it may be worth considering if we want to use this functionality with a
    large number of variants.
  """

  def __init__(self, variants, header=None):
    """Creates a VCFReader backed by a collection of variants.

    Args:
      variants: list of third_party.nucleus.protos.Variant protos we will "read"
        from.
      header: a VCFHeader object to provide as a result to calls to self.header,
        or None, indicating that we don't have a header associated with this
        reader.
    """
    super(InMemoryVcfReader, self).__init__()
    self.variants = list(variants)
    self.header = header

  def iterate(self):
    return iter(self.variants)

  def query(self, region):
    return iter(
        variant for variant in self.variants
        if ranges.ranges_overlap(variant_utils.variant_range(variant), region)
    )
