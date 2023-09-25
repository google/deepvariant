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
"""Classes for reading and writing VCF files.

The VCF format is described at
https://samtools.github.io/hts-specs/VCFv4.3.pdf

API for reading:

```python
from third_party.nucleus.io import vcf

with vcf.VcfReader(input_path) as reader:
  for variant in reader:
    print(variant)
```

API for writing:

```python
from third_party.nucleus.io import vcf

# variants is an iterable of nucleus.genomics.v1.Variant protocol buffers.
variants = ...

with vcf.VcfWriter(output_path, header=header) as writer:
  for variant in variants:
    writer.write(variant)
```

The class attempts to infer the file format (`TFRecord` vs VCF) from the file
path provided to the constructor.

1. For files that end with '.tfrecord' and '.tfrecord.gz' (a gzipped version),
  a `TFRecord` file is assumed and is attempted to be read or written.

2. For all other files, the VCF format will be used.

  VCF format used in writing is inferred from file paths:
    - ending in '.bcf.gz': BGZF compressed BCF format will be written;
    - ending in '.bcf': uncompressed BCF format will be written;
    - ending in '.gz' and not in '.bcf.gz': BGZP compressed VCF format will be
        written;
    - all other suffixes: uncompressed VCF format will be written.

  VCF format used in reading is inferred from the contents of the file.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import genomics_writer
from third_party.nucleus.io.python import vcf_reader
from third_party.nucleus.io.python import vcf_writer
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

  NOTE: Users should not need to interact with this class at all. It is used
  by the variant_utils.{get,set}_info and variantcall_utils.{get,set}_format
  functions for interacting with the INFO and FORMAT fields in a Variant proto.
  """

  def __init__(self, header):
    """Initializer.

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
               excluded_info_fields=None,
               excluded_format_fields=None,
               store_gl_and_pl_in_info_map=False,
               header=None):
    """Initializer for NativeVcfReader.

    Args:
      input_path: str. The path to the VCF file to read.
      excluded_info_fields: list(str). A list of INFO field IDs that should not
        be parsed into the Variants. If None, all INFO fields are included.
      excluded_format_fields: list(str). A list of FORMAT field IDs that should
        not be parsed into the Variants. If None, all FORMAT fields are
        included.
      store_gl_and_pl_in_info_map: bool. If True, the "GL" and "PL" FORMAT
        fields are stored in the VariantCall.info map rather than as top-level
        values in the VariantCall.genotype_likelihood field.
      header: If not None, specifies the variants_pb2.VcfHeader. The file at
        input_path must not contain any header information.
    """
    super(NativeVcfReader, self).__init__()

    options = variants_pb2.VcfReaderOptions(
        excluded_info_fields=excluded_info_fields,
        excluded_format_fields=excluded_format_fields,
        store_gl_and_pl_in_info_map=store_gl_and_pl_in_info_map)
    if header is not None:
      self._reader = vcf_reader.VcfReader.from_file_with_header(
          input_path.encode('utf8'), options, header)
    else:
      self._reader = vcf_reader.VcfReader.from_file(
          input_path.encode('utf8'), options)

    self.header = self._reader.header
    self.field_access_cache = VcfHeaderCache(self.header)

  def iterate(self):
    """Returns an iterable of Variant protos in the file."""
    return self._reader.iterate()

  def query(self, region):
    """Returns an iterator for going through variants in the region."""
    return self._reader.query(region)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)

  @property
  def c_reader(self):
    """Returns the underlying C++ reader."""
    return self._reader


class VcfReader(genomics_reader.DispatchingGenomicsReader):
  """Class for reading Variant protos from VCF or TFRecord files."""

  def _native_reader(self, input_path, **kwargs):
    return NativeVcfReader(input_path, **kwargs)

  def _record_proto(self):
    return variants_pb2.Variant

  def _post_init_hook(self):
    # Initialize field_access_cache.  If we are dispatching to a
    # NativeVcfReader, we use its field_access_cache. Otherwise, we need to
    # create a new one.
    self.field_access_cache = getattr(
        self._reader, 'field_access_cache', VcfHeaderCache(self.header))

  @property
  def c_reader(self):
    """Returns the underlying C++ reader.

    Note that the C++ reader might be a VcfReader or it might be a
    TFRecordReader, depending on the input_path's extension.
    """
    return self._reader.c_reader


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
               excluded_format_fields=None,
               retrieve_gl_and_pl_from_info_map=False,
               exclude_header=False):
    """Initializer for NativeVcfWriter.

    Args:
      output_path: str. The path to which to write the VCF file.
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
      retrieve_gl_and_pl_from_info_map: bool. If True, the "GL" and "PL" FORMAT
        fields are retrieved from the VariantCall.info map rather than from the
        top-level value in the VariantCall.genotype_likelihood field.
      exclude_header: bool. If True, write a headerless VCF.
    """
    super(NativeVcfWriter, self).__init__()

    if header is None:
      header = variants_pb2.VcfHeader()
    writer_options = variants_pb2.VcfWriterOptions(
        round_qual_values=round_qualities,
        excluded_info_fields=excluded_info_fields,
        excluded_format_fields=excluded_format_fields,
        retrieve_gl_and_pl_from_info_map=retrieve_gl_and_pl_from_info_map,
        exclude_header=exclude_header,
    )
    self._writer = vcf_writer.VcfWriter.to_file(output_path, header,
                                                writer_options)
    self.field_access_cache = VcfHeaderCache(header)

  def write(self, proto):
    self._writer.write(proto)

  def write_somatic(self, proto):
    self._writer.write_somatic(proto)

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._writer.__exit__(exit_type, exit_value, exit_traceback)


class VcfWriter(genomics_writer.DispatchingGenomicsWriter):
  """Class for writing Variant protos to VCF or TFRecord files."""

  def _native_writer(self,
                     output_path,
                     header,
                     round_qualities=False,
                     excluded_info_fields=None,
                     excluded_format_fields=None,
                     retrieve_gl_and_pl_from_info_map=False,
                     exclude_header=False):
    return NativeVcfWriter(
        output_path,
        header=header,
        round_qualities=round_qualities,
        excluded_info_fields=excluded_info_fields,
        excluded_format_fields=excluded_format_fields,
        retrieve_gl_and_pl_from_info_map=retrieve_gl_and_pl_from_info_map,
        exclude_header=exclude_header)

  def _post_init_hook(self):
    # Initialize field_access_cache.  If we are dispatching to a
    # NativeVcfWriter, we use its field_access_cache.  Otherwise, we
    # need to create a new one.
    self.field_access_cache = getattr(
        self._writer, 'field_access_cache', VcfHeaderCache(self.header))


class InMemoryVcfReader(genomics_reader.GenomicsReader):
  """Class for "reading" Variant protos from an in-memory cache of variants.

  ```python
  from third_party.nucleus.io import vcf
  from third_party.nucleus.protos import variants_pb2

  variants = [... Variant protos ...]
  header = variants_pb2.VcfHeader()

  with vcf.InMemoryVcfReader(variants, header) as reader:
    for variant in reader:
      print(variant)
  ```

  This class accepts a collection of variants and optionally a header and
  provides all of the standard API functions of VcfReader but instead of
  fetching variants from a file the variants are queried from an in-memory cache
  of variant protos.

  Note that the input variants provided to this class aren't checked in any way,
  and their ordering determines the order of variants emitted by this class for
  the iterate() and query() operations. This is intentional, to make this class
  easy to use for testing where you often want to use less-than-perfectly formed
  inputs. In order to fully meet the contract of a standard VcfReader, variants
  should be sorted by their contig ordering and then by their start and finally
  by their ends.

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
      variants: list of nucleus.genomics.v1.Variant protos we will "read"
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
