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
"""Class for "reading" Variant protos from an in-memory store of variants.

API for reading:
  variants = [... Variant protos ...]
  header = ... vcf_header ...
  with InMemoryVcfReader(variants, header) as reader:
    for variant in reader:
      process(reader.header, variant)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from deepvariant.util.io import genomics_reader
from deepvariant.util import ranges
from deepvariant.util import variant_utils


class InMemoryVcfReader(genomics_reader.GenomicsReader):
  """Class for "reading" Variant protos from an in-memory cache of variants.

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
    self.variants = list(variants)
    self.header = header
    super(InMemoryVcfReader, self).__init__()

  def iterate(self):
    return iter(self.variants)

  def query(self, region):
    return iter(
        variant for variant in self.variants
        if ranges.ranges_overlap(variant_utils.variant_range(variant), region)
    )

