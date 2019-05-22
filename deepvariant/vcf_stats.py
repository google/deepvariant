# Copyright 2019 Google LLC.
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
r"""Creates a JSON summary of variants from a VCF file."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import json

import tensorflow as tf

from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

_VARIANT_STATS_COLUMNS = [
    'reference_name', 'position', 'reference_bases', 'alternate_bases',
    'variant_type', 'is_transition', 'is_transversion', 'genotype_quality'
]

VariantStats = collections.namedtuple('VariantStats', _VARIANT_STATS_COLUMNS)


def _get_variant_type(variant):
  """Returns the type of variant as a string."""
  if variant_utils.is_biallelic(variant):
    if variant_utils.is_snp(variant):
      return 'SNP'
    elif variant_utils.is_insertion(variant.reference_bases,
                                    variant.alternate_bases[0]):
      return 'Insertion'
    elif variant_utils.is_deletion(variant.reference_bases,
                                   variant.alternate_bases[0]):
      return 'Deletion'

  return 'Other'


def _tstv(variant, vtype):
  """Returns a pair of bools indicating Transition, Transversion status."""
  if vtype == 'SNP':
    is_transition = variant_utils.is_transition(variant.reference_bases,
                                                variant.alternate_bases[0])
    is_transversion = not is_transition
  else:
    is_transition = is_transversion = False

  return is_transition, is_transversion


def get_variant_stats(variant):
  """Returns a VariantStats object corresponding to the input variant."""
  vtype = _get_variant_type(variant)
  is_transition, is_transversion = _tstv(variant, vtype)

  return VariantStats(
      reference_name=variant.reference_name,
      position=(variant.start + 1),
      reference_bases=variant.reference_bases,
      alternate_bases=list(variant.alternate_bases),
      variant_type=vtype,
      is_transition=is_transition,
      is_transversion=is_transversion,
      genotype_quality=variantcall_utils.get_gq(
          variant_utils.only_call(variant)))


def variants_to_stats_json(variants):
  """Computes variant statistics of each variant.

  Args:
    variants: iterable(Variant).

  Returns:
    A JSON representation of statistics for all variants.
  """
  records = [get_variant_stats(v) for v in variants]
  transposed_records = zip(*records)
  transposed_dict = dict(zip(_VARIANT_STATS_COLUMNS, transposed_records))
  stats_json = json.dumps(
      transposed_dict, indent=4, sort_keys=True, separators=(',', ': '))
  return stats_json


def write_json(json_of_variant_stats, outfile):
  """Writes the JSON representation of variant stats to the output file."""

  with tf.io.gfile.GFile(outfile, 'w') as writer:
    writer.write(json_of_variant_stats)
