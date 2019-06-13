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
r"""Library to produce variant statistics from a VCF file."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import json
import numpy as np

import tensorflow as tf

from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

_VARIANT_STATS_COLUMNS = [
    'reference_name', 'position', 'reference_bases', 'alternate_bases',
    'variant_type', 'is_variant', 'is_transition', 'is_transversion', 'depth',
    'genotype_quality'
]

VariantStats = collections.namedtuple('VariantStats', _VARIANT_STATS_COLUMNS)

_VARIANT_SUMMARY_STATS_COLUMNS = [
    'record_count', 'variant_count', 'snv_count', 'insertion_count',
    'deletion_count', 'mnp_count', 'complex_count', 'depth_mean', 'depth_stdev',
    'gq_mean', 'gq_stdev'
]

BIALLELIC_SNP = 'Biallelic_SNP'
BIALLELIC_INSERTION = 'Biallelic_Insertion'
BIALLELIC_DELETION = 'Biallelic_Deletion'
BIALLELIC_MNP = 'Biallelic_MNP'
MULTIALLELIC_SNP = 'Multiallelic_SNP'
MULTIALLELIC_INSERTION = 'Multiallelic_Insertion'
MULTIALLELIC_DELETION = 'Multiallelic_Deletion'
MULTIALLELIC_COMPLEX = 'Multiallelic_Complex'
REFCALL = 'RefCall'


class VCFSummaryStats(
    collections.namedtuple('VCFSummaryStats', _VARIANT_SUMMARY_STATS_COLUMNS)):
  __slots__ = ()

  def asdict(self):
    return {k: getattr(self, k) for k in _VARIANT_SUMMARY_STATS_COLUMNS}


def _get_variant_type(variant):
  """Returns the type of variant as a string."""
  if variant_utils.is_variant_call(variant):
    biallelic = variant_utils.is_biallelic(variant)
    snp = variant_utils.is_snp(variant)
    insertion = all(
        variant_utils.is_insertion(variant.reference_bases, alt)
        for alt in variant.alternate_bases)
    deletion = all(
        variant_utils.is_deletion(variant.reference_bases, alt)
        for alt in variant.alternate_bases)

    if biallelic:
      if snp:
        return BIALLELIC_SNP
      elif insertion:
        return BIALLELIC_INSERTION
      elif deletion:
        return BIALLELIC_DELETION
      else:
        return BIALLELIC_MNP
    else:
      if snp:
        return MULTIALLELIC_SNP
      elif insertion:
        return MULTIALLELIC_INSERTION
      elif deletion:
        return MULTIALLELIC_DELETION
      else:
        return MULTIALLELIC_COMPLEX
  else:
    return REFCALL


def _tstv(variant, vtype):
  """Returns a pair of bools indicating Transition, Transversion status."""
  if vtype == BIALLELIC_SNP:
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
      is_variant=variant_utils.is_variant_call(variant),
      depth=variantcall_utils.get_format(
          variant_utils.only_call(variant), 'DP'),
      genotype_quality=variantcall_utils.get_gq(
          variant_utils.only_call(variant)))


def single_variant_stats(variants):
  return [get_variant_stats(v) for v in variants]


def summary_stats(single_stats):
  """Computes summary statistics for a set of variants.

  Args:
    single_stats: list of VariantStats objects.

  Returns:
    A VCFSummaryStats object.
  """
  transposed_records = zip(*single_stats)
  transposed_dict = dict(zip(_VARIANT_STATS_COLUMNS, transposed_records))

  return VCFSummaryStats(
      record_count=len(single_stats),
      variant_count=sum(transposed_dict['is_variant']),
      snv_count=sum([
          t == BIALLELIC_SNP or t == MULTIALLELIC_SNP
          for t in transposed_dict['variant_type']
      ]),
      insertion_count=sum([
          t == BIALLELIC_INSERTION or t == MULTIALLELIC_INSERTION
          for t in transposed_dict['variant_type']
      ]),
      deletion_count=sum([
          t == BIALLELIC_DELETION or t == MULTIALLELIC_DELETION
          for t in transposed_dict['variant_type']
      ]),
      mnp_count=sum(
          [t == BIALLELIC_MNP for t in transposed_dict['variant_type']]),
      complex_count=sum(
          [t == MULTIALLELIC_COMPLEX for t in transposed_dict['variant_type']]),
      depth_mean=np.mean(transposed_dict['depth']),
      depth_stdev=np.std(transposed_dict['depth']),
      gq_mean=np.mean(transposed_dict['genotype_quality']),
      gq_stdev=np.std(transposed_dict['genotype_quality']))


def variants_to_stats_json(variants):
  """Computes variant statistics of each variant.

  Args:
    variants: iterable(Variant).

  Returns:
    A tuple of (stats_json, summary_json), where state_json is a JSON of single
    variant statistics, and summary_json is a JSON of single sample statistics.
  """
  single_var_stats = single_variant_stats(variants)
  transposed_records = zip(*single_var_stats)
  transposed_dict = dict(zip(_VARIANT_STATS_COLUMNS, transposed_records))
  stats_json = json.dumps(
      transposed_dict, indent=4, sort_keys=True, separators=(',', ': '))

  summ_stats = summary_stats(single_var_stats)
  summary_json = json.dumps(
      summ_stats.asdict(), indent=4, sort_keys=True, separators=(',', ': '))

  return stats_json, summary_json


def write(stats, outfile):
  """Writes stats to the output file."""

  with tf.io.gfile.GFile(outfile, 'w') as writer:
    writer.write(stats)
