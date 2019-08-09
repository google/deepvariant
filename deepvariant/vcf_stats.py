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
import itertools
import json
import numpy as np

import tensorflow as tf

from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from deepvariant import vcf_stats_vis

_VARIANT_STATS_COLUMNS = [
    'reference_name', 'position', 'reference_bases', 'alternate_bases',
    'variant_type', 'is_variant', 'is_transition', 'is_transversion', 'depth',
    'genotype_quality', 'genotype', 'vaf', 'qual'
]

VariantStats = collections.namedtuple('VariantStats', _VARIANT_STATS_COLUMNS)

_VARIANT_SUMMARY_STATS_COLUMNS = [
    'record_count', 'variant_count', 'snv_count', 'insertion_count',
    'deletion_count', 'mnp_count', 'complex_count', 'depth_mean', 'depth_stdev',
    'gq_mean', 'gq_stdev', 'transition_count', 'transversion_count'
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


def _get_vaf(variant, vcf_reader):
  """Gets the VAF (variant allele frequency)."""
  vafs = variantcall_utils.get_format(
      variant_utils.only_call(variant), 'VAF', vcf_reader)
  # Simple sum of alleles for multi-allelic variants
  return sum(vafs)


def _get_variant_stats(variant, vcf_reader=None):
  """Returns a VariantStats object corresponding to the input variant."""
  vtype = _get_variant_type(variant)
  is_transition, is_transversion = _tstv(variant, vtype)
  vaf = None
  if vcf_reader is not None:
    vaf = _get_vaf(variant, vcf_reader)

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
          variant_utils.only_call(variant)),
      genotype=str(
          sorted(variantcall_utils.get_gt(variant_utils.only_call(variant)))),
      vaf=vaf,
      qual=variant.quality)


def _single_variant_stats(variants, vcf_reader=None):
  return [_get_variant_stats(v, vcf_reader=vcf_reader) for v in variants]


def _format_histogram_for_vega(counts, bins):
  """Format histogram counts and bins for vega.

  Args:
    counts: list of bin counts from np.histogram
    bins: list of bins from np.histogram

  Returns:
    A list of objects with bin_start, bin_end, and count for each bin in the
    histogram.
  """
  # Avoid floats becoming 0.6000000000000001 to save space in output json file
  rounded_bins = [round(x, 10) for x in bins]
  # pylint: disable=g-complex-comprehension
  vega_formatted_hist = [{
      'bin_start': rounded_bins[idx],
      'bin_end': rounded_bins[idx + 1],
      'count': count
  } for idx, count in enumerate(counts)]
  # pylint: enable=g-complex-comprehension
  return vega_formatted_hist


def _vaf_histograms_by_genotype(single_stats, number_of_bins=10):
  """Computes histograms of allele frequency for each genotype.

  Args:
    single_stats: list of VariantStats objects.
    number_of_bins: integer, number of bins in allele frequency histogram.

  Returns:
    A dictionary keyed by genotype where each value is a list of bins.
  """

  # Group by genotype
  sorted_by_genotype = sorted(single_stats, key=lambda x: x.genotype)
  grouped_by_genotype = itertools.groupby(sorted_by_genotype,
                                          lambda x: x.genotype)
  stats_by_genotype = {}
  for genotype, group in grouped_by_genotype:
    # Get VAF for each variant where it is defined
    vafs = [x.vaf for x in group if x.vaf is not None]
    counts, bins = np.histogram(vafs, bins=number_of_bins, range=(0, 1))
    stats_by_genotype[genotype] = _format_histogram_for_vega(counts, bins)

  return stats_by_genotype


def _count_base_changes_and_indel_sizes(single_stats):
  """Count each base change, such as A->G or C->T, and count the number of indels of each size.

  Args:
    single_stats: list of VariantStats objects.

  Returns:
    base_changes: {(ref, alt): count, ...}
    indel_sizes: {size: count, ...}
  """

  base_changes = collections.defaultdict(int)
  indel_sizes = collections.defaultdict(int)
  for v in single_stats:
    ref = v.reference_bases
    alts = v.alternate_bases
    # RefCalls are ignored
    if v.is_variant:
      # Multiallelic variants ignored here because they have different indel
      # sizes and/or base changes
      if v.variant_type == BIALLELIC_SNP:
        # SNV: get base change
        base_changes[(ref, alts[0])] += 1
      elif v.variant_type in [BIALLELIC_INSERTION, BIALLELIC_DELETION]:
        # indel: get size
        # + = insertion
        # - = deletion
        size = len(alts[0]) - len(ref)
        indel_sizes[size] += 1

  base_changes_for_json = []
  for key in base_changes:
    ref, alt = key
    base_changes_for_json.append([ref, alt, base_changes[key]])

  indel_sizes_for_json = []
  for key in indel_sizes:
    indel_sizes_for_json.append([int(key), indel_sizes[key]])

  return base_changes_for_json, indel_sizes_for_json


def _compute_qual_histogram(single_var_stats):
  """Compute a histogram over variant quality (QUAL column in VCF).

  Args:
    single_var_stats: list of VariantStats objects.

  Returns:
    histogram of variant quality scores.
  """
  quals = [round(v.qual, 4) for v in single_var_stats]

  counts, bins = np.histogram(quals, bins=100)
  return _format_histogram_for_vega(counts, bins)


def _compute_gq_histogram(single_var_stats):
  """Compute a histogram over genotype quality (GQ sub-column under FORMAT in VCF).

  Args:
    single_var_stats: list of VariantStats objects.

  Returns:
    histogram of genotype quality scores.
  """
  quals = [v.genotype_quality for v in single_var_stats]
  counts, bins = np.histogram(quals, bins=100)
  return _format_histogram_for_vega(counts, bins)


def _compute_summary_stats(single_stats):
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
      gq_stdev=np.std(transposed_dict['genotype_quality']),
      transition_count=sum(transposed_dict['is_transition']),
      transversion_count=sum(transposed_dict['is_transversion']))


def _variants_to_stats_json(variants, vcf_reader=None, histogram_bins=10):
  """Computes variant statistics of each variant.

  Args:
    variants: iterable(Variant)
    vcf_reader: VcfReader
    histogram_bins: integer, number of bins for histogram

  Returns:
    A tuple of (stats_json, summary_json), where state_json is a JSON of single
    variant statistics, and summary_json is a JSON of single sample statistics.
  """

  single_var_stats = _single_variant_stats(variants, vcf_reader=vcf_reader)
  transposed_records = zip(*single_var_stats)
  stats = dict(
      zip(_VARIANT_STATS_COLUMNS, [list(x) for x in transposed_records]))

  summ_stats = _compute_summary_stats(single_var_stats).asdict()

  base_changes, indel_sizes = _count_base_changes_and_indel_sizes(
      single_var_stats)

  histograms = _vaf_histograms_by_genotype(
      single_var_stats, number_of_bins=histogram_bins)

  qual_histogram = _compute_qual_histogram(single_var_stats)
  gq_hist = _compute_gq_histogram(single_var_stats)

  vis_data = {
      'vaf_histograms_by_genotype': histograms,
      'indel_sizes': indel_sizes,
      'base_changes': base_changes,
      'qual_histogram': qual_histogram,
      'gq_histogram': gq_hist
  }

  return stats, summ_stats, vis_data


def _write_json(stats, outfile):
  """Writes stats to the output file."""
  stats_json_string = json.dumps(stats, sort_keys=True, separators=(',', ':'))

  with tf.io.gfile.GFile(outfile, 'w') as writer:
    writer.write(stats_json_string)


def create_vcf_report(variants,
                      output_basename,
                      sample_name,
                      vcf_reader=None,
                      histogram_bins=10,
                      include_individual_variant_stats=True):
  """Calculate VCF stats, save JSON files, and create a visual report."""
  stats, summary_stats, vis_data = _variants_to_stats_json(
      variants, vcf_reader=vcf_reader, histogram_bins=histogram_bins)

  if include_individual_variant_stats:
    _write_json(stats, output_basename + '.per_record.json')
  _write_json(summary_stats, output_basename + '.summary.json')
  _write_json(vis_data, output_basename + '.vis_data.json')
  vcf_stats_vis.create_visual_report(output_basename, summary_stats, vis_data,
                                     sample_name)
