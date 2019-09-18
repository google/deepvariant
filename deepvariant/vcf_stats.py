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
import math
import numpy as np

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
  # Avoid floats becoming 0.6000000000000001 to save space in output json
  rounded_bins = [round(x, 10) for x in bins]
  # pylint: disable=g-complex-comprehension
  vega_formatted_hist = [{
      'bin_start': rounded_bins[idx],
      'bin_end': rounded_bins[idx + 1],
      'count': count
  } for idx, count in enumerate(counts)]
  # pylint: enable=g-complex-comprehension
  return vega_formatted_hist


def _fraction_histogram(values, number_of_bins=10):
  counts, bins = np.histogram(values, bins=number_of_bins, range=(0, 1))
  return _format_histogram_for_vega(counts, bins)


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
  # Fill in empty placeholders for genotypes to populate all five charts
  stats_by_genotype = {}
  required_genotypes = ['[0, 0]', '[0, 1]', '[1, 1]', '[-1, -1]', '[1, 2]']
  for genotype in required_genotypes:
    stats_by_genotype[genotype] = _fraction_histogram([], number_of_bins)
  # Count vafs from variants (replacing placeholders)
  for genotype, group in grouped_by_genotype:
    # Get VAF for each variant where it is defined
    vafs = [x.vaf for x in group if x.vaf is not None]
    stats_by_genotype[genotype] = _fraction_histogram(vafs, number_of_bins)

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


def _default_limits_histogram(numbers, default_min=0, default_max=100):
  """Create a histogram with default but flexible x-axis limits.

  The purpose is to standardize the x-axis of the plots without the risk of
  hiding some real data in rare cases that don't conform to expectations.

  Args:
    numbers: a list of numbers to serve as the data for the histogram.
    default_min: an integer.
    default_max: an integer.

  Returns:
    histogram data as a list of bins.
  """
  # Set range from e.g. 0 to 100 by default but allow it to expand to include
  # all the numbers if any fall outside of the range. E.g. a range of 0 to 100
  # will expand to 0 to 190 if one of the numbers is 190.

  bin_range = (int(math.floor(min(min(numbers), default_min))),
               int(math.ceil(max(max(numbers), default_max))))
  counts, bins = np.histogram(
      numbers, range=bin_range, bins=bin_range[1] - bin_range[0])
  return _format_histogram_for_vega(counts, bins)


def _compute_qual_histogram(single_var_stats):
  """Compute a histogram over variant quality (QUAL column in VCF).

  Args:
    single_var_stats: list of VariantStats objects.

  Returns:
    histogram of variant quality scores.
  """
  quals = [round(v.qual, 4) for v in single_var_stats]
  # set range from 0 to 150 by default but allow it to expand
  return _default_limits_histogram(quals, default_min=0, default_max=150)


def _compute_gq_histogram(single_var_stats):
  """Compute a histogram over genotype quality (GQ sub-column under FORMAT in VCF).

  Args:
    single_var_stats: list of VariantStats objects.

  Returns:
    histogram of genotype quality scores.
  """
  quals = [v.genotype_quality for v in single_var_stats]
  return _default_limits_histogram(quals, default_min=0, default_max=150)


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


def _count_variant_types(single_stats):
  count_all_variant_types = collections.defaultdict(int)
  for v in single_stats:
    count_all_variant_types[v.variant_type] += 1

  return count_all_variant_types


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

  summ_stats = _compute_summary_stats(single_var_stats).asdict()

  base_changes, indel_sizes = _count_base_changes_and_indel_sizes(
      single_var_stats)

  histograms = _vaf_histograms_by_genotype(
      single_var_stats, number_of_bins=histogram_bins)

  qual_histogram = _compute_qual_histogram(single_var_stats)
  gq_hist = _compute_gq_histogram(single_var_stats)

  variant_type_counts = _count_variant_types(single_var_stats)

  vis_data = {
      'vaf_histograms_by_genotype': histograms,
      'indel_sizes': indel_sizes,
      'base_changes': base_changes,
      'qual_histogram': qual_histogram,
      'gq_histogram': gq_hist,
      'variant_type_counts': variant_type_counts
  }

  return summ_stats, vis_data


def create_vcf_report(variants,
                      output_basename,
                      sample_name,
                      vcf_reader=None,
                      histogram_bins=10):
  """Calculate VCF stats and create a visual report."""
  summary_stats, vis_data = _variants_to_stats_json(
      variants, vcf_reader=vcf_reader, histogram_bins=histogram_bins)

  vcf_stats_vis.create_visual_report(output_basename, summary_stats, vis_data,
                                     sample_name)
