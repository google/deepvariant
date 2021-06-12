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

BIALLELIC_SNP = 'Biallelic_SNP'
BIALLELIC_INSERTION = 'Biallelic_Insertion'
BIALLELIC_DELETION = 'Biallelic_Deletion'
BIALLELIC_MNP = 'Biallelic_MNP'
MULTIALLELIC_SNP = 'Multiallelic_SNP'
MULTIALLELIC_INSERTION = 'Multiallelic_Insertion'
MULTIALLELIC_DELETION = 'Multiallelic_Deletion'
MULTIALLELIC_COMPLEX = 'Multiallelic_Complex'
REFCALL = 'RefCall'


def _get_variant_type(variant):
  """Returns the type of variant as a string."""
  if variant_utils.is_variant_call(variant):
    biallelic = variant_utils.is_biallelic(variant)
    snp = variant_utils.is_snp(variant)
    insertion = variant_utils.variant_is_insertion(variant)
    deletion = variant_utils.variant_is_deletion(variant)

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
  return sum(vafs)


def _get_variant_stats(variant, vaf_available=False, vcf_reader=None):
  """Returns a VariantStats object corresponding to the input variant."""
  vtype = _get_variant_type(variant)
  is_transition, is_transversion = _tstv(variant, vtype)
  vaf = None
  if vaf_available:
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


def _single_variant_stats(variants, vaf_available=False, vcf_reader=None):
  return [
      _get_variant_stats(v, vaf_available=vaf_available, vcf_reader=vcf_reader)
      for v in variants
  ]


def _format_histogram_for_vega(counts, bins):
  """Format histogram counts and bins for vega.

  Args:
    counts: list of bin counts from np.histogram
    bins: list of bins from np.histogram

  Returns:
    A list of objects with s (bin start), e (bin end), and c (bin count) for
    each bin in the histogram.
  """
  # Avoid floats becoming 0.6000000000000001 to save space in output json
  rounded_bins = [round(x, 10) for x in bins]
  # pylint: disable=g-complex-comprehension
  vega_formatted_hist = [{
      's': rounded_bins[idx],
      'e': rounded_bins[idx + 1],
      'c': count
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
    # Create a few placeholder bins
    stats_by_genotype[genotype] = _fraction_histogram([], 2)
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


def _round_down(num):
  return int(math.floor(num))


def _round_up(num):
  return int(math.ceil(num))


def _compute_qual_histogram(single_var_stats):
  """Compute a histogram over variant quality (QUAL column in VCF).

  Args:
    single_var_stats: list of VariantStats objects.

  Returns:
    histogram of variant quality scores.
  """
  quals = [round(v.qual, 4) for v in single_var_stats]

  if quals:
    bin_range = (_round_down(min(quals)), _round_up(max(quals) + 1))
    counts, bins = np.histogram(
        quals, range=bin_range, bins=bin_range[1] - bin_range[0])
    hist = _format_histogram_for_vega(counts, bins)
    return [x for x in hist if x['c'] > 0]
  else:
    return []


def _get_integer_counts(nums):
  """Turn a list of integers into a list of counts of those integers.

  Args:
    nums: a list of numbers (e.g. [1,2,2,4])

  Returns:
    a list of [num, count] (e.g. [[1,1],[2,2],[4,1]]) for all integers with
    non-zero counts
  """
  bin_counts = np.bincount(nums)
  non_zero_counts = [[i, x] for i, x in enumerate(bin_counts) if x > 0]
  return non_zero_counts


def _compute_gq_histogram(single_var_stats):
  """Compute a histogram over genotype quality (GQ sub-column under FORMAT in VCF).

  Args:
    single_var_stats: list of VariantStats objects.

  Returns:
    histogram of genotype quality scores.
  """
  quals = [
      v.genotype_quality
      for v in single_var_stats
      if not isinstance(v.genotype_quality, list)
  ]
  return _get_integer_counts(quals)


def _compute_depth_histogram(single_var_stats):
  """Compute a histogram on the depth, with larger bins as depth increases."""
  depths = [v.depth for v in single_var_stats if not isinstance(v.depth, list)]
  return _get_integer_counts(depths)


def _count_variant_types(single_stats):
  count_all_variant_types = collections.defaultdict(int)
  for v in single_stats:
    count_all_variant_types[v.variant_type] += 1

  return count_all_variant_types


def _count_titv(single_stats):
  titv_counts = {'Transition': 0, 'Transversion': 0}
  titv_counts['Transition'] = sum([v.is_transition for v in single_stats])
  titv_counts['Transversion'] = sum([v.is_transversion for v in single_stats])
  return titv_counts


def _compute_variant_stats_for_charts(variants, vcf_reader=None):
  """Computes variant statistics of each variant.

  Args:
    variants: iterable(Variant).
    vcf_reader: VcfReader.

  Returns:
    A dict with summarized data prepared for charts.
  """
  vaf_available = False
  if vcf_reader:
    vcf_columns = [col.id for col in vcf_reader.header.formats]
    vaf_available = 'VAF' in vcf_columns

  single_var_stats = _single_variant_stats(
      variants, vaf_available=vaf_available, vcf_reader=vcf_reader)

  titv_counts = _count_titv(single_var_stats)
  variant_type_counts = _count_variant_types(single_var_stats)

  base_changes, indel_sizes = _count_base_changes_and_indel_sizes(
      single_var_stats)

  histograms = _vaf_histograms_by_genotype(single_var_stats, number_of_bins=50)

  qual_histogram = _compute_qual_histogram(single_var_stats)
  gq_hist = _compute_gq_histogram(single_var_stats)
  depth_histogram = _compute_depth_histogram(single_var_stats)

  vis_data = {
      'vaf_histograms_by_genotype': histograms,
      'indel_sizes': indel_sizes,
      'base_changes': base_changes,
      'qual_histogram': qual_histogram,
      'gq_histogram': gq_hist,
      'variant_type_counts': variant_type_counts,
      'depth_histogram': depth_histogram,
      'titv_counts': titv_counts
  }

  return vis_data


def create_vcf_report(variants, output_basename, sample_name, vcf_reader=None):
  """Calculate VCF stats and create a visual report."""
  vis_data = _compute_variant_stats_for_charts(variants, vcf_reader=vcf_reader)

  vcf_stats_vis.create_visual_report(output_basename, vis_data, sample_name)
