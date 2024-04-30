# Copyright 2017 Google LLC.
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
"""Haplotype-based labeling algorithm for DeepVariant.

This module provides a haplotype-aware labeling algorithm. This is a more
sophisticated approach to labeling that allows for slight representational
differences between candidate and truth variant sets. See:

https://github.com/ga4gh/benchmarking-tools
https://www.biorxiv.org/content/early/2018/03/15/270157

for an introduction to the concepts and why this is important.


The module is implemented in two big pieces of functionality:

find_best_matching_haplotypes(candidates, truths) provides an function that
accepts a list of candidate variants and a list of truth variants with known
genotypes and finds an assignment of genotypes for candidates and truth that
results in the same two haplotype sequences in the region. Since the truth
variants have known genotypes, the search there is constrained to those
genotypes and their potential set of false negatives (e.g., if truth is (0, 1)
we may have missed the variant so we consider both (0, 1) and (0, 0)). The
returned value is a HaplotypeMatch object describing the genotype assignments
for candidates and truth.

HaplotypeLabeler implements the variant_labeler.VariantLabeler API by calling
our find_best_matching_haplotypes function to get the HaplotypeMatch objects and
returning variant_labeler.VariantLabel objects for each candidate variant.
"""

import collections
import copy
import enum
import heapq
import itertools
import operator

from absl import logging

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io import fasta
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

VariantAndGenotypes = collections.namedtuple(
    'VariantAndGenotype', ['variant', 'genotypes']
)

# The default maximum size of a variant group we'll try to label. See
# the HaplotypeLabeler class for more information.
_MAX_GROUP_SIZE = 8

# The default maximum distance between subsequent variants within a group. See
# the HaplotypeLabeler class for more information.
_MAX_SEPARATION_WITHIN_VARIANT_GROUP = 30

# The default maximum product of possible genotypes combinations.
_MAX_GT_OPTIONS_PRODUCT = 100000

# The variants that are within this value will be forcefully grouped together.
# This is to ensure that we don't decouple a candidate with its truth variant
# in variant dense regions.
# The value of this is set to 0bp because 1bp creates a lot of overhead.
# Context about the value: internal#comment3
_FORCE_GROUP_WITHIN_BP = 0

# True we will generate enough information into our logs to help debug bad
# regions.
_DEBUG_PRINTING_IS_ENABLED = False


class HaplotypeLabeler(variant_labeler.VariantLabeler):
  """Haplotype-based variant labeler."""

  def __init__(
      self,
      truth_vcf_reader,
      ref_reader,
      confident_regions,
      max_group_size=_MAX_GROUP_SIZE,
      max_separation=_MAX_SEPARATION_WITHIN_VARIANT_GROUP,
      max_gt_options_product=_MAX_GT_OPTIONS_PRODUCT,
  ):
    """Creates a new HaplotypeVariantLabeler.

    Args:
      truth_vcf_reader: a VcfReader object that points to our truth variant set.
      ref_reader: A FastaReader object we can use to get reference bases.
      confident_regions: A RangeSet containing all of the confidently called
        regions. A variant that falls outside of one of these regions will be
        receive a special not-confident marker.
      max_group_size: int >= 1. The maximum number of variants we'll attempt to
        label together. Larger values increase the runtime of the algorithm.
      max_separation: int >= 0. The maximum distance between variants within a
        group. Sequential variants separated by more than this value will be
        placed in separate groups for labeling.
      max_gt_options_product: int >= 0. The maximum number of combinations of
        genotypes (product of all genotypes in the group).

    Raises:
      ValueError: if vcf_reader is None.
    """
    super().__init__(
        truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions
    )
    if confident_regions is None:
      raise ValueError('confident_regions cannot be None for HaplotypeLabeler.')
    self._ref_reader = ref_reader
    self.max_group_size = max_group_size
    self.max_separation = max_separation
    self.max_gt_options_product = max_gt_options_product
    self._metrics = deepvariant_pb2.LabelingMetrics()

  def label_variants(self, variants, region):
    # Grab our truth variants and group up variants + truth into small enough
    # chunks that we can safely send them into our find_best_matching_haplotypes
    # function.
    truths = list(self._get_truth_variants(region))
    if truths:
      # Filter out homozygous reference labels.
      truths = [
          y
          for x, y in zip(
              map(lambda x: sum(x) > 0, _variant_genotypes(truths)), truths
          )
          if x
      ]
    grouped = group_variants(
        candidates=list(variants),
        truths=truths,
        max_group_size=self.max_group_size,
        max_separation=self.max_separation,
        max_gt_options_product=self.max_gt_options_product,
    )

    # Now loop over our grouped variants, labeling them, and yielding
    # VariantLabel objects.
    for candidates_group, truth_group in grouped:
      ref = self.make_labeler_ref(candidates_group, truth_group)
      labeling = find_best_matching_haplotypes(
          candidates_group, truth_group, ref
      )
      if labeling is None:
        # Note this test must be 'is None' since label_variants can return an
        # empty list.
        raise ValueError(
            'Failed to assign labels for variants',
            candidates_group,
            truth_group,
            ref,
        )

      self._update_metrics(labeling)
      for labeled in labeling.candidates_with_assigned_genotypes():
        # This logic doesn't make a huge amount of sense when you are doing
        # haplotype-based labeling. Currently we only say a variant is confident
        # if it overlaps the confident regions, which is the baseline behavior.
        # However, it may be useful to rethink how we establish a variant is
        # confident, as the "event" may be within the confident regions but
        # shifted outside due to differences in representational choices. Seems
        # like another approach would be to assign confidence if it has a
        # non-ref genotype (as we only consider confident truth variants) or if
        # it overlaps the confident regions.
        yield variant_labeler.VariantLabel(
            is_confident=self._confident_regions.variant_overlaps(labeled),
            genotype=tuple(labeled.calls[0].genotype),
            variant=labeled,
        )

  @property
  def metrics(self):
    """Gets the LabelingMetrics proto tracking metrics for this labeler."""
    return self._metrics

  def _update_metrics(self, labeling):
    """Update self._metrics with the HaplotypeMatch labeling results.

    This function updates the LabelingMetrics information in self._metrics using
    the labeling results in labeling.

    Args:
      labeling: HaplotypeMatch. The labeling information to use to update our
        LabelingMetrics.
    """

    def _n_alts_by_genotype(gt):
      """Returns the number of distinct alt alleles with non-zero genotype."""
      return len({g for g in gt if g > 0})

    def _is_hom_ref(gt):
      """Are all genotypes in gt the reference alleles (i.e., == 0)?"""
      return all(g == 0 for g in gt)

    def _has_alt_genotypes(gt):
      """Is any genotype in gt a non-ref (> 0) genotype?"""
      return any(g > 0 for g in gt)

    # Iterate over the truth variant and its associated original genotypes
    # (those provided by the input VCF) and the assigned genotypes (i.e., the
    # genotypes assigned to truth to make candidates and truth match haplotypes)
    # and compute a few metric values.
    for truth, original_gt, assigned_gt in zip(
        labeling.truths,
        labeling.original_truth_genotypes,
        labeling.truth_genotypes,
    ):
      n_alts_original = _n_alts_by_genotype(original_gt)

      self._metrics.n_truth_variant_sites += 1
      self._metrics.n_truth_variant_alleles += n_alts_original
      self._metrics.n_true_positive_sites += _has_alt_genotypes(assigned_gt)
      self._metrics.n_false_negative_sites += _is_hom_ref(assigned_gt)

      # If we have more than one alt allele in the original genotypes and the
      # assigned genotypes imply more or more are missing then we've got a
      # multi-allelic truth variant with some missing alleles.
      if (
          n_alts_original > 1
          and _n_alts_by_genotype(assigned_gt) < n_alts_original
      ):
        self._metrics.n_truth_multiallelics_sites_with_missed_alleles += 1

      # Iterate over the original and assigned genotypes for the truth variants
      # and count up the number of true positive alleles (i.e. original and
      # assigned genotypes are non-ref) and false negative alleles (i.e.,
      # original is non-ref but assigned is ref).
      for og, ag in zip(original_gt, assigned_gt):
        if og > 0:
          if ag > 0:
            self._metrics.n_true_positive_alleles += 1
          elif ag == 0:
            self._metrics.n_false_negative_alleles += 1

    # Create a dict from the start of truth to the truth variant itself and its
    # assigned genotypes. This is needed to compute site match counts below.
    truth_by_pos = {
        truth.start: (truth, gt)
        for truth, gt in zip(labeling.truths, labeling.truth_genotypes)
    }

    # Iterate over the candidates and their assigned genotypes to compute
    # the remaining metrics.
    #
    # Note that this counts all candidates, not just the ones in the confident
    # regions of the genome. This seems like a reasonable first approach, but it
    # may be necessary to restrict ourselves to only those overlapping the
    # confident regions.
    for candidate, genotype in zip(
        labeling.candidates, labeling.candidate_genotypes
    ):
      # If candidate isn't confident, add it to the non_confident count and
      # continue as the other metrics are only computed over confident
      # candidates.
      if not self._confident_regions.variant_overlaps(candidate):
        self._metrics.n_non_confident_candidate_variant_sites += 1
        continue

      n_alt_alleles = len(candidate.alternate_bases)
      self._metrics.n_candidate_variant_sites += 1
      self._metrics.n_candidate_variant_alleles += n_alt_alleles
      self._metrics.n_false_positive_sites += _is_hom_ref(genotype)
      self._metrics.n_false_positive_alleles += (
          n_alt_alleles - _n_alts_by_genotype(genotype)
      )

      # Use the truth_by_pos dict to determine which candidates occur at the
      # same position as a truth variant. If there is one, grab it and its
      # genotypes so we can compute metrics on exact position, allele, genotype
      # matches. If not, update the number of inexact matches if our candidate
      # is non-reference itself.
      truth, assigned_gt = truth_by_pos.get(candidate.start, (None, None))
      if truth:
        self._metrics.n_exact_position_matches += 1
        if sorted(candidate.alternate_bases) == sorted(truth.alternate_bases):
          self._metrics.n_exact_position_and_allele_matches += 1
          if sorted(genotype) == sorted(assigned_gt):
            self._metrics.n_exact_position_and_allele_and_genotype_matches += 1
      elif _has_alt_genotypes(genotype):
        self._metrics.n_inexact_position_matches += 1

  def make_labeler_ref(self, candidates, true_variants, bufsize=20):
    all_variants = candidates + true_variants
    contig = all_variants[0].reference_name
    start = min(x.start for x in all_variants)
    end = max(x.end for x in all_variants)
    contig_nbp = self._ref_reader.contig(contig).n_bases
    region = ranges.make_range(
        contig, max(start - 1, 0), min(end + bufsize, contig_nbp)
    )
    ref_bases = self._ref_reader.query(region)
    return ReferenceRegion(ref_bases, start=region.start)


class ReferenceRegion(fasta.InMemoryFastaReader):
  """Allows us to get bases from a cached reference interval."""

  # We don't want to worry about the chromosome we are working on for code
  # clarity, so we create an InMemoryFastaReader that has a single chromosome
  # named _PLACEHOLDER_CHROM_NAME which allows us to provide a bases(start, end)
  # function for convenient reading of bases.
  _PLACEHOLDER_CHROM_NAME = '*'

  def __init__(self, bases, start):
    super(ReferenceRegion, self).__init__(
        [(self._PLACEHOLDER_CHROM_NAME, start, bases)]
    )
    self.start = start
    self.end = start + len(bases)

  def bases(self, start, end):
    return self.query(
        ranges.make_range(self._PLACEHOLDER_CHROM_NAME, start, end)
    )


_CANDIDATE_MARKER = 'candidate'
_TRUTH_MARKER = 'truth'
_VariantToGroup = collections.namedtuple(
    '_VariantToGroup', ['start', 'type', 'variant']
)


def _raise_if_not_sorted_or_not_on_same_chromosome(variants):
  """Raises a ValueError if variants isn't sorted on the same chromosome."""
  if not variant_utils.variants_are_sorted(variants):
    raise ValueError('Variants must be sorted', variants)
  for v in variants[1:]:
    if variants[0].reference_name != v.reference_name:
      raise ValueError(
          'Variants (v1={}, v2={}) not on the same chromosome'.format(
              v.reference_name, variants[0].reference_name
          )
      )


def group_variants(
    candidates,
    truths,
    max_group_size=_MAX_GROUP_SIZE,
    max_separation=_MAX_SEPARATION_WITHIN_VARIANT_GROUP,
    max_gt_options_product=_MAX_GT_OPTIONS_PRODUCT,
    force_group_within_bp=_FORCE_GROUP_WITHIN_BP,
):
  """Splits candidate and truth variants into smaller groups if necessary.

  This function takes in a list of candidate and truth variants and splits up
  those lists into groups that respect the requirements of the max_group_size
  and max_separation arguments. This is necessary because the labeling algorithm
  is very expensive as a function of the number of input variants, so to avoid
  excessive runtime we break up our potentially large list of candidate and
  truth variants into smaller groups (max number controlled by max_group_size)
  based on a maximum distance allowed between the closest variants within the
  group.

  The current algorithm is a simple greedy one; we effectively merge the two
  variant lists together, make groups greedily on that list until either the
  maximum number of elements of a specific type (i.e., max_group_size of 2
  implies we can have up to two candidate variants or truth variants within a
  group) or we encounter a variant further away from the closest variant within
  the current group than allowed by max_separation.

  Args:
    candidates: list[nucleus.proto.Variant]. A sorted list of candidate variants
      on the same chromosome.
    truths: list[nucleus.proto.Variant]. A sorted list of truth variants on the
      same chromosome.
    max_group_size: int >= 0. The maximum number of variants of a specific type
      allowed within a group.
    max_separation: int >= 0. The maximum distance, in basepairs, allowed
      between the closest variants within a group.
    max_gt_options_product: int >= 0. The maximum number of combinations of
      genotypes (product of all genotypes in the group).
    force_group_within_bp: int >= 0. Variants within this many bps will be
      forced to be put in the same group. This is to ensure that we do not
      decouple candidates and truths in variant-dense regions. This value can be
      set to -1 for unit-test purposes. Setting -1 will not force any grouping
      of variants.

  Returns:
    A list of grouped variants in 2-tuples, such as:

      [(candidates1, truth_variants1), ...]

    where each tuple contains the candidate and truth variants for that group.

  Raises:
    ValueError: if any of the inputs are malformed.
  """
  if max_group_size < 0:
    raise ValueError('max_group_size={} must be >= 0'.format(max_group_size))
  if max_separation < 0:
    raise ValueError('max_separation={} must be >= 0'.format(max_separation))
  if max_gt_options_product < 0:
    raise ValueError(
        'max_gt_options_product={} must be >= 0'.format(max_gt_options_product)
    )
  _raise_if_not_sorted_or_not_on_same_chromosome(candidates)
  _raise_if_not_sorted_or_not_on_same_chromosome(truths)

  def _num_genotypes(variant):
    """For example, if there is 2 alts, there are 6 (==4*3/2) possible GTs."""
    num_ref_and_alts = len(variant.alternate_bases) + 1
    return (num_ref_and_alts + 1) * num_ref_and_alts / 2

  def to_grouped_variants(variants, candidate_type):
    """Converts a Variant proto to a _VariantToGroup tuple."""
    return [_VariantToGroup(v.start, candidate_type, v) for v in variants]

  def _of_type(group, required_type):
    """Selects a list of Variant protos from list[_VariantToGroup] of type."""
    return [gv.variant for gv in group if gv.type == required_type]

  def _split_grouped_variants(group):
    """Splits a list of _VariantToGroup into candidate and truth variants."""
    return _of_type(group, _CANDIDATE_MARKER), _of_type(group, _TRUTH_MARKER)

  def _include_in_variant_group(group, group_variant, new_gt_options_product):
    if not group:
      return True
    if new_gt_options_product >= max_gt_options_product:
      logging.info(
          (
              'Not including more because genotype_options_product will '
              'be %s, which exceeds max(=%s)'
          ),
          new_gt_options_product,
          max_gt_options_product,
      )
      return False
    n_of_type = sum(1 for g in group if g.type == group_variant.type)
    if n_of_type >= max_group_size:
      return False
    else:
      return any(
          group_variant.variant.start - g.variant.end + 1 <= max_separation
          for g in group
      )

  def _include_group_by_end_in_variant_group(
      group, group_by_end, new_gt_options_product
  ):
    for variant in group_by_end:
      if not _include_in_variant_group(group, variant, new_gt_options_product):
        return False
    return True

  def _regroup_by_end(current_group, force_group_within_bp):
    """Regroups variants by end positions, if force_group_within_bp >= 0."""
    # If force_group_within_bp is negative, we don't need to group variants
    # by .end. Each variant should be in its own group.
    if force_group_within_bp < 0:
      return [[v] for v in current_group]

    # Group the elements based on the 'variant.end' attribute.
    # Here, I didn't sort by variant.end. The list should be sorted by
    # variant.start already, and I'm using that order.
    new_groups = []
    for _, group_items in itertools.groupby(
        current_group, key=operator.attrgetter('variant.end')
    ):
      new_groups.append(list(group_items))
    return new_groups

  # Convert our lists of variant protos into _VariantToGroup tuples compatible
  # with the heapq API (sorts on tuples), so we get a single iterable of
  # variants with marked types and sorted by start position (first element of
  # each tuple).
  groupable_variants = heapq.merge(
      to_grouped_variants(candidates, _CANDIDATE_MARKER),
      to_grouped_variants(truths, _TRUTH_MARKER),
  )

  # Go through groupable_variants:
  # First, group by variant.end if force_group_within_bp >= 0. Because any
  # variant that has the same end has to be in the same group anyway.
  # Then, consider each group, and merge them into groups according to
  # the predicate _include_group_by_end_in_variant_group.
  groups = []
  current_group = []
  current_gt_options_product = 1
  previous_pos_end = 0
  variants_groups_by_end = _regroup_by_end(
      groupable_variants, force_group_within_bp
  )
  for group_by_end in variants_groups_by_end:
    new_gt_options_product = current_gt_options_product
    for group_variant in group_by_end:
      new_gt_options_product *= _num_genotypes(group_variant.variant)
    distance_from_prev_variant = group_by_end[0].variant.end - previous_pos_end
    if (
        _include_group_by_end_in_variant_group(
            current_group, group_by_end, new_gt_options_product
        )
        or distance_from_prev_variant <= force_group_within_bp
    ):
      current_group.extend(group_by_end)
      current_gt_options_product = new_gt_options_product
    else:
      groups.append(current_group)
      current_group = group_by_end
      current_gt_options_product = 1
      for v in group_by_end:
        current_gt_options_product *= _num_genotypes(v.variant)
    previous_pos_end = group_by_end[0].variant.end
  if current_group:
    groups.append(current_group)

  # Finally split up each group into candidates and truths.
  return [_split_grouped_variants(g) for g in groups]


def with_false_negative_genotypes(gt):
  """Returns a set of genotypes that includes false negatives.

  This function takes a concrete genotype for a Variant, such as (0, 1), and
  returns a set of genotypes that includes gt as well as all possible genotypes
  consistent with some of the alleles in gt being missed. For example, here are
  a few outputs for help understand what this means:

    input genotype (gt) => returned set of genotypes
    ------------------------------------------------

    # A hom-ref genotype doesn't have any alleles to miss.
    (0, 0)  => {(0, 0)}

    # Might miss the 1 allele, or not.
    (0, 1)  => {(0, 0), (0, 1)}

    # We could miss one, or both of the 1 alleles.
    (1, 1)  => {(0, 0), (0, 1), (1, 1)}

    # Multi-allelics are more complex, in that we could miss either the 1 or the
    # 2 allele, or both.
    (1, 2)  => {(0, 0), (0, 1), (0, 2), (1, 2)}

  Args:
    gt: iterable[int]: A genotype for a Variant, such as [0, 1] or [1, 1].

  Returns:
    A sorted list of tuples containing diploid genotypes.
  """
  alts = set(gt) - {0}
  return sorted(list({(0, 0), tuple(gt)} | {(0, alt) for alt in alts}))


class ImpossibleHaplotype(Exception):
  """Indicates that an impossible haplotype configuration has been observed."""

  pass


def enumerate_all_possible_haplotypes(variants, ref, enumeration_type):
  """Returns all possible haplotype/genotype combinations for variants.

  Args:
    variants: list[nucleus.protos.Variant]. A list of candidate variants, in
      coordinate-sorted order, all on the same chromosome.
    ref: ReferenceRegion. Used to get reference bases for variants. Must cover
      at least the span of the variants.
    enumeration_type: EnumerationType enum value. What kind of enumeration do we
      want to do? Can be either CANDIDATES or TRUTH.

  Returns:
    Dict[Haplotypes, List[Genotypes]]
      where
      Genotypes = List[Genotype]
      Haplotypes = FrozenSet[str]

    Haplotypes is a set of either one or two strings, where each string is a
    haplotype (i.e., a series of bases) generated by the genotypes assigned to
    each variant. genotypes is a list of genotype tuples, in the same order as
    variants, indicating the genotype assignment for each variant. These
    genotypes are phased, so [(0, 1), (0, 1)] is not the same as
    [(0, 1), (1, 0)].
  """

  def create_haplotypes_recursive(variants_and_genotypes, last_pos):
    """Recursive driver to enumerate all haplotypes."""
    if not variants_and_genotypes:
      yield {ref.bases(last_pos, ref.end)} if last_pos != ref.end else {''}
    else:
      group, remaining = split_independent_variants(variants_and_genotypes)
      group_haplotypes, next_pos = phased_genotypes_to_haplotypes(
          group, last_pos, ref
      )
      prefix_haplotypes = list(all_diploid_haplotypes(group, group_haplotypes))

      if not prefix_haplotypes:
        # prefix_haplotypes can be empty when group contains incompatible
        # variants making it impossible to construct any haplotypes for group.
        # For example, if group is:
        #   variant(start=6, alleles=("AAT", "A"), genotype=(0, 1))
        #   variant(start=7, alleles=("AT", "T"), genotype=(1, 1))
        # prefix_haplotypes will be empty because there's no way to construct
        # the haplotype where the variant@6 has a 1 genotype since it is
        # deleting away bases that overlap the variant@7 which has a genotype of
        # (1, 1), meaning it *has* to be present in some haplotype. In this
        # situation we raise a ImpossibleHaplotype exception, which is caught
        # in the outer loop, allowing us to bail out of the search ASAP.
        raise ImpossibleHaplotype

      for haplotypes in create_haplotypes_recursive(remaining, next_pos):
        for result in extend_haplotypes(prefix_haplotypes, haplotypes):
          yield result

  def create_haplotypes(variants_and_genotypes, last_pos):
    try:
      for r in create_haplotypes_recursive(variants_and_genotypes, last_pos):
        yield r
    except ImpossibleHaplotype:
      # See comment in create_haplotypes_recursive for more information, but in
      # this case we simply `pass`, as we cannot construct any valid haplotypes.
      pass

  genotype_options = genotype_options_for_variants(variants, enumeration_type)
  haplotypes_to_genotypes_dict = collections.OrderedDict()
  for genotypes in itertools.product(*genotype_options):
    paired = [VariantAndGenotypes(v, g) for v, g in zip(variants, genotypes)]
    for haplotypes in create_haplotypes(paired, ref.start):
      key = frozenset(haplotypes)
      if key not in haplotypes_to_genotypes_dict:
        haplotypes_to_genotypes_dict[key] = []
      haplotypes_to_genotypes_dict[key].append(genotypes)
  return haplotypes_to_genotypes_dict


def all_diploid_haplotypes(variants_and_genotypes, genotypes2haplotype):
  """Returns all diploid haplotypes for variants given their genotypes."""

  def complement_haploid_genotype(haploid_genotype, genotypes):
    assert len(haploid_genotype) == len(genotypes)
    return tuple(
        g1[1] if hg1 == g1[0] and len(g1) == 2 else g1[0]
        for hg1, g1 in zip(haploid_genotype, genotypes)
    )

  genotypes = [vg.genotypes for vg in variants_and_genotypes]
  generated_already = set()
  for haploid_genotype, haplotype in genotypes2haplotype.items():
    complement = complement_haploid_genotype(haploid_genotype, genotypes)
    complement_haplotype = genotypes2haplotype.get(complement, None)
    if complement_haplotype is not None and complement not in generated_already:
      generated_already.add(haploid_genotype)
      yield {haplotype, complement_haplotype}


class EnumerationType(enum.Enum):
  """Enumeration type indicating how we should explore genotype configurations.

  See genotype_options_for_variants for more information.
  """

  # This enumeration produces all possible genotype combinations for a variant.
  CANDIDATES = 1
  # This enumeration type will produce all combinations of the provided genotype
  # for the variant with genotypes that allow for one or more of the
  # non-reference genotypes to be missed.
  TRUTH = 2
  # This enumeration type will only produce a single (0, 0) genotype option for
  # each variant.
  ONLY_HOM_REF = 3


def genotype_options_for_variants(variants, enumeration_type):
  """Returns a list of sets of possible genotypes for each variant in variants.

  This function takes a list of variants and enumeration_type and produces a
  list of possible genotypes for each variant in order.

  If enumeration_type is ONLY_HOM_REF, then we return a singleton set for each
  variant containing only the hom-ref genotype (0, 0). If enumeration_type is
  TRUTH, then each variant must have an associated genotype field values, say
  (A, B), and we return the set genotype as well as all possible false negative
  genotypes. In our example, this means we'd return {(A, B), (0, A), (0, B),
  (0, 0)} as we could miss either the A, the B, or both alleles. If the
  enumeration_type is CANDIDATES, we don't require the Variant protos to have
  existing genotype field values and instead enumerate all possible unphased
  genotypes for each variant given its alternative alleles of each variant. For
  example, if we have a Variant with alleles = 'A' and 'C', we would return the
  three possible diploid genotypes {(0, 0), (0, 1), (1, 1)}.

  Args:
    variants: List[nucleus.protos.Variant]. A list of Variant protos to provide
      genotype options for. Some enumeration types may require the protos to
      have existing genotypes in their calls[] subfield.
    enumeration_type: EnumerationType. The kind of genotypes we want to explore
      for each variant.

  Returns:
    A list of sets with the same length and "order" as variants. Each set
    contains one or more diploid genotype tuples [e.g., (0, 1)] that
    collectively represent the possible genotypes we need to explore.

  Raises:
    ValueError: if enumeration_type isn't one of the valid options.
  """
  if enumeration_type == EnumerationType.TRUTH:
    return [
        with_false_negative_genotypes(x) for x in _variant_genotypes(variants)
    ]
  elif enumeration_type == EnumerationType.CANDIDATES:
    return [
        {
            (i, j)
            for i, j, _, _ in variant_utils.genotype_ordering_in_likelihoods(v)
        }
        for v in variants
    ]
  elif enumeration_type == EnumerationType.ONLY_HOM_REF:
    return [{(0, 0)}] * len(variants)
  else:
    raise ValueError('Unexpected EnumerationType', enumeration_type)


def split_independent_variants(variants_and_genotypes):
  """Splits variants_and_genotypes into an overlapping group and remaining."""
  if not variants_and_genotypes:
    raise ValueError('Expected at least one value in variants_and_genotypes')

  overlaps = [variants_and_genotypes[0]]
  for i in range(1, len(variants_and_genotypes)):
    vgi = variants_and_genotypes[i].variant
    if any(variant_utils.variants_overlap(vg.variant, vgi) for vg in overlaps):
      overlaps.append(variants_and_genotypes[i])
    else:
      return overlaps, variants_and_genotypes[i:]
  return overlaps, []


def extend_haplotypes(prefix_haplotypes_list, haplotypes):
  """Yields all diploid combinations of prefix_haplotypes_list x haplotypes.

  Args:
    prefix_haplotypes_list: list[set[string]]: prefix_haplotypes_list contains a
      list of set[string], which are just like haplotypes (i.e., contains 1 or 2
      strings), that collectively represent all possible prefixes of haplotypes.
    haplotypes: set[string]. A set containing 1 or 2 haplotype strings. So it
      looks like {h} or {h1, h2}.

  Yields:
    A series of set[string], each containing 1 or 2 haplotype strings.

  Raises:
    ValueError: if any of the arguments are invalid.
  """
  if not prefix_haplotypes_list:
    raise ValueError('prefix_haplotypes_list cannot be empty')
  if len(haplotypes) not in {1, 2}:
    raise ValueError('haplotypes must have exactly 1 or 2 elements', haplotypes)

  for prefix_haplotypes in prefix_haplotypes_list:
    if len(prefix_haplotypes) == 1:
      (f,) = prefix_haplotypes
      yield {f + h for h in haplotypes}
    else:
      f1, f2 = prefix_haplotypes
      if len(haplotypes) == 1:
        (h,) = haplotypes
        yield {f1 + h, f2 + h}
      else:
        h1, h2 = haplotypes
        yield {f1 + h1, f2 + h2}
        yield {f1 + h2, f2 + h1}


def phased_genotypes_to_haplotypes(variants_and_genotypes, start, ref):
  """Returns a map from phased genotypes => haplotype sequences.

  This function creates a map from all possible haploid genotypes of the
  genotypes in variants_and_genotypes to their corresponding haplotype sequences
  implied by the variants, ref, start, and their genotypes. This map can be used
  to efficiently look up the haplotype sequence for any haploid genotype.

  Args:
    variants_and_genotypes: list[VariantAndGenotypes]. The variants and
      associated genotypes to use to build the dictionary.
    start: int >= 0. The position on the genome to start constructing our
      haplotypes at.
    ref: ReferenceRegion. Object containing the reference genome bases we use to
      construct our haplotypes.

  Returns:
    A 2-tuple. The first element is a dictionary[tuple, string], where each key
    is a phased haploid genotype and its value is the haplotype sequence implied
    by that genotype given the variants and the reference genome. The second
    position is the ending position of the haplotype on the reference genome.
  """
  genotypes_to_haplotypes = {}
  genotypes = [vg.genotypes for vg in variants_and_genotypes]
  variants = [vg.variant for vg in variants_and_genotypes]
  all_haploid_genotypes = sorted(set(itertools.product(*genotypes)))
  end = max(v.end for v in variants)
  for phased in all_haploid_genotypes:
    haplotype = build_haplotype(variants, phased, ref, start, end)
    if haplotype:
      genotypes_to_haplotypes[phased] = haplotype
  return genotypes_to_haplotypes, end


def build_haplotype(variants, allele_indices, ref, ref_start, ref_end):
  """Builds the haplotype string from variants and its phased gneotypes.

  This function takes a list of variants and associated phased genotypes and
  constructs the haplotype sequence implied by variants and its genotypes. For
  example, suppose we have two variants:

    ref: CAGC where the first base (C) is at position 10.
    chr20:10 A=>C
    chr20:11 G=>T
    allele_indices: [0, 1]

  We would receive arguments here of a list of two variants and a list of the
  allele_indices [0, 1]. We now look up which base is implied for the variant
  (e.g., [0, 1] takes the reference bases from variant 1 at 10 and then the
  first alternate allele of variant 2 at 11). If ref_start is 9, we would then
  construct the haplotype as:

  haplotype is 'CATC' derived as follows:
    'C' [ref_prefix, since ref_start=]
    +
    'A' (variant1 has a genotype of 0)
    +
    'T' (variant2 has a genotype of 1)
    +
    'C' [ref_postfix if ref_end == 13]
    = 'CATC'

  Args:
    variants: list[nucleus.protos.Variant]: The variants to use in constructing
      the haplotype.
    allele_indices: list[int]: The list of allele_indexes (where 0 means
      reference_bases and > 0 implies alternative bases, following the
      VariantCall genotypes semantics).
    ref: ReferenceRegion. Used to get the reference bases.
    ref_start: int >= 0. The first position (zero-indexed, inclusive) in the
      genome where we want to start constructing our haplotype.
    ref_end: int >= 0 and > ref_start. The last position (zero-indexed,
      exclusive) in the genome where we want to end constructing our haplotype.

  Returns:
    A string containing the haplotype bases, or None, if any variant starts
    before ref_start and has a non-reference genotype.

  Raises:
    ValueError: If any of the input arguments are malformed or otherwise violate
    the assumptions of this algorithm.
  """
  if len(variants) != len(allele_indices):
    raise ValueError(
        'Expected the same number of variants {} as allele_indices {}'.format(
            len(variants), len(allele_indices)
        )
    )
  if ref_start < 0 or ref_start >= ref_end:
    raise ValueError(
        'expected ref_start {} < ref_end {}'.format(ref_start, ref_end)
    )

  parts = []
  position = ref_start
  for variant, allele_index in zip(variants, allele_indices):
    if variant.start < position:
      if allele_index != 0:
        return None
    else:
      ref_prefix = ref.bases(position, variant.start)
      allele = _allele_from_index(variant, allele_index)
      if allele_index == 0:
        # Update our position variable to be the next reference base we want to
        # use when further constructing our haplotype string. If we are using
        # the reference base, we start our position at the base after variant
        # start, whereas if we are using a non-reference base we use the
        # variant.end.
        #
        # This special-case is needed to handle deletion alleles properly. If we
        # have a deletion (e.g., AA => A with start = 10 and end = 12) then we
        # only want to skip to position 12 for the next reference bases if we
        # have have the deletion, otherwise we'd miss the second 'A' base which
        # is really there (the variant isn't present, after all). Another
        # consequence of this choice we only want to add the first base of the
        # reference allele, not the whole string, since this would append all of
        # deletion bases inappropriately to our haplotype.
        allele = allele[0]
        position = variant.start + 1
      else:
        position = variant.end
      parts.append(ref_prefix + allele)

  # We have some bases left to add between the position of our last variant
  # and the ref_end, so append those now.
  if position < ref_end:
    parts.append(ref.bases(position, ref_end))

  return ''.join(parts)


class HaplotypeMatch:
  """DataClass holding information about a matching of variants.

  The haplotype labeling algorithm, at its core, searches for an assignment of
  genotypes to the candidate variants and the truth variants that result in the
  same diploid haplotype sequences, which we call a match. All of the
  information in that previous sentence is captured here as class attributes:

  Attributes:
    haplotypes: list[str]. The sorted list of haplotypes produced by this match.
    candidates: list[nucleus.proto.Variant]: The list of candidate variants.
    truths: list[nucleus.proto.Variant]: The list of true variants.
    candidate_genotypes: list[tuple]: The genotypes that, when assigned to the
      candidate variants, give rise to haplotypes.
    truth_genotypes: list[tuple]: The genotypes that, when assigned to the known
      variants, give rise to haplotypes.
  """

  def __init__(
      self, haplotypes, candidates, candidate_genotypes, truths, truth_genotypes
  ):
    if len(haplotypes) not in {1, 2}:
      raise ValueError('Expected 1 or 2 haplotypes but got', haplotypes)
    if len(candidates) != len(candidate_genotypes):
      raise ValueError(
          'candidates and candidate_genotypes should have the same length'
      )
    if len(truths) != len(truth_genotypes):
      raise ValueError('truths and truth_genotypes should have the same length')

    self.haplotypes = sorted(haplotypes)
    self.candidates = candidates
    self.truths = truths
    self.candidate_genotypes = candidate_genotypes
    self.truth_genotypes = truth_genotypes

    # Computed on-demand.
    self._n_false_positives = None
    self._n_false_negatives = None

  def __str__(self):
    return (
        'HaplotypeMatch(haplotypes={}, false_negatives={}, '
        'false_positives={} true_positives={} match_metrics={}, '
        'variant_gts={}, true_gts={})'
    ).format(
        self.haplotypes,
        self.n_false_negatives,
        self.n_false_positives,
        self.n_true_positives,
        self.match_metrics,
        self.candidate_genotypes,
        self.truth_genotypes,
    )

  __repr__ = __str__

  @property
  def original_truth_genotypes(self):
    return _variant_genotypes(self.truths)

  @property
  def match_metrics(self):
    """Quality of this match.

    Lower scores are better.

    Returns:
      tuple[int] where all elements are >= 0: The tuple is suitable for sorting
      matches, so that sorted(matches, key=lambda x: x.match_metrics) will rank
      matches so that the best option is first.
    """
    return (
        self.n_false_negatives,
        self.n_false_positives,
        self.n_true_positives,
    )

  @property
  def n_true_positives(self):
    """Gets the number of candidates whose matched genotype is not (0, 0).

    Since the candidates don't have expected genotypes, we can only count each
    site instead of each genotype. So this is the number of candidates whose
    matched genotype is not (0, 0).

    Returns:
      int >= 0.
    """
    return len(self.candidate_genotypes) - self.n_false_positives

  @property
  def n_false_positives(self):
    """Gets the number of candidates whose matched genotype is (0, 0).

    Since the candidates don't have expected genotypes, we can only count each
    site instead of each genotype. So this is the number of candidates whose
    matched genotype is (0, 0).

    Returns:
      int >= 0.
    """
    if self._n_false_positives is None:
      self._n_false_positives = sum(
          sum(gt) == 0 for gt in self.candidate_genotypes
      )
    return self._n_false_positives

  @property
  def n_false_negatives(self):
    """Gets the number of missed true genotypes.

    This is the sum of missed non-ref genotypes over all truth variants. So if
    we have a matched truth genotype of (0, 1) and the true genotype is (1, 1),
    then we have 1 FN. If the matched genotype were (0, 0), we'd have 2 FNs.

    Returns:
      int >= 0.
    """
    if self._n_false_negatives is None:
      self._n_false_negatives = sum(
          n_zeroes(assigned_gt) - n_zeroes(original_gt)
          for original_gt, assigned_gt in zip(
              self.original_truth_genotypes, self.truth_genotypes
          )
      )
    return self._n_false_negatives

  def candidates_with_assigned_genotypes(self):
    """Gets a copy of our candidates with their matched genotypes.

    Returns:
      list[Variant protobuf]: Returns a copy of self.candidates in order, with
      genotypes corresponding to their matched genotypes. Any previous genotypes
      in these Variants will be overwrite. If no VariantCall is present one will
      be added.
    """
    with_gts = [copy.deepcopy(v) for v in self.candidates]
    for variant, gt in zip(with_gts, self.candidate_genotypes):
      call = variant.calls[0] if variant.calls else variant.calls.add()
      variantcall_utils.set_gt(call, gt)
    return with_gts


def deduplicate_haplotypes(haplotypes_to_genotypes_dict):
  """Returns a new dictionary with deduplicated value.

  Type description:
    Genotype = Tuple[int, int]
    Genotypes = List[Genotype]
    Haplotypes = FrozenSet[str]

  The type of the input `haplotypes_to_genotypes_dict` is:
    Dict[Haplotypes, List[Genotypes]]

  whereas the return type of deduplicate_haplotypes (this function) is:

  Dict[Haplotypes, Genotypes].

  This function goes through the values in `haplotypes_to_genotypes_dict` and
  keeps only a single example of Genotypes if there are multiple elements of
  the list that have the same haplotypes. Duplicates are expected
  in the list because different genotype configurations can sometimes produce
  the same set of haplotypes, and analyzing a dict of possible
  haplotypes/genotypes combinations with duplicates is much harder and less
  efficient than the deduplicated dict.

  Args:
    haplotypes_to_genotypes_dict: Dict[Haplotypes, List[Genotypes]].

  Returns:
    Dict[Haplotypes, Genotypes].
  """
  retval = {}
  for haplotypes, genotypes in haplotypes_to_genotypes_dict.items():
    # Keep the last element in the dedup process. The reason is for the
    # behavior to be the same as the previous implementation using list.
    retval[haplotypes] = genotypes[-1]
  return retval


# TODO: Create a comparison engine that accepts an iterable of
# variants and truths, and yields information about each variant and
# truth variant sequentially. This should be the primary API. Refactor
# label_examples to use this new API. Then create a new implementation that does
# the fast version.
def find_best_matching_haplotypes(candidates, truths, ref):
  """Assigns genotypes to each variant to best match truths.

  See the module-level documentation for general information on how this
  algorithm works.

  Args:
    candidates: list[nucleus.protos.Variant]. A list of candidate variants, in
      coordinate-sorted order, all on the same chromosome.
    truths: list[nucleus.protos.Variant]. A list of truth variants, in
      coordinate-sorted order, for the same interval on the genome as variants.
    ref: ReferenceRegion. Used to get reference bases for variants. Must cover
      at least the span of the variants.

  Returns:
    A HaplotypeMatch object describing the best assignment of genotypes between
    the candidates and truth_variants, or None, if no consistent assignment can
    be found.

  Raises:
    ValueError: If any inputs are malformed.
  """
  candidates = list(candidates)
  truths = list(truths)

  if _DEBUG_PRINTING_IS_ENABLED:
    _log_variants('candidates', candidates)
    _log_variants('truth', truths)

  if not variant_utils.variants_are_sorted(candidates):
    raise ValueError('candidates are not sorted', candidates)
  if not variant_utils.variants_are_sorted(truths):
    raise ValueError('truths are not sorted', truths)

  def _hom_ref_enum_if_empty(list_of_variants, non_empty_enum):
    """If list_of_variants is empty, use a ONLY_HOM_REF enum for speed."""
    return non_empty_enum if list_of_variants else EnumerationType.ONLY_HOM_REF

  truth_haplotypes = deduplicate_haplotypes(
      enumerate_all_possible_haplotypes(
          truths, ref, _hom_ref_enum_if_empty(candidates, EnumerationType.TRUTH)
      )
  )

  # Note, it may be worth deduplicating these haplotypes as well.
  variant_haplotypes = enumerate_all_possible_haplotypes(
      candidates,
      ref,
      _hom_ref_enum_if_empty(truths, EnumerationType.CANDIDATES),
  )

  found = []
  for vh, vgt_list in variant_haplotypes.items():
    tgt = truth_haplotypes.get(vh)
    # Explicitly check tgt is  None, because () is a valid value for tgt
    # here to continue.
    if tgt is None:
      continue
    for vgt in vgt_list:
      found.append(
          HaplotypeMatch(
              haplotypes=vh,
              candidates=candidates,
              candidate_genotypes=vgt,
              truths=truths,
              truth_genotypes=tgt,
          )
      )

  if not found:
    return None
  else:
    return select_best_haplotype_match(found)


def select_best_haplotype_match(all_matches):
  """Returns the best HaplotypeMatch among all_matches.

  The best matching HaplotypeMatch is the one with the lowest match_metrics
  score.

  Args:
    all_matches: iterable[HaplotypeMatch]. An iterable of HaplotypeMatch objects
      we want to select the best match from.

  Returns:
    The best matching HaplotypeMatch object.
  """
  sorted_matches = sorted(all_matches, key=lambda x: x.match_metrics)
  best = sorted_matches[0]
  equivalents = [
      f for f in all_matches if f.match_metrics == best.match_metrics
  ]

  # TODO: Why is this triggering so much?
  if len(equivalents) > 1:
    for i, f in enumerate(equivalents):
      extra_info = 'best' if i == 0 else i
      logging.warning('Equivalent match to best: %s [%s]', f, extra_info)

  return equivalents[0]


# -----------------------------------------------------------------------------
#
# Private utility functions
#
# -----------------------------------------------------------------------------


def _variant_genotypes(variants, missing_genotypes_default=(-1, -1)):
  """Returns the genotypes of variants as a list of tuples.

  Args:
    variants: iterable[nucleus.protos.Variant]. The variants whose genotypes we
      want to get.
    missing_genotypes_default: tuple. If a variant in variants doesn't have
      genotypes, this value is returned. The default value is (-1, -1), the
      standard representation for "missing" diploid genotypes.

  Returns:
    list[nucleus.protos.Variant] protos in the same order as variants.
  """
  return [
      tuple(v.calls[0].genotype) if v.calls else missing_genotypes_default
      for v in variants
  ]


def _log_haplotypes(name, haplotypes):
  """Write basic information about haplotypes to logging.info."""
  logging.info('haplotypes: %s', name)
  for haplotype, genotypes in haplotypes:
    logging.info('  %s with %s', haplotype, genotypes)


def _log_variants(name, variants):
  """Write basic information about variants to logging.info."""
  logging.info('variants: %s [%d]', name, len(variants))
  for v in variants:
    logging.info(
        '  %s gt=%s', variant_utils.variant_key(v), _variant_genotypes([v])[0]
    )


def n_zeroes(l):
  """Returns the number of elements of l that are 0."""
  return sum(1 for x in l if x == 0)


def _allele_from_index(variant, allele_index):
  """Gets the reference_bases or alternative_bases for allele_index."""
  if allele_index == 0:
    return variant.reference_bases
  else:
    return variant.alternate_bases[allele_index - 1]
