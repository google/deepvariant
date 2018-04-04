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
"""Haplotype-based labeling algorithm for DeepVariant.

redacted
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import copy
import itertools

from absl import logging
import enum

from third_party.nucleus.io import fasta
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from deepvariant.labeler import variant_labeler


VariantAndGenotypes = collections.namedtuple('VariantAndGenotype',
                                             ['variant', 'genotypes'])

# The default maximum size of a variant group we'll try to label. See
# the HaplotypeLabeler class for more information.
_MAX_GROUP_SIZE = 6
# The default maximum distance between subsequent variants within a group. See
# the HaplotypeLabeler class for more information.
_MAX_DISTANCE_WITHIN_VARIANT_GROUP = 30
# When querying for the truth variants within our span, we extend the query
# region by this many basepairs to capture nearby truth variants.
_TRUTH_VARIANTS_QUERY_REGION_EXPANSION_IN_BP = 10

# redacted
# True we will generate enough information into our logs to help debug bad
# regions.
_DEBUG_PRINTING_IS_ENABLED = True


class HaplotypeLabeler(variant_labeler.VariantLabeler):
  """Haplotype-based variant labeler."""

  def __init__(
      self,
      truth_vcf_reader,
      ref_reader,
      confident_regions,
      max_group_size=_MAX_GROUP_SIZE,
      max_distance_within_grouped_variants=_MAX_DISTANCE_WITHIN_VARIANT_GROUP):
    """Creates a new HaplotypeVariantLabeler.

    Args:
      truth_vcf_reader: a VcfReader object that points to our truth variant set.
      ref_reader: A FastaReader object we can use to get reference bases.
      confident_regions: A RangeSet containing all of the confidently called
        regions. A variant that falls outside of one of these regions will be
        receive a special not-confident marker.
      max_group_size: int >= 1. The maximum number of variants we'll attempt to
        label together. Larger values increase the runtime of the algorithm.
      max_distance_within_grouped_variants: int >= 0. The maximum distance
        between variants within a group. Sequential variants separated by more
        than this value will be placed in separate groups for labeling.

    Raises:
      ValueError: if vcf_reader is None.
    """
    super(HaplotypeLabeler, self).__init__(
        truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions)
    self._ref_reader = ref_reader
    self.max_group_size = max_group_size
    self.max_distance_within_grouped_variants = (
        max_distance_within_grouped_variants)

  def label_variants(self, variants):
    # redacted
    for variant_group in self.group_variants(variants):
      for label in self._label_grouped_variants(variant_group):
        yield label

  def group_variants(self, variants):
    # redacted
    # redacted
    # we don't want to miss any truths (for FN counts) or double count truths.
    groups = []
    current_group = []
    for variant in variants:
      if self._include_in_variant_group(current_group, variant):
        current_group.append(variant)
      else:
        groups.append(current_group)
        current_group = [variant]

    if current_group:
      groups.append(current_group)
    return groups

  def _include_in_variant_group(self, group, variant):
    if not group:
      return True
    elif len(group) >= self.max_group_size:
      return False
    else:
      last_variant = group[-1]
      assert variant.reference_name == last_variant.reference_name
      return (variant.start - last_variant.start <=
              self.max_distance_within_grouped_variants)

  def _label_grouped_variants(self, variants):
    # redacted

    # redacted
    # they should be computed in the grouping.
    span = ranges.span([variant_utils.variant_range(v) for v in variants])
    truths = list(
        self._get_truth_variants(
            ranges.expand(span, _TRUTH_VARIANTS_QUERY_REGION_EXPANSION_IN_BP)))

    if len(truths) > self.max_group_size:
      logging.warning(
          ('Found a large number of variants to label (n_candidates=%d, '
           'n_truth=%d) relative to candidate cap of %d. This may make the '
           'algorithm very slow.'), len(variants), len(truths),
          self.max_group_size)
      # redacted
      logging.warning('Returning all variants with not-confident markers.')
      for variant in variants:
        yield variant_labeler.VariantLabel(
            is_confident=False, genotype=(-1, -1), variant=variant)
      return
    ref = self.make_labeler_ref(variants, truths)
    labeled_variants = label_variants(variants, truths, ref)

    if not labeled_variants:
      raise ValueError('Failed to assign labels for variants', variants)
    else:
      for labeled in labeled_variants:
        yield variant_labeler.VariantLabel(
            # redacted
            # now. Rethink how we establish a variant is confident. Seems like
            # it'd be confident if it has a non-ref genotype (as we only
            # consider confident truth variants) or if it overlaps the confident
            # regions.
            is_confident=self._confident_regions.variant_overlaps(labeled),
            genotype=tuple(labeled.calls[0].genotype),
            variant=labeled)

  def make_labeler_ref(self, candidate_variants, true_variants, bufsize=20):
    all_variants = candidate_variants + true_variants
    contig = all_variants[0].reference_name
    start = min(x.start for x in all_variants)
    end = max(x.end for x in all_variants)
    region = ranges.make_range(contig, start - 1, end + bufsize)
    ref_bases = self._ref_reader.query(region)
    return ReferenceRegion(ref_bases, start=region.start)


class EnumerationType(enum.Enum):
  CANDIDATES = 1
  TRUTH = 2


class ReferenceRegion(fasta.InMemoryRefReader):
  """Allows us to get bases from a cached reference interval."""

  # We don't want to worry about the chromosome we are working on for code
  # clarity, so we create a InMemoryRefReader that has a single chromosome named
  # _DUMMY_CHROM_NAME which allows us to provide a bases(start, end) function
  # for convenient reading of bases.
  _DUMMY_CHROM_NAME = '*'

  def __init__(self, bases, start):
    super(ReferenceRegion, self).__init__([(self._DUMMY_CHROM_NAME, start,
                                            bases)])
    self.start = start
    self.end = start + len(bases)

  def bases(self, start, end):
    return self.query(ranges.make_range(self._DUMMY_CHROM_NAME, start, end))


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
    A set of tuples containing diploid genotypes.
  """
  alts = set(gt) - {0}
  return {(0, 0), tuple(gt)} | {(0, alt) for alt in alts}


def enumerate_all_possible_haplotypes(variants, ref, enumeration_type):
  """Yields all possible haplotype/genotype combinations for variants.

  Args:
    variants: list[nucleus.protos.Variant]. A list of candidate variants, in
      coordinate-sorted order, all on the same chromosome.
    ref: ReferenceRegion. Used to get reference bases for variants. Must cover
      at least the span of the variants.
    enumeration_type: EnumerationType enum value. What kind of enumeration do we
      want to do? Can be either CANDIDATES or TRUTH.

  Yields:
    2-tuple of haplotypes and genotypes. Haplotypes is a set of either one or
    two strings, where each string is a haplotype (i.e., a series of bases)
    generated by the genotypes assigned to each variant. genotypes is a list of
    genotype tuples, in the same order as variants, indicating the genotype
    assignment for each variant. These genotypes are phased, so [(0, 1), (0, 1)]
    is not the same as [(0, 1), (1, 0)].
  """
  def create_haplotypes(variants_and_genotypes, last_pos, depth=0):

    if not variants_and_genotypes:
      yield {ref.bases(last_pos, ref.end)}
    else:
      # redacted
      group, remaining = split_independent_variants(variants_and_genotypes)
      group_haplotypes, next_pos = build_all_haplotypes(group, last_pos, ref)
      frags = list(all_diploid_haplotypes(group, group_haplotypes))

      for haplotypes in create_haplotypes(remaining, next_pos, depth + 1):
        for result in extend_haplotypes(frags, haplotypes):
          yield result
  genotype_options = genotype_options_for_variants(variants, enumeration_type)
  for genotypes in itertools.product(*genotype_options):
    paired = [VariantAndGenotypes(v, g) for v, g in zip(variants, genotypes)]
    for haplotypes in create_haplotypes(paired, ref.start):
      yield haplotypes, genotypes


# redacted
# review of all of the function and variable names and replace with more
# meaningful options.
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
  for haploid_genotype, haplotype in genotypes2haplotype.iteritems():
    complement = complement_haploid_genotype(haploid_genotype, genotypes)
    complement_haplotype = genotypes2haplotype.get(complement, None)
    if complement_haplotype is not None and complement not in generated_already:
      generated_already.add(haploid_genotype)
      yield {haplotype, complement_haplotype}


def genotype_options_for_variants(variants, enumeration_type):
  if enumeration_type == EnumerationType.TRUTH:
    return [
        with_false_negative_genotypes(x) for x in _variant_genotypes(variants)
    ]
  elif enumeration_type == EnumerationType.CANDIDATES:
    return [
        [(i, j)
         for i, j, _, _ in variant_utils.genotype_ordering_in_likelihoods(v)]
        for v in variants
    ]
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


def extend_haplotypes(fragments_list, haplotypes):
  """Yields all diploid combinations of frags x haplotypes.

  Args:
    fragments_list: list[set[string]]: fragments_list contains a list of
      set[string], which are just like haplotypes (i.e., contains 1 or 2
      strings), that collectively represent all possible prefixes of haplotypes.
    haplotypes: set[string]. A set containing 1 or 2 haplotype strings. So it
      looks like {h} or {h1, h2}.

  Yields:
    A series of set[string], each containing 1 or 2 haplotype strings.
  """
  # redacted
  for frags in fragments_list:
    if len(frags) == 1:
      f = next(iter(frags))
      yield {f + h for h in haplotypes}
    else:
      f1, f2 = list(frags)
      if len(haplotypes) == 1:
        h = next(iter(haplotypes))
        yield {f1 + h, f2 + h}
      else:
        h1, h2 = list(haplotypes)
        yield {f1 + h1, f2 + h2}
        yield {f1 + h2, f2 + h1}


# redacted
# redacted
def build_all_haplotypes(variants_and_genotypes, last_pos, ref):
  # redacted
  if len(variants_and_genotypes) == 1:
    variant = variants_and_genotypes[0].variant
    genotypes = variants_and_genotypes[0].genotypes
    ref_prefix = ref.bases(last_pos, variant.start)
    base_choices = [variant.reference_bases] + list(variant.alternate_bases)
    frags = {(alt,): ref_prefix + base_choices[alt] for alt in set(genotypes)}
    return frags, variant.end
  else:
    frags = {}
    genotypes = [vg.genotypes for vg in variants_and_genotypes]
    variants = [vg.variant for vg in variants_and_genotypes]
    all_haploid_genotypes = sorted(set(itertools.product(*genotypes)))
    end = max(v.end for v in variants)
    for phased in all_haploid_genotypes:
      haplotype = build_haplotype(variants, phased, last_pos, end, ref)
      if haplotype:
        frags[phased] = haplotype
    return frags, end


# redacted
# redacted
# redacted
# redacted
def build_haplotype(variants, phased_genotypes, last_pos, end, ref):
  """Builds the haplotype string from variants and its phased gneotypes."""
  parts = []
  for variant, genotype in zip(variants, phased_genotypes):
    if variant.start < last_pos:
      if genotype != 0:
        return None
    else:
      ref_prefix = ref.bases(last_pos, variant.start)
      base_choices = [variant.reference_bases[0]] + list(
          variant.alternate_bases)
      last_pos = variant.start + 1 if genotype == 0 else variant.end
      parts.append(ref_prefix + base_choices[genotype])

  if last_pos < end:
    parts.append(ref.bases(last_pos, end))

  return ''.join(parts)


def print_haplotypes(name, haplotypes):
  logging.info('haplotypes: %s', name)
  for haplotype, genotypes in haplotypes:
    logging.info('  %s with %s', haplotype, genotypes)


def print_variants(name, variants):
  logging.info('variants: %s [%d]', name, len(variants))
  for v in variants:
    logging.info('  %s gt=%s', variant_utils.variant_key(v),
                 _variant_genotypes([v])[0])


# redacted
# redacted
class LabelerMatch(object):
  """DataClass holding information about a labeling of variants."""

  def __init__(self, haplotypes, variants, matched_variant_genotypes,
               truth_variants, matched_truth_genotypes):
    self.haplotypes = sorted(haplotypes)
    self.variants = variants
    self.truth_variants = truth_variants
    self.matched_variant_genotypes = matched_variant_genotypes
    self.matched_truth_genotypes = matched_truth_genotypes

    # Computed on-demand.
    self._n_false_positives = None
    self._n_false_negatives = None

    if any(sum(gt) == 0 for gt in _variant_genotypes(self.truth_variants)):
      raise ValueError('No truth genotypes should be hom-ref')
    assert len(self.matched_variant_genotypes) == len(self.variants)
    assert len(self.matched_truth_genotypes) == len(self.truth_variants)

  def __str__(self):
    return ('LabelerMatch(haplotypes={}, false_negatives={}, '
            'false_positives={} true_positives={} match_quality={}, '
            'variant_gts={}, true_gts={})').format(
                self.haplotypes, self.n_false_negatives, self.n_false_positives,
                self.n_true_positives, self.match_quality,
                self.matched_variant_genotypes, self.truth_genotypes)

  __repr__ = __str__

  @property
  def truth_genotypes(self):
    return _variant_genotypes(self.truth_variants)

  # redacted
  @property
  def match_quality(self):
    """Quality of this match. Lower scores are better.

    Returns:
      tuple[int] where all elements are >= 0: The tuple is suitable for sorting
      matches, so that sorted(matches, key=lambda x: x.match_quality) will rank
      matches so that the best option is first.
    """
    return (self.n_false_negatives, self.n_false_positives,
            self.n_true_positives)

  # redacted
  @property
  def n_true_positives(self):
    """Gets the number of variants whose matched genotype is not (0, 0).

    Since the variants don't have expected genotypes, we can only count each
    site instead of each genotype. So this is the number of variants whose
    matched genotype is not (0, 0).

    Returns:
      int >= 0.
    """
    return len(self.matched_variant_genotypes) - self.n_false_positives

  # redacted
  @property
  def n_false_positives(self):
    """Gets the number of variants whose matched genotype is (0, 0).

    Since the variants don't have expected genotypes, we can only count each
    site instead of each genotype. So this is the number of variants whose
    matched genotype is (0, 0).

    Returns:
      int >= 0.
    """

    if self._n_false_positives is None:
      self._n_false_positives = sum(
          sum(gt) == 0 for gt in self.matched_variant_genotypes)
    return self._n_false_positives

  # redacted
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

      def n_zeroes(l):
        return sum(1 for x in l if x == 0)

      # redacted
      self._n_false_negatives = sum(
          n_zeroes(assigned_gt) - n_zeroes(true_gt)
          for true_gt, assigned_gt in zip(self.truth_genotypes,
                                          self.matched_truth_genotypes))

    return self._n_false_negatives

  def variants_with_assigned_genotypes(self):
    """Gets a copy of our variants with their matched genotypes.

    Returns:
      list[Variant protobuf]: Returns a copy of self.variants in order, with
      genotypes corresponding to their matched genotypes. Any previous genotypes
      in these Variants will be overwrite. If no VariantCall is present one will
      be added.
    """
    with_gts = [copy.deepcopy(v) for v in self.variants]
    for variant, gt in zip(with_gts, self.matched_variant_genotypes):
      if variant.calls:
        variant.calls[0].genotype[:] = gt
      else:
        variant.calls.add(genotype=gt)
    return with_gts


def deduplicate_haplotypes(haplotypes):
  haplotypes = list(haplotypes)
  return [(vh1, vg1)
          for i, (vh1, vg1) in enumerate(haplotypes)
          if not any(vh1 == vh2 for (vh2, _) in haplotypes[i + 1:])]


# redacted
# variants and truth_variants, and yields information about each variant and
# truth variant sequentially. This should be the primary API. Refactor
# label_examples to use this new API. Then create a new implementation that does
# the fast version.
def label_variants(variants, truth_variants, ref):
  """Assigns genotypes to each variant to best match truth_variants.

  See the module-level documentation for general information on how this
  algorithm works.

  Args:
    variants: list[nucleus.protos.Variant]. A list of candidate variants, in
      coordinate-sorted order, all on the same chromosome.
    truth_variants: list[nucleus.protos.Variant]. A list of truth variants, in
      coordinate-sorted order, for the same interval on the genome as variants.
    ref: ReferenceRegion. Used to get reference bases for variants. Must cover
      at least the span of the variants.

  Returns:
    A list of new variants, copied from variants, but with their
    call[0].genotype field assigned values for the optimal labeling of variants.

  Raises:
    ValueError: If any inputs are malformed.
  """
  variants = list(variants)
  truth_variants = list(truth_variants)

  if _DEBUG_PRINTING_IS_ENABLED:
    print_variants('candidates', variants)
    print_variants('truth', truth_variants)

  if not variant_utils.variants_are_sorted(variants):
    raise ValueError('Variants are not sorted', variants)
  if not variant_utils.variants_are_sorted(truth_variants):
    raise ValueError('truth_variants are not sorted', truth_variants)

  truth_haplotypes = deduplicate_haplotypes(
      enumerate_all_possible_haplotypes(truth_variants, ref,
                                        EnumerationType.TRUTH))

  # redacted
  variant_haplotypes = enumerate_all_possible_haplotypes(
      variants, ref, EnumerationType.CANDIDATES)

  found = []
  for vh, vgt in variant_haplotypes:
    for th, tgt in truth_haplotypes:
      if th == vh:
        found.append(
            LabelerMatch(
                haplotypes=th,
                variants=variants,
                matched_variant_genotypes=vgt,
                truth_variants=truth_variants,
                matched_truth_genotypes=tgt))

  if not found:
    return None
  else:
    best = select_best_match(found)
    return best.variants_with_assigned_genotypes()


def select_best_match(all_matches):
  """Returns the best LabelerMatch among all_matches.

  The best matching LabelerMatch is the one with the lowest match_quality score.

  Args:
    all_matches: iterable[LabelerMatch]. An iterable of LabelerMatch objects we
      want to select the best match from.

  Returns:
    The best matching LabelerMatch object.
  """
  sorted_matches = sorted(all_matches, key=lambda x: x.match_quality)
  best = sorted_matches[0]
  equivalents = [
      f for f in all_matches if f.match_quality == best.match_quality
  ]

  # redacted
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
  return [
      tuple(v.calls[0].genotype) if v.calls else missing_genotypes_default
      for v in variants
  ]
