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
"""Library for resolving variants into consistent haplotypes.

The convolutional neural network that evaluates the probability of a candidate
variant being non-reference evaluates each candidate variant independently.
This can lead to overlapping variant calls that cannot actually exist in an
organism: for example, a diploid human cannot have overlapping variants for
which one is homozygous alternate and the other is heterozygous alternate, since
that implies three total alternate alleles.

This library tries to resolve overlapping variant calls into consistent
haplotypes by using the most likely configuration based on individual call
probabilities that is a valid set of two haplotypes. In rare cases where this
is not possible, the haplotypes are left unmodified.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy
import itertools
from tensorflow import flags
import numpy as np

from absl import logging
from third_party.nucleus.util import genomics_math
from third_party.nucleus.util import variant_utils

FLAGS = flags.FLAGS

flags.DEFINE_bool(
    'disable_haplotype_resolution', False,
    'If True, makes `maybe_resolve_conflicting_variants` a no-op.')

# The maximum number of overlapping variants to try to resolve into compatible
# haplotypes. This corresponds to generating 3^12 (= 531,441) possible variant
# configurations for diploid individuals.
_MAX_OVERLAPPING_VARIANTS_TO_RESOLVE = 12


def maybe_resolve_conflicting_variants(sorted_variants):
  """Yields Variant protos in sorted order after fixing conflicting haplotypes.

  The input is an iterable of Variants in chromosome and position sorted order,
  with potential incompatibilies as described in this module's docstring. This
  function tries to resolve variants into valid haplotypes, though is not
  guaranteed to do so if the variant composition is not amenable to this or it
  would be computationally intractable.

  Args:
    sorted_variants: Iterable of Variant protos. Sorted in coordinate order, but
      with potentially incompatible haplotypes.

  Yields:
    Variant protos in coordinate-sorted order with no incompatible haplotypes.
  """
  if FLAGS.disable_haplotype_resolution:
    logging.info('disable_haplotype_resolution is True. '
                 '`maybe_resolve_conflicting_variants` has no effect.')
    for v in sorted_variants:
      yield v
  else:
    for overlapping_candidates in _group_overlapping_variants(sorted_variants):
      for resolved_candidate in _maybe_resolve_mixed_calls(
          overlapping_candidates):
        yield resolved_candidate


def _group_overlapping_variants(sorted_variants):
  """Yields lists of Variant protos that overlap on the reference sequence.

  Args:
    sorted_variants: Iterable of Variant protos, sorted in coordinate order.

  Yields:
    Lists of variants within `sorted_variants` that overlap with each other on
    the reference sequence.
  """
  curr_variants = []
  prev_chrom = None
  prev_max_end = -1
  for variant in sorted_variants:
    if variant.reference_name != prev_chrom or variant.start >= prev_max_end:
      if curr_variants:
        yield curr_variants
      curr_variants = [variant]
      prev_chrom = variant.reference_name
      prev_max_end = variant.end
    else:
      curr_variants.append(variant)
      prev_max_end = max(prev_max_end, variant.end)
  # Fencepost.
  if curr_variants:
    yield curr_variants


def _maybe_resolve_mixed_calls(overlapping_candidates):
  """Yields variants with compatible genotype calls in order.

  This function differs from `_resolve_overlapping_variants` below in that the
  input here is a block of all candidate calls that overlap in a region, which
  may contain candidates that are deemed to be most likely reference calls.
  We often tune DeepVariant to be highly sensitive. Consequently, there can be
  many candidate calls that are predicted as reference. Since those do not
  contribute to potential incompatibilities, we split them out from variants
  predicted to contain non-reference genotypes since the computation of
  compatible haplotypes is exponential in the number of inputs.

  Args:
    overlapping_candidates: list(Variant). A non-empty list of Variant protos in
      coordinate-sorted order that overlap on the reference genome.

  Yields:
    Variant protos in coordinate-sorted order that try to resolve incompatible
    haplotypes.
  """
  # Short circuit the simplest case: A single variant in a region is compatible
  # with itself by definition.
  if len(overlapping_candidates) == 1:
    yield overlapping_candidates[0]
    return

  def has_variation(candidate):
    return _nonref_genotype_count(candidate) > 0

  reference_calls = [c for c in overlapping_candidates if not has_variation(c)]
  variant_calls = [v for v in overlapping_candidates if has_variation(v)]

  resolved_variant_calls = []
  for variant_group in _group_overlapping_variants(variant_calls):
    resolved_variant_calls.extend(_resolve_overlapping_variants(variant_group))

  # Merge the reference and resolved variants back together in sorted order.
  # Note: This could be done in an interleaving fashion, but since the total
  # number of variants in the input is nearly always < 20 this is not an issue.
  for variant in sorted(
      reference_calls + resolved_variant_calls,
      key=variant_utils.variant_range_tuple):
    yield variant


class _VariantCompatibilityCalculator(object):
  """Represents the reference genome spanned by overlapping Variants.

  Each Variant affects a portion of the reference genome that is determined by
  its start and end coordinates. For a given set of Variants, they are deemed
  compatible if the total area along the reference genome that is called as
  non-reference genotypes never exceeds the ploidy of the organism.
  """

  def __init__(self, overlapping_variants):
    """Constructor.

    Args:
      overlapping_variants: list(Variant). The Variant protos of interest.
    """
    min_start = min(v.start for v in overlapping_variants)
    self.variant_indices = [
        (v.start - min_start, v.end - min_start) for v in overlapping_variants
    ]
    self.size = max(v.end - min_start for v in overlapping_variants)

  def all_variants_compatible(self, nonref_genotype_counts, ploidy=2):
    """Returns True if and only if all variants are compatible.

    Args:
      nonref_genotype_counts: list of ints in [0, ploidy]. Element i in this
        list represents the number of non-reference genotypes for the i'th
        variant.
      ploidy: int. The ploidy of the individual.

    Returns:
      True if and only if the variants are compatible.

    Raises:
      ValueError: nonref_genotype_counts is not the same length as
        self.variant_indices.
      ValueError: nonref_genotype_counts does not contain elements in [0,
      ploidy].
    """
    if len(nonref_genotype_counts) != len(self.variant_indices):
      raise ValueError(
          'Variant counts must have same length as variant indices.')
    if not all(0 <= cnt <= ploidy for cnt in nonref_genotype_counts):
      raise ValueError('Invalid variant allele count for ploidy {}: {}'.format(
          ploidy, nonref_genotype_counts))

    alts_in_span = np.zeros(self.size, dtype=int)
    for cnt, (start, end) in zip(nonref_genotype_counts, self.variant_indices):
      alts_in_span[start:end] += cnt
    return np.all(alts_in_span <= ploidy)


class _LikelihoodAggregator(object):
  """Container class for genotype likelihoods of allele configurations.

  When evaluating valid genotype configurations across multiple variants, we
  calculate the likelihood of each configuration. To then calculate the marginal
  likelihoods for each variant's genotypes, for each genotype we need to sum the
  probabilities of all configurations that include that genotype.

  For numerical stability we do this by storing the genotype likelihoods
  = log10(p) and then aggregate using the log-sum-exp trick.
  """

  def __init__(self, num_alts):
    """Constructor.

    Args:
      num_alts: int. The number of alternate alleles in the variant.
    """
    self._num_likelihoods = variant_utils.genotype_likelihood_index(
        (num_alts, num_alts)) + 1

    # At each GL index, we keep a list that will include the joint GL across all
    # variants that include that particular set of allele indices for this
    # variant.
    self._genotype_likelihood_containers = []
    for _ in xrange(self._num_likelihoods):
      self._genotype_likelihood_containers.append([])

  def add(self, allele_indices, likelihood):
    """Add some likelihood to a particular allele configuration.

    Args:
      allele_indices: Pair of (g1, g2) ints representing the genotype.
      likelihood: float. log10(probability of this genotype configuration).
    """
    ix = variant_utils.genotype_likelihood_index(allele_indices)
    self._genotype_likelihood_containers[ix].append(likelihood)

  def scaled_likelihoods(self):
    """Returns the scaled likelihood of each genotype."""
    if not all(bool(x) for x in self._genotype_likelihood_containers):
      raise ValueError(
          'All genotypes must have some probability mass: {}'.format(
              self._genotype_likelihood_containers))

    return genomics_math.normalize_log10_probs([
        genomics_math.log10sumexp(unscaled)
        for unscaled in self._genotype_likelihood_containers
    ])

  def most_likely_allele_indices(self):
    """Returns allele indices for the genotype with the largest likelihood."""
    ix = np.argmax(self.scaled_likelihoods())
    return variant_utils.allele_indices_for_genotype_likelihood_index(
        ix, ploidy=2)


def _resolve_overlapping_variants(overlapping_variants):
  """Yields variants with compatible haplotypes, if possible.

  Args:
    overlapping_variants: list(Variant). A non-empty list of Variant protos in
      coordinate-sorted order that overlap on the reference genome and are
      predicted to contain alternate allele genotypes.

  Yields:
    Variant protos in coordinate-sorted order that try to resolve incompatible
    haplotypes.
  """
  # Short circuit the simplest case: A single variant in a region is compatible
  # with itself by definition.
  if len(overlapping_variants) == 1:
    yield overlapping_variants[0]
    return

  # If the actual genotype calls are compatible, we can safely return those
  # since they would be the most likely configuration also when restricting to
  # only valid configurations of genotype calls.
  calculator = _VariantCompatibilityCalculator(overlapping_variants)
  nonref_counts = [_nonref_genotype_count(v) for v in overlapping_variants]
  if calculator.all_variants_compatible(nonref_counts):
    logging.info('Overlapping variants are naturally compatible: %s',
                 overlapping_variants)
    for variant in overlapping_variants:
      yield variant
    return

  # The actual genotype calls produce an inconsistent haplotype. If the number
  # of affected variants is "too large", avoid processing since this is an
  # exponential process.
  if len(overlapping_variants) > _MAX_OVERLAPPING_VARIANTS_TO_RESOLVE:
    logging.warning(
        'Overlapping variants are not naturally compatible, and there are too '
        'many to exhaustively search (%s). Returning variants without '
        'modification, beginning with %s.', len(overlapping_variants),
        overlapping_variants[0])
    for variant in overlapping_variants:
      yield variant
    return

  # Otherwise, the actual genotype calls are incompatible. Since the genotype
  # likelihoods are generally well-calibrated, we examine all configurations of
  # genotypes that create compatible haplotypes and retain the single
  # configuration with the highest joint likelihood across all variants as the
  # proposed genotype assignment. Separately, we rescale the likelihood of each
  # individual variant using only the valid genotype configurations. If the
  # results are concordant (i.e., the genotype predicted by the marginal
  # likelihood for each variant is the same as the genotype predicted when
  # maximizing the joint likelihood across all variants), we return variants
  # with those calls and the rescaled likelihoods. Otherwise, we log a warning
  # and emit the original (incompatible) variants.
  #
  # For example, a biallelic deletion with probabilities of homref, het, homalt
  # = 0.01, 0.9, 0.09 and inside it a biallelic SNP with probs 0.02, 0.48, 0.5.
  # Naively this would be called as a heterozygous indel and a homozygous SNP,
  # which is impossible as there are three total alternate genotypes. The
  # algorithm does the following:
  #
  #   Indel    SNP    Joint prob
  #   0/0      0/0    0.01 * 0.02 = 0.0002
  #   0/0      0/1    0.01 * 0.48 = 0.0048
  #   0/0      1/1    0.01 * 0.50 = 0.0050
  #   0/1      0/0    0.90 * 0.02 = 0.0180
  #   0/1      0/1    0.90 * 0.48 = 0.4320*
  #   0/1      1/1    <invalid>   = 0
  #   1/1      0/0    0.09 * 0.02 = 0.0018
  #   1/1      0/1    <invalid>   = 0
  #   1/1      1/1    <invalid>   = 0
  #
  #   So using the highest joint likelihood, we predict het indel and het SNP.
  #
  #   The marginal probability of each genotype for the indel is:
  #   0/0:  0.0002 + 0.0048 + 0.0050 = 0.01
  #   0/1:  0.0180 + 0.4320          = 0.45
  #   1/1:  0.0018                   = 0.0018
  #
  #   which after normalizing to sum to 1 is roughly 0.022, 0.974, 0.004.
  #   The marginal probability for the SNP, after performing similar
  #   calculations, is 0.043, 0.946, 0.011. So the marginals also predict a het
  #   indel and a het SNP. Since the two calculations agree, we use this
  #   genotype call and modified likelihoods.
  #
  # First, we find all non-reference count configurations that are compatible.
  # This represents each variant solely based on its number of non-reference
  # genotypes, and assumes that variants are compatible if the total number of
  # non-reference genotypes at a single position is at most two. By using
  # non-reference counts, we avoid testing multiple allele configurations that
  # will return the same result (e.g. a variant with two possible alternate
  # alleles has three allele configurations that are homozygous alternate
  # [1/1, 1/2, 2/2] and either all or none of them will be valid depending on
  # the variants it interacts with).
  valid_nonref_count_configurations = [
      conf for conf in itertools.product(
          [0, 1, 2], repeat=len(overlapping_variants))
      if calculator.all_variants_compatible(conf)
  ]

  # Next, we find the single compatible variant assignment with the individually
  # highest likelihood and track the total likelihood distributed to all variant
  # genotypes.
  likelihood_aggregators = [
      _LikelihoodAggregator(len(v.alternate_bases))
      for v in overlapping_variants
  ]
  most_likely_allele_indices_config = None
  most_likely_likelihood = None
  for nonref_count_config in valid_nonref_count_configurations:
    for allele_indices_config in _get_all_allele_indices_configurations(
        overlapping_variants, nonref_count_config):
      config_likelihood = _allele_indices_configuration_likelihood(
          overlapping_variants, allele_indices_config)
      if (most_likely_likelihood is None or
          config_likelihood > most_likely_likelihood):
        most_likely_likelihood = config_likelihood
        most_likely_allele_indices_config = allele_indices_config
      for aggregator, allele_indices in zip(likelihood_aggregators,
                                            allele_indices_config):
        aggregator.add(allele_indices, config_likelihood)

  marginal_allele_indices_config = tuple(
      agg.most_likely_allele_indices() for agg in likelihood_aggregators)
  if marginal_allele_indices_config == most_likely_allele_indices_config:
    logging.info(
        'Overlapping variants are not naturally compatible, but the genotype '
        'configuration with the most likely joint likelihood is the same as '
        'that from the scaled marginal likelihoods: %s',
        overlapping_variants[0])
    # Collapse the probabilities of all configurations to a single GL for each
    # allele, independently for each variant.
    scaled_gls = [agg.scaled_likelihoods() for agg in likelihood_aggregators]

    for variant, allele_indices, gls in zip(
        overlapping_variants, most_likely_allele_indices_config, scaled_gls):
      newvariant = copy.deepcopy(variant)
      call = variant_utils.only_call(newvariant)
      call.genotype[:] = allele_indices
      call.genotype_likelihood[:] = gls
      yield newvariant
  else:
    logging.warning(
        'Overlapping variants are not naturally compatible, and the genotype '
        'configuration with the most likely joint likelihood is different from '
        'that using the scaled marginal likelihoods: %s',
        overlapping_variants[0])
    # redacted
    for variant in overlapping_variants:
      yield variant


def _get_all_allele_indices_configurations(variants,
                                           nonref_count_configuration):
  """Returns an iterable of allele configurations that satisfy the genotype.

  Args:
    variants: list(Variant). The list of variants for which to generate
      configurations of valid allele_indices.
    nonref_count_configuration: list(int). The list of numbers of non-reference
      genotypes that should be generated for each variant.

  Returns:
    Iterable of lists of allele indices to assign to each Variant to satisfy the
    desired configuration of number of non-reference genotypes for each variant.

  Raises:
    ValueError: variants and nonref_count_configuration do not have the same
    length.
  """
  if len(variants) != len(nonref_count_configuration):
    raise ValueError(
        'len(variants) must equal len(nonref_count_configuration): {} vs {}'.
        format(len(variants), len(nonref_count_configuration)))

  allele_indices_configs = [
      variant_utils.allele_indices_with_num_alts(variant, num_alts, ploidy=2)
      for variant, num_alts in zip(variants, nonref_count_configuration)
  ]
  return itertools.product(*allele_indices_configs)


def _allele_indices_configuration_likelihood(variants, allele_indices_config):
  """Returns the joint likelihood of the alleles given to the variants.

  Args:
    variants: list(Variant). The variants with associated likelihoods.
    allele_indices_config: list((int, int)). The allele indices to assign to
      each variant.

  Returns:
    The joint likelihood of the particular allele configuration.

  Raises:
    ValueError: variants and allele_indices_config do not have the same length.
  """
  if len(variants) != len(allele_indices_config):
    raise ValueError(
        'len(variants) must equal len(allele_indices_config): {} vs {}'.format(
            len(variants), len(allele_indices_config)))

  retval = 0
  for variant, alleles in zip(variants, allele_indices_config):
    retval += variant_utils.genotype_likelihood(
        variant_utils.only_call(variant), alleles)
  return retval


def _nonref_genotype_count(variant):
  """Returns the number of non-reference alleles in the called genotype."""
  return sum(g > 0 for g in variant_utils.only_call(variant).genotype)
