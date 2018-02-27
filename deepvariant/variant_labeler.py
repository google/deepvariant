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
"""variant_labeler for DeepVariant."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import abc



from absl import logging

from deepvariant.util.genomics import variants_pb2
from deepvariant.util import variant_utils


def make_labeler(truth_variants_reader, ref_reader, confident_regions, options):
  # redacted
  del options, ref_reader

  return PositionalVariantLabeler(
      vcf_reader=truth_variants_reader, confident_regions=confident_regions)

  # redacted
  # elif options.version == 'v2':
  #   return HaplotypeAwareLabeler(
  #       vcf_reader=truth_variants_reader,
  #       ref_reader=ref_reader,
  #       confident_regions=confident_regions,
  #       options=options)


# ---------------------------------------------------------------------------
# Variant counters
#
# This code was originally in make_examples and is now moved to the variant
# labeler itself because we will be refactoring the variant labeler to keep
# track of stats in the future. This code will be deleted when b/73783601 is
# completed.
# ---------------------------------------------------------------------------


class _Counter(object):

  def __init__(self, name, selectp):
    self.name = name
    self.selectp = selectp
    self.n_selected = 0


class VariantCounters(object):
  """Provides stats about the number of variants satisfying pfuncs."""

  def __init__(self, names_and_selectors):
    self.counters = []
    self.n_total = 0
    for name, selector in names_and_selectors:
      self.counters.append(_Counter(name, selector))

  def update(self, variant):
    self.n_total += 1
    for counter in self.counters:
      if counter.selectp(variant):
        counter.n_selected += 1

  def log(self):
    logging.info('----- VariantCounts -----')
    for counter in self.counters:
      percent = (100.0 * counter.n_selected) / (max(self.n_total, 1.0))
      logging.info('%s: %s/%s (%.2f%%)', counter.name, counter.n_selected,
                   self.n_total, percent)


def make_counters():
  """Creates all of the VariantCounters we want to track."""

  def _gt_selector(*gt_types):
    return lambda v: variant_utils.genotype_type(v) in gt_types

  return VariantCounters([
      ('All', lambda v: True),
      ('SNPs', variant_utils.is_snp),
      ('Indels', variant_utils.is_indel),
      ('BiAllelic', variant_utils.is_biallelic),
      ('MultiAllelic', variant_utils.is_multiallelic),
      ('HomRef', _gt_selector(variant_utils.GenotypeType.hom_ref)),
      ('Het', _gt_selector(variant_utils.GenotypeType.het)),
      ('HomAlt', _gt_selector(variant_utils.GenotypeType.hom_var)),
      ('NonRef',
       _gt_selector(variant_utils.GenotypeType.het,
                    variant_utils.GenotypeType.hom_var)),
  ])


# ---------------------------------------------------------------------------
# VariantLabel
#


# redacted
class VariantLabel(object):
  """Dataclass containing information about a label assigned to a variant.

  Attributes:
    is_confident: bool. True if we could confidently assign a label to this
      variant, False otherwise.
    variant: nucleus.protos.Variant proto that we assigned a label for.
    truth_variant: nucleus.protos.Variant proto containing the truth variant we
      will use to assign the label for this variant.
  """

  def __init__(self, is_confident, variant, truth_variant):
    self.is_confident = is_confident
    self.variant = variant
    self.truth_variant = truth_variant

  def label_for_alt_alleles(self, alt_alleles_indices):
    """Computes the label value for an example using alt_alleles_indices.

    Args:
      alt_alleles_indices: list[int]. A list of the alt_allele_indices
        used to compute the tf.Example for this candidate.

    Returns:
      int >= 0. The number of copies of alt_allele_indices we'd expect to be
      called for this example.
    """
    alt_alleles = [self.variant.alternate_bases[i] for i in alt_alleles_indices]
    return _match_to_alt_count(self.variant, self.truth_variant, alt_alleles)


# ---------------------------------------------------------------------------
# _VariantLabeler base class
#


class _VariantLabeler(object):
  """BaseClass for systems that want to provide training labels for examples.

  A VariantLabeler provides methods to assign a genotype label to each of a
  series of candidate variants using data from a truth set of variants
  accessible with vcf_reader and an optional RangeSet of confident regions.

  The basic logic of this class is something like:

  candidates = [third_party.nucleus.protos.Variant(...), ...]
  labeler = ConcreteSubclassOfVariantLabeler(vcf_reader, confident_regions)
  for label in labeler.label_variants(candidates):
    if label.is_confident:
      for i in range(len(label.variant.alternate_bases)
        genotype_label_value = label.label_for_alt_alleles([i])

  See the docs on each individual function to get a better understanding of what
  each function does and the meaning of the return values.

  """

  __metaclass__ = abc.ABCMeta

  # redacted
  def __init__(self, vcf_reader, confident_regions):
    if vcf_reader is None:
      raise ValueError('vcf_reader cannot be None')
    self._vcf_reader = vcf_reader
    self._confident_regions = confident_regions
    self.counters = make_counters()

  def log(self):
    """Logs information about the types and counts of labeled variants.

    This logs some basic summary information, like the number of labeled
    variants and their types, to logging.info.
    """
    self.counters.log()

  @abc.abstractmethod
  def label_variants(self, variants):
    """Gets label information for each variant in variants.

    This is the primary API for assigning labels to variants. This function
    takes and iterable of variants and yield a VariantLabel object for each
    variant. The VariantLabel can be used to determine the genotype label for
    each variant suitable for training a DeepVariant model. The API accepts
    an iterable of Variants because, in the general case, the labeling of
    variants aren't independent, in that the label assigned to one variant may
    impact the label we assign to a nearby variant.

    Args:
      variants: iterable[nucleus.protos.Variant]: An iterable of variants to
        label. The variants should be in coordinate-sorted order and all on the
        same chromosome.

    Yields:
      A VariantLabel object for each variant in variants, in order.
    """
    raise NotImplementedError


class PositionalVariantLabeler(_VariantLabeler):
  """Finds matching "truth" variants using a position-specific labeler.

  This is the original variant labeler from DeepVariant used up until v0.5,
  which assigns labels to variant calls by matching the chrom:position of a
  candidate variant with ones in truth, and then if one exists, assigns the
  label based on the genotype of the matched truth variant. This method works
  reasonably well but cannot handle complex representational differences between
  the candidate variants and the truth variants.
  """

  def __init__(self, vcf_reader, confident_regions=None):
    """Creates a new VariantLabeler.

    Args:
      vcf_reader: a VcfReader object that points to our truth variant set.
      confident_regions: A RangeSet containing all of the confidently called
        regions. A variant that falls outside of one of these regions will be
        receive a special not-confident marker.

    Raises:
      ValueError: if vcf_reader is None.
    """
    super(PositionalVariantLabeler, self).__init__(
        vcf_reader, confident_regions=confident_regions)

  def label_variants(self, variants):
    for variant in variants:
      is_confident, truth_variant = self._match(variant)
      if is_confident:
        self.counters.update(truth_variant)
      yield VariantLabel(
          is_confident=is_confident,
          variant=variant,
          truth_variant=truth_variant)

  def _match(self, variant):
    """Get a truth variant matching variant.

    A matching variant is defined here as one that starts at the same position
    on the genome as variant, regardless of the alleles of variant. This allows
    the client to make decisions on how to translate a matched between variant
    and truth_variant into a label (e.g. by comparing the alleles).

    This code will emit a logging.warning() if it detects multiple
    variants with the same chrom/start as variant provided by the
    vcf_reader and simply return the first variant. Though technically
    correct - VCF allows this - most files in practice merge variants
    that occur at the same location into a single multi-allelic variant
    record. So this assumption is reasonable. A future extension of the
    code could attempt to choose among the competing options, refuse to
    provide any answer, or even except out (presumably controlled by a
    configuration option).

    Args:
      variant: Our
        candidate third_party.nucleus.protos.Variant
        variant.

    Returns:
      A tuple of (match_status, truth_variant) where match_status is True if
      we are confident in our truth_variant call or False if not. truth_variant
      is a third_party.nucleus.protos.Variant object of
      the truth variant that matched
      variant, or None if none was found and we aren't confident in being
      hom-ref here, or a synthetic variant with the same position and alleles as
      variant but with a hom-ref genotype.
    """
    matched_variant = self._find_matching_variant_in_reader(variant)
    if self._confident_regions is None:
      confident = matched_variant is not None
    else:
      confident = self._confident_regions.variant_overlaps(
          variant, empty_set_return_value=False)
      if matched_variant is None and confident:
        matched_variant = self._make_synthetic_hom_ref(variant)
    return confident, matched_variant

  # redacted
  def _make_synthetic_hom_ref(self, variant):
    """Creates a version of variant with a hom-ref genotype.

    Args:
      variant: Our
        candidate third_party.nucleus.protos.Variant
        variant.

    Returns:
      A new Variant with the same position and alleles as variant but with a
      hom-ref genotype.
    """
    return variants_pb2.Variant(
        reference_name=variant.reference_name,
        start=variant.start,
        end=variant.end,
        reference_bases=variant.reference_bases,
        alternate_bases=variant.alternate_bases,
        calls=[variants_pb2.VariantCall(genotype=[0, 0])])

  def _find_matching_variant_in_reader(self, variant):
    """Finds a variant in vcf_reader compatible with variant, if one exists."""

    # redacted
    def _usable_truth(truth_variant):
      return (variant.start == truth_variant.start and
              not variant_utils.is_filtered(truth_variant))

    region = variant_utils.variant_position(variant)
    matches = [m for m in self._vcf_reader.query(region) if _usable_truth(m)]
    if not matches:
      return None
    elif len(matches) > 1:
      logging.warning(
          'Multiple matches detected, keeping first, for variant %s: %s',
          variant, matches)
    return matches[0]


# redacted
def _match_to_alt_count(candidate_variant, truth_variant, alt_alleles):
  """Returns the number of copies of alt_alleles that occur in truth_variant.

  This method figures out how many alternate copies of alt_alleles occur
  in truth_variant. For example, if candidate is A/C and truth is A/C with
  a 0/1 genotype, then this function would return 1 indicating there's one
  copy of the C allele in truth. If the true genotype is 1/1, then this
  routine would return 2.

  The routine allows candidate_variant and truth_variant to differ in both
  the number of alternate alleles, and even in the representation of the same
  alleles due to those differences. For example, candidate could be:

      AGT/A/AGTGT => 2 bp deletion and 2 bp insertion

  and truth could have:

      A/AGT => just the simplified 2 bp insertion

  And this routine will correctly equate the AGT/AGTGT allele in candidate
  with the A/AGT in truth and use the number of copies of AGT in truth to
  compute the number of copies of AGTGT.

  alt_alleles has a bit of a special meaning here. All of the alleles in
  alt_alleles are counted up, so if alt_alleles contains two alleles we'll
  compute the number of copies of each allele and return the sum. This allows
  us to compute the genotype of a synthetic combined allele that includes
  both alt1 and alt2.

  Args:
    candidate_variant: Our candidate third_party.nucleus.protos.Variant variant.
    truth_variant: Our third_party.nucleus.protos.Variant truth variant as
      returned by match().
    alt_alleles: An iterable of strings. Each element should be an alternate
      allele in variant that is considered "alt" when labeling variant.

  Returns:
    Number of copies of alt_alleles in the true genotype.

  Raises:
    ValueError: If candidate_variant or truth_variant is None, truth_variant
      doesn't have genotypes, or any of alt_alleles aren't found in
      candidate_variant.
  """
  if candidate_variant is None:
    raise ValueError('candidate_variant cannot be None')
  if truth_variant is None:
    raise ValueError('truth_variant cannot be None')
  if not variant_utils.has_genotypes(truth_variant):
    raise ValueError('truth_variant needs genotypes to be used for labeling',
                     truth_variant)
  if any(alt not in candidate_variant.alternate_bases for alt in alt_alleles):
    raise ValueError('All alt_alleles must be present in variant', alt_alleles,
                     candidate_variant)

  # If our candidate_variant is a reference call, return a label of 0.
  if variant_utils.is_ref(candidate_variant):
    return 0

  def _simplify_alleles(variant, alleles):
    return [
        variant_utils.simplify_alleles(variant.reference_bases, allele)
        for allele in alleles
        if allele != variant.reference_bases
    ]

  simplified_alt_alleles = _simplify_alleles(candidate_variant, alt_alleles)
  return sum(
      true_alt in simplified_alt_alleles for true_alt in _simplify_alleles(
          truth_variant, variant_utils.genotype_as_alleles(truth_variant)))
