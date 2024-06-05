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
"""variant_labeler for DeepVariant."""

import abc



from absl import logging

from deepvariant import dv_utils
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

# ---------------------------------------------------------------------------
# VariantLabel
#


class VariantLabel(object):
  """Dataclass containing information about a label assigned to a variant.

  Attributes:
    is_confident: bool. True if we could confidently assign a label to this
      variant, False otherwise.
    variant: nucleus.protos.Variant proto that we assigned a label for.
    genotype: tuple of ints. The labeled genotype (e.g., (0, 1) for a het) in
      the standard nucleus.proto.VariantCall style. Genotype can be None if the
      labeler doesn't have any genotype to assign. If Genotype is not None, the
      genotype of variant will be set to genotype.
  """

  def __init__(self, is_confident, variant, genotype=None):
    if genotype is not None:
      if not variant.calls:
        variant.calls.add(genotype=genotype)
      else:
        variant.calls[0].genotype[:] = genotype

    self.is_confident = is_confident
    self.variant = variant
    self.genotype = genotype

  def label_for_alt_alleles(self, alt_alleles_indices):
    """Computes the label value for an example using alt_alleles_indices.

    This function computes the TensorFlow label value (0, 1, 2) we train
    DeepVariant to predict. The label value is an int >= which is the number of
    copies of the alt allele present, which is computed from the true genotypes
    (self.genotypes) and the alt_allele_indices ([0] for the first alt, [1] for
    the second, [0, 1] to combine the first and second). For example, suppose we
    have a variant with alts A and C, and a true genotype of (0, 1), indicating
    that we have 1 copy of the A allele. We'd expect:

      label_for_alt_alleles([0]) => 1 since there's 1 copy of the first alt.
      label_for_alt_alleles([1]) => 0 since there's 0 copies of the second alt.
      label_for_alt_alleles([0, 1]) => 1 since there's 1 copy of the first or
        second allele.

    Args:
      alt_alleles_indices: list[int]. A list of the alt_allele_indices used to
        compute the tf.Example for this candidate.

    Returns:
      int >= 0. The number of copies of alt_allele_indices we'd expect to be
      called for this example.
    """
    return sum(gt - 1 in alt_alleles_indices for gt in self.genotype if gt != 0)

  def convert_to_class(self, example):
    """Convert label to the class id."""
    alt_alleles_indices = dv_utils.example_alt_alleles_indices(example)
    # Set the label of the example to the # alts given our alt_alleles_indices.
    return self.label_for_alt_alleles(alt_alleles_indices)

  def get_class(self):
    """Returns the class of a given label.

    This function does not consider multi-allelic variants. It is used for the
    purpose of downsampling.
    """
    raise NotImplementedError

  def set_variant_genotype(self, variant):
    if not variant.calls:
      variant.calls.add(genotype=self.genotype)
    else:
      variant.calls[0].genotype[:] = self.genotype


# ---------------------------------------------------------------------------
# _VariantLabeler base class
#


class VariantLabeler(object):
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

  def __init__(self, truth_vcf_reader, confident_regions=None):
    if truth_vcf_reader is None:
      raise ValueError('truth_vcf_reader cannot be None')
    if confident_regions is None:
      logging.warning(
          'Note: confident_regions for VariantLabeler is None. '
          'It is possible that this is not allowed for some '
          'subtype of VariantLabelers.'
      )
    self._truth_vcf_reader = truth_vcf_reader
    self._confident_regions = confident_regions

  @property
  def metrics(self):
    """Gets the LabelingMetrics proto tracking metrics for this labeler.

    A variant labeler may provide information information about the labeling of
    variants via the metrics property. If the labeler provides metrics, a
    filled in LabelingMetrics protobuf will be returned by this property. If
    metrics aren't supported, a None value will be returned.
    """
    return None

  @abc.abstractmethod
  def label_variants(self, variants, region):
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
      region: A nucleus.genomics.v1.Range object specifying the region over
        which we are labeling variants. This should span at least the span of
        variants, but may be larger. Statistics about the labeling will be
        computed over region.

    Yields:
      A VariantLabel object for each variant in variants, in order.
    """
    raise NotImplementedError

  def _get_truth_variants(self, region):
    """Gets truth variants within region to use in labeling calculations.

    This function queries _truth_vcf_reader in region to get a complete list of
    truth variants that overlap region, and then filters them down by removing
    filtered truth variants and ones that aren't contained in the truth
    intervals.

    Args:
      region: nucleus.Range proto describing the region on the genome where we
        want to get our truth variants.

    Yields:
      nucleus.Variant proto.
    """
    for variant in self._truth_vcf_reader.query(region):
      if not variant_utils.is_filtered(variant) and (
          self._confident_regions is None
          or self._confident_regions.variant_overlaps(
              variant, empty_set_return_value=False
          )
      ):
        yield variant


def _genotype_from_matched_truth(candidate_variant, truth_variant):
  """Gets the diploid genotype for candidate_variant from matched truth_variant.

  This method figures out the genotype for candidate_variant by matching alleles
  in candidate_variant with those used by the genotype assigned to
  truth_variant. For example, if candidate is A/C and truth is A/C with a 0/1
  genotype, then this function would return (0, 1) indicating that there's one
  copy of the A allele and one of C in truth. If the true genotype is 1/1, then
  this routine would return (1, 1).

  The routine allows candidate_variant and truth_variant to differ in both
  the number of alternate alleles, and even in the representation of the same
  alleles due to those differences. For example, candidate could be:

      AGT/A/AGTGT => 2 bp deletion and 2 bp insertion

  and truth could have:

      A/AGT => just the simplified 2 bp insertion

  And this routine will correctly equate the AGT/AGTGT allele in candidate
  with the A/AGT in truth and use the number of copies of AGT in truth to
  compute the number of copies of AGTGT when determining the returned genotype.

  Args:
    candidate_variant: Our candidate third_party.nucleus.protos.Variant variant.
    truth_variant: Our third_party.nucleus.protos.Variant truth variant
      containing true alleles and genotypes.

  Returns:
    A tuple genotypes with the same semantics at the genotype field of the
    VariantCall proto.

  Raises:
    ValueError: If candidate_variant is None, truth_variant is None, or
      truth_variant doesn't have genotypes.
  """
  if candidate_variant is None:
    raise ValueError('candidate_variant cannot be None')
  if truth_variant is None:
    raise ValueError('truth_variant cannot be None')
  if not variantcall_utils.has_genotypes(
      variant_utils.only_call(truth_variant)
  ):
    raise ValueError(
        'truth_variant needs genotypes to be used for labeling', truth_variant
    )

  def _match_one_allele(true_allele):
    if true_allele == truth_variant.reference_bases:
      return 0
    else:
      simplifed_true_allele = variant_utils.simplify_alleles(
          truth_variant.reference_bases, true_allele
      )
      for alt_index, alt_allele in enumerate(candidate_variant.alternate_bases):
        simplifed_alt_allele = variant_utils.simplify_alleles(
            candidate_variant.reference_bases, alt_allele
        )
        if simplifed_true_allele == simplifed_alt_allele:
          return alt_index + 1
      # If nothing matched, we don't have this alt, so the alt allele index for
      # should be 0 (i.e., not any alt).
      return 0

  # If our candidate_variant is a reference call, return a (0, 0) genotype.
  if variant_utils.is_ref(candidate_variant):
    return (0, 0)
  else:
    return tuple(
        sorted(
            _match_one_allele(true_allele)
            for true_allele in variant_utils.genotype_as_alleles(truth_variant)
        )
    )
