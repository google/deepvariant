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



from absl import logging

from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from deepvariant.labeler import variant_labeler


class PositionalVariantLabeler(variant_labeler.VariantLabeler):
  """Finds matching "truth" variants using a position-specific labeler.

  This is the original variant labeler from DeepVariant used up until v0.5,
  which assigns labels to variant calls by matching the chrom:position of a
  candidate variant with ones in truth, and then if one exists, assigns the
  label based on the genotype of the matched truth variant. This method works
  reasonably well but cannot handle complex representational differences between
  the candidate variants and the truth variants.
  """

  def __init__(self, truth_vcf_reader, confident_regions):
    """Creates a new VariantLabeler.

    Args:
      truth_vcf_reader: a VcfReader object that points to our truth variant set.
      confident_regions: A RangeSet containing all of the confidently called
        regions. A variant that falls outside of one of these regions will be
        receive a special not-confident marker.

    Raises:
      ValueError: if vcf_reader is None.
    """
    super(PositionalVariantLabeler, self).__init__(
        truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions)

  def label_variants(self, variants):
    for variant in variants:
      is_confident, truth_variant = self._match(variant)

      genotype = None
      if truth_variant is not None:
        genotype = _genotype_from_matched_truth(variant, truth_variant)

      yield variant_labeler.VariantLabel(
          is_confident=is_confident, variant=variant, genotype=genotype)

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
    confident = self._confident_regions.variant_overlaps(
        variant, empty_set_return_value=False)
    if matched_variant is None and confident:
      matched_variant = self._make_synthetic_hom_ref(variant)
    return confident, matched_variant

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
    region = variant_utils.variant_position(variant)
    matches = [
        truth_variant for truth_variant in self._get_truth_variants(region)
        if variant.start == truth_variant.start
    ]

    if not matches:
      return None
    elif len(matches) > 1:
      logging.warning(
          'Multiple matches detected, keeping first, for variant %s: %s',
          variant, matches)
    return matches[0]


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
      variant_utils.only_call(truth_variant)):
    raise ValueError('truth_variant needs genotypes to be used for labeling',
                     truth_variant)

  def _match_one_allele(true_allele):
    if true_allele == truth_variant.reference_bases:
      return 0
    else:
      simplified_true_allele = variant_utils.simplify_alleles(
          truth_variant.reference_bases, true_allele)
      for alt_index, alt_allele in enumerate(candidate_variant.alternate_bases):
        simplified_alt_allele = variant_utils.simplify_alleles(
            candidate_variant.reference_bases, alt_allele)
        if simplified_true_allele == simplified_alt_allele:
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
            _match_one_allele(true_allele) for true_allele in
            variant_utils.genotype_as_alleles(truth_variant)))
