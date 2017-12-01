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

from deepvariant.core import variantutils
from deepvariant.core.genomics import variants_pb2


class VariantLabeler(object):
  """Finds matching "truth" variants.

  A VariantLabeler provides methods to find a compatible "truth" variant for a
  given candidate variant using data from a vcf_reader and an optional RangeSet
  of confident regions. The basic logic of this class is something like:

  candidate = learning.genomics.deepvariant.core.genomics.Variant(...)
  labeler = VariantLabeler(vcf_reader, confident_regions)
  is_confident, truth_variant = labeler.match(candidate)
  if is_confident:
    alts = candidate.alternate_bases  # logic here assumes biallelic call.
    label = labeler.match_to_alt_count(self, candidate, truth_variant, alts)

  See the docs on each individual function to get a better understanding of what
  each function does and the meaning of the return values.
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
    if vcf_reader is None:
      raise ValueError('vcf_reader cannot be None')
    self._vcf_reader = vcf_reader
    self._confident_regions = confident_regions

  def match(self, variant):
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
        candidate learning.genomics.deepvariant.core.genomics.Variant
        variant.

    Returns:
      A tuple of (match_status, truth_variant) where match_status is True if
      we are confident in our truth_variant call or False if not. truth_variant
      is a learning.genomics.deepvariant.core.genomics.Variant object of
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

  def _make_synthetic_hom_ref(self, variant):
    """Creates a version of variant with a hom-ref genotype.

    Args:
      variant: Our
        candidate learning.genomics.deepvariant.core.genomics.Variant
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

    def _usable_truth(truth_variant):
      return (variant.start == truth_variant.start and
              not variantutils.is_filtered(truth_variant))

    region = variantutils.variant_position(variant)
    matches = [m for m in self._vcf_reader.query(region) if _usable_truth(m)]
    if not matches:
      return None
    elif len(matches) > 1:
      logging.warning(
          'Multiple matches detected, keeping first, for variant %s: %s',
          variant, matches)
    return matches[0]

  def match_to_alt_count(self, candidate_variant, truth_variant, alt_alleles):
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
      candidate_variant: Our
        candidate learning.genomics.deepvariant.core.genomics.Variant
        variant.
      truth_variant: Our
        learning.genomics.deepvariant.core.genomics.Variant
        truth variant as returned by match().
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
    if variantutils.is_ref(candidate_variant):
      raise ValueError(
          'candidate_variant must have at least one alternate allele',
          candidate_variant)
    if not variantutils.has_genotypes(truth_variant):
      raise ValueError('truth_variant needs genotypes to be used for labeling',
                       truth_variant)
    if any(alt not in candidate_variant.alternate_bases for alt in alt_alleles):
      raise ValueError('All alt_alleles must be present in variant',
                       alt_alleles, candidate_variant)

    def _simplify_alleles(variant, alleles):
      return [
          variantutils.simplify_alleles(variant.reference_bases, allele)
          for allele in alleles
          if allele != variant.reference_bases
      ]

    simplified_alt_alleles = _simplify_alleles(candidate_variant, alt_alleles)
    return sum(
        true_alt in simplified_alt_alleles
        for true_alt in _simplify_alleles(
            truth_variant, variantutils.genotype_as_alleles(truth_variant)))
