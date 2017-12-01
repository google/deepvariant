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
"""Tests for deepvariant .variant_labeler."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized
import mock

from deepvariant import test_utils
from deepvariant import variant_labeler
from deepvariant.core import ranges
from deepvariant.core import variantutils


def setUpModule():
  test_utils.init()


def mock_vcf_reader(variants):
  """Creates a Mock vcf_reader returning variants from variants.

  This function creates a mock vcf_reader with a query method that searches for
  overlapping variant protos from variants.

  Args:
    variants: list of
      learning.genomics.deepvariant.core.genomics.Variant protos.
      These variants will be used to return results from this vcf_reader.

  Returns:
    A Mock.
  """

  def _fake_query(region):
    return [
        variant for variant in variants
        if ranges.ranges_overlap(variantutils.variant_range(variant), region)
    ]

  reader = mock.MagicMock()
  reader.query.side_effect = _fake_query
  return reader


class VariantLabelerTest(parameterized.TestCase):
  # Confident variants: SNP, deletion, and multi-allelic.
  snp = test_utils.make_variant(start=10, alleles=['A', 'C'], gt=[0, 1])
  deletion = test_utils.make_variant(start=20, alleles=['ACG', 'A'])
  multiallelic = test_utils.make_variant(
      start=30, alleles=['ACT', 'ACTGT', 'A'])
  # Outside our confident regions.
  non_confident = test_utils.make_variant(start=200, alleles=['A', 'C'])
  filtered = test_utils.make_variant(
      start=40, alleles=['A', 'C'], filters='FAILED')
  filtered_match = test_utils.make_variant(
      start=40, alleles=['A', 'C'], gt=[0, 0])

  variants = [snp, deletion, multiallelic, non_confident, filtered]

  def setUp(self):
    self.labeler = variant_labeler.VariantLabeler(
        vcf_reader=mock_vcf_reader(self.variants),
        confident_regions=ranges.RangeSet(
            [ranges.make_range(self.snp.reference_name, 10, 100)]))

  @parameterized.parameters(
      # Simple tests: we get back our matching variants in the confident regions
      (snp, True, snp),
      (deletion, True, deletion),
      (multiallelic, True, multiallelic),

      # Test the behavior outside of our confident regions.
      # We get back non_confident since it matches but we're not confident.
      (non_confident, False, non_confident),
      # No matching variant, so we get a None as well as False.
      (test_utils.make_variant(start=300, alleles=['A', 'C']), False, None),

      # This variant doesn't have any match but we're confident in it.
      (test_utils.make_variant(start=15, alleles=['C', 'A']), True,
       test_utils.make_variant(start=15, alleles=['C', 'A'], gt=[0, 0])),

      # These variant start at our SNP but has a different allele. We are
      # confident and we get back the true snp variant, despite having the
      # different alleles.
      (test_utils.make_variant(start=snp.start, alleles=['A', 'G']), True, snp),
      (test_utils.make_variant(start=snp.start, alleles=['AC', 'C']), True,
       snp),
      (test_utils.make_variant(start=snp.start, alleles=['A', 'CA']), True,
       snp),

      # We don't match filtered variants.
      (filtered, True, filtered_match),
  )
  def test_match(self, candidate, expected_confident, expected_variant):
    actual_confident, actual_variant = self.labeler.match(candidate)
    self.assertEqual(expected_confident, actual_confident)
    self.assertEqual(expected_variant, actual_variant)

  def test_match_selects_variant_by_start(self):
    # Tests that match() selects the variant at the same start even if that
    # variant doesn't have the same alleles at candidate and there's an
    # overlapping with the same alleles.
    overlapping = [
        test_utils.make_variant(start=20, alleles=['CC', 'A']),
        test_utils.make_variant(start=21, alleles=['AAA', 'A']),
        test_utils.make_variant(start=22, alleles=['AA', 'A']),
    ]
    self.labeler = variant_labeler.VariantLabeler(
        vcf_reader=mock_vcf_reader(overlapping))
    candidate = test_utils.make_variant(start=21, alleles=['CC', 'A'])
    self.assertEqual(self.labeler.match(candidate)[1], overlapping[1])

  @parameterized.parameters(
      # Make sure we get the right alt counts for all diploid genotypes.
      (['A', 'C'], ['C'], ['A', 'C'], [0, 0], 0),
      (['A', 'C'], ['C'], ['A', 'C'], [0, 1], 1),
      (['A', 'C'], ['C'], ['A', 'C'], [1, 0], 1),
      (['A', 'C'], ['C'], ['A', 'C'], [1, 1], 2),

      # Basic multi-allelic tests, without having to deal with simplifying
      # alleles as all of the alleles are SNPs. Our candidates have an extra
      # allele, but the true GT is A/C.
      (['A', 'C', 'G'], ['C'], ['A', 'C'], [0, 1], 1),
      (['A', 'C', 'G'], ['C'], ['A', 'C'], [1, 1], 2),

      # When considering A/G our answer should be 0 as we have no copies
      # of the G allele.
      (['A', 'C', 'G'], ['G'], ['A', 'C'], [0, 1], 0),
      (['A', 'C', 'G'], ['G'], ['A', 'C'], [1, 1], 0),

      # We are considering the het-alt configuration here of A vs. C+G. We've
      # got one copy of the C allele so our true genotype is het. If truth is
      # hom-var for the C, though, we again label the composite as hom_var as
      # we have two copies of the C/G alt.
      (['A', 'C', 'G'], ['C', 'G'], ['A', 'C'], [0, 1], 1),
      (['A', 'C', 'G'], ['C', 'G'], ['A', 'C'], [1, 1], 2),

      # Here we have an extra allele in truth, while candidate is bi-allelic.
      # This example 'G' is unused in truth, so we are simply the normal
      # bi-allelic result.
      (['A', 'C'], ['C'], ['A', 'C', 'G'], [0, 0], 0),
      (['A', 'C'], ['C'], ['A', 'C', 'G'], [0, 1], 1),
      (['A', 'C'], ['C'], ['A', 'C', 'G'], [1, 1], 2),

      # We check here that we get the bi-allelic result even when the extra
      # allele is in position 1 not 2.
      (['A', 'G'], ['G'], ['A', 'C', 'G'], [0, 0], 0),
      (['A', 'G'], ['G'], ['A', 'C', 'G'], [0, 2], 1),
      (['A', 'G'], ['G'], ['A', 'C', 'G'], [2, 2], 2),

      # Now for a real het-alt. We've got three alleles in both, and the true
      # genotype is 1/2.
      (['A', 'C', 'G'], ['C'], ['A', 'C', 'G'], [1, 2], 1),
      (['A', 'C', 'G'], ['G'], ['A', 'C', 'G'], [1, 2], 1),
      (['A', 'C', 'G'], ['C', 'G'], ['A', 'C', 'G'], [1, 2], 2),

      # Test ll possible values in candidate against het-alt:
      (['A', 'C', 'G', 'T'], ['C'], ['A', 'C', 'G'], [1, 2], 1),
      (['A', 'C', 'G', 'T'], ['G'], ['A', 'C', 'G'], [1, 2], 1),
      (['A', 'C', 'G', 'T'], ['T'], ['A', 'C', 'G'], [1, 2], 0),
      (['A', 'C', 'G', 'T'], ['C', 'G'], ['A', 'C', 'G'], [1, 2], 2),
      (['A', 'C', 'G', 'T'], ['C', 'T'], ['A', 'C', 'G'], [1, 2], 1),
      (['A', 'C', 'G', 'T'], ['G', 'T'], ['A', 'C', 'G'], [1, 2], 1),

      # Simple start for indel alleles => exact matching works here.
      (['A', 'AC'], ['AC'], ['A', 'AC'], [0, 0], 0),
      (['A', 'AC'], ['AC'], ['A', 'AC'], [0, 1], 1),
      (['A', 'AC'], ['AC'], ['A', 'AC'], [1, 1], 2),

      # We've got a multi-allelic truth, but again exact matching is enough.
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [0, 0], 0),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [0, 1], 1),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [1, 1], 2),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [0, 2], 0),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [1, 2], 1),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [2, 2], 0),

      # This case has an extra allele (A) in truth but the true genotype
      # corresponds to our candidate alleles exactly.
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [0, 2], 1),
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [2, 2], 2),
      # If the true genotype involved just the deletion (A) allele, we don't
      # have that allele in our candidate so we always get 0 copies.
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [0, 1], 0),
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [1, 1], 0),
      # If the truth is het-alt, we can't match the deletion A allele but we do
      # in fact have the A => AC allele as this matches the AC => ACC allele in
      # truth set.
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [1, 2], 1),

      # We have a multi-allelic candidate but a simple bi-allelic truth. Make
      # sure we match correctly. This is an key case, as we should expect that
      # our candidates frequently have extra alleles changing the represention
      # relative to our truth candidates.
      (['ACT', 'A', 'AACT'], ['A'], ['A', 'AA'], [0, 1], 0),
      (['ACT', 'A', 'AACT'], ['A'], ['A', 'AA'], [1, 1], 0),
      (['ACT', 'A', 'AACT'], ['AACT'], ['A', 'AA'], [0, 1], 1),
      (['ACT', 'A', 'AACT'], ['AACT'], ['A', 'AA'], [1, 1], 2),
      (['ACT', 'A', 'AACT'], ['A', 'AACT'], ['A', 'AA'], [0, 1], 1),
      (['ACT', 'A', 'AACT'], ['A', 'AACT'], ['A', 'AA'], [1, 1], 2),

      # The whole complexity: multi-allelic candidate and truth, all with
      # different allele representations.
      # True genotype here is A/AGTGT where ref is AGT [common
      # dinucleotide expansion]. Both candidate and truth have this but each
      # as a different ref so none of the alleles exactly match.
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['A'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], 0),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['AGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['AGTGTGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['A', 'AGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['A', 'AGTGTGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['AGT', 'AGTGTGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], 2),

      # Misc. checks with block substititions.
      (['AT', 'A', 'GC'], ['A'], ['ATT', 'AT', 'A'], [0, 1], 1),
      (['AT', 'A', 'GT'], ['A'], ['A', 'G'], [0, 1], 0),
      (['AT', 'A', 'GT'], ['GT'], ['A', 'G'], [0, 1], 1),
  )
  def test_match_to_genotype_label(self, variant_alleles, alt_alleles,
                                   truth_alleles, truth_gt, expected_n_alts):
    variant = test_utils.make_variant(start=10, alleles=variant_alleles)
    truth_variant = test_utils.make_variant(
        start=10, alleles=truth_alleles, gt=truth_gt)
    self.assertEqual(expected_n_alts,
                     self.labeler.match_to_alt_count(variant, truth_variant,
                                                     alt_alleles))

  def test_match_to_genotype_label_none_truth_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'truth_variant cannot be None'):
      self.labeler.match_to_alt_count(self.snp, None, self.snp.alternate_bases)

  def test_match_to_genotype_label_no_gt_truth_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'truth_variant needs genotypes'):
      self.labeler.match_to_alt_count(self.snp,
                                      test_utils.make_variant(
                                          start=10, alleles=['A', 'C']),
                                      self.snp.alternate_bases)

  def test_match_to_genotype_label_none_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'variant cannot be None'):
      self.labeler.match_to_alt_count(None, self.snp, self.snp.alternate_bases)

  def test_match_to_genotype_label_ref_variant_raises(self):
    with self.assertRaisesRegexp(
        ValueError, 'variant must have at least one alternate allele'):
      self.labeler.match_to_alt_count(
          test_utils.make_variant(start=10, alleles=['A']), self.snp,
          self.snp.alternate_bases)


if __name__ == '__main__':
  absltest.main()
