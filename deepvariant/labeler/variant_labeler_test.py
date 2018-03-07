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

from deepvariant.util import ranges
from deepvariant.util import variant_utils
from deepvariant import test_utils
from deepvariant.labeler import variant_labeler


def setUpModule():
  test_utils.init()


def mock_vcf_reader(variants):
  """Creates a Mock vcf_reader returning variants from variants.

  This function creates a mock vcf_reader with a query method that searches for
  overlapping variant protos from variants.

  Args:
    variants: list of
      third_party.nucleus.protos.Variant protos.
      These variants will be used to return results from this vcf_reader.

  Returns:
    A Mock.
  """

  def _fake_query(region):
    return [
        variant for variant in variants
        if ranges.ranges_overlap(variant_utils.variant_range(variant), region)
    ]

  reader = mock.MagicMock()
  reader.query.side_effect = _fake_query
  return reader


class PositionalVariantLabelerTest(parameterized.TestCase):
  # Confident variants: SNP, deletion, and multi-allelic.
  snp = test_utils.make_variant(start=10, alleles=['A', 'C'], gt=[0, 1])
  deletion = test_utils.make_variant(start=20, alleles=['ACG', 'A'], gt=[1, 1])
  multiallelic = test_utils.make_variant(
      start=30, alleles=['ACT', 'ACTGT', 'A'], gt=[1, 2])
  # Outside our confident regions.
  non_confident = test_utils.make_variant(
      start=200, alleles=['A', 'C'], gt=[0, 1])
  filtered = test_utils.make_variant(
      start=40, alleles=['A', 'C'], filters='FAILED', gt=[0, 1])
  filtered_match = test_utils.make_variant(
      start=40, alleles=['A', 'C'], gt=[0, 0])

  variants = [snp, deletion, multiallelic, non_confident, filtered]

  def _make_labeler(self, variants, confident_regions):
    return variant_labeler.make_labeler(
        truth_variants_reader=mock_vcf_reader(variants),
        confident_regions=confident_regions,
        ref_reader=None,
        options=None)

  @parameterized.parameters(
      # Simple tests: we get back our matching variants in the confident regions
      dict(candidate=snp, expected_confident=True, expected_truth=snp),
      dict(
          candidate=deletion, expected_confident=True, expected_truth=deletion),
      dict(
          candidate=multiallelic,
          expected_confident=True,
          expected_truth=multiallelic),

      # Test the behavior outside of our confident regions.
      # We get back non_confident since it matches but we're not confident.
      dict(
          candidate=non_confident,
          expected_confident=False,
          expected_truth=non_confident),
      # No matching variant, so we get a None as well as False.
      dict(
          candidate=test_utils.make_variant(start=300, alleles=['A', 'C']),
          expected_confident=False,
          expected_truth=None),

      # This variant doesn't have any match but we're confident in it.
      dict(
          candidate=test_utils.make_variant(start=15, alleles=['C', 'A']),
          expected_confident=True,
          expected_genotype=(0, 0),
          expected_truth=test_utils.make_variant(
              start=15, alleles=['C', 'A'], gt=[0, 0])),

      # These variant start at our SNP but has a different allele. We are
      # confident and we get back the true snp variant, despite having the
      # different alleles. snp has alleles=['A', 'C'] and gt=[0, 1].
      dict(
          candidate=test_utils.make_variant(
              start=snp.start, alleles=['A', 'G']),
          expected_confident=True,
          expected_genotype=(0, 0),
          expected_truth=snp),
      dict(
          candidate=test_utils.make_variant(
              start=snp.start, alleles=['AC', 'C']),
          expected_confident=True,
          expected_genotype=(0, 0),
          expected_truth=snp),
      dict(
          candidate=test_utils.make_variant(
              start=snp.start, alleles=['A', 'CA']),
          expected_confident=True,
          expected_genotype=(0, 0),
          expected_truth=snp),

      # We don't match filtered variants.
      dict(
          candidate=filtered,
          expected_confident=True,
          expected_genotype=(0, 0),
          expected_truth=filtered_match),
  )
  def test_label_variants(self,
                          candidate,
                          expected_confident,
                          expected_truth,
                          expected_genotype=None):
    labeler = self._make_labeler(
        self.variants,
        ranges.RangeSet([ranges.make_range(self.snp.reference_name, 10, 100)]))

    # Call _match so we can compare our expected truth with the actual one.
    is_confident, truth_variant = labeler._match(candidate)
    self.assertEqual(expected_truth, truth_variant)
    self.assertEqual(is_confident, expected_confident)

    # Now call label_variants to exercise the higher-level API.
    if expected_genotype is None and expected_truth is not None:
      expected_genotype = tuple(expected_truth.calls[0].genotype)
    labels = list(labeler.label_variants([candidate]))
    self.assertEqual(len(labels), 1)
    self.assertEqual(candidate, labels[0].variant)
    self.assertEqual(expected_confident, labels[0].is_confident)
    self.assertEqual(expected_genotype, labels[0].genotype)

  def test_match_selects_variant_by_start(self):
    # Tests that match() selects the variant at the same start even if that
    # variant doesn't have the same alleles at candidate and there's an
    # overlapping with the same alleles.
    overlapping = [
        test_utils.make_variant(start=20, alleles=['CC', 'A'], gt=[1, 1]),
        test_utils.make_variant(start=21, alleles=['AAA', 'A'], gt=[0, 1]),
        test_utils.make_variant(start=22, alleles=['AA', 'A'], gt=[1, 1]),
    ]
    candidate = test_utils.make_variant(start=21, alleles=['CC', 'A'])

    labeler = self._make_labeler(overlapping, None)
    is_confident, truth_variant = labeler._match(candidate)
    self.assertEqual(is_confident, True)
    self.assertEqual(truth_variant, overlapping[1])


class VariantLabelerTest(parameterized.TestCase):
  snp = test_utils.make_variant(start=10, alleles=['A', 'C'], gt=[0, 1])

  @parameterized.parameters(
      # Make sure we get the right alt counts for all diploid genotypes.
      (['A', 'C'], ['C'], ['A', 'C'], [0, 0], (0, 0), 0),
      (['A', 'C'], ['C'], ['A', 'C'], [0, 1], (0, 1), 1),
      (['A', 'C'], ['C'], ['A', 'C'], [1, 0], (0, 1), 1),
      (['A', 'C'], ['C'], ['A', 'C'], [1, 1], (1, 1), 2),

      # Make sure get back a zero alt count for a reference variant.
      (['A'], [], ['A'], [0, 0], (0, 0), 0),

      # Basic multi-allelic tests, without having to deal with simplifying
      # alleles as all of the alleles are SNPs. Our candidates have an extra
      # allele, but the true GT is A/C.
      (['A', 'C', 'G'], ['C'], ['A', 'C'], [0, 1], (0, 1), 1),
      (['A', 'C', 'G'], ['C'], ['A', 'C'], [1, 1], (1, 1), 2),

      # When considering A/G our answer should be 0 as we have no copies
      # of the G allele.
      (['A', 'C', 'G'], ['G'], ['A', 'C'], [0, 1], (0, 1), 0),
      (['A', 'C', 'G'], ['G'], ['A', 'C'], [1, 1], (1, 1), 0),

      # We are considering the het-alt configuration here of A vs. C+G. We've
      # got one copy of the C allele so our true genotype is het. If truth is
      # hom-var for the C, though, we again label the composite as hom_var as
      # we have two copies of the C/G alt.
      (['A', 'C', 'G'], ['C', 'G'], ['A', 'C'], [0, 1], (0, 1), 1),
      (['A', 'C', 'G'], ['C', 'G'], ['A', 'C'], [1, 1], (1, 1), 2),

      # Here we have an extra allele in truth, while candidate is bi-allelic.
      # This example 'G' is unused in truth, so we are simply the normal
      # bi-allelic result.
      (['A', 'C'], ['C'], ['A', 'C', 'G'], [0, 0], (0, 0), 0),
      (['A', 'C'], ['C'], ['A', 'C', 'G'], [0, 1], (0, 1), 1),
      (['A', 'C'], ['C'], ['A', 'C', 'G'], [1, 1], (1, 1), 2),

      # We check here that we get the bi-allelic result even when the extra
      # allele is in position 1 not 2.
      (['A', 'G'], ['G'], ['A', 'C', 'G'], [0, 0], (0, 0), 0),
      (['A', 'G'], ['G'], ['A', 'C', 'G'], [0, 2], (0, 1), 1),
      (['A', 'G'], ['G'], ['A', 'C', 'G'], [2, 2], (1, 1), 2),

      # Now for a real het-alt. We've got three alleles in both, and the true
      # genotype is 1/2.
      (['A', 'C', 'G'], ['C'], ['A', 'C', 'G'], [1, 2], (1, 2), 1),
      (['A', 'C', 'G'], ['G'], ['A', 'C', 'G'], [1, 2], (1, 2), 1),
      (['A', 'C', 'G'], ['C', 'G'], ['A', 'C', 'G'], [1, 2], (1, 2), 2),

      # Test all possible values in candidate against het-alt:
      (['A', 'C', 'G', 'T'], ['C'], ['A', 'C', 'G'], [1, 2], (1, 2), 1),
      (['A', 'C', 'G', 'T'], ['G'], ['A', 'C', 'G'], [1, 2], (1, 2), 1),
      (['A', 'C', 'G', 'T'], ['T'], ['A', 'C', 'G'], [1, 2], (1, 2), 0),
      (['A', 'C', 'G', 'T'], ['C', 'G'], ['A', 'C', 'G'], [1, 2], (1, 2), 2),
      (['A', 'C', 'G', 'T'], ['C', 'T'], ['A', 'C', 'G'], [1, 2], (1, 2), 1),
      (['A', 'C', 'G', 'T'], ['G', 'T'], ['A', 'C', 'G'], [1, 2], (1, 2), 1),

      # Simple start for indel alleles => exact matching works here.
      (['A', 'AC'], ['AC'], ['A', 'AC'], [0, 0], (0, 0), 0),
      (['A', 'AC'], ['AC'], ['A', 'AC'], [0, 1], (0, 1), 1),
      (['A', 'AC'], ['AC'], ['A', 'AC'], [1, 1], (1, 1), 2),

      # We've got a multi-allelic truth, but again exact matching is enough.
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [0, 0], (0, 0), 0),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [0, 1], (0, 1), 1),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [1, 1], (1, 1), 2),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [0, 2], (0, 0), 0),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [1, 2], (0, 1), 1),
      (['A', 'AC'], ['AC'], ['A', 'AC', 'ACC'], [2, 2], (0, 0), 0),

      # This case has an extra allele (A) in truth but the true genotype
      # corresponds to our candidate alleles exactly.
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [0, 2], (0, 1), 1),
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [2, 2], (1, 1), 2),
      # If the true genotype involved just the deletion (A) allele, we don't
      # have that allele in our candidate so we always get 0 copies.
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [0, 1], (0, 0), 0),
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [1, 1], (0, 0), 0),
      # If the truth is het-alt, we can't match the deletion A allele but we do
      # in fact have the A => AC allele as this matches the AC => ACC allele in
      # truth set.
      (['A', 'AC'], ['AC'], ['AC', 'A', 'ACC'], [1, 2], (0, 1), 1),

      # We have a multi-allelic candidate but a simple bi-allelic truth. Make
      # sure we match correctly. This is an key case, as we should expect that
      # our candidates frequently have extra alleles changing the represention
      # relative to our truth candidates.
      (['ACT', 'A', 'AACT'], ['A'], ['A', 'AA'], [0, 1], (0, 2), 0),
      (['ACT', 'A', 'AACT'], ['A'], ['A', 'AA'], [1, 1], (2, 2), 0),
      (['ACT', 'A', 'AACT'], ['AACT'], ['A', 'AA'], [0, 1], (0, 2), 1),
      (['ACT', 'A', 'AACT'], ['AACT'], ['A', 'AA'], [1, 1], (2, 2), 2),
      (['ACT', 'A', 'AACT'], ['A', 'AACT'], ['A', 'AA'], [0, 1], (0, 2), 1),
      (['ACT', 'A', 'AACT'], ['A', 'AACT'], ['A', 'AA'], [1, 1], (2, 2), 2),

      # The whole complexity: multi-allelic candidate and truth, all with
      # different allele representations.
      # True genotype here is A/AGTGT where ref is AGT [common
      # dinucleotide expansion]. Both candidate and truth have this but each
      # as a different ref so none of the alleles exactly match.
      #
      # Truth     : AGT   => A [1] + AGTGT [2]
      # Candidate : AGTGT => AGT [2] + AGTGTGT [3]
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['A'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], (2, 3), 0),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['AGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], (2, 3), 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['AGTGTGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], (2, 3), 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['A', 'AGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], (2, 3), 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['A', 'AGTGTGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], (2, 3), 1),
      (['AGTGT', 'A', 'AGT', 'AGTGTGT'], ['AGT', 'AGTGTGT'],
       ['AGT', 'A', 'AGTGT', 'AGTGTGT'], [1, 2], (2, 3), 2),

      # Misc. checks with block substititions.
      (['AT', 'A', 'GC'], ['A'], ['ATT', 'AT', 'A'], [0, 1], (0, 1), 1),
      (['AT', 'A', 'GT'], ['A'], ['A', 'G'], [0, 1], (0, 2), 0),
      (['AT', 'A', 'GT'], ['GT'], ['A', 'G'], [0, 1], (0, 2), 1),
  )
  def test_genotype_from_matched_truth(self, variant_alleles, alt_alleles,
                                       truth_alleles, truth_gt,
                                       expected_genotype, expected_label):
    variant = test_utils.make_variant(start=10, alleles=variant_alleles)
    truth_variant = test_utils.make_variant(
        start=10, alleles=truth_alleles, gt=truth_gt)
    self.assertEqual(expected_genotype,
                     variant_labeler._genotype_from_matched_truth(
                         variant, truth_variant))
    labeled = variant_labeler.VariantLabel(
        is_confident=True,
        variant=variant,
        genotype=expected_genotype)
    indices = [variant_alleles.index(alt) - 1 for alt in alt_alleles]
    self.assertEqual(labeled.label_for_alt_alleles(indices), expected_label)

  def test_genotype_from_matched_truth_none_truth_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'truth_variant cannot be None'):
      variant_labeler._genotype_from_matched_truth(self.snp, None)

  def test_genotype_from_matched_truth_no_gt_truth_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'truth_variant needs genotypes'):
      variant_labeler._genotype_from_matched_truth(self.snp,
                                                   test_utils.make_variant(
                                                       start=10,
                                                       alleles=['A', 'C']))

  def test_genotype_from_matched_truth_none_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'variant cannot be None'):
      variant_labeler._genotype_from_matched_truth(None, self.snp)


if __name__ == '__main__':
  absltest.main()
