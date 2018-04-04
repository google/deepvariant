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

from third_party.nucleus.io import vcf
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from deepvariant import testdata
from deepvariant.labeler import variant_labeler


def setUpModule():
  testdata.init()


class DummyVariantLabeler(variant_labeler.VariantLabeler):
  """A dummy VariantLabeler.

  This class provides a label_variants implementation and so allows the base
  class to be instantiated and its methods tested.
  """

  def __init__(self, *pos, **kwargs):
    super(DummyVariantLabeler, self).__init__(*pos, **kwargs)

  def label_variants(self, variants):
    raise NotImplementedError


class VariantLabelerTest(parameterized.TestCase):
  snp = test_utils.make_variant(start=10, alleles=['A', 'C'], gt=[0, 1])

  def test_get_truth_variants(self):
    v1 = test_utils.make_variant(chrom='1', start=10)
    v2 = test_utils.make_variant(chrom='1', start=20)
    v3_filtered = test_utils.make_variant(chrom='1', start=30, filters=['FAIL'])
    v4_del = test_utils.make_variant(chrom='1', start=40, alleles=['AAAA', 'A'])
    v5_non_confident = test_utils.make_variant(chrom='1', start=150)

    variants = [v1, v2, v3_filtered, v4_del, v5_non_confident]
    reader = vcf.InMemoryVcfReader(variants=variants)
    confident_regions = ranges.RangeSet([ranges.make_range('1', 1, 100)])
    labeler = DummyVariantLabeler(
        truth_vcf_reader=reader, confident_regions=confident_regions)

    # Check that we get v1 and v2 specifically when only they are covered by the
    # query.
    self.assertEqual(
        list(labeler._get_truth_variants(ranges.parse_literal('1:1-15'))), [v1])
    self.assertEqual(
        list(labeler._get_truth_variants(ranges.parse_literal('1:15-25'))),
        [v2])

    # We don't include filtered variants.
    self.assertEqual(
        list(labeler._get_truth_variants(ranges.parse_literal('1:25-35'))), [])

    # Check that we get all overlapping variants of our query.
    for del_query in ['1:35-45', '1:42-43', '1:38-42', '1:42-50']:
      self.assertEqual(
          list(labeler._get_truth_variants(ranges.parse_literal(del_query))),
          [v4_del])

    # Checks that a simple query gets all our non-filtered variants.
    self.assertEqual(
        list(labeler._get_truth_variants(ranges.parse_literal('1:1-100'))),
        [v1, v2, v4_del])
    # Even through our query covers v5, it's not confident, so we don't get it.
    self.assertEqual(
        list(labeler._get_truth_variants(ranges.parse_literal('1:1-1000'))),
        [v1, v2, v4_del])

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

  def test_genotype_from_matched_truth_no_call_truth_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'Expected exactly one VariantCal'):
      variant_labeler._genotype_from_matched_truth(self.snp,
                                                   test_utils.make_variant(
                                                       start=10,
                                                       alleles=['A', 'C'],
                                                   ))

  def test_genotype_from_matched_truth_no_gt_truth_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'truth_variant needs genotypes'):
      variant_labeler._genotype_from_matched_truth(self.snp,
                                                   test_utils.make_variant(
                                                       start=10,
                                                       alleles=['A', 'C'],
                                                       gt=[-1, -1],
                                                   ))

  def test_genotype_from_matched_truth_none_variant_raises(self):
    with self.assertRaisesRegexp(ValueError, 'variant cannot be None'):
      variant_labeler._genotype_from_matched_truth(None, self.snp)


if __name__ == '__main__':
  absltest.main()
