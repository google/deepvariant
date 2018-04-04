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
from deepvariant.labeler import positional_labeler


def setUpModule():
  testdata.init()


class PositionalVariantLabelerTest(parameterized.TestCase):
  # Confident variants: SNP, deletion, and multi-allelic.
  snp = test_utils.make_variant(start=10, alleles=['A', 'C'], gt=[0, 1])
  deletion = test_utils.make_variant(start=20, alleles=['ACG', 'A'], gt=[1, 1])
  multiallelic = test_utils.make_variant(
      start=30, alleles=['ACT', 'ACTGT', 'A'], gt=[1, 2])
  # Outside our confident regions.
  non_confident = test_utils.make_variant(
      start=200, alleles=['A', 'C'], gt=[0, 1])
  filtered = test_utils.make_variant(start=40, filters='FAILED', gt=[0, 1])

  variants = [snp, deletion, multiallelic, non_confident, filtered]

  def _make_labeler(self, variants, confident_regions):
    return positional_labeler.PositionalVariantLabeler(
        truth_vcf_reader=vcf.InMemoryVcfReader(variants),
        confident_regions=confident_regions)

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
      # If we provide a variant outside the confident regions (non_confident) we
      # don't get back any expected_truth variants.
      dict(
          candidate=non_confident,
          expected_confident=False,
          expected_truth=None),
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
      # Checks that we don't match against the filtered truth variant in our
      # database. This means that we return not the filtered variant but one
      # with a (0, 0) genotype.
      dict(
          candidate=test_utils.make_variant(start=filtered.start),
          expected_confident=True,
          expected_genotype=(0, 0),
          expected_truth=test_utils.make_variant(
              start=filtered.start, gt=(0, 0))),
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

    labeler = self._make_labeler(
        overlapping,
        ranges.RangeSet(
            [ranges.make_range(overlapping[0].reference_name, 0, 100)]))
    is_confident, truth_variant = labeler._match(candidate)
    self.assertEqual(is_confident, True)
    self.assertEqual(truth_variant, overlapping[1])


if __name__ == '__main__':
  absltest.main()
