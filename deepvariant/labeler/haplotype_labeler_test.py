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
"""Tests for deepvariant.haplotype_labeler."""

import itertools
from unittest import mock

from absl.testing import absltest
from absl.testing import parameterized
from third_party.nucleus.io import fasta
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils

from deepvariant.labeler import haplotype_labeler
from deepvariant.protos import deepvariant_pb2


def _test_variant(start=10, alleles=('A', 'C'), gt=None):
  variant = variants_pb2.Variant(
      reference_name='20',
      start=start,
      end=start + len(alleles[0]),
      reference_bases=alleles[0],
      alternate_bases=alleles[1:],
  )

  if gt:
    variant.calls.add(genotype=gt)

  return variant


def _variants_from_grouped_positions(grouped_positions):
  groups = [
      [_test_variant(start=s) for s in starts] for starts in grouped_positions
  ]
  variants = [v for subgroup in groups for v in subgroup]
  return variants, [(g, []) for g in groups]


def _make_labeler(
    truths=None, confident_regions=None, ref_reader=None, **kwargs
):
  if ref_reader is None:
    ref_reader = mock.MagicMock()

  if confident_regions is None:
    # Use the reference of the truth variants if possible, otherwise just use
    # a placeholder placeholder value for the contig name and make the confident
    # region a giant span.
    contig = truths[0].reference_name if truths else 'placeholder'
    confident_regions = ranges.RangeSet(
        [ranges.make_range(contig, 0, 1000000000)]
    )

  return haplotype_labeler.HaplotypeLabeler(
      truth_vcf_reader=vcf.InMemoryVcfReader(truths or []),
      ref_reader=ref_reader,
      confident_regions=confident_regions,
      **kwargs,
  )


class HaplotypeLabelerClassUnitTest(parameterized.TestCase):

  @parameterized.parameters(
      # If there are no variants, we don't get any groups.
      dict(grouped_positions=[]),
      # A single variant at 10 gets put in a single group.
      dict(grouped_positions=[[10]]),
      # A two variants (at 10 and 11) get grouped because 1 < the max_dist=10.
      dict(grouped_positions=[[10, 11]]),
      # A two variants at 10 and 50 get grouped separately because their
      # distance > max_dist=10.
      dict(grouped_positions=[[10], [50]]),
      # Check the behavior right around max_dist.
      dict(grouped_positions=[[10, 19]]),
      dict(grouped_positions=[[10, 20]]),
      dict(grouped_positions=[[10], [21]]),
      # A few complex examples with multiple variants getting grouped.
      dict(grouped_positions=[[10, 15], [45, 50], [65]]),
      dict(grouped_positions=[[10, 12, 15], [45, 49, 50], [65, 70]]),
      # The distance calculation compares not first variant in the group but the
      # closest one, so we can have a group with collective distance > max_dist
      # as long as the distance between successive variants <= max_dist
      dict(grouped_positions=[[10, 20, 30, 40]]),
  )
  def test_group_variants_examples(self, grouped_positions):
    variants, groups = _variants_from_grouped_positions(grouped_positions)
    self.assertEqual(
        groups,
        haplotype_labeler.group_variants(variants, [], max_separation=10),
    )

  @parameterized.parameters(
      dict(separation=s, max_separation=d)
      for d in range(5)
      for s in range(d + 1)
  )
  def test_group_variants_respects_max_dist(self, separation, max_separation):
    self.assertLessEqual(separation, max_separation)
    variants = [
        _test_variant(start=10),
        _test_variant(start=10 + separation),
    ]
    # Because separation <= max_separation, all variants should be in a
    # single group.
    self.assertEqual(
        [(variants, [])],
        haplotype_labeler.group_variants(
            variants, [], max_separation=max_separation
        ),
    )

  @parameterized.parameters(range(1, 10))
  def test_group_variants_works_with_any_number_of_variants(self, n_variants):
    variants = [_test_variant(start=10 + i) for i in range(n_variants)]
    self.assertEqual(
        [(variants, [])],
        haplotype_labeler.group_variants(
            variants, [], max_group_size=n_variants, max_separation=n_variants
        ),
    )

  @parameterized.parameters(
      dict(
          grouped_positions=[[10, 11, 12, 13]],
          max_group_size=4,
          force_group_within_bp=-1,
      ),
      dict(
          grouped_positions=[[10, 11, 12], [13]],
          max_group_size=3,
          force_group_within_bp=-1,
      ),
      dict(
          grouped_positions=[[10, 11], [12, 13]],
          max_group_size=2,
          force_group_within_bp=-1,
      ),
      dict(
          grouped_positions=[[10], [11], [12], [13]],
          max_group_size=1,
          force_group_within_bp=-1,
      ),
  )
  def test_group_variants_group_size_works(
      self, grouped_positions, max_group_size, force_group_within_bp
  ):
    variants, groups = _variants_from_grouped_positions(grouped_positions)
    self.assertEqual(
        groups,
        haplotype_labeler.group_variants(
            variants,
            [],
            max_group_size=max_group_size,
            force_group_within_bp=force_group_within_bp,
        ),
    )

  @parameterized.parameters(
      # A few basic tests of functionality to start:
      # We group a single variant without associated truth variants.
      dict(
          variant_positions=[10],
          truth_positions=[],
          expected_positions=[([10], [])],
          force_group_within_bp=-1,
      ),
      # We group an isolated truth variant in isolation without any candidates.
      dict(
          variant_positions=[],
          truth_positions=[10],
          expected_positions=[([], [10])],
          force_group_within_bp=-1,
      ),
      # A single candidate + truth at the same position are grouped together.
      dict(
          variant_positions=[10],
          truth_positions=[10],
          expected_positions=[([10], [10])],
          force_group_within_bp=-1,
      ),
      # We respect our distance between variants across positional boundaries.
      dict(
          variant_positions=[10],
          truth_positions=[11],
          expected_positions=[([10], [11])],
          force_group_within_bp=-1,
      ),
      # Another test of grouping across candidates and truth.
      dict(
          variant_positions=[9, 11],
          truth_positions=[10, 12],
          expected_positions=[([9, 11], [10, 12])],
          force_group_within_bp=-1,
      ),
      # These tests actually result in broken up groups, with isolated truth
      # and candidates as well as grouped ones.
      dict(
          variant_positions=[1, 9, 11, 20, 25, 50],
          truth_positions=[10, 12, 23, 45, 100],
          expected_positions=[
              ([1], []),
              ([9, 11], [10, 12]),
              ([20, 25], [23]),
              ([50], [45]),
              ([], [100]),
          ],
          force_group_within_bp=-1,
      ),
      # Now some tests to exercise the max group size with both candidates and
      # truth variants. We vary the max group size to make sure the grouping
      # algorithm splits correctly.
      dict(
          variant_positions=[1, 2, 3, 4, 5],
          truth_positions=[1, 2, 3, 4, 5],
          expected_positions=[
              ([1, 2], [1, 2]),
              ([3, 4], [3, 4]),
              ([5], [5]),
          ],
          force_group_within_bp=-1,
          max_group_size=2,
      ),
      dict(
          variant_positions=[1, 2, 3, 4, 5],
          truth_positions=[1, 2, 3, 4, 5],
          expected_positions=[
              ([1, 2, 3], [1, 2, 3]),
              ([4, 5], [4, 5]),
          ],
          force_group_within_bp=-1,
          max_group_size=3,
      ),
      # Force group variants within 1bp distance
      dict(
          variant_positions=[10, 11],
          truth_positions=[10, 11],
          expected_positions=[([10, 11], [10, 11])],
          force_group_within_bp=1,
          max_group_size=1,
      ),
      # Force group variants within 1bp distance with 8 variants with
      # max_group_size of 1
      dict(
          variant_positions=[1, 2, 3, 4, 5, 6, 7, 8],
          truth_positions=[1, 2, 3, 4, 5, 6, 7, 8],
          expected_positions=[(
              [1, 2, 3, 4, 5, 6, 7, 8],
              [1, 2, 3, 4, 5, 6, 7, 8],
          )],
          force_group_within_bp=1,
          max_group_size=1,
      ),
  )
  def test_group_variants_with_truth(
      self,
      variant_positions,
      truth_positions,
      expected_positions,
      force_group_within_bp,
      max_group_size=10,
  ):
    variants = [_test_variant(start=pos) for pos in variant_positions]
    truths = [_test_variant(start=pos) for pos in truth_positions]
    actual = haplotype_labeler.group_variants(
        variants,
        truths,
        max_group_size=max_group_size,
        max_separation=5,
        force_group_within_bp=force_group_within_bp,
    )
    self.assertEqual(
        [
            ([v.start for v in variant_group], [v.start for v in truth_group])
            for variant_group, truth_group in actual
        ],
        expected_positions,
    )

  @parameterized.parameters(
      dict(truth_position=7, expected_together=False),
      dict(truth_position=8, expected_together=False),
      dict(truth_position=9, expected_together=True),
      dict(truth_position=10, expected_together=True),
      dict(truth_position=11, expected_together=True),
      dict(truth_position=12, expected_together=True),
      dict(truth_position=13, expected_together=True),
      dict(truth_position=14, expected_together=True),
      dict(truth_position=15, expected_together=True),
      dict(truth_position=16, expected_together=False),
      dict(truth_position=17, expected_together=False),
  )
  def test_group_variants_respects_end(self, truth_position, expected_together):
    # Deletion spans from 10-15 as it's a deletion.
    deletion = _test_variant(start=10, alleles=('AAAAA', 'A'))
    truth = _test_variant(start=truth_position)

    # The actual result returned by group_variants is a list of tuples
    # containing the grouped candidates and truth. The order they appear depends
    # on the truth_position, since our deletion starts at 10.
    if expected_together:
      expected = [([deletion], [truth])]
    elif truth_position < 10:
      expected = [([], [truth]), ([deletion], [])]
    else:
      expected = [([deletion], []), ([], [truth])]

    self.assertEqual(
        haplotype_labeler.group_variants([deletion], [truth], max_separation=1),
        expected,
    )

  @parameterized.parameters(
      # The 8 variants each have #GT: 10, 6, 6, 3, 15, 15, 6, 6.
      # Combining any of them will be 10*6*6*3*15*15*6*6 = 8748000.
      # Here we set a max_gt_options_product larger than that. So they are all
      # grouped together.
      dict(
          max_gt_options_product=10000000,
          force_group_within_bp=-1,
          expected_group_split=[8],
      ),
      # The 8 variants each have #GT: 10, 6, 6, 3, 15, 15, 6, 6.
      # Because _MAX_GT_OPTIONS_PRODUCT is 100000, the split is:
      # 5 variants in one group: 10*6*6*3*15 = 16200
      # 3 variants in the next group: 15*6*6 = 540
      dict(
          max_gt_options_product=haplotype_labeler._MAX_GT_OPTIONS_PRODUCT,
          force_group_within_bp=-1,
          expected_group_split=[5, 3],
      ),
      # The 8 variants each have #GT: 10, 6, 6, 3, 15, 15, 6, 6.
      # Combining any of them will be more than max_gt_options_product=10.
      # So each individual is split into their own group.
      dict(
          max_gt_options_product=10,
          force_group_within_bp=-1,
          expected_group_split=[1] * 8,
      ),
      # The 8 variants each have #GT: 10, 6, 6, 3, 15, 15, 6, 6.
      # Combining any of them will be more than max_gt_options_product=10.
      # However, because two of the candidates have the same variant.end,
      # When I set force_group_within_bp to 0, they are grouped together.
      # So two of the variants are grouped together, and the rest are grouped
      # separately.
      dict(
          max_gt_options_product=10,
          force_group_within_bp=0,
          expected_group_split=[1, 2, 1, 1, 1, 1, 1],
      ),
  )
  def test_group_variants_exceeds_max_gt_options_product(
      self, max_gt_options_product, force_group_within_bp, expected_group_split
  ):
    # This is an extreme case from a DeepTrio exome training run. internal.
    candidates = [
        _test_variant(36206921, alleles=('A', 'C', 'G', 'T')),
        _test_variant(36206925, alleles=('GA', 'G', 'TA')),
        _test_variant(36206926, alleles=('A', 'C', 'T')),
        _test_variant(36206927, alleles=('T', 'G')),
        _test_variant(36206930, alleles=('G', 'C', 'GC', 'GT', 'T')),
        _test_variant(36206931, alleles=('G', 'A', 'C', 'GT', 'T')),
        _test_variant(36206934, alleles=('A', 'C', 'T')),
        _test_variant(36206937, alleles=('G', 'C', 'T')),
    ]
    grouped = haplotype_labeler.group_variants(
        candidates,
        [],
        max_group_size=haplotype_labeler._MAX_GROUP_SIZE,
        max_separation=haplotype_labeler._MAX_SEPARATION_WITHIN_VARIANT_GROUP,
        max_gt_options_product=max_gt_options_product,
        force_group_within_bp=force_group_within_bp,
    )
    # Check the end of each variant:
    self.assertEqual(
        [c.end for c in candidates],
        [
            36206922,
            36206927,
            36206927,
            36206928,
            36206931,
            36206932,
            36206935,
            36206938,
        ],
    )

    self.assertLen(grouped, len(expected_group_split))
    self.assertEqual([len(g[0]) for g in grouped], expected_group_split)

  @parameterized.parameters(
      # Check a simple case of two SNPs.
      dict(
          candidates=[
              _test_variant(start=10, alleles=['A', 'C']),
              _test_variant(start=20, alleles=['A', 'C']),
          ],
          truths=[],
          expected_start=9,
          expected_end=21,
          bufsize=0,
      ),
      # Check that we respect the deletion's span in the last candidate.
      dict(
          candidates=[
              _test_variant(start=10, alleles=['A', 'C']),
              _test_variant(start=20, alleles=['AAA', 'C']),
          ],
          truths=[],
          expected_start=9,
          expected_end=23,
          bufsize=0,
      ),
      # Check that we respect handle truth variants, as the interval is entirely
      # determined by truth variants here.
      dict(
          candidates=[
              _test_variant(start=15, alleles=['A', 'C']),
          ],
          truths=[
              _test_variant(start=10, alleles=['A', 'C']),
              _test_variant(start=20, alleles=['AAA', 'C']),
          ],
          expected_start=9,
          expected_end=23,
          bufsize=0,
      ),
      # Check that bufsize is respected.
      dict(
          candidates=[
              _test_variant(start=10, alleles=['A', 'C']),
              _test_variant(start=20, alleles=['AAA', 'C']),
          ],
          truths=[],
          expected_start=9,
          expected_end=23 + 10,  # 10 is the bufsize.
          bufsize=10,
      ),
  )
  def test_make_labeler_ref(
      self, candidates, truths, expected_start, expected_end, bufsize
  ):
    expected_bases = 'A' * (expected_end - expected_start)

    labeler = _make_labeler()
    labeler._ref_reader.query.return_value = expected_bases
    labeler._ref_reader.contig.return_value = reference_pb2.ContigInfo(
        name='20', n_bases=50
    )

    labeler_ref = labeler.make_labeler_ref(candidates, truths, bufsize=bufsize)

    labeler._ref_reader.query.assert_called_once_with(
        ranges.make_range('20', expected_start, expected_end)
    )
    self.assertEqual(labeler_ref.start, expected_start)
    self.assertEqual(labeler_ref.end, expected_end)
    self.assertEqual(
        labeler_ref.bases(expected_start, expected_end), expected_bases
    )

  # Check that we don't issue a query with a bad start if the variant is at
  # position 0 on the genome. See internal.
  def test_make_labeler_ref_handles_variant_at_pos_zero(self):
    labeler = _make_labeler(
        ref_reader=fasta.InMemoryFastaReader([('20', 0, 'ACGT')])
    )
    labeler_ref = labeler.make_labeler_ref(
        candidates=[
            _test_variant(start=0, alleles=['A', 'C']),
        ],
        true_variants=[],
        bufsize=2,
    )
    expected_start, expected_end = 0, 3
    self.assertEqual(labeler_ref.start, expected_start)
    self.assertEqual(labeler_ref.end, expected_end)
    self.assertEqual(labeler_ref.bases(expected_start, expected_end), 'ACG')

  # Check that we don't issue a query with a bad end if the variant is at
  # the last position on the genome. See internal.
  def test_make_labeler_ref_handles_variant_at_end_of_chrom(self):
    labeler = _make_labeler(
        ref_reader=fasta.InMemoryFastaReader([('20', 0, 'ACGT')])
    )
    labeler_ref = labeler.make_labeler_ref(
        candidates=[
            _test_variant(start=3, alleles=['T', 'A']),
        ],
        true_variants=[],
        bufsize=2,
    )
    expected_start, expected_end = 2, 4
    self.assertEqual(labeler_ref.start, expected_start)
    self.assertEqual(labeler_ref.end, expected_end)
    self.assertEqual(labeler_ref.bases(expected_start, expected_end), 'GT')

  def test_label_variants(self):
    variants = [
        _test_variant(start=10),
        _test_variant(start=11),
        _test_variant(start=20),
        _test_variant(start=30),
        _test_variant(start=31),
        _test_variant(start=42),
    ]
    truths = [
        _test_variant(start=10, gt=(0, 1)),
        _test_variant(start=30, gt=(1, 1)),
        _test_variant(start=42, gt=(0, 0)),
    ]

    labeler = _make_labeler(
        truths=truths,
        max_separation=5,
        ref_reader=fasta.InMemoryFastaReader([('20', 0, 'A' * 100)]),
    )
    region = ranges.make_range('20', 1, 50)
    result = list(labeler.label_variants(variants, region))

    expected_genotypes_by_pos = {
        10: (0, 1),
        11: (0, 0),
        20: (0, 0),
        30: (1, 1),
        31: (0, 0),
        42: (0, 0),
    }
    self.assertEqual(len(result), len(variants))
    for variant, label in zip(variants, result):
      self.assertEqual(variant.start, label.variant.start)
      self.assertTrue(label.is_confident)
      self.assertEqual(
          tuple(label.variant.calls[0].genotype),
          expected_genotypes_by_pos[variant.start],
          'Bad genotype for ' + str(label.variant),
      )

  def test_label_variants_bug1(self):
    # Test for a bug encountered in make_examples.
    #
    # variants: candidates [0]
    # variants: truth [1]
    #    20:6299587:C->T gt=(1, 1)
    # Top-level exception: ('Failed to assign labels for variants', [])
    labeler = _make_labeler(
        truths=[_test_variant(6299586, alleles=('C', 'T'), gt=(1, 1))],
        ref_reader=fasta.InMemoryFastaReader(
            [('20', 6299585, 'TCCTGCTTTCTCTTGTGGGCAT')]
        ),
    )
    region = ranges.make_range('20', 6299000, 6299999)
    self.assertIsNotNone(labeler.label_variants([], region))

  @parameterized.parameters(
      # A single TP bi-allelic variant.
      dict(
          candidates=[
              _test_variant(start=2, alleles=('A', 'C')),
          ],
          truths=[
              _test_variant(start=2, alleles=('A', 'C'), gt=(0, 1)),
          ],
          n_truth_variant_sites=1,
          n_truth_variant_alleles=1,
          n_candidate_variant_sites=1,
          n_candidate_variant_alleles=1,
          n_true_positive_sites=1,
          n_true_positive_alleles=1,
          n_false_negative_sites=0,
          n_false_negative_alleles=0,
          n_false_positive_sites=0,
          n_false_positive_alleles=0,
          n_inexact_position_matches=0,
          n_exact_position_matches=1,
          n_exact_position_and_allele_matches=1,
          n_exact_position_and_allele_and_genotype_matches=1,
          n_truth_multiallelics_sites_with_missed_alleles=0,
      ),
      # A single TP tri-allelic variant.
      dict(
          candidates=[
              _test_variant(start=2, alleles=('A', 'C', 'G')),
          ],
          truths=[
              _test_variant(start=2, alleles=('A', 'C', 'G'), gt=(1, 2)),
          ],
          n_truth_variant_sites=1,
          n_truth_variant_alleles=2,
          n_candidate_variant_sites=1,
          n_candidate_variant_alleles=2,
          n_true_positive_sites=1,
          n_true_positive_alleles=2,
          n_false_negative_sites=0,
          n_false_negative_alleles=0,
          n_false_positive_sites=0,
          n_false_positive_alleles=0,
          n_inexact_position_matches=0,
          n_exact_position_matches=1,
          n_exact_position_and_allele_matches=1,
          n_exact_position_and_allele_and_genotype_matches=1,
          n_truth_multiallelics_sites_with_missed_alleles=0,
      ),
      # Here we have an extra alt in our candidate and one tri-allelic truth.
      # Because of this we have a FP allele and our exact matching counts are
      # different than the above example.
      dict(
          candidates=[
              _test_variant(start=2, alleles=('A', 'C', 'G', 'T')),
          ],
          truths=[
              _test_variant(start=2, alleles=('A', 'C', 'G'), gt=(1, 2)),
          ],
          n_truth_variant_sites=1,
          n_truth_variant_alleles=2,
          n_candidate_variant_sites=1,
          n_candidate_variant_alleles=3,
          n_true_positive_sites=1,
          n_true_positive_alleles=2,
          n_false_negative_sites=0,
          n_false_negative_alleles=0,
          n_false_positive_sites=0,
          n_false_positive_alleles=1,
          n_inexact_position_matches=0,
          n_exact_position_matches=1,
          n_exact_position_and_allele_matches=0,
          n_exact_position_and_allele_and_genotype_matches=0,
          n_truth_multiallelics_sites_with_missed_alleles=0,
      ),
      # A single FP variant in candidates with 3 alt alleles.
      dict(
          candidates=[
              _test_variant(start=2, alleles=('A', 'C', 'G')),
          ],
          truths=[],
          n_truth_variant_sites=0,
          n_truth_variant_alleles=0,
          n_candidate_variant_sites=1,
          n_candidate_variant_alleles=2,
          n_true_positive_sites=0,
          n_true_positive_alleles=0,
          n_false_negative_sites=0,
          n_false_negative_alleles=0,
          n_false_positive_sites=1,
          n_false_positive_alleles=2,
          n_inexact_position_matches=0,
          n_exact_position_matches=0,
          n_exact_position_and_allele_matches=0,
          n_exact_position_and_allele_and_genotype_matches=0,
          n_truth_multiallelics_sites_with_missed_alleles=0,
      ),
      # A single FN variant in truth with 4 alt alleles.
      dict(
          candidates=[],
          truths=[
              _test_variant(start=2, alleles=('A', 'C', 'G', 'T'), gt=(1, 2)),
          ],
          n_truth_variant_sites=1,
          # We only have 2 non-ref alleles in truth, so this is 2 not 3.
          n_truth_variant_alleles=2,
          n_candidate_variant_sites=0,
          n_candidate_variant_alleles=0,
          n_true_positive_sites=0,
          n_true_positive_alleles=0,
          n_false_negative_sites=1,
          n_false_negative_alleles=2,
          n_false_positive_sites=0,
          n_false_positive_alleles=0,
          n_inexact_position_matches=0,
          n_exact_position_matches=0,
          n_exact_position_and_allele_matches=0,
          n_exact_position_and_allele_and_genotype_matches=0,
          n_truth_multiallelics_sites_with_missed_alleles=1,
      ),
  )
  def test_metrics(self, candidates, truths, **kwargs):
    self.assertMetricsEqual(candidates=candidates, truths=truths, **kwargs)

  @parameterized.parameters(range(1, 5))
  def test_metrics_multiple_variants(self, max_separation):
    # This is parameterized over the max_separation so we can test that the
    # metrics are properly calculated no matter the grouping. The candidates and
    # truth variants below should give the same metrics regardless of grouping.
    self.assertMetricsEqual(
        candidates=[
            # Example one from our narrative in LabelerMetrics proto.
            _test_variant(start=2, alleles=('A', 'C')),
            # Example two from our narrative in LabelerMetrics proto.
            _test_variant(start=3, alleles=('A', 'C', 'T')),
            # A genuine false positive => no corresponding truth variant.
            _test_variant(start=4, alleles=('A', 'G', 'C')),
        ],
        truths=[
            _test_variant(start=2, alleles=('A', 'C'), gt=(0, 1)),
            _test_variant(start=3, alleles=('A', 'C', 'G'), gt=(1, 2)),
            # A genuine false negative => no corresponding variant in truth.
            _test_variant(start=5, alleles=('A', 'C'), gt=(0, 1)),
        ],
        max_separation=max_separation,
        n_truth_variant_sites=3,
        n_truth_variant_alleles=4,
        n_candidate_variant_sites=3,
        n_candidate_variant_alleles=5,
        n_true_positive_sites=2,
        n_true_positive_alleles=2,
        n_false_negative_sites=1,
        n_false_negative_alleles=2,
        n_false_positive_sites=1,
        n_false_positive_alleles=3,
        n_inexact_position_matches=0,
        n_exact_position_matches=2,
        n_exact_position_and_allele_matches=1,
        n_exact_position_and_allele_and_genotype_matches=1,
        n_truth_multiallelics_sites_with_missed_alleles=1,
    )

  def test_metrics_inexact_matches(self):
    # ref looks like AACTG. Truth is just a single SNP turning the C into a G.
    # Candidates do the same but via an insertion + deletion. This test ensures
    # that the metrics work even in the case where we have different
    # representations for the same haplotype.
    self.assertMetricsEqual(
        candidates=[
            _test_variant(start=1, alleles=('A', 'AG')),
            _test_variant(start=2, alleles=('CT', 'T')),
        ],
        truths=[
            _test_variant(start=2, alleles=('C', 'G'), gt=(0, 1)),
        ],
        ref_prefix='AACTG',
        n_truth_variant_sites=1,
        n_truth_variant_alleles=1,
        n_candidate_variant_sites=2,
        n_candidate_variant_alleles=2,
        n_true_positive_sites=1,
        n_true_positive_alleles=1,
        n_false_negative_sites=0,
        n_false_negative_alleles=0,
        n_false_positive_sites=0,
        n_false_positive_alleles=0,
        n_inexact_position_matches=1,
        n_exact_position_matches=1,
        n_exact_position_and_allele_matches=0,
        n_exact_position_and_allele_and_genotype_matches=0,
        n_truth_multiallelics_sites_with_missed_alleles=0,
    )

  def assertMetricsEqual(
      self,
      candidates,
      truths,
      ref_prefix='',
      max_separation=10,
      **metric_kwargs,
  ):
    labeler = _make_labeler(
        truths=truths,
        max_separation=max_separation,
        confident_regions=ranges.RangeSet([ranges.make_range('20', 0, 1000)]),
        ref_reader=fasta.InMemoryFastaReader(
            [('20', 0, ref_prefix + 'A' * 50)]
        ),
    )
    region = ranges.make_range('20', 1, 10)
    _ = list(labeler.label_variants(candidates, region))
    self.assertEqual(
        labeler.metrics, deepvariant_pb2.LabelingMetrics(**metric_kwargs)
    )

  def test_metrics_respects_confident_regions(self):
    # The confident region is 2-4, so we should only count the variant starting
    # at 3.
    candidates = [
        _test_variant(start=1),
        _test_variant(start=3),
        _test_variant(start=5),
    ]
    labeler = _make_labeler(
        truths=[],
        confident_regions=ranges.RangeSet(
            [ranges.make_range(candidates[0].reference_name, 2, 4)]
        ),
        ref_reader=fasta.InMemoryFastaReader([('20', 0, 'A' * 50)]),
    )
    _ = list(labeler.label_variants(candidates, ranges.make_range('20', 1, 10)))
    self.assertEqual(labeler.metrics.n_candidate_variant_sites, 1)
    self.assertEqual(labeler.metrics.n_candidate_variant_alleles, 1)
    self.assertEqual(labeler.metrics.n_non_confident_candidate_variant_sites, 2)
    self.assertEqual(labeler.metrics.n_false_positive_sites, 1)
    self.assertEqual(labeler.metrics.n_false_positive_alleles, 1)


class HaplotypeMatchTests(parameterized.TestCase):

  def setUp(self):
    self.haplotypes = ['AC', 'GT']
    self.variants = [
        _test_variant(42, ['A', 'G']),
        _test_variant(43, ['G', 'A']),
        _test_variant(44, ['C', 'T']),
    ]
    self.truths = [
        _test_variant(42, ['A', 'G'], [0, 1]),
        _test_variant(44, ['C', 'T'], [0, 1]),
        _test_variant(45, ['G', 'A'], [1, 1]),
    ]
    self.candidate_genotypes = [(0, 1), (0, 0), (0, 1)]
    self.truth_genotypes = [(0, 1), (0, 1), (0, 0)]
    self.match = haplotype_labeler.HaplotypeMatch(
        self.haplotypes,
        self.variants,
        self.candidate_genotypes,
        self.truths,
        self.truth_genotypes,
    )

  def test_fields_are_expected(self):
    self.assertEqual(self.match.haplotypes, self.haplotypes)
    self.assertEqual(self.match.candidates, self.variants)
    self.assertEqual(self.match.truths, self.truths)
    self.assertEqual(self.match.candidate_genotypes, self.candidate_genotypes)
    self.assertEqual(self.match.truth_genotypes, self.truth_genotypes)

    # Computed fields.
    self.assertEqual(
        self.match.original_truth_genotypes,
        haplotype_labeler._variant_genotypes(self.truths),
    )
    self.assertEqual(self.match.n_false_positives, 1)
    self.assertEqual(self.match.n_false_negatives, 2)
    self.assertEqual(self.match.n_true_positives, 2)
    self.assertEqual(self.match.match_metrics, (2, 1, 2))

  def test_str(self):
    self.assertIn('HaplotypeMatch(', str(self.match))

  def test_candidates_with_assigned_genotypes(self):
    self.assertEqual(
        self.match.candidates_with_assigned_genotypes(),
        [
            _test_variant(42, ['A', 'G'], [0, 1]),
            _test_variant(43, ['G', 'A'], [0, 0]),
            _test_variant(44, ['C', 'T'], [0, 1]),
        ],
    )
    # Assert that we have no genotypes in self.variants to check that
    # candidates_with_assigned_genotypes isn't modifying our variants.
    for v in self.match.candidates:
      self.assertFalse(v.calls, 'Variant genotypes modified')

  @parameterized.parameters(
      # All genotypes match, so we have no FNs and no FPs.
      dict(
          vgenotypes=[(0, 1), (0, 1)],
          matched_tgenotypes=[(0, 1), (0, 1)],
          tgenotypes=[(0, 1), (0, 1)],
          expected_fns=0,
          expected_fps=0,
          expected_tps=2,
      ),
      #
      # Here are a few cases where we false negatives.
      #
      dict(
          vgenotypes=[(0, 1), (0, 1)],
          # True het is 0, 0 in second true variant.
          matched_tgenotypes=[(0, 1), (0, 0)],
          tgenotypes=[(0, 1), (0, 1)],
          expected_fns=1,
          expected_fps=0,
          expected_tps=2,
      ),
      dict(
          vgenotypes=[(0, 1), (0, 1)],
          # True hom-alt is het in second true variant.
          matched_tgenotypes=[(0, 1), (0, 1)],
          tgenotypes=[(0, 1), (1, 1)],
          expected_fns=1,
          expected_fps=0,
          expected_tps=2,
      ),
      dict(
          vgenotypes=[(0, 1), (0, 1)],
          # True hom-alts are het in first and second true variants.
          matched_tgenotypes=[(0, 1), (0, 1)],
          tgenotypes=[(1, 1), (1, 1)],
          expected_fns=2,
          expected_fps=0,
          expected_tps=2,
      ),
      dict(
          vgenotypes=[(0, 1), (0, 1)],
          # True hom-alts are hom-ref in first and second true variants.
          matched_tgenotypes=[(0, 0), (0, 0)],
          tgenotypes=[(1, 1), (1, 1)],
          expected_fns=4,
          expected_fps=0,
          expected_tps=2,
      ),
      #
      # Here are a few cases where we false positives.
      #
      dict(
          # The first variant genotype is hom-ref so we have one FP.
          vgenotypes=[(0, 0), (0, 1)],
          matched_tgenotypes=[(0, 1), (0, 1)],
          tgenotypes=[(0, 1), (0, 1)],
          expected_fns=0,
          expected_fps=1,
          expected_tps=1,
      ),
      dict(
          # The second variant genotype is hom-ref so we have one FP.
          vgenotypes=[(0, 1), (0, 0)],
          matched_tgenotypes=[(0, 1), (0, 1)],
          tgenotypes=[(0, 1), (0, 1)],
          expected_fns=0,
          expected_fps=1,
          expected_tps=1,
      ),
      dict(
          # Both variant genotypes are hom-ref, so we have two FPs.
          vgenotypes=[(0, 0), (0, 0)],
          matched_tgenotypes=[(0, 1), (0, 1)],
          tgenotypes=[(0, 1), (0, 1)],
          expected_fns=0,
          expected_fps=2,
          expected_tps=0,
      ),
  )
  def test_fns_fps(
      self,
      vgenotypes,
      matched_tgenotypes,
      tgenotypes,
      expected_fns,
      expected_fps,
      expected_tps,
  ):
    match = haplotype_labeler.HaplotypeMatch(
        haplotypes=self.haplotypes,
        candidates=[
            _test_variant(42, ['A', 'G']),
            _test_variant(43, ['G', 'A']),
        ],
        candidate_genotypes=vgenotypes,
        truths=[
            _test_variant(42, ['A', 'G'], tgenotypes[0]),
            _test_variant(43, ['G', 'A'], tgenotypes[1]),
        ],
        truth_genotypes=matched_tgenotypes,
    )
    self.assertEqual(match.n_false_negatives, expected_fns)
    self.assertEqual(match.n_false_positives, expected_fps)
    self.assertEqual(match.n_true_positives, expected_tps)
    self.assertEqual(match.n_true_positives + match.n_false_positives, 2)
    self.assertEqual(
        match.match_metrics, (expected_fns, expected_fps, expected_tps)
    )


class LabelExamplesTest(parameterized.TestCase):
  # Many of these tests are cases from our labeler analysis doc:
  # https://docs.google.com/document/d/1V89IIT0YM3P0gH_tQb-ahodf8Jvnz0alXEnjCf6JVNo

  def assertGetsCorrectLabels(
      self,
      candidates,
      true_variants,
      ref,
      expected_genotypes,
      start=None,
      end=None,
  ):
    start = start or ref.start
    end = end or ref.end
    labeling = haplotype_labeler.find_best_matching_haplotypes(
        candidates, true_variants, ref
    )
    self.assertIsNotNone(labeling)

    # Check that the genotypes of our labeled variants are the ones we expect.
    labeled_variants = labeling.candidates_with_assigned_genotypes()
    self.assertEqual(
        haplotype_labeler._variant_genotypes(labeled_variants),
        [tuple(x) for x in expected_genotypes],
    )

  @parameterized.parameters(
      dict(genotype=[0, 0], expected=[(0, 0)]),
      dict(genotype=[0, 1], expected=[(0, 0), (0, 1)]),
      dict(genotype=[1, 1], expected=[(0, 0), (0, 1), (1, 1)]),
      dict(genotype=[0, 2], expected=[(0, 0), (0, 2)]),
      dict(genotype=[2, 2], expected=[(0, 0), (0, 2), (2, 2)]),
      dict(genotype=[1, 2], expected=[(0, 0), (0, 1), (0, 2), (1, 2)]),
  )
  def test_with_false_negative_genotypes(self, genotype, expected):
    self.assertEqual(
        haplotype_labeler.with_false_negative_genotypes(genotype), expected
    )

  @parameterized.parameters(
      dict(
          prefix_haplotypes_list=[{'a'}],
          haplotypes={'A'},
          expected=[{'aA'}],
      ),
      dict(
          prefix_haplotypes_list=[{'a', 'b'}],
          haplotypes={'A'},
          expected=[{'aA', 'bA'}],
      ),
      dict(
          prefix_haplotypes_list=[{'a'}],
          haplotypes={'A', 'B'},
          expected=[{'aA', 'aB'}],
      ),
      dict(
          prefix_haplotypes_list=[{'a', 'b'}],
          haplotypes={'A', 'B'},
          expected=[{'aA', 'bB'}, {'aB', 'bA'}],
      ),
      dict(
          prefix_haplotypes_list=[{'a', 'b'}, {'c'}, {'d', 'e'}],
          haplotypes={'A'},
          expected=[{'aA', 'bA'}, {'cA'}, {'dA', 'eA'}],
      ),
      dict(
          prefix_haplotypes_list=[{'a', 'b'}, {'c'}, {'d', 'e'}],
          haplotypes={'A', 'B'},
          expected=[
              {'aA', 'bB'},
              {'aB', 'bA'},
              {'cA', 'cB'},
              {'dA', 'eB'},
              {'dB', 'eA'},
          ],
      ),
  )
  def test_extend_haplotypes(
      self, prefix_haplotypes_list, haplotypes, expected
  ):
    self.assertCountEqual(
        expected,
        list(
            haplotype_labeler.extend_haplotypes(
                prefix_haplotypes_list, haplotypes
            )
        ),
    )

  def test_extend_haplotypes_raises_on_empty_prefix_list(self):
    with self.assertRaisesRegex(ValueError, 'prefix_haplotypes_list.*empty'):
      list(haplotype_labeler.extend_haplotypes([], {'A'}))

  def test_extend_haplotypes_raises_on_empty_haplotypes(self):
    with self.assertRaisesRegex(ValueError, 'haplotypes'):
      list(haplotype_labeler.extend_haplotypes([{'A'}], set()))

  def test_extend_haplotypes_raises_on_too_many_haplotypes(self):
    with self.assertRaisesRegex(ValueError, 'haplotypes'):
      list(haplotype_labeler.extend_haplotypes([{'A'}], {'a', 'b', 'c'}))

  # The reference sequence is xAAAAAy.
  @parameterized.parameters(
      # Test a SNP at a few positions.
      dict(
          variants=[
              _test_variant(1, alleles=('A', 'C')),
          ],
          allele_indices_and_expected={
              (0,): 'xAAAAAy',
              (1,): 'xCAAAAy',
          },
      ),
      dict(
          variants=[
              _test_variant(3, alleles=('A', 'C')),
          ],
          allele_indices_and_expected={
              (0,): 'xAAAAAy',
              (1,): 'xAACAAy',
          },
      ),
      # Test an insertion at a few positions.
      dict(
          variants=[
              _test_variant(1, alleles=('A', 'CC')),
          ],
          allele_indices_and_expected={
              (0,): 'xAAAAAy',
              (1,): 'xCCAAAAy',
          },
      ),
      dict(
          variants=[
              _test_variant(3, alleles=('A', 'CC')),
          ],
          allele_indices_and_expected={
              (0,): 'xAAAAAy',
              (1,): 'xAACCAAy',
          },
      ),
      # Test a deletion at a few positions.
      dict(
          variants=[
              _test_variant(1, alleles=('AAA', 'A')),
          ],
          allele_indices_and_expected={
              (0,): 'xAAAAAy',
              (1,): 'xAAAy',
          },
      ),
      dict(
          variants=[
              _test_variant(3, alleles=('AAA', 'A')),
          ],
          allele_indices_and_expected={
              (0,): 'xAAAAAy',
              (1,): 'xAAAy',
          },
      ),
      # A complete example with multiple variants.
      dict(
          variants=[
              _test_variant(1, alleles=('A', 'C')),
              _test_variant(2, alleles=('A', 'CC')),
              _test_variant(4, alleles=('AA', 'A')),
          ],
          allele_indices_and_expected={
              (0, 0, 0): 'xAAAAAy',
              (0, 0, 1): 'xAAAAy',
              (0, 1, 0): 'xACCAAAy',
              (0, 1, 1): 'xACCAAy',
              (1, 0, 0): 'xCAAAAy',
              (1, 0, 1): 'xCAAAy',
              (1, 1, 0): 'xCCCAAAy',
              (1, 1, 1): 'xCCCAAy',
          },
      ),
  )
  def test_build_haplotype(self, variants, allele_indices_and_expected):
    refseq = 'xAAAAAy'
    for allele_indices, expected in allele_indices_and_expected.items():
      self.assertEqual(
          expected,
          haplotype_labeler.build_haplotype(
              variants,
              allele_indices,
              ref=haplotype_labeler.ReferenceRegion(refseq, 0),
              ref_start=0,
              ref_end=len(refseq),
          ),
      )

  @parameterized.parameters(
      # All possible genotypes for a simple tri-allelic case.
      (
          dict(
              variants=[
                  _test_variant(11, ['TG', 'A', 'TGC'], gt),
              ],
              ref=haplotype_labeler.ReferenceRegion('TG', 11),
              expected_frags=expected,
              expected_next_pos=13,
          )
          for gt, expected in {
              # Simple bi-allelic configurations:
              (0, 0): {(0,): 'TG'},
              (0, 1): {(0,): 'TG', (1,): 'A'},
              (1, 0): {(0,): 'TG', (1,): 'A'},
              (1, 1): {(1,): 'A'},
              # Multi-allelic configurations:
              (0, 2): {(0,): 'TG', (2,): 'TGC'},
              (1, 2): {(1,): 'A', (2,): 'TGC'},
              (2, 2): {(2,): 'TGC'},
          }.items()
      ),
  )
  def test_phased_genotypes_to_haplotypes_single_variant(
      self, variants, ref, expected_frags, expected_next_pos
  ):
    variants_and_genotypes = [
        haplotype_labeler.VariantAndGenotypes(v, tuple(v.calls[0].genotype))
        for v in variants
    ]
    frags, next_pos = haplotype_labeler.phased_genotypes_to_haplotypes(
        variants_and_genotypes, ref.start, ref
    )
    self.assertEqual(frags, expected_frags)
    self.assertEqual(next_pos, expected_next_pos)

  @parameterized.parameters(
      ('G', 'A'),
      ('GG', 'A'),
      ('GGG', 'A'),
      ('GGGG', 'A'),
      ('A', 'G'),
      ('A', 'GG'),
      ('A', 'GGG'),
      ('A', 'GGGG'),
  )
  def test_phased_genotypes_to_haplotypes_next_pos_is_correct(self, ref, alt):
    # Check that the next_pos calculation is working.
    pos = 10
    for gt in [(0, 0), (0, 1), (1, 1)]:
      _, next_pos = haplotype_labeler.phased_genotypes_to_haplotypes(
          [
              haplotype_labeler.VariantAndGenotypes(
                  _test_variant(pos, [ref, alt]), gt
              )
          ],
          start=pos,
          ref=haplotype_labeler.ReferenceRegion(ref, pos),
      )
      self.assertEqual(next_pos, pos + len(ref))

  @parameterized.parameters(
      # A single deletion overlapping a SNP:
      # ref: xTG
      # v1:   A-
      # v2:    C
      dict(
          variants=[
              _test_variant(11, ['TG', 'A'], (0, 1)),
              _test_variant(12, ['G', 'C'], (0, 1)),
          ],
          ref=haplotype_labeler.ReferenceRegion('xTG', 10),
          expected_frags={
              (0, 0): 'xTG',  # haplotype 0|0.
              (0, 1): 'xTC',  # haplotype 0|1.
              (1, 0): 'xA',  # haplotype 1|0.
              (1, 1): None,  # haplotype 1|1 => invalid.
          },
          expected_next_pos=13,
      ),
      # Deletion overlapping two downstream events (SNP and insertion):
      # ref: xTGC
      # v1:   A--
      # v2:    C
      # v3:     TTT
      dict(
          variants=[
              _test_variant(11, ['TGC', 'A'], (0, 1)),
              _test_variant(12, ['G', 'C'], (0, 1)),
              _test_variant(13, ['C', 'TTT'], (0, 1)),
          ],
          ref=haplotype_labeler.ReferenceRegion('xTGC', 10),
          expected_frags={
              (0, 0, 0): 'xTGC',  # haplotype 0|0|0.
              (0, 0, 1): 'xTGTTT',  # haplotype 0|0|1.
              (0, 1, 0): 'xTCC',  # haplotype 0|1|0.
              (0, 1, 1): 'xTCTTT',  # haplotype 0|1|1.
              (1, 0, 0): 'xA',  # haplotype 1|0|0.
              (1, 0, 1): None,  # haplotype 1|0|1 => invalid.
              (1, 1, 0): None,  # haplotype 1|1|0 => invalid.
              (1, 1, 1): None,  # haplotype 1|1|1 => invalid.
          },
          expected_next_pos=14,
      ),
      # Two incompatible deletions to check that the end extension is working:
      # pos: 01234
      # ref: xTGCA
      # v1:   T-
      # v2:    G-
      dict(
          variants=[
              _test_variant(11, ['TG', 'T'], (0, 1)),
              _test_variant(12, ['GC', 'G'], (0, 1)),
          ],
          ref=haplotype_labeler.ReferenceRegion('xTGCA', 10),
          expected_frags={
              (0, 0): 'xTGC',  # haplotype 0|0.
              (0, 1): 'xTG',  # haplotype 0|1.
              (1, 0): 'xTC',  # haplotype 1|0.
              (1, 1): None,  # haplotype 1|1 => invalid.
          },
          expected_next_pos=14,
      ),
      # Multiple overlapping deletions with complex incompatibilities:
      # ref: xTGCGA
      # v1:   A--
      # v2:    G---  [conflicts with v1]
      # v3:     C-   [conflicts with v1 and v2]
      # v4:      G-  [conflicts with v2 and v3, ok with v1]
      # v5:       C  [conflicts with v2 and v4, ok with v1, v3]
      dict(
          variants=[
              _test_variant(11, ['TGC', 'A'], (0, 1)),
              _test_variant(12, ['GCGA', 'G'], (0, 1)),
              _test_variant(13, ['CG', 'C'], (0, 1)),
              _test_variant(14, ['GA', 'G'], (0, 1)),
              _test_variant(15, ['A', 'C'], (0, 1)),
          ],
          ref=haplotype_labeler.ReferenceRegion('xTGCGA', 10),
          expected_frags={
              (0, 0, 0, 0, 0): 'xTGCGA',  # haplotype 0|0|0|0|0.
              (0, 0, 0, 0, 1): 'xTGCGC',  # haplotype 0|0|0|0|1.
              (0, 0, 0, 1, 0): 'xTGCG',  # haplotype 0|0|0|1|0.
              (0, 0, 0, 1, 1): None,  # haplotype 0|0|0|1|1.
              (0, 0, 1, 0, 0): 'xTGCA',  # haplotype 0|0|1|0|0.
              (0, 0, 1, 0, 1): 'xTGCC',  # haplotype 0|0|1|0|1.
              (0, 0, 1, 1, 0): None,  # haplotype 0|0|1|1|0.
              (0, 0, 1, 1, 1): None,  # haplotype 0|0|1|1|1.
              (0, 1, 0, 0, 0): 'xTG',  # haplotype 0|1|0|0|0.
              (0, 1, 0, 0, 1): None,  # haplotype 0|1|0|0|1.
              (0, 1, 0, 1, 0): None,  # haplotype 0|1|0|1|0.
              (0, 1, 0, 1, 1): None,  # haplotype 0|1|0|1|1.
              (0, 1, 1, 0, 0): None,  # haplotype 0|1|1|0|0.
              (0, 1, 1, 0, 1): None,  # haplotype 0|1|1|0|1.
              (0, 1, 1, 1, 0): None,  # haplotype 0|1|1|1|0.
              (0, 1, 1, 1, 1): None,  # haplotype 0|1|1|1|1.
              (1, 0, 0, 0, 0): 'xAGA',  # haplotype 1|0|0|0|0.
              (1, 0, 0, 0, 1): 'xAGC',  # haplotype 1|0|0|0|1.
              (1, 0, 0, 1, 0): 'xAG',  # haplotype 1|0|0|1|0.
              (1, 0, 0, 1, 1): None,  # haplotype 1|0|0|1|1.
              (1, 0, 1, 0, 0): None,  # haplotype 1|0|1|0|0.
              (1, 0, 1, 0, 1): None,  # haplotype 1|0|1|0|1.
              (1, 0, 1, 1, 0): None,  # haplotype 1|0|1|1|0.
              (1, 0, 1, 1, 1): None,  # haplotype 1|0|1|1|1.
              (1, 1, 0, 0, 0): None,  # haplotype 1|1|0|0|0.
              (1, 1, 0, 0, 1): None,  # haplotype 1|1|0|0|1.
              (1, 1, 0, 1, 0): None,  # haplotype 1|1|0|1|0.
              (1, 1, 0, 1, 1): None,  # haplotype 1|1|0|1|1.
              (1, 1, 1, 0, 0): None,  # haplotype 1|1|1|0|0.
              (1, 1, 1, 0, 1): None,  # haplotype 1|1|1|0|1.
              (1, 1, 1, 1, 0): None,  # haplotype 1|1|1|1|0.
              (1, 1, 1, 1, 1): None,  # haplotype 1|1|1|1|1.
          },
          expected_next_pos=16,
      ),
  )
  def test_phased_genotypes_to_haplotypes_overlapping(
      self, variants, ref, expected_frags, expected_next_pos
  ):
    variants_and_genotypes = [
        haplotype_labeler.VariantAndGenotypes(v, tuple(v.calls[0].genotype))
        for v in variants
    ]
    frags, next_pos = haplotype_labeler.phased_genotypes_to_haplotypes(
        variants_and_genotypes, ref.start, ref
    )
    self.assertEqual(
        frags, {k: v for k, v in expected_frags.items() if v is not None}
    )
    self.assertEqual(next_pos, expected_next_pos)

  @parameterized.parameters(
      # Check that simple bi-allelic matching works for all possible possible
      # genotypes and a variety of types of alleles.
      (
          dict(
              candidate_alleles=alleles,
              truth_alleles=alleles,
              truth_genotype=gt,
              # Returns [0, 1] even if truth is [1, 0], so sort the genotypes for
              # the expected value.
              expected_genotype=sorted(gt),
          )
          for gt in [[0, 1], [1, 0], [1, 1]]
          for alleles in [['A', 'C'], ['ACC', 'A'], ['A', 'ATG'], ['AC', 'GT']]
      ),
  )
  def test_single_variants(
      self, candidate_alleles, truth_alleles, truth_genotype, expected_genotype
  ):
    candidate = _test_variant(42, candidate_alleles)
    truth = _test_variant(42, truth_alleles, truth_genotype)
    ref_allele = sorted([candidate_alleles[0], truth_alleles[0]], key=len)[0]
    self.assertGetsCorrectLabels(
        candidates=[candidate],
        true_variants=[truth],
        ref=haplotype_labeler.ReferenceRegion('x' + ref_allele + 'y', 41),
        expected_genotypes=[expected_genotype],
    )

  @parameterized.parameters(
      dict(true_genotype=(0, 1)),
      dict(true_genotype=(1, 1)),
  )
  def test_no_candidates_only_truth_variants(self, true_genotype):
    labeling = haplotype_labeler.find_best_matching_haplotypes(
        candidates=[],
        truths=[_test_variant(42, gt=true_genotype)],
        ref=haplotype_labeler.ReferenceRegion('xAy', 41),
    )
    self.assertIsNotNone(labeling)

    # Since we don't have any candidates, our relabeled variants should be [].
    self.assertEqual(labeling.candidates_with_assigned_genotypes(), [])

  @parameterized.parameters(
      dict(empty='variants'),
      dict(empty='truth'),
  )
  def test_no_variants_or_truth_is_fast(self, empty):
    # This test will time out if we aren't able to efficiently handle the case
    # where we have a lot of candidate or truth variants but none of the other.
    many_variants = [_test_variant(i, gt=(0, 1)) for i in range(10, 50)]
    if empty == 'truth':
      variants, truth = many_variants, []
    else:
      variants, truth = [], many_variants

    labeling = haplotype_labeler.find_best_matching_haplotypes(
        candidates=variants,
        truths=truth,
        ref=haplotype_labeler.ReferenceRegion('A' * 50, 10),
    )

    # Since we don't have any truth variants, all of the variants should get a
    # (0, 0) [i.e., hom-ref] genotype assigned.
    self.assertEqual(
        haplotype_labeler._variant_genotypes(
            labeling.candidates_with_assigned_genotypes()
        ),
        [(0, 0)] * len(variants),
    )

  def test_genotype_options_for_variants_truth_enum(self):
    # Check all configurations for the TRUTH enumeration:
    enum_type = haplotype_labeler.EnumerationType.TRUTH

    # Bi-allelic cases.
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, gt=(0, 0))], enum_type
        ),
        [[(0, 0)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, gt=(0, 1))], enum_type
        ),
        [[(0, 0), (0, 1)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, gt=(1, 1))], enum_type
        ),
        [[(0, 0), (0, 1), (1, 1)]],
    )

    # Multi-allelic cases.
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'), gt=(0, 0))], enum_type
        ),
        [[(0, 0)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'), gt=(0, 1))], enum_type
        ),
        [[(0, 0), (0, 1)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'), gt=(1, 1))], enum_type
        ),
        [[(0, 0), (0, 1), (1, 1)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'), gt=(0, 2))], enum_type
        ),
        [[(0, 0), (0, 2)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'), gt=(2, 2))], enum_type
        ),
        [[(0, 0), (0, 2), (2, 2)]],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'), gt=(1, 2))], enum_type
        ),
        [[(0, 0), (0, 1), (0, 2), (1, 2)]],
    )

  def test_genotype_options_for_variants_candidates_enum(self):
    # Check all configurations for the CANDIDATES enumeration:
    enum_type = haplotype_labeler.EnumerationType.CANDIDATES
    # Note we don't need to provide a genotype for the candidate enumeration.
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1)], enum_type
        ),
        [{(0, 0), (0, 1), (1, 1)}],
    )
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(
            [_test_variant(1, alleles=('A', 'C', 'G'))], enum_type
        ),
        [{(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)}],
    )

  def test_genotype_options_for_variants_only_hom_ref(self):
    # Check all configurations for the ONLY_HOM_REF enumeration:
    enum_type = haplotype_labeler.EnumerationType.ONLY_HOM_REF
    variants = [
        _test_variant(1),
        _test_variant(1, gt=(0, 1)),
        _test_variant(1, gt=(1, 1)),
        _test_variant(1, alleles=('A', 'C', 'G')),
        _test_variant(1, alleles=('A', 'C', 'G'), gt=(1, 2)),
    ]
    self.assertEqual(
        haplotype_labeler.genotype_options_for_variants(variants, enum_type),
        [{(0, 0)}] * len(variants),
    )

  def test_genotype_options_for_variants_no_variants(self):
    # All enumeration types return [] if not provided with any variants.
    for enum_type in haplotype_labeler.EnumerationType:
      self.assertEqual(
          haplotype_labeler.genotype_options_for_variants([], enum_type), []
      )

  @parameterized.parameters(
      dict(
          candidate_alleles=['A', 'C'],
          truth_alleles=['A', 'G', 'C'],
          truth_genotypes_and_expected={
              (0, 2): (0, 1),  # A/C => 0/C
              (1, 2): (0, 1),  # G/C => 0/C
              (1, 1): (0, 0),  # G/G => 0/0
              (2, 2): (1, 1),  # C/C => C/C
          },
      ),
      dict(
          candidate_alleles=['A', 'C'],
          truth_alleles=['A', 'C', 'G'],
          truth_genotypes_and_expected={
              (0, 1): (0, 1),  # A/C => 0/C
              (2, 1): (0, 1),  # G/C => 0/C
              (2, 2): (0, 0),  # G/G => 0/0
              (1, 1): (1, 1),  # C/C => C/C
          },
      ),
      dict(
          candidate_alleles=['A', 'TT', 'TTT'],
          truth_alleles=['A', 'C', 'G'],
          truth_genotypes_and_expected={
              (i, j): (0, 0) for i, j in itertools.combinations([0, 1, 2], 2)
          },
      ),
      # Here the candidate is also multi-allelic
      dict(
          candidate_alleles=['A', 'G', 'C'],
          truth_alleles=['A', 'C', 'G'],
          truth_genotypes_and_expected={
              (0, 1): (0, 2),
              (0, 2): (0, 1),
              (1, 1): (2, 2),
              (1, 2): (1, 2),
              (2, 2): (1, 1),
          },
      ),
  )
  def test_multi_allelic(
      self, candidate_alleles, truth_alleles, truth_genotypes_and_expected
  ):
    candidate = _test_variant(42, candidate_alleles)
    for true_gt, expected_gt in truth_genotypes_and_expected.items():
      truth = _test_variant(42, truth_alleles, true_gt)
      ref_allele = sorted([candidate_alleles[0], truth_alleles[0]], key=len)[0]
      self.assertGetsCorrectLabels(
          candidates=[candidate],
          true_variants=[truth],
          ref=haplotype_labeler.ReferenceRegion('x' + ref_allele + 'y', 41),
          expected_genotypes=[expected_gt],
      )

  def test_false_variants_get_homref_genotype(self):
    ref = haplotype_labeler.ReferenceRegion('xACGTAy', 10)
    v1 = _test_variant(11, ['A', 'T'], [0, 1])
    v2 = _test_variant(13, ['G', 'GG'], [1, 1])
    all_fps = [
        _test_variant(12, ['C', 'G'], [0, 0]),
        _test_variant(14, ['T', 'A'], [0, 0]),
        _test_variant(15, ['A', 'AA'], [0, 0]),
    ]
    for n_fps in range(1, len(all_fps) + 1):
      for fps in itertools.combinations(all_fps, n_fps):
        candidates = variant_utils.sorted_variants([v1, v2] + list(fps))
        self.assertGetsCorrectLabels(
            candidates=candidates,
            true_variants=[v1, v2],
            ref=ref,
            expected_genotypes=haplotype_labeler._variant_genotypes(candidates),
        )

  def test_false_negatives(self):
    ref = haplotype_labeler.ReferenceRegion('xACGTAy', 10)
    v1 = _test_variant(11, ['A', 'T'], [0, 1])
    v2 = _test_variant(13, ['G', 'GG'], [1, 1])
    all_fns = [
        _test_variant(12, ['C', 'G'], [0, 1]),
        _test_variant(14, ['T', 'A', 'G'], [1, 2]),
        _test_variant(15, ['A', 'AA'], [1, 1]),
    ]
    for n_fns in [1]:
      for fns in itertools.combinations(all_fns, n_fns):
        candidates = [v1, v2]
        self.assertGetsCorrectLabels(
            candidates=candidates,
            true_variants=variant_utils.sorted_variants([v1, v2] + list(fns)),
            ref=ref,
            expected_genotypes=haplotype_labeler._variant_genotypes(candidates),
        )

  # example 20:3528533 and 20:3528534
  def test_example1(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(3528531, ['ATAG', 'A']),
            _test_variant(3528537, ['A', 'ATT']),
        ],
        true_variants=[
            _test_variant(3528533, ['A', 'T'], [1, 1]),
            _test_variant(3528534, ['G', 'A'], [1, 1]),
            _test_variant(3528536, ['TA', 'T'], [1, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('xATAGTTATC', 3528530),
        expected_genotypes=[
            [1, 1],
            [1, 1],
        ],
    )

  # example 20:4030071
  def test_example2(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(4030067, ['TC', 'T']),
            _test_variant(4030072, ['C', 'G']),
        ],
        true_variants=[
            _test_variant(4030071, ['CC', 'G'], [1, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('xTCCCCCA', 4030066),
        expected_genotypes=[
            [1, 1],
            [1, 1],
        ],
    )

  # example 20:4568152
  def test_example3(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(4568151, ['AC', 'A']),
            _test_variant(4568154, ['TG', 'T']),
            _test_variant(4568156, ['G', 'T']),
            _test_variant(4568157, ['A', 'ATACCCTTT']),
        ],
        true_variants=[
            _test_variant(4568152, ['C', 'A'], [1, 1]),
            _test_variant(4568153, ['A', 'T'], [1, 1]),
            _test_variant(4568155, ['G', 'A'], [1, 1]),
            _test_variant(4568156, ['G', 'T'], [1, 1]),
            _test_variant(4568157, ['A', 'ACCCTTT'], [1, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('xACATGGATGGA', 4568150),
        expected_genotypes=[
            [1, 1],
            [1, 1],
            [1, 1],
            [1, 1],
        ],
    )

  # example 20:1689636, 20:1689639, 20:1689640, 20:1689641
  def test_example4(self):
    # CTGTAAACAGAA [phased alts] + CGTGAATGAAA [phased ref]
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(1689633, ['C', 'CT']),
            _test_variant(1689635, ['TG', 'T']),
            _test_variant(1689638, ['ATG', 'A']),
            _test_variant(1689641, ['A', 'ACAG']),
        ],
        true_variants=[
            _test_variant(1689633, ['C', 'CT'], [1, 0]),
            _test_variant(1689636, ['G', 'A'], [1, 0]),
            _test_variant(1689639, ['T', 'C'], [1, 0]),
            _test_variant(1689640, ['G', 'A'], [1, 0]),
            _test_variant(1689641, ['A', 'G'], [1, 0]),
        ],
        ref=haplotype_labeler.ReferenceRegion('xCGTGAATGAAA', 1689632),
        expected_genotypes=[
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
        ],
    )

  # 20:2401511
  def test_example5(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(2401510, ['ATGT', 'A']),
            _test_variant(2401515, ['C', 'T']),
        ],
        true_variants=[
            _test_variant(2401511, ['TG', 'A'], [1, 1]),
            _test_variant(2401513, ['TAC', 'T'], [1, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('xATGTACACAG', 2401509),
        expected_genotypes=[
            [1, 1],
            [1, 1],
        ],
    )

  # 20:2525695: genotype assign was incorrect in a previous run. This is because
  # the candidate variants overlap:
  #
  # ref: AAATT
  #  v1:  A--
  #  v2:   A-
  #
  # And this is causing us to construct incorrect haplotypes.
  def test_example6(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(2525696, ['AAT', 'A']),
            _test_variant(2525697, ['AT', 'T']),
        ],
        true_variants=[
            _test_variant(2525696, ['AAT', 'A'], [0, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('xAATT', 2525695),
        expected_genotypes=[
            [0, 1],
            [0, 0],
        ],
    )

  # Variants were getting incorrect genotypes due to complex region.
  #
  # variants: candidates
  #   20:279768:G->C gt=(-1, -1)
  #   20:279773:ATA->C/CTA gt=(-1, -1)
  # variants: truth
  #   20:279773:A->C gt=(1, 0)
  #
  # pos    : 789012345678901
  # ref    : CGCCCCATACCTTTT
  # truth  :       C          => CGCCCCCTACCTTTT
  # DV 1   :  C               => CCCCCCATACCTTTT [bad]
  # DV 2.a :       C--        => CGCCCCCCTTTT    [bad]
  # DV 2.b :       CTA        => CGCCCCCTACCTTTT [match]
  #
  def test_example7(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(279768, ['G', 'C']),
            _test_variant(279773, ['ATA', 'C', 'CTA']),
        ],
        true_variants=[
            _test_variant(279773, ['A', 'C'], [0, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('CGCCCCATACCTTTT', 279767),
        expected_genotypes=[
            [0, 0],
            [0, 2],
        ],
    )

  # TODO: retarget this test to a higher-level version of the API
  # that accepts a whole region of variants so we make sure it divides up the
  # problem into more fine-grained pieces that run quickly. The current call is
  # to a lower-level API that doesn't do variant chunking.
  # Commented out because this remains super slow.
  # def test_super_slow_example(self):
  #   self.assertGetsCorrectLabels(
  #       candidates=[
  #           _test_variant(32274452, ['C', 'G']),
  #           _test_variant(32274453, ['T', 'G']),
  #           _test_variant(32274456, ['A', 'G']),
  #           _test_variant(32274459, ['C', 'G']),
  #           _test_variant(32274461, ['T', 'G']),
  #           _test_variant(32274465, ['GACA', 'G']),
  #           _test_variant(32274467, ['CA', 'C']),
  #           _test_variant(32274470, ['C', 'G']),
  #           _test_variant(32274473, ['A', 'G']),
  #           _test_variant(32274474, ['AC', 'A']),
  #           _test_variant(32274475, ['C', 'A']),
  #           _test_variant(32274477, ['T', 'A']),
  #           _test_variant(32274480, ['G', 'C']),
  #       ],
  #       true_variants=[
  #           _test_variant(32274470, ['C', 'G'], (1, 1)),
  #       ],
  #       ref=haplotype_labeler.ReferenceRegion(
  #           'GCTGGAGGCGTGGGGACACCGGAACATAGGCCCCGCCCCGCCCCGACGC', 32274451),
  #       expected_genotypes=[
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [1, 1],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #           [0, 0],
  #       ])

  # Variants were getting incorrect genotypes in an exome callset.
  #
  # ref: AGACACACACACACAAAAAAAAATCATAAAATGAAG, start=214012389
  # candidates 2:214012390:G->GAC
  # candidates 2:214012402:CAA->C
  # candidates 2:214012404:A->C
  # true_variants 2:214012404:A->C
  #
  # 2:214012390:G->GAC => gt=(1, 1) new_label=2 old_label=0 alts=[0]
  # 2:214012402:CAA->C => gt=(1, 1) new_label=2 old_label=0 alts=[0]
  # 2:214012404:A->C => gt=(0, 0) new_label=0 old_label=2 alts=[0]
  #
  #           90--------- 0---------10--------20---
  # pos    : 90  1234567890123456789012345678901234
  # ref    : AG  ACACACACACACAAAAAAAAATCATAAAATGAAG
  # truth  :                  C => AGACACACACACACACAAAAAAATCATAAAATGAAG
  # DV 1   :  GAC               => [doesn't match]
  # DV 2   :                C-- => [doesn't match]
  # DV 1+2 : AGACACACACACACAC  AAAAAAATCATAAAATGAAG
  # DV 1+2 :                    => AGACACACACACACACAAAAAAATCATAAAATGAAG [match]
  # DV 3   :                  C => AGACACACACACACACAAAAAAATCATAAAATGAAG [match]
  #
  # So this is an interesting case. G->GAC + CAA->C matches the true haplotype,
  # and the SNP itself gets assigned a FP status since we can have either two
  # FPs (dv1 and dv2 candidates) or have just one (dv3). What's annoying here is
  # that DV3 exactly matches the variant as described in the truth set. It's
  # also strange that we've generated multiple equivalent potential variants
  # here.
  #
  # This test ensures that we are picking the most parsimonous genotype
  # assignment (e.g., fewest number of TPs) needed to explain the truth, after
  # accounting for minimizing the number of FNs and FPs.
  def test_exome_variants_multiple_equivalent_representations(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(214012390, ['G', 'GAC']),
            _test_variant(214012402, ['CAA', 'C']),
            _test_variant(214012404, ['A', 'C']),
        ],
        true_variants=[
            _test_variant(214012404, ['A', 'C'], [1, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion(
            'AGACACACACACACAAAAAAAAATCAT', 214012389
        ),
        expected_genotypes=[
            [0, 1],
            [0, 1],
            [0, 1],
            # This configuration makes the most sense but we cannot choose it
            # if we want to minimize the number of FNs, FPs, and then TPs.
            # [0, 0],
            # [0, 0],
            # [1, 1],
        ],
    )

  # Variant group: 5 candidates 2 truth variants
  # ref: ReferenceRegion(bases=TGTTTTTTTTTAAAAAAATTATTTCTTCTTT, start=167012239)
  #   candidates 4:167012240:GT->G
  #   candidates 4:167012246:TTTT->A
  #   candidates 4:167012247:T->A
  #   candidates 4:167012248:T->A
  #   candidates 4:167012249:T->A
  #   true_variants 4:167012240:GTTT->G/GTT [2, 1]
  #   true_variants 4:167012249:T->A/TAA [2, 1]
  #   4:167012240:GT->G => gt=(1, 1) new_label=2 old_label=1 alts=[0]
  #   4:167012246:TTTT->A => gt=(0, 0) new_label=0 old_label=0 alts=[0]
  #   4:167012247:T->A => gt=(0, 0) new_label=0 old_label=0 alts=[0]
  #   4:167012248:T->A => gt=(0, 1) new_label=1 old_label=0 alts=[0]
  #   4:167012249:T->A => gt=(1, 1) new_label=2 old_label=1 alts=[0]
  def test_exome_complex_example(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(167012240, ['GT', 'G']),
            _test_variant(167012246, ['TTTT', 'A']),
            _test_variant(167012247, ['T', 'A']),
            _test_variant(167012248, ['T', 'A']),
            _test_variant(167012249, ['T', 'A']),
        ],
        true_variants=[
            _test_variant(167012240, ['GTTT', 'G', 'GTT'], [1, 2]),
            _test_variant(167012249, ['T', 'A', 'TAA'], [1, 2]),
        ],
        ref=haplotype_labeler.ReferenceRegion(
            'TGTTTTTTTTTAAAAAAATTATTTCTTCTTT', 167012239
        ),
        expected_genotypes=[
            [1, 1],
            [0, 0],
            [0, 0],
            [0, 1],
            [1, 1],
        ],
    )

  # ref: GGGTGTGTGTGTGTGTGTGTGTGTGTGCGTGTGTGTGTTTGTGTTG, start=9508942
  #   candidates 20:9508943:GGT->G
  #   candidates 20:9508967:T->C/TGC
  #   candidates 20:9508967:T->C/TGC
  #   candidates 20:9508967:T->C/TGC
  #   true_variants 20:9508943:GGT->G [0, 1]
  #   true_variants 20:9508967:T->C/TGC [1, 2]
  #   20:9508943:GGT->G => gt=(0, 0) new_label=0 old_label=1 alts=[0
  #   20:9508967:T->C/TGC => gt=(1, 1) new_label=2 old_label=1 alts=[0]
  #   20:9508967:T->C/TGC => gt=(1, 1) new_label=0 old_label=1 alts=[1]
  #   20:9508967:T->C/TGC => gt=(1, 1) new_label=2 old_label=2 alts=[0, 1]
  #
  # This test fixes a bug where we weren't scoring our matches properly.
  # Previously we were not accounting for FPs in our score, so we were taking
  # a match with 0 FN, 1 FP, 1 TP over one with 0 FN, 0 FP, and 2 TP!
  #
  #      40------50--------60---------
  # pos: 2345678901234567890123456789012345678901234567
  # ref: GGGTGTGTGTGTGTGTGTGTGTGTGTGCGTGTGTGTGTTTGTGTTG
  # t1:   G--
  # t2a:                          C
  # t2b:                          Tgc
  #
  def test_bad_scoring_bug(self):
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(9508943, ['GGT', 'G']),
            _test_variant(9508967, ['T', 'C', 'TGC']),
        ],
        true_variants=[
            _test_variant(9508943, ['GGT', 'G'], [0, 1]),
            _test_variant(9508967, ['T', 'C', 'TGC'], [1, 2]),
        ],
        ref=haplotype_labeler.ReferenceRegion(
            'GGGTGTGTGTGTGTGTGTGTGTGTGTGCGTGTGTGTGTTTGTGTTG', 9508942
        ),
        expected_genotypes=[
            [0, 1],
            [1, 2],
        ],
    )

  def test_variants_at_edge_of_contig_work_end_to_end(self):
    # This test checks that we can label end-to-end variants at occur at the
    # start and at the end of a chromosome. This is unlikely in humans but can
    # occur in bacterial genomes. See internal for a motivating example.
    self.assertGetsCorrectLabels(
        candidates=[
            # At chrom start.
            _test_variant(0, ['A', 'G']),
            # At chrom end. I've included an insertion here because that's a
            # common representation when there are bases at the start.
            _test_variant(3, ['T', 'C', 'TCCC']),
        ],
        true_variants=[
            _test_variant(0, ['A', 'G'], [1, 1]),
            _test_variant(3, ['T', 'C'], [0, 1]),
        ],
        ref=haplotype_labeler.ReferenceRegion('ACGT', 0),
        expected_genotypes=[
            [1, 1],
            [0, 1],
        ],
    )


if __name__ == '__main__':
  absltest.main()
