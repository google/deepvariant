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
"""Tests for deepvariant .haplotypes."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import types


from tensorflow import flags
from absl.testing import absltest
from absl.testing import parameterized
import numpy as np

from third_party.nucleus.testing import test_utils

from deepvariant import haplotypes
from deepvariant.testing import flagsaver

FLAGS = flags.FLAGS


def _var(chrom='1',
         start=5,
         end=None,
         ref=None,
         alt=None,
         qual=50,
         genotype=None,
         likelihoods=None,
         sample_name='NA12878'):
  """Creates a Variant record for testing.

  Args:
    chrom: reference name for this variant
    start: start position on the contig
    end: end position on the contig
    ref: reference base(s)
    alt: list(str). alternate base(s)
    qual: PHRED scaled detection probability
    genotype: list of integers corresponding to the called genotype
    likelihoods: genotype likelihoods for this variant
    sample_name: sample name for the single call in the variant

  Returns:
    A Variant record created with the specified arguments.

  Raises:
    ValueError: Both ref and end are specified, and are inconsistent.
  """
  if ref is None and end is None:
    ref = 'A'
  elif ref is None:
    ref = 'A' * (end - start)
  elif ref is not None and end is not None and end != start + len(ref):
    raise ValueError('Inconsistent end and reference allele.')

  if alt is None:
    alt = ['C']
  if genotype is None:
    genotype = [0, 1]
  if likelihoods is None:
    likelihoods = [-1.0, -0.0506099933550872, -2.0]
  return test_utils.make_variant(
      chrom=chrom,
      start=start,
      alleles=[ref] + alt,
      qual=qual,
      filters=None,
      gt=genotype,
      gls=likelihoods,
      sample_name=sample_name)


def _resolvable_incompatible_inputs():
  """Returns a list of variants that are incompatible but can be resolved."""
  return [
      _var(
          start=20,
          ref='ACCCCC',
          alt=['A'],
          genotype=[0, 1],
          likelihoods=[-2., -0.0506099933550872, -1.]),
      _var(
          start=23,
          ref='C',
          alt=['T'],
          genotype=[1, 1],
          likelihoods=[-2., -0.3098039199714863, -0.3010299956639812])
  ]


def _resolved_compatible_outputs():
  """Returns a list of het variants that are correctly resolved."""
  return [
      _var(
          start=20,
          ref='ACCCCC',
          alt=['A'],
          genotype=[0, 1],
          likelihoods=[
              -1.658964842664435, -0.010604831683503404, -2.6589648426644352
          ]),
      _var(
          start=23,
          ref='C',
          alt=['T'],
          genotype=[0, 1],
          likelihoods=[
              -1.658964842664435, -0.014526253196596468, -1.9599948383284163
          ])
  ]


class ResolveOverlappingVariantsTest(parameterized.TestCase):

  @flagsaver.FlagSaver
  def test_maybe_resolve_conflicting_variants(self):
    FLAGS.disable_haplotype_resolution = False
    # Note: Most of the resolution code is tested below in the
    # test_resolve_overlapping_variants function. This test mostly just ensures
    # that the interaction with RefCall variants is properly handled.
    ref_call_deletion = _var(
        start=1,
        end=30,
        genotype=[0, 0],
        # Not a real likelihood -- if we weren't just punting
        # this would get rescaled to sum to 1.
        likelihoods=[-1, -2, -3])
    independent_hom_alts = [
        _var(start=i, genotype=[1, 1], likelihoods=[-3, -2, -1])
        for i in range(3, 20)
    ]

    variants = ([ref_call_deletion] + independent_hom_alts +
                _resolvable_incompatible_inputs())
    expected = ([ref_call_deletion] + independent_hom_alts +
                _resolved_compatible_outputs())
    actual = haplotypes.maybe_resolve_conflicting_variants(variants)
    self._assert_generator_of_variants_equals_expected(actual, expected)

  @parameterized.parameters(
      dict(
          disable_haplotype_resolution=False,
          variants=_resolvable_incompatible_inputs(),
          expected=_resolved_compatible_outputs()),
      dict(
          disable_haplotype_resolution=True,
          variants=_resolvable_incompatible_inputs(),
          expected=_resolvable_incompatible_inputs()),
  )
  @flagsaver.FlagSaver
  def test_can_disable_haplotype_resolution(self, disable_haplotype_resolution,
                                            variants, expected):
    FLAGS.disable_haplotype_resolution = disable_haplotype_resolution
    actual = haplotypes.maybe_resolve_conflicting_variants(variants)
    self._assert_generator_of_variants_equals_expected(actual, expected)

  @parameterized.parameters(
      # The simple case where there is a single variant.
      dict(
          variants=[
              _var(
                  start=10,
                  ref='A',
                  alt=['C'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.])
          ],
          expected=[
              _var(
                  start=10,
                  ref='A',
                  alt=['C'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.])
          ]),
      # Cases where the actual genotype calls are compatible.
      dict(
          variants=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.])
          ],
          expected=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.])
          ]),
      dict(
          variants=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=21,
                  ref='C',
                  alt=['G'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.])
          ],
          expected=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=21,
                  ref='C',
                  alt=['G'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.])
          ]),
      dict(
          variants=[
              _var(
                  start=20,
                  ref='ACCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=22,
                  ref='CCCGAGAGAG',
                  alt=['C'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.]),
              _var(
                  start=25,
                  ref='G',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.])
          ],
          expected=[
              _var(
                  start=20,
                  ref='ACCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=22,
                  ref='CCCGAGAGAG',
                  alt=['C'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.]),
              _var(
                  start=25,
                  ref='G',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[-3., -0.004803708402820599, -2.])
          ]),
      # Cases where the actual genotype calls are incompatible.
      dict(
          variants=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T'],
                  genotype=[1, 1],
                  likelihoods=[-2., -0.3098039199714863, -0.3010299956639812])
          ],
          expected=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[
                      -1.658964842664435, -0.010604831683503404,
                      -2.6589648426644352
                  ]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T'],
                  genotype=[0, 1],
                  likelihoods=[
                      -1.658964842664435, -0.014526253196596468,
                      -1.9599948383284163
                  ])
          ]),
      dict(
          variants=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[-2., -0.0506099933550872, -1.]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T', 'G'],
                  genotype=[1, 2],
                  likelihoods=[
                      -2.0, -1.0, -0.6989700043360187, -0.958607314841775,
                      -0.4814860601221125, -0.6020599913279624
                  ])
          ],
          expected=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[0, 1],
                  likelihoods=[
                      -1.315550534421905, -0.02373784695478589,
                      -2.315550534421905
                  ]),
              _var(
                  start=23,
                  ref='C',
                  alt=['T', 'G'],
                  genotype=[0, 2],
                  likelihoods=[
                      -1.315550534421905, -0.36130802498257997,
                      -2.0145205387579237, -0.319915339824355,
                      -1.7970365945440174, -1.9176105257498672
                  ])
          ]),
      # Issues we can't currently resolve.
      dict(
          variants=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[1, 1],
                  likelihoods=[
                      -1.5228787452803376, -0.09691001300805639,
                      -0.7695510786217261
                  ]),
              _var(
                  start=23,
                  ref='CCCGATGAT',
                  alt=['C'],
                  genotype=[1, 1],
                  likelihoods=[
                      -1.3979400086720375, -0.1366771398795441,
                      -0.638272163982407
                  ]),
              _var(
                  start=24,
                  ref='C',
                  alt=['G'],
                  genotype=[1, 1],
                  likelihoods=[
                      -1.5228787452803376, -0.13076828026902382,
                      -0.638272163982407
                  ])
          ],
          expected=[
              _var(
                  start=20,
                  ref='ACCCCC',
                  alt=['A'],
                  genotype=[1, 1],
                  likelihoods=[
                      -1.5228787452803376, -0.09691001300805639,
                      -0.7695510786217261
                  ]),
              _var(
                  start=23,
                  ref='CCCGATGAT',
                  alt=['C'],
                  genotype=[1, 1],
                  likelihoods=[
                      -1.3979400086720375, -0.1366771398795441,
                      -0.638272163982407
                  ]),
              _var(
                  start=24,
                  ref='C',
                  alt=['G'],
                  genotype=[1, 1],
                  likelihoods=[
                      -1.5228787452803376, -0.13076828026902382,
                      -0.638272163982407
                  ])
          ]),
      # Too many variants to resolve.
      dict(
          variants=[
              _var(
                  start=1,
                  end=30,
                  genotype=[0, 1],
                  # Not a real likelihood -- if we weren't just punting
                  # this would get rescaled to sum to 1.
                  likelihoods=[-2, -1, -3])
          ] + [
              _var(start=i, genotype=[1, 1], likelihoods=[-3, -2, -1])
              for i in range(3, 25)
          ],
          expected=[
              _var(
                  start=1,
                  end=30,
                  genotype=[0, 1],
                  # Not a real likelihood -- if we weren't just punting
                  # this would get rescaled to sum to 1.
                  likelihoods=[-2, -1, -3])
          ] + [
              _var(start=i, genotype=[1, 1], likelihoods=[-3, -2, -1])
              for i in range(3, 25)
          ]),
  )
  def test_resolve_overlapping_variants(self, variants, expected):
    actual = haplotypes._resolve_overlapping_variants(variants)
    self._assert_generator_of_variants_equals_expected(actual, expected)

  @parameterized.parameters(
      dict(variants=[], expected=[]),
      dict(
          variants=[_var(chrom='1', start=5, end=6)],
          expected=[[_var(chrom='1', start=5, end=6)]]),
      # Test finding overlaps and different chromosomes not overlapping.
      dict(
          variants=[
              _var(chrom='1', start=1, end=5),
              _var(chrom='1', start=5, end=7),
              _var(chrom='1', start=6, end=8),
              _var(chrom='2', start=6, end=10)
          ],
          expected=[[_var(chrom='1', start=1, end=5)], [
              _var(chrom='1', start=5, end=7),
              _var(chrom='1', start=6, end=8)
          ], [_var(chrom='2', start=6, end=10)]]),
      # Test one large variant spanning multiple others.
      dict(
          variants=[
              _var(chrom='1', start=1, end=25),
              _var(chrom='1', start=3, end=5),
              _var(chrom='1', start=7, end=8),
              _var(chrom='1', start=14, end=20),
              _var(chrom='1', start=24, end=27)
          ],
          expected=[[
              _var(chrom='1', start=1, end=25),
              _var(chrom='1', start=3, end=5),
              _var(chrom='1', start=7, end=8),
              _var(chrom='1', start=14, end=20),
              _var(chrom='1', start=24, end=27)
          ]]),
      # Test mix of non-overlapping and ending with an overlap.
      dict(
          variants=[
              _var(chrom='1', start=1, end=5),
              _var(chrom='2', start=3, end=5),
              _var(chrom='2', start=7, end=10),
              _var(chrom='2', start=9, end=10)
          ],
          expected=[[_var(chrom='1', start=1, end=5)],
                    [_var(chrom='2', start=3, end=5)], [
                        _var(chrom='2', start=7, end=10),
                        _var(chrom='2', start=9, end=10)
                    ]]),
      # Same as prior test but using a generator as input.
      dict(
          variants=(_var(chrom='1', start=1, end=5),
                    _var(chrom='2', start=3, end=5),
                    _var(chrom='2', start=7, end=10),
                    _var(chrom='2', start=9, end=10)),
          expected=[[_var(chrom='1', start=1, end=5)],
                    [_var(chrom='2', start=3, end=5)], [
                        _var(chrom='2', start=7, end=10),
                        _var(chrom='2', start=9, end=10)
                    ]]),
  )
  def test_group_overlapping_variants(self, variants, expected):
    actual = haplotypes._group_overlapping_variants(variants)
    self.assertIsInstance(actual, types.GeneratorType)
    actual_list = list(actual)
    self.assertEqual(actual_list, expected)

  @parameterized.parameters(
      dict(variants=[_var(start=10, end=20)], nonref_counts=[2], expected=True),
      dict(variants=[_var(start=10, end=20)], nonref_counts=[1], expected=True),
      dict(
          variants=[_var(start=10, end=20),
                    _var(start=15, end=20)],
          nonref_counts=[1, 1],
          expected=True),
      dict(
          variants=[_var(start=10, end=20),
                    _var(start=15, end=20)],
          nonref_counts=[1, 2],
          expected=False),
      dict(
          variants=[
              _var(start=10, end=20),
              _var(start=15, end=25),
              _var(start=20, end=25)
          ],
          nonref_counts=[1, 1, 1],
          expected=True),
      dict(
          variants=[
              _var(start=10, end=20),
              _var(start=15, end=25),
              _var(start=20, end=25)
          ],
          nonref_counts=[1, 2, 1],
          expected=False),
      dict(
          variants=[
              _var(start=10, end=20),
              _var(start=15, end=25),
              _var(start=20, end=25)
          ],
          nonref_counts=[1, 1, 2],
          expected=False),
      dict(
          variants=[
              _var(start=10, end=20),
              _var(start=15, end=25),
              _var(start=19, end=25)
          ],
          nonref_counts=[2, 0, 2],
          expected=False),
      dict(
          variants=[
              _var(start=10, end=20),
              _var(start=15, end=25),
              _var(start=20, end=25)
          ],
          nonref_counts=[2, 0, 2],
          expected=True),
  )
  def test_all_variants_compatible(self, variants, nonref_counts, expected):
    calculator = haplotypes._VariantCompatibilityCalculator(variants)
    actual = calculator.all_variants_compatible(nonref_counts)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(variants=[_var(start=10, end=20)], nonref_counts=[3]),
      dict(variants=[_var(start=10, end=20)], nonref_counts=[1, 2]),
      dict(variants=[_var(start=10, end=20)], nonref_counts=[]),
      dict(
          variants=[_var(start=10, end=20),
                    _var(start=15, end=20)],
          nonref_counts=[1]),
      dict(
          variants=[_var(start=10, end=20),
                    _var(start=15, end=20)],
          nonref_counts=[1, 2, 1]),
  )
  def test_invalid_all_variants_compatible(self, variants, nonref_counts):
    calculator = haplotypes._VariantCompatibilityCalculator(variants)
    with self.assertRaisesRegexp(ValueError, 'variant'):
      calculator.all_variants_compatible(nonref_counts)

  @parameterized.parameters(
      dict(num_alts_list=[1, 1], config=[0, 1], expected=[
          ((0, 0), (0, 1)),
      ]),
      dict(
          num_alts_list=[1, 2],
          config=[0, 1],
          expected=[((0, 0), (0, 1)), ((0, 0), (0, 2))]),
      dict(
          num_alts_list=[2, 2],
          config=[2, 1],
          expected=[((1, 1), (0, 1)), ((1, 1), (0, 2)), ((1, 2), (0, 1)),
                    ((1, 2), (0, 2)), ((2, 2), (0, 1)), ((2, 2), (0, 2))]),
  )
  def test_get_all_allele_indices_configurations(self, num_alts_list, config,
                                                 expected):

    def get_alt_bases(num_alts):
      return ['C' + 'A' * i for i in xrange(1, num_alts + 1)]

    variants = [
        _var(ref='C', alt=get_alt_bases(num_alts)) for num_alts in num_alts_list
    ]
    actual = haplotypes._get_all_allele_indices_configurations(variants, config)
    self.assertEqual(list(actual), expected)

  def test_invalid_get_all_allele_indices_configurations(self):
    with self.assertRaisesRegexp(ValueError, r'len\(variants\) must equal'):
      haplotypes._get_all_allele_indices_configurations(
          _resolved_compatible_outputs(), [1])

  def test_invalid_allele_indices_configuration_likelihood(self):
    with self.assertRaisesRegexp(ValueError, r'equal len\(allele_indices_conf'):
      haplotypes._allele_indices_configuration_likelihood(
          _resolved_compatible_outputs(), [(1, 1)])

  @parameterized.parameters(
      dict(genotype=[-1, -1], expected=0),
      dict(genotype=[0, 0], expected=0),
      dict(genotype=[0, 1], expected=1),
      dict(genotype=[1, 0], expected=1),
      dict(genotype=[1, 1], expected=2),
      dict(genotype=[1, 2], expected=2),
      dict(genotype=[2, 2], expected=2),
  )
  def test_nonref_genotype_count(self, genotype, expected):
    variant = _var(ref='AC', alt=['A', 'ACC'], genotype=genotype)
    actual = haplotypes._nonref_genotype_count(variant)
    self.assertEqual(actual, expected)

  def test_invalid_nonref_genotype_count(self):
    zero_calls_variant = test_utils.make_variant()
    with self.assertRaisesRegexp(ValueError, 'Expected exactly one VariantCal'):
      haplotypes._nonref_genotype_count(zero_calls_variant)

  def _assert_generator_of_variants_equals_expected(self, actual, expected):
    """Helper method to compare a generator of Variant protos to a list."""
    self.assertIsInstance(actual, types.GeneratorType)
    actual_list = list(actual)
    self.assertEqual(len(actual_list), len(expected))
    for actual_variant, expected_variant in zip(actual_list, expected):
      self._assert_variants_equal_with_likelihood_tolerance(
          actual_variant, expected_variant)

  def _assert_variants_equal_with_likelihood_tolerance(self,
                                                       v1,
                                                       v2,
                                                       tolerance=1e-10):
    """Asserts variant equality allowing numerical differences in GLs."""
    gl1 = list(v1.calls[0].genotype_likelihood)
    gl2 = list(v2.calls[0].genotype_likelihood)
    np.testing.assert_allclose(gl1, gl2, rtol=tolerance)
    v1.calls[0].genotype_likelihood[:] = []
    v2.calls[0].genotype_likelihood[:] = []
    self.assertEqual(v1, v2)


if __name__ == '__main__':
  absltest.main()
