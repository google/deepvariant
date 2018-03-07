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
"""Tests for deepvariant .postprocess_variants."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import errno
import itertools
import sys


from tensorflow import flags
from absl.testing import absltest
from absl.testing import parameterized
import mock
import numpy as np
import tensorflow as tf

from absl import logging
from deepvariant.core.genomics import struct_pb2
from deepvariant.core.genomics import variants_pb2
from deepvariant import postprocess_variants
from deepvariant import test_utils
from deepvariant.core import genomics_math
from deepvariant.core import io_utils
from deepvariant.core import variantutils
from deepvariant.core.protos import core_pb2
from deepvariant.protos import deepvariant_pb2
from deepvariant.testing import flagsaver

FLAGS = flags.FLAGS

_DEFAULT_SAMPLE_NAME = 'NA12878'

# Test contigs for gVCF merging code.
_CONTIGS = [
    core_pb2.ContigInfo(name='1', n_bases=100),
    core_pb2.ContigInfo(name='2', n_bases=200),
    core_pb2.ContigInfo(name='10', n_bases=300),
]


class DummyReferenceReader(object):

  def __init__(self):
    self._bases = {
        '1': 'AACCGGTTACGTTCGATTTTAAAACCCCGGGG',
        '2': 'GCAGTGACGTAGCGATGACGTAGACGCTTACG'
    }

  def bases(self, region):
    chrom = region.reference_name
    if chrom not in self._bases or region.end > len(self._bases[chrom]):
      raise ValueError('Invalid region for dummy reader: {}'.format(region))
    return self._bases[chrom][region.start:region.end]


def setUpModule():
  test_utils.init()


def _create_variant(ref_name, start, ref_base, alt_bases, qual, filter_field,
                    genotype, gq, likelihoods):
  """Creates a Variant record for testing.

  Args:
    ref_name: reference name for this variant
    start: start position on the contig
    ref_base: reference base(s)
    alt_bases: list(str). alternate base(s)
    qual: PHRED scaled detection probability
    filter_field: filter string for this variant
    genotype: list of integers corresponding to the called genotype
    gq: PHRED scaled genotype quality
    likelihoods: genotype likelihoods for this variant

  Returns:
    A Variant record created with the specified arguments.
  """
  return test_utils.make_variant(
      chrom=ref_name,
      start=start,
      alleles=[ref_base] + alt_bases,
      qual=qual,
      filters=filter_field,
      gt=genotype,
      gq=gq,
      gls=likelihoods,
      sample_name=_DEFAULT_SAMPLE_NAME)


def _create_variant_with_alleles(ref=None, alts=None, start=0):
  """Creates a Variant record with specified alternate_bases."""
  return variants_pb2.Variant(
      reference_bases=ref,
      alternate_bases=alts,
      start=start,
      calls=[variants_pb2.VariantCall(call_set_name=_DEFAULT_SAMPLE_NAME)])


def _create_call_variants_output(indices,
                                 probabilities,
                                 ref=None,
                                 alts=None,
                                 variant=None):
  if alts is None != variant is None:
    raise ValueError('Exactly one of either `alts` or `variant` should be set.')
  if not variant:
    variant = _create_variant_with_alleles(ref=ref, alts=alts)
  return deepvariant_pb2.CallVariantsOutput(
      genotype_probabilities=probabilities,
      alt_allele_indices=deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
          indices=indices),
      variant=variant)


def _simple_variant(ref_name, start, ref_base):
  """Creates a Variant record for testing variant and non-variant merge.

  Args:
    ref_name: str. Reference name for this variant.
    start: int. start position on the contig [0-based, half open).
    ref_base: str. reference base(s).

  Returns:
    A Variant record created with the specified arguments.
  """
  return test_utils.make_variant(
      chrom=ref_name,
      start=start,
      end=start + len(ref_base),
      alleles=[ref_base, 'A' if ref_base != 'A' else 'C'])


def _create_nonvariant(ref_name, start, end, ref_base):
  """Creates a non-variant Variant record for testing.

  Args:
    ref_name: str. Reference name for this variant.
    start: int. start position on the contig [0-based, half open).
    end: int. end position on the contig [0-based, half open).
    ref_base: str. reference base at the start position.

  Returns:
    A non-variant Variant record created with the specified arguments.
  """
  return test_utils.make_variant(
      chrom=ref_name,
      start=start,
      end=end,
      alleles=[ref_base, variantutils.GVCF_ALT_ALLELE])


def make_golden_dataset(compressed_inputs=False):
  if compressed_inputs:
    source_path = test_utils.test_tmpfile(
        'golden.postprocess_single_site_input.tfrecord.gz')
    io_utils.write_tfrecords(
        io_utils.read_tfrecords(
            test_utils.GOLDEN_POSTPROCESS_INPUT,
            proto=deepvariant_pb2.CallVariantsOutput), source_path)
  else:
    source_path = test_utils.GOLDEN_POSTPROCESS_INPUT
  return source_path


class AlleleRemapperTest(parameterized.TestCase):

  @parameterized.parameters(
      (list('C'), [], [True]),
      (list('C'), ['C'], [False]),
      (list('CT'), [], [True, True]),
      (list('CT'), ['C'], [False, True]),
      (list('CT'), ['T'], [True, False]),
      (list('CT'), ['C', 'T'], [False, False]),
      (list('CTG'), ['C'], [False, True, True]),
      (list('CTG'), ['T'], [True, False, True]),
      (list('CTG'), ['G'], [True, True, False]),
      (list('CTG'), ['C', 'G'], [False, True, False]),
  )
  def test_basic(self, alt_alleles, remove, keep_index_expected):
    remapper = postprocess_variants.AlleleRemapper(alt_alleles, remove)
    self.assertEqual(remapper.original_alts, alt_alleles)
    self.assertEqual(remapper.alleles_to_remove, set(remove))
    self.assertEqual(keep_index_expected,
                     [remapper.keep_index(i) for i in range(len(alt_alleles))])
    # When our i is 1 for the first alt allele, we expect that we will get back
    # our keep_index_expected but also that keep_index(i==0) is True for the
    # reference allele.
    self.assertEqual([True] + keep_index_expected, [
        remapper.keep_index(i, ref_is_zero=True)
        for i in range(len(alt_alleles) + 1)
    ])

  def test_makes_copy_of_inputs(self):
    alt_alleles = ['A', 'B']
    removes = {'B'}
    remapper = postprocess_variants.AlleleRemapper(alt_alleles, removes)
    del alt_alleles[0]
    removes -= {'B'}
    self.assertEqual(remapper.original_alts, ['A', 'B'])
    self.assertEqual(remapper.alleles_to_remove, {'B'})


class PostprocessVariantsTest(parameterized.TestCase):

  @parameterized.parameters(False, True)
  @flagsaver.FlagSaver
  def test_call_end2end(self, compressed_inputs):
    FLAGS.infile = make_golden_dataset(compressed_inputs)
    FLAGS.ref = test_utils.CHR20_FASTA
    FLAGS.outfile = test_utils.test_tmpfile('calls.vcf')
    FLAGS.nonvariant_site_tfrecord_path = (
        test_utils.GOLDEN_POSTPROCESS_GVCF_INPUT)
    FLAGS.gvcf_outfile = test_utils.test_tmpfile('gvcf_calls.vcf')

    postprocess_variants.main(['postprocess_variants.py'])

    self.assertEqual(
        tf.gfile.FastGFile(FLAGS.outfile).readlines(),
        tf.gfile.FastGFile(test_utils.GOLDEN_POSTPROCESS_OUTPUT).readlines())
    self.assertEqual(
        tf.gfile.FastGFile(FLAGS.gvcf_outfile).readlines(),
        tf.gfile.FastGFile(
            test_utils.GOLDEN_POSTPROCESS_GVCF_OUTPUT).readlines())

  def test_extract_single_variant_name(self):
    record = _create_call_variants_output(
        indices=[0], probabilities=[0.19, 0.75, 0.06], ref='A', alts=['C', 'T'])
    expected = _DEFAULT_SAMPLE_NAME
    actual = postprocess_variants._extract_single_sample_name(record)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      ([],),
      (['multiple', 'names'],),
  )
  def test_exception_extract_single_variant_name(self, names):
    variant_calls = [
        variants_pb2.VariantCall(call_set_name=name) for name in names
    ]
    variant = variants_pb2.Variant(calls=variant_calls)
    record = deepvariant_pb2.CallVariantsOutput(variant=variant)
    with self.assertRaisesRegexp(ValueError, 'Error extracting name:'):
      postprocess_variants._extract_single_sample_name(record)

  @parameterized.parameters(
      ([
          _create_call_variants_output(
              indices=[0],
              probabilities=[0.19, 0.75, 0.06],
              ref='A',
              alts=['C', 'T']),
          _create_call_variants_output(
              indices=[1],
              probabilities=[0.03, 0.93, 0.04],
              ref='A',
              alts=['C', 'T']),
          _create_call_variants_output(
              indices=[0, 1],
              probabilities=[0.03, 0.92, 0.05],
              ref='A',
              alts=['C', 'T'])
      ], set(), {
          ('A', 'A'): [0.19, 0.03, 0.03],
          ('A', 'C'): [0.75, 0.92],
          ('C', 'C'): [0.06, 0.05],
          ('A', 'T'): [0.93, 0.92],
          ('T', 'T'): [0.04, 0.05],
          ('C', 'T'): [0.05],
          ('T', 'C'): [0.05],
      }),
      # Example where all alt alleles are below qual_filter, but we keep one
      # where the qual is highest among the ones filtered out ('T')
      ([
          _create_call_variants_output(
              indices=[0],
              probabilities=[0.19, 0.75, 0.06],
              ref='A',
              alts=['C', 'T']),
          _create_call_variants_output(
              indices=[1],
              probabilities=[0.03, 0.93, 0.04],
              ref='A',
              alts=['C', 'T']),
          _create_call_variants_output(
              indices=[0, 1],
              probabilities=[0.03, 0.92, 0.05],
              ref='A',
              alts=['C', 'T'])
      ], set(['C']), {
          ('A', 'A'): [0.03],
          ('A', 'T'): [0.93],
          ('T', 'T'): [0.04],
      }),
  )
  def test_convert_call_variants_outputs_to_probs_dict(
      self, call_variants_outputs, alt_alleles_to_remove, expected_probs_dict):
    # In the current code, all call_variants_outputs have the same variant
    # field.
    canonical_variant = call_variants_outputs[0].variant
    self.assertEqual(
        postprocess_variants.convert_call_variants_outputs_to_probs_dict(
            canonical_variant,
            call_variants_outputs,
            alt_alleles_to_remove=alt_alleles_to_remove), expected_probs_dict)

  @parameterized.parameters(
      # Example with 2 alternate_bases:
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
      ([
          _create_call_variants_output(
              indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['C', 'T']),
          _create_call_variants_output(
              indices=[1], probabilities=[0.03, 0.93, 0.04], alts=['C', 'T']),
          _create_call_variants_output(
              indices=[0, 1], probabilities=[0.03, 0.92, 0.05], alts=['C', 'T'])
      ], [0.03, 0.75, 0.05, 0.92, 0.05, 0.04]),
      # One more example with 2 alternate_bases:
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
      ([
          _create_call_variants_output(
              indices=[1], probabilities=[0.978, 0.03, 0.002], alts=['C', 'T']),
          _create_call_variants_output(
              indices=[0, 1],
              probabilities=[0.992, 0.007, 0.001],
              alts=['C', 'T']),
          _create_call_variants_output(
              indices=[0],
              probabilities=[0.99997, 0.00002, 0.00001],
              alts=['C', 'T']),
      ], [0.978, 0.00002, 0.00001, 0.007, 0.001, 0.001]),
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1.
      ([
          _create_call_variants_output(
              indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['A']),
      ], [0.19, 0.75, 0.06]),
      # expected_unnormalized_probs is min of
      # 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, 0/3, 1/3, 2/3, 3/3.
      ([
          _create_call_variants_output(
              indices=[0],
              probabilities=[0.999, 0.001, 0],
              alts=['C', 'G', 'T']),
          _create_call_variants_output(
              indices=[0, 1], probabilities=[0, 1, 0], alts=['C', 'G', 'T']),
          _create_call_variants_output(
              indices=[0, 2],
              probabilities=[0.0001, 0.9996, 0.0003],
              alts=['C', 'G', 'T']),
          _create_call_variants_output(
              indices=[1], probabilities=[0, 1, 0], alts=['C', 'G', 'T']),
          _create_call_variants_output(
              indices=[1, 2],
              probabilities=[0.0001, 0.0002, 0.9997],
              alts=['C', 'G', 'T']),
          _create_call_variants_output(
              indices=[2],
              probabilities=[0.00004, 0.9999, 0.00006],
              alts=['C', 'G', 'T']),
      ], [0, 0.001, 0, 0.0002, 0, 0, 0.0002, 0.0003, 0.9997, 0.00006]),
  )
  def test_merge_predictions(self, inputs, expected_unnormalized_probs):
    denominator = sum(expected_unnormalized_probs)
    for permuted_inputs in itertools.permutations(inputs):
      _, predictions = postprocess_variants.merge_predictions(permuted_inputs)
      np.testing.assert_almost_equal(
          predictions, [x / denominator for x in expected_unnormalized_probs])

  @parameterized.parameters(
      # With 1 alt allele, we expect to see 1 alt_allele_indices: [0].
      (
          [
              _create_call_variants_output(
                  indices=[1], probabilities=[0.19, 0.75, 0.06], alts=['A']),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # With 2 alt alleles, we expect to see 3 alt_allele_indices.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.19, 0.75, 0.06],
                  alts=['G', 'T']),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.03, 0.93, 0.04],
                  alts=['G', 'T']),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # With 2 alt alleles, we expect to see 3 alt_allele_indices:
      # [0], [1], [0, 1].
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.19, 0.75, 0.06],
                  alts=['G', 'T']),
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.03, 0.93, 0.04],
                  alts=['G', 'T']),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.93, 0.04],
                  alts=['G', 'T']),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.19, 0.75, 0.06],
                  alts=['G', 'T']),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      ([], 'Expected 1 or more call_variants_outputs.'),
      # With 3 alt alleles, we expect to see 6 alt_allele_indices.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  alts=['AA', 'T', 'AAA']),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0, 1, 0],
                  alts=['AA', 'T', 'AAA']),
              _create_call_variants_output(
                  indices=[0, 2],
                  probabilities=[0.0001, 0.9996, 0.0003],
                  alts=['AA', 'T', 'AAA'])
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # reference_bases have to be exactly the same.
      ([
          _create_call_variants_output(
              indices=[0],
              probabilities=[0.999, 0.001, 0],
              variant=_create_variant_with_alleles(ref='A', alts=['T', 'C'])),
          _create_call_variants_output(
              indices=[1],
              probabilities=[0.2, 0.8, 0],
              variant=_create_variant_with_alleles(ref='A', alts=['T', 'C'])),
          _create_call_variants_output(
              indices=[0, 1],
              probabilities=[0.2, 0.8, 0],
              variant=_create_variant_with_alleles(ref='G', alts=['T', 'C'])),
      ], '`call_variants_outputs` did not pass sanity check.'),
      # alternate_bases have to be exactly the same. Different orders are
      # not acceptable either.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  variant=_create_variant_with_alleles(alts=['T', 'C'])),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.2, 0.8, 0],
                  variant=_create_variant_with_alleles(alts=['T', 'C'])),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.2, 0.8, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T'])),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
  )
  def test_exception_merge_predictions(self, inputs, text):
    with self.assertRaisesRegexp(ValueError, text):
      postprocess_variants.merge_predictions(inputs)

  @parameterized.parameters(
      ([0.01, 0.0, 0.99],
       _create_variant('GL000220.1', 1, 'A', ['.'], 20.0,
                       postprocess_variants.DEEP_VARIANT_PASS, [1, 1], 20,
                       [-2.0, -15.0003472607, -0.0043648054])),
      ([0.01, 0.0, 0.99],
       _create_variant('GL000220.1', 10000210, 'C', ['T'], 20.0,
                       postprocess_variants.DEEP_VARIANT_PASS, [1, 1], 20,
                       [-2.0, -15.0003472607, -0.0043648054])),
      ([0.001, 0.999, 0.0],
       _create_variant('20', 10000210, 'C', ['CT'], 30.0,
                       postprocess_variants.DEEP_VARIANT_PASS, [0, 1], 30,
                       [-3.0, -0.00043451177, -15.0003472607])),
      ([0.0001, 0.0, 0.9999],
       _create_variant('1', 1, 'C', ['T'], 40.0,
                       postprocess_variants.DEEP_VARIANT_PASS, [1, 1], 40,
                       [-4.0, -15.0003472607, -0.00004343161])),
      ([0.1, 0.90, 0.0],
       _create_variant('20', 10000210, 'A', ['T'], 10.0,
                       postprocess_variants.DEEP_VARIANT_PASS, [0, 1], 10,
                       [-1.0, -0.04575749056, -15.0003472607])),
      ([0.99, 0.005, 0.005],
       _create_variant('X', 10000210, 'CACA', ['C'], 0.04364805402,
                       postprocess_variants.DEEP_VARIANT_REF_FILTER, [0, 0], 20,
                       [-0.0043648054, -2.30102999566, -2.30102999566])),
      ([0.9999, 0.0001, 0.0],
       _create_variant('chrY', 10000210, 'C', ['T'], 0.00043431619,
                       postprocess_variants.DEEP_VARIANT_REF_FILTER, [0, 0], 40,
                       [-0.00004343161, -4.0, -15.0003472607])),
      # Multi-allelic test examples.
      ([0.995, 0.001, 0.001, 0.001, 0.001, 0.001],
       _create_variant('X', 10000210, 'CACA', ['C', 'A'], 0.0217691925,
                       postprocess_variants.DEEP_VARIANT_REF_FILTER, [0, 0], 23,
                       [-0.00217691925, -3, -3, -3, -3, -3])),
      ([0.001, 0.001, 0.001, 0.995, 0.001, 0.001],
       _create_variant('X', 10000210, 'CACA', ['C', 'A'], 30,
                       postprocess_variants.DEEP_VARIANT_PASS, [0, 2], 23,
                       [-3, -3, -3, -0.00217691925, -3, -3])),
  )
  def test_add_call_to_variant(self, probs, expected):
    raw_variant = variants_pb2.Variant(
        reference_name=expected.reference_name,
        reference_bases=expected.reference_bases,
        alternate_bases=expected.alternate_bases,
        start=expected.start,
        end=expected.end,
        calls=[variants_pb2.VariantCall(call_set_name=_DEFAULT_SAMPLE_NAME)])
    variant = postprocess_variants.add_call_to_variant(
        variant=raw_variant,
        predictions=probs,
        sample_name=_DEFAULT_SAMPLE_NAME)
    self.assertEqual(variant.reference_bases, expected.reference_bases)
    self.assertEqual(variant.alternate_bases, expected.alternate_bases)
    self.assertEqual(variant.reference_name, expected.reference_name)
    self.assertEqual(variant.start, expected.start)
    self.assertEqual(variant.end, expected.end)
    self.assertAlmostEquals(variant.quality, expected.quality, places=6)
    self.assertEqual(variant.filter, expected.filter)
    self.assertEqual(len(variant.calls), 1)
    self.assertEqual(len(expected.calls), 1)
    self.assertEqual(variant.calls[0].genotype, expected.calls[0].genotype)
    self.assertEqual(variant.calls[0].info['GQ'], expected.calls[0].info['GQ'])
    for gl, expected_gl in zip(variant.calls[0].genotype_likelihood,
                               expected.calls[0].genotype_likelihood):
      self.assertAlmostEquals(gl, expected_gl, places=6)

  @parameterized.parameters(
      (
          0,
          [0, 0],
      ),
      (
          1,
          [0, 1],
      ),
      (
          2,
          [1, 1],
      ),
      (
          3,
          [0, 2],
      ),
      (
          4,
          [1, 2],
      ),
      (
          5,
          [2, 2],
      ),
  )
  def test_triallelic_genotype_in_add_call_to_variant(
      self, highest_prob_position, expected_best_genotype):
    """Ensures the order of GLs are interpreted correctly for triallelics."""
    raw_variant = _create_variant_with_alleles(ref='CACA', alts=['C', 'A'])
    # Create a probability with 6 elements, one of them 0.995 (best genotype),
    # and the rest 0.001.
    probs = [0.001] * 6
    assert 0 <= highest_prob_position <= len(probs)
    probs[highest_prob_position] = 0.995
    variant = postprocess_variants.add_call_to_variant(
        variant=raw_variant,
        predictions=probs,
        sample_name=_DEFAULT_SAMPLE_NAME)
    self.assertEqual(variant.calls[0].genotype, expected_best_genotype)

  @parameterized.parameters(
      # Q20 tests
      ([0.01, 0.0, 0.99], 0, 0, 20.0),
      ([0.01, 0.0, 0.99], 1, 0, 20.0),
      ([0.01, 0.0, 0.99], 2, 20, 20.0),
      # Q30 tests
      ([0.001, 0.0, 0.999], 0, 0, 30.0),
      ([0.001, 0.0, 0.999], 1, 0, 30.0),
      ([0.001, 0.0, 0.999], 2, 30, 30.0),
      # Q40 tests
      ([0.0001, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.9999], 1, 0, 40.0),
      ([0.0001, 0.0, 0.9999], 2, 40, 40.0),
      # Test that this works with any sized genotype vector.
      ([0.0001, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9999], 0, 0, 40.0),
  )
  def test_compute_quals(self, probs, call, expected_gq, expected_qual):
    gq, qual = postprocess_variants.compute_quals(probs, call)
    self.assertEquals(gq, expected_gq)
    self.assertAlmostEquals(qual, expected_qual, places=6)

  @parameterized.parameters(
      # Make sure code is robust to minor numerical issues where the sum of
      # the vector isn't exactly 1.0.
      # This vector sums to ~1.0, minus any parsing uncertainty, and will
      # return a GQ of 40 but a qual of MAX_QUAL.
      ([0.0, 0.0001, 0.9999], 2, 40),
      # This vector sums to >1.0, but we should get back a max qual value
      # instead of throwing an exception.
      ([0.0, 0.00011, 0.9999], 2, 40),
  )
  def test_compute_quals_numerical_stability(self, probs, call, expected_gq):
    max_qual = round(
        genomics_math.ptrue_to_bounded_phred(1.0),
        postprocess_variants._QUAL_PRECISION)
    gq, qual = postprocess_variants.compute_quals(probs, call)
    self.assertEquals(expected_gq, gq)
    self.assertEquals(max_qual, qual)

  @parameterized.parameters(
      # Standard diploid case with 1 alt allele.
      ([1, 0, 0], 1, (0, [0, 0])),
      ([0, 1, 0], 1, (1, [0, 1])),
      ([0, 0, 1], 1, (2, [1, 1])),
      # Diploid case with 2 alt alleles.
      ([1, 0, 0, 0, 0, 0], 2, (0, [0, 0])),
      ([0, 1, 0, 0, 0, 0], 2, (1, [0, 1])),
      ([0, 0, 1, 0, 0, 0], 2, (2, [1, 1])),
      ([0, 0, 0, 1, 0, 0], 2, (3, [0, 2])),
      ([0, 0, 0, 0, 1, 0], 2, (4, [1, 2])),
      ([0, 0, 0, 0, 0, 1], 2, (5, [2, 2])),
  )
  def test_most_likely_genotype(self, probs, n_alts, expected):
    self.assertEqual(expected, postprocess_variants.most_likely_genotype(probs))

  @parameterized.parameters(1, 3, 4)
  def test_most_likely_genotype_raises_with_non2_ploidy(self, bad_ploidy):
    with self.assertRaises(NotImplementedError):
      postprocess_variants.most_likely_genotype([1, 0, 0], ploidy=bad_ploidy)

  def test_compute_filter_fields(self):
    # This generates too many tests as a parameterized test.
    for qual, min_qual in itertools.product(range(100), range(100)):
      # First test with no call and filter threshold
      variant = variants_pb2.Variant()
      variant.quality = qual
      expected = []
      expected.append(postprocess_variants.DEEP_VARIANT_PASS if qual >= min_qual
                      else postprocess_variants.DEEP_VARIANT_QUAL_FILTER)
      self.assertEqual(
          postprocess_variants.compute_filter_fields(variant, min_qual),
          expected)
      # Now add hom ref genotype --> qual shouldn't affect filter field
      del variant.filter[:]
      variant.calls.add(genotype=[0, 0])
      expected = []
      expected.append(postprocess_variants.DEEP_VARIANT_REF_FILTER)
      self.assertEqual(
          postprocess_variants.compute_filter_fields(variant, min_qual),
          expected)
      # Now add variant genotype --> qual filter should matter again
      del variant.filter[:]
      del variant.calls[:]
      variant.calls.add(genotype=[0, 1])
      expected = []
      expected.append(postprocess_variants.DEEP_VARIANT_PASS if qual >= min_qual
                      else postprocess_variants.DEEP_VARIANT_QUAL_FILTER)
      self.assertEqual(
          postprocess_variants.compute_filter_fields(variant, min_qual),
          expected)

  @parameterized.parameters(
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.01, 0.98, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT'])),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT'])),
              _create_call_variants_output(
                  indices=[2],
                  probabilities=[0.01, 0.97, 0.02],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT'])),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.01, 0.98, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT'])),
              _create_call_variants_output(
                  indices=[0, 2],
                  probabilities=[0.04, 0.95, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT'])),
              _create_call_variants_output(
                  indices=[1, 2],
                  probabilities=[0.01, 0.98, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT'])),
          ],
          6,
          set(['T']),
      ),
      # Example where all alt alleles are below qual_filter, but we keep one
      # where the qual is highest among the ones filtered out.
      (
          [
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T'])),
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.99, 0.01, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T'])),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T'])),
          ],
          6,
          set(['T']),
      ),
  )
  def test_get_alt_alleles_to_remove(self, call_variants_outputs, qual_filter,
                                     expected_output):
    self.assertEqual(
        postprocess_variants.get_alt_alleles_to_remove(
            call_variants_outputs, qual_filter), expected_output)

  @parameterized.parameters(
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set([]),
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['C']),
          _create_variant_with_alleles(alts=['T', 'TT']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['T']),
          _create_variant_with_alleles(alts=['C', 'TT']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['TT']),
          _create_variant_with_alleles(alts=['C', 'T']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['C', 'TT']),
          _create_variant_with_alleles(alts=['T']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['C', 'T', 'TT']),
          _create_variant_with_alleles(alts=[]),
      ),
  )
  def test_prune_alleles(self, canonical_variant, alt_alleles_to_remove,
                         expected_variant):
    self.assertEqual(
        postprocess_variants.prune_alleles(
            canonical_variant, alt_alleles_to_remove), expected_variant)

  @parameterized.parameters(
      # Check that we are simplifying alleles and that the simplification deps
      # on the alleles we've removed.
      dict(
          alleles=['CAA', 'CA', 'C'],
          start=5,
          alt_alleles_to_remove=[],
          expected_alleles=['CAA', 'CA', 'C'],
          expected_end=8),
      # Removing the C allele allows us to simplify CAA + CA => CA + C.
      dict(
          alleles=['CAA', 'CA', 'C'],
          start=4,
          alt_alleles_to_remove=['C'],
          expected_alleles=['CA', 'C'],
          expected_end=6),
      # Removing the CA allele doens't allow any simplification.
      dict(
          alleles=['CAA', 'CA', 'C'],
          start=3,
          alt_alleles_to_remove=['CA'],
          expected_alleles=['CAA', 'C'],
          expected_end=6),
      # Make sure we keep at least one anchor base when pruning.
      dict(
          alleles=['CCA', 'CA', 'T'],
          start=2,
          alt_alleles_to_remove=['T'],
          expected_alleles=['CC', 'C'],
          expected_end=4),
  )
  def test_simplify_alleles(self, alleles, start, alt_alleles_to_remove,
                            expected_alleles, expected_end):
    """Test that prune_alleles + simplify_alleles works as expected."""
    variant = _create_variant_with_alleles(
        ref=alleles[0], alts=alleles[1:], start=start)
    pruned = postprocess_variants.prune_alleles(variant, alt_alleles_to_remove)
    simplified = postprocess_variants.simplify_alleles(pruned)
    self.assertEqual(simplified.reference_bases, expected_alleles[0])
    self.assertEqual(simplified.alternate_bases, expected_alleles[1:])
    self.assertEqual(simplified.start, start)
    self.assertEqual(simplified.end, expected_end)

  def test_merge_predictions_simplifies_alleles(self):
    """Checks that merge_predictions simplifies alleles."""
    ref, alts = 'CCA', ['CA', 'C']
    inputs = [
        _create_call_variants_output(
            ref=ref, indices=[0], probabilities=[0.0, 1.0, 0.0], alts=alts),
        _create_call_variants_output(
            ref=ref, indices=[1], probabilities=[1.0, 0.0, 0.0], alts=alts),
        _create_call_variants_output(
            ref=ref, indices=[0, 1], probabilities=[0.0, 1.0, 0.0], alts=alts),
    ]

    for permuted_inputs in itertools.permutations(inputs):
      # qual_filter=2 is needed so we remove our middle 'C' allele.
      variant, probs = postprocess_variants.merge_predictions(
          permuted_inputs, qual_filter=2)
      np.testing.assert_almost_equal(probs, [0.0, 1.0, 0.0])
      self.assertEqual(variant.reference_bases, 'CC')
      self.assertEqual(variant.alternate_bases, ['C'])

  @parameterized.parameters(
      (['A'], {}, [1, 2], [1, 2]),
      (['A'], {'A'}, [1, 2], [1]),
      (['A', 'C'], {}, [1, 2, 3], [1, 2, 3]),
      (['A', 'C'], {'A'}, [1, 2, 3], [1, 3]),
      (['A', 'C'], {'C'}, [1, 2, 3], [1, 2]),
      (['A', 'C'], {'A', 'C'}, [1, 2, 3], [1]),
  )
  def test_prune_alleles_handles_format_fields(self, alts, to_remove, orig_ad,
                                               expected_ad):
    variant = _create_variant_with_alleles(alts=alts)
    test_utils.set_list_values(variant.calls[0].info['AD'], orig_ad)
    actual = postprocess_variants.prune_alleles(variant, to_remove)
    self.assertEqual(
        [v.number_value for v in actual.calls[0].info['AD'].values],
        expected_ad)

  @parameterized.parameters(
      (1, [[0]]),
      (2, [[0], [0, 1], [1]]),
      (3, [[0], [0, 1], [0, 2], [1], [1, 2], [2]]),
      (4, [[0], [0, 1], [0, 2], [0, 3], [1], [1, 2], [1, 3], [2], [2, 3], [3]]),
      (8, [[0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6], [0, 7], [1],
           [1, 2], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [2], [2, 3], [2, 4],
           [2, 5], [2, 6], [2, 7], [3], [3, 4], [3, 5], [3, 6], [3, 7], [4],
           [4, 5], [4, 6], [4, 7], [5], [5, 6], [5, 7], [6], [6, 7], [7]]),
  )
  def test_expected_alt_allele_indices(self, num_alternate_bases,
                                       expected_indices):
    self.assertEqual(
        postprocess_variants.expected_alt_allele_indices(num_alternate_bases),
        expected_indices)

  @flagsaver.FlagSaver
  def test_catches_bad_argv(self):
    # Define valid flags to ensure raise occurs due to argv issues.
    FLAGS.infile = make_golden_dataset(False)
    FLAGS.ref = test_utils.CHR20_FASTA
    FLAGS.outfile = test_utils.test_tmpfile('nonempty_outfile.vcf')
    with mock.patch.object(logging, 'error') as mock_logging,\
        mock.patch.object(sys, 'exit') as mock_exit:
      postprocess_variants.main(['postprocess_variants.py', 'extra_arg'])
    mock_logging.assert_called_once_with(
        'Command line parsing failure: postprocess_variants does not accept '
        'positional arguments but some are present on the command line: '
        '"[\'postprocess_variants.py\', \'extra_arg\']".')
    mock_exit.assert_called_once_with(errno.ENOENT)

  @flagsaver.FlagSaver
  def test_catches_bad_flags(self):
    FLAGS.infile = make_golden_dataset(False)
    FLAGS.ref = test_utils.CHR20_FASTA
    FLAGS.outfile = 'nonempty_outfile.vcf'
    FLAGS.nonvariant_site_tfrecord_path = (
        test_utils.GOLDEN_POSTPROCESS_GVCF_INPUT)
    # This is the bad flag.
    FLAGS.gvcf_outfile = ''
    with mock.patch.object(logging, 'error') as mock_logging, \
        mock.patch.object(sys, 'exit') as mock_exit:
      postprocess_variants.main(['postprocess_variants.py'])
    mock_logging.assert_called_once_with(
        'gVCF creation requires both nonvariant_site_tfrecord_path and '
        'gvcf_outfile flags to be set.')
    mock_exit.assert_called_once_with(errno.ENOENT)


class MergeVcfAndGvcfTest(parameterized.TestCase):

  @parameterized.parameters(
      # Smaller chromosome should return less than.
      (('1', 10, 20), ('2', 1, 2), True),
      (('2', 10, 20), ('10', 1, 2), True),
      # Larger chromosome should not return less than.
      (('2', 1, 2), ('1', 10, 20), False),
      (('10', 1, 2), ('2', 10, 20), False),
      # Same chromosome, smaller.
      (('1', 1, 2), ('1', 3, 4), True),
      (('1', 1, 3), ('1', 3, 4), True),
      # Same chromosome, overlapping.
      (('1', 1, 4), ('1', 3, 6), False),
      (('1', 1, 9), ('1', 3, 6), False),
      (('1', 4, 9), ('1', 3, 6), False),
      # Same chromosome, larger.
      (('1', 6, 9), ('1', 3, 4), False),
  )
  def test_lessthan_comparison(self, v1, v2, expected):
    variant1 = _create_nonvariant(*v1, ref_base='A')
    variant2 = _create_nonvariant(*v2, ref_base='C')
    lessthan = postprocess_variants._get_contig_based_lessthan(_CONTIGS)
    self.assertEqual(lessthan(variant1, variant2), expected)
    self.assertEqual(lessthan(variant1, None), True)
    self.assertEqual(lessthan(None, variant1), False)

  @parameterized.parameters(
      (_create_nonvariant('1', 10, 15, 'G'), 10, 12,
       _create_nonvariant('1', 10, 12, 'G')),
      (_create_nonvariant('2', 1, 9, 'C'), 3, 4,
       _create_nonvariant('2', 3, 4, 'G')),
  )
  def test_create_record_from_template(self, template, start, end, expected):
    reader = DummyReferenceReader()
    actual = postprocess_variants._create_record_from_template(
        template, start, end, reader)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      # One alt, with het GLs.
      (['C'], [-2.0457574905606752, -0.004364805402450088, -3.0], [0.75]),
      # Multi alts.
      (['G', 'C'], [
          -1.1368906918484387, -0.5279124552610386, -0.5923808731731073,
          -0.8155431286425007, -0.8415961054266092, -1.108308924501657
      ], [0.5, 0.1]),
      (['G', 'C', 'T'], [
          -0.7956722868920258, -0.663917423732382, -1.493986734511771,
          -0.8202531343562444, -0.9377869397242453, -1.0415699718993066,
          -1.4176189291054515, -1.5795151893394743, -1.8101482990393198,
          -0.8139951558313916
      ], [0.5, 0.1, 0.05]),
  )
  def test_transform_to_gvcf_add_allele(self, prior_alts, prior_gls, prior_vaf):
    variant = _create_variant(
        ref_name='chr1',
        start=10,
        ref_base='A',
        alt_bases=prior_alts,
        qual=40,
        filter_field='PASS',
        genotype=[0, 1],
        gq=None,
        likelihoods=prior_gls)
    prior_vaf_values = [struct_pb2.Value(number_value=v) for v in prior_vaf]
    variant.calls[0].info['VAF'].values.extend(prior_vaf_values)
    expected = _create_variant(
        ref_name='chr1',
        start=10,
        ref_base='A',
        alt_bases=prior_alts + [variantutils.GVCF_ALT_ALLELE],
        qual=40,
        filter_field='PASS',
        genotype=[0, 1],
        gq=None,
        likelihoods=prior_gls + ([postprocess_variants._GVCF_ALT_ALLELE_GL] *
                                 (len(prior_alts) + 2)))
    expected.calls[0].info['VAF'].values.extend(
        prior_vaf_values + [struct_pb2.Value(number_value=0)])
    actual = postprocess_variants._transform_to_gvcf_record(variant)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      # One alt, with het GLs.
      ([variantutils.GVCF_ALT_ALLELE],
       [-2.0457574905606752, -0.004364805402450088, -3.0], [0.5]),
      # Multi alts.
      (['G', variantutils.GVCF_ALT_ALLELE], [
          -1.1368906918484387, -0.5279124552610386, -0.5923808731731073,
          -0.8155431286425007, -0.8415961054266092, -1.108308924501657
      ], [0.5, 0.1]),
      (['G', 'C', variantutils.GVCF_ALT_ALLELE], [
          -0.7956722868920258, -0.663917423732382, -1.493986734511771,
          -0.8202531343562444, -0.9377869397242453, -1.0415699718993066,
          -1.4176189291054515, -1.5795151893394743, -1.8101482990393198,
          -0.8139951558313916
      ], [0, 0.5, 0.1]),
  )
  def test_transform_to_gvcf_no_allele_addition(self, alts, gls, vaf):
    variant = _create_variant(
        ref_name='chr1',
        start=10,
        ref_base='A',
        alt_bases=alts,
        qual=40,
        filter_field='PASS',
        genotype=[0, 1],
        gq=None,
        likelihoods=gls)
    vaf_values = [struct_pb2.Value(number_value=v) for v in vaf]
    variant.calls[0].info['VAF'].values.extend(vaf_values)
    expected = variants_pb2.Variant()
    expected.CopyFrom(variant)
    actual = postprocess_variants._transform_to_gvcf_record(variant)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      # Simple inputs.
      # Empty.
      ([], [], []),
      # Non-overlapping records.
      ([('1', 1, 'A')], [], [_simple_variant('1', 1, 'A')]),
      ([('1', 3, 'C'), ('1', 7, 'T'),
        ('2', 6, 'A')], [('2', 3, 6, 'G'), ('2', 7, 9, 'C')], [
            _simple_variant('1', 3, 'C'),
            _simple_variant('1', 7, 'T'),
            _create_nonvariant('2', 3, 6, 'G'),
            _simple_variant('2', 6, 'A'),
            _create_nonvariant('2', 7, 9, 'C')
        ]),
      # Non-variant record overlaps a variant from the left.
      ([('1', 5, 'GTTACG')], [('1', 2, 8, 'C')],
       [_create_nonvariant('1', 2, 5, 'C'),
        _simple_variant('1', 5, 'GTTACG')]),
      # Non-variant record overlaps a variant from the right.
      ([('1', 5, 'GTTACG')], [('1', 8, 15, 'A')], [
          _simple_variant('1', 5, 'GTTACG'),
          _create_nonvariant('1', 11, 15, 'T')
      ]),
      # Non-variant record is subsumed by a variant.
      ([('1', 5, 'GTTACG')], [('1', 5, 11, 'G')],
       [_simple_variant('1', 5, 'GTTACG')]),
      # Non-variant record subsumes a variant.
      ([('1', 5, 'GTTACG')], [('1', 4, 12, 'G')], [
          _create_nonvariant('1', 4, 5, 'G'),
          _simple_variant('1', 5, 'GTTACG'),
          _create_nonvariant('1', 11, 12, 'T')
      ]),
      # Non-variant record subsumes multiple overlapping variants.
      ([('1', 3, 'CGGTTAC'), ('1', 5, 'G')], [('1', 1, 15, 'A')], [
          _create_nonvariant('1', 1, 3, 'A'),
          _simple_variant('1', 3, 'CGGTTAC'),
          _simple_variant('1', 5, 'G'),
          _create_nonvariant('1', 10, 15, 'G')
      ]),
      ([('1', 3, 'CGGTTAC'), ('1', 5, 'GTTACGT')], [('1', 1, 15, 'A')], [
          _create_nonvariant('1', 1, 3, 'A'),
          _simple_variant('1', 3, 'CGGTTAC'),
          _simple_variant('1', 5, 'GTTACGT'),
          _create_nonvariant('1', 12, 15, 'T')
      ]),
  )
  def test_merge_variants_and_nonvariants(self, variants, nonvariants,
                                          expected):
    viter = (_simple_variant(*v) for v in variants)
    nonviter = (_create_nonvariant(*nv) for nv in nonvariants)
    lessthan = postprocess_variants._get_contig_based_lessthan(_CONTIGS)
    reader = DummyReferenceReader()
    actual = postprocess_variants.merge_variants_and_nonvariants(
        viter, nonviter, lessthan, reader)
    self.assertEqual(list(actual), expected)


if __name__ == '__main__':
  absltest.main()
