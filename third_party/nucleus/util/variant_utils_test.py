# Copyright 2018 Google Inc.
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

"""Tests for variant_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import functools
import itertools

from absl.testing import absltest
from absl.testing import parameterized
import mock

from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import struct_utils
from third_party.nucleus.util import variant_utils

NO_MISMATCH = set()
EVAL_DUP = variant_utils.AlleleMismatchType.duplicate_eval_alleles
TRUE_DUP = variant_utils.AlleleMismatchType.duplicate_true_alleles
TRUE_MISS = variant_utils.AlleleMismatchType.unmatched_true_alleles
EVAL_MISS = variant_utils.AlleleMismatchType.unmatched_eval_alleles


class VariantUtilsTests(parameterized.TestCase):

  def test_only_call(self):
    expected = variants_pb2.VariantCall(call_set_name='name', genotype=[0, 1])
    variant = variants_pb2.Variant(calls=[expected])
    actual = variant_utils.only_call(variant)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      0,
      2,
      3,
  )
  def test_invalid_only_call(self, num_calls):
    calls = [
        variants_pb2.VariantCall(call_set_name=str(x)) for x in range(num_calls)
    ]
    variant = variants_pb2.Variant(calls=calls)
    with self.assertRaisesRegexp(ValueError,
                                 'Expected exactly one VariantCall'):
      variant_utils.only_call(variant)

  def test_modify_only_call(self):
    variant = variants_pb2.Variant(calls=[variants_pb2.VariantCall()])
    call = variant_utils.only_call(variant)
    call.call_set_name = 'name'
    call.genotype[:] = [0, 1]
    self.assertLen(variant.calls, 1)
    self.assertEqual(variant.calls[0].call_set_name, 'name')
    self.assertEqual(variant.calls[0].genotype, [0, 1])

  def test_decode_variants(self):
    variants = [
        test_utils.make_variant(start=1),
        test_utils.make_variant(start=2)
    ]
    encoded = [variant.SerializeToString() for variant in variants]
    actual = variant_utils.decode_variants(encoded)
    # We have an iterable, so actual isn't equal to variants.
    self.assertNotEqual(actual, variants)
    # Making actual a list now makes it equal.
    self.assertEqual(list(actual), variants)

  def test_variant_position_and_range(self):
    v1 = test_utils.make_variant(chrom='1', alleles=['A', 'C'], start=10)
    v2 = test_utils.make_variant(chrom='1', alleles=['AGCT', 'C'], start=10)
    pos = ranges.make_range('1', 10, 11)
    range_ = ranges.make_range('1', 10, 14)
    v1_range_tuple = ('1', 10, 11)
    v2_range_tuple = ('1', 10, 14)
    self.assertEqual(pos, variant_utils.variant_position(v1))
    self.assertEqual(pos, variant_utils.variant_position(v2))
    self.assertEqual(pos, variant_utils.variant_range(v1))
    self.assertEqual(range_, variant_utils.variant_range(v2))
    self.assertEqual(v1_range_tuple, variant_utils.variant_range_tuple(v1))
    self.assertEqual(v2_range_tuple, variant_utils.variant_range_tuple(v2))

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), 'A/C'),
      (test_utils.make_variant(alleles=['A', 'C', 'T']), 'A/C,T'),
      (test_utils.make_variant(alleles=['A', 'AT']), 'A/AT'),
      (test_utils.make_variant(alleles=['AT', 'A']), 'AT/A'),
      (test_utils.make_variant(alleles=['AT', 'A', 'CT']), 'AT/A,CT'),
  )
  def test_format_alleles(self, variant, expected):
    self.assertEqual(variant_utils.format_alleles(variant), expected)

  @parameterized.parameters(
      (None, '.'),
      (['.'], '.'),
      (['PASS'], 'PASS'),
      (['FILTER1', 'FILTER2'], 'FILTER1,FILTER2'),
      (['FILTER1', 'FILTER2', 'FILTER3'], 'FILTER1,FILTER2,FILTER3'),
  )
  def test_format_filters(self, filters, expected):
    variant = test_utils.make_variant(filters=filters)
    if filters is None:
      variant.ClearField('filter')
    self.assertEqual(variant_utils.format_filters(variant), expected)

  @parameterized.parameters(
      # variant => status if we require non_ref genotype / status if we don't.
      (test_utils.make_variant(alleles=['A', 'C']), True, True),
      (test_utils.make_variant(alleles=['A', 'C'], gt=None), True, True),
      (test_utils.make_variant(alleles=['A', 'C', 'AT']), True, True),
      (test_utils.make_variant(alleles=['A']), False, False),
      (test_utils.make_variant(filters=['FAIL']), False, False),
      (test_utils.make_variant(gt=[-1, -1]), False, True),
      (test_utils.make_variant(gt=[0, 0]), False, True),
      (test_utils.make_variant(gt=[0, 1]), True, True),
      (test_utils.make_variant(gt=[1, 1]), True, True),
  )
  def test_is_variant_call(self, variant, expected_req_non_ref,
                           expected_any_genotype):
    # Check that default call checks for genotypes.
    self.assertEqual(
        variant_utils.is_variant_call(variant), expected_req_non_ref)
    # Ask explicitly for genotypes to be included.
    self.assertEqual(
        variant_utils.is_variant_call(variant, require_non_ref_genotype=True),
        expected_req_non_ref)
    # Don't require non_ref genotypes.
    self.assertEqual(
        variant_utils.is_variant_call(variant, require_non_ref_genotype=False),
        expected_any_genotype)

    with self.assertRaises(Exception):
      variant_utils.is_variant_call(None)

  def test_is_variant_call_no_calls_are_variant(self):

    def check_is_variant(variant, expected, **kwargs):
      self.assertEqual(
          variant_utils.is_variant_call(variant, **kwargs), expected)

    no_call = test_utils.make_variant(gt=[-1, -1])
    hom_ref = test_utils.make_variant(gt=[0, 0])
    het = test_utils.make_variant(gt=[0, 1])
    hom_var = test_utils.make_variant(gt=[1, 1])

    check_is_variant(no_call, False, no_calls_are_variant=False)
    check_is_variant(no_call, True, no_calls_are_variant=True)
    check_is_variant(hom_ref, False, no_calls_are_variant=False)
    check_is_variant(hom_ref, False, no_calls_are_variant=True)
    check_is_variant(het, True, no_calls_are_variant=False)
    check_is_variant(het, True, no_calls_are_variant=True)
    check_is_variant(hom_var, True, no_calls_are_variant=False)
    check_is_variant(hom_var, True, no_calls_are_variant=True)

  @parameterized.parameters(
      (test_utils.make_variant(filters=None), False),
      (test_utils.make_variant(filters=['.']), False),
      (test_utils.make_variant(filters=['PASS']), False),
      (test_utils.make_variant(filters=['FAIL']), True),
      (test_utils.make_variant(filters=['FAIL1', 'FAIL2']), True),
      # These two are not allowed in VCF, but worth testing our
      # code's behavior
      (test_utils.make_variant(filters=['FAIL1', 'PASS']), True),
      (test_utils.make_variant(filters=['FAIL1', '.']), True),
  )
  def test_is_filtered(self, variant, expected):
    self.assertEqual(variant_utils.is_filtered(variant), expected)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']),
       variant_utils.VariantType.snp),
      (test_utils.make_variant(alleles=['A', 'C', 'T']),
       variant_utils.VariantType.snp),
      (test_utils.make_variant(alleles=['A']), variant_utils.VariantType.ref),
      (test_utils.make_variant(alleles=['A', '.']),
       variant_utils.VariantType.ref),
      (test_utils.make_variant(alleles=['A', 'AC']),
       variant_utils.VariantType.indel),
      (test_utils.make_variant(alleles=['AC', 'A']),
       variant_utils.VariantType.indel),
      (test_utils.make_variant(alleles=['A', 'AC', 'ACC']),
       variant_utils.VariantType.indel),
      (test_utils.make_variant(alleles=['ACC', 'AC', 'A']),
       variant_utils.VariantType.indel),
  )
  def test_variant_type(self, variant, expected):
    self.assertEqual(variant_utils.variant_type(variant), expected)

  @parameterized.parameters(
      (test_utils.make_variant('chr1', 10), 'chr1:11'),
      (test_utils.make_variant('chr2', 100), 'chr2:101'),
  )
  def test_format_position(self, variant, expected):
    self.assertEqual(variant_utils.format_position(variant), expected)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), True),
      (test_utils.make_variant(alleles=['A', 'C', 'T']), True),
      (test_utils.make_variant(alleles=['A', 'AT']), False),
      (test_utils.make_variant(alleles=['AT', 'A']), False),
      (test_utils.make_variant(alleles=['AT', 'A', 'CT']), False),
      (test_utils.make_variant(alleles=['A', 'C', 'AT']), False),
      (test_utils.make_variant(alleles=['A']), False),
      (test_utils.make_variant(alleles=['A', '.']), False),
  )
  def test_is_snp(self, variant, expected):
    self.assertEqual(variant_utils.is_snp(variant), expected)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), False),
      (test_utils.make_variant(alleles=['A', 'C', 'T']), False),
      (test_utils.make_variant(alleles=['A', 'AT']), True),
      (test_utils.make_variant(alleles=['AT', 'A']), True),
      (test_utils.make_variant(alleles=['AT', 'A', 'CT']), True),
      (test_utils.make_variant(alleles=['A', 'C', 'AT']), True),
      (test_utils.make_variant(alleles=['A']), False),
      (test_utils.make_variant(alleles=['A', '.']), False),
  )
  def test_is_indel(self, variant, expected):
    self.assertEqual(variant_utils.is_indel(variant), expected)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), False),
      (test_utils.make_variant(alleles=['A', 'C', 'T']), True),
      (test_utils.make_variant(alleles=['A', 'AT']), False),
      (test_utils.make_variant(alleles=['AT', 'A']), False),
      (test_utils.make_variant(alleles=['AT', 'A', 'CT']), True),
      (test_utils.make_variant(alleles=['A', 'C', 'AT']), True),
  )
  def test_is_multiallelic(self, variant, expected):
    self.assertEqual(variant_utils.is_multiallelic(variant), expected)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), True),
      (test_utils.make_variant(alleles=['A', 'C', 'T']), False),
      (test_utils.make_variant(alleles=['A', 'AT']), True),
      (test_utils.make_variant(alleles=['AT', 'A']), True),
      (test_utils.make_variant(alleles=['AT', 'A', 'CT']), False),
      (test_utils.make_variant(alleles=['AT']), False),
  )
  def test_is_biallelic(self, variant, expected):
    self.assertEqual(variant_utils.is_biallelic(variant), expected)

  @parameterized.parameters(
      (['A', 'C'], ['A', 'C']),
      (['AA', 'CA'], ['A', 'C']),
      (['AAG', 'CAG'], ['A', 'C']),
      (['AAGAG', 'CAGAG'], ['A', 'C']),
      (['AACAG', 'CAGAG'], ['AAC', 'CAG']),
      (['AACAC', 'CAGAG'], ['AACAC', 'CAGAG']),
      (['ACT', 'A'], ['ACT', 'A']),
      (['ACTCT', 'ACT'], ['ACT', 'A']),
      (['ACTCT', 'A'], ['ACTCT', 'A']),
      (['CAG', 'GAG'], ['C', 'G']),
      # Make sure we don't reduce an allele to nothing.
      (['AT', 'ATAT'], ['A', 'ATA']),
      # Tests for multi-allelics.
      # There's one extra T here.
      (['ATT', 'AT', 'ATTT'], ['AT', 'A', 'ATT']),
      # Another single base postfix where we can remove a 'G'.
      (['CAG', 'GAG', 'TCG'], ['CA', 'GA', 'TC']),
      # There are two extra Ts to remove.
      (['ATTT', 'ATT', 'ATTTT'], ['AT', 'A', 'ATT']),
      # One pair can simplify, but not the other, so nothing can reduce.
      (['CAG', 'GAG', 'TCA'], ['CAG', 'GAG', 'TCA']),
      # Example from b/64022627.
      (['CGGCGG', 'CGG', 'CAACGG'], ['CGGC', 'C', 'CAAC']),
  )
  def test_simplify_alleles(self, alleles, expected):
    self.assertEqual(variant_utils.simplify_alleles(*alleles), tuple(expected))
    self.assertEqual(
        variant_utils.simplify_alleles(*reversed(alleles)),
        tuple(reversed(expected)))

  @parameterized.parameters(
      (['A', 'C'], ['A', 'C'], NO_MISMATCH),
      (['A', 'AC'], ['A', 'AC'], NO_MISMATCH),
      (['AC', 'A'], ['AC', 'A'], NO_MISMATCH),
      (['AC', 'A', 'ACT'], ['AC', 'A', 'ACT'], NO_MISMATCH),
      (['AC', 'A', 'ACT'], ['AC', 'ACT', 'A'], NO_MISMATCH),
      # Alleles are incompatible, so we have mismatches in both directions.
      (['A', 'C'], ['A', 'T'], {TRUE_MISS, EVAL_MISS}),
      (['A', 'C'], ['G', 'C'], {TRUE_MISS, EVAL_MISS}),
      # Missing alts specific to eval and truth.
      (['A', 'C', 'G'], ['A', 'C'], {EVAL_MISS}),
      (['A', 'C'], ['A', 'C', 'G'], {TRUE_MISS}),
      # Duplicate alleles.
      (['A', 'C', 'C'], ['A', 'C'], {EVAL_DUP}),
      (['A', 'C'], ['A', 'C', 'C'], {TRUE_DUP}),
      (['A', 'C', 'C'], ['A', 'C', 'C'], {EVAL_DUP, TRUE_DUP}),
      # Dups in truth, discordant alleles.
      (['A', 'C'], ['A', 'G', 'G'], {TRUE_DUP, EVAL_MISS, TRUE_MISS}),
      # Simplification of alleles does the right matching.
      (['A', 'C'], ['AA', 'CA'], NO_MISMATCH),  # trailing A.
      # preceding A, doesn't simplify so it's a mismatch.
      (['A', 'C'], ['AA', 'AC'], {EVAL_MISS, TRUE_MISS}),
      # both training preceding A, doesn't simplify, so mismatches
      (['A', 'C'], ['AAA', 'ACA'], {EVAL_MISS, TRUE_MISS}),
      # # Eval has 1 of the two alt alleles, so no eval mismatch.
      (['ACT', 'A'], ['ACTCT', 'ACT', 'A'], {TRUE_MISS}),
      # Eval has extra unmatched alleles, so it's got a mismatch.
      (['ACTCT', 'ACT', 'A'], ['ACT', 'A'], {EVAL_MISS}),
  )
  def test_allele_mismatch(self, a1, a2, expected):
    v1 = test_utils.make_variant(alleles=a1)
    v2 = test_utils.make_variant(alleles=a2)
    self.assertEqual(variant_utils.allele_mismatches(v1, v2), expected)

  @parameterized.parameters(
      (['A', 'C'], False),
      (['A', 'G'], True),
      (['A', 'T'], False),
      (['C', 'G'], False),
      (['C', 'T'], True),
      (['G', 'T'], False),
  )
  def test_is_transition(self, ordered_alleles, expected):
    for alleles in [ordered_alleles, reversed(ordered_alleles)]:
      self.assertEqual(variant_utils.is_transition(*alleles), expected)

  def test_is_transition_raises_with_bad_args(self):
    with self.assertRaises(ValueError):
      variant_utils.is_transition('A', 'A')
    with self.assertRaises(ValueError):
      variant_utils.is_transition('A', 'AA')
    with self.assertRaises(ValueError):
      variant_utils.is_transition('AA', 'A')

  @parameterized.parameters(
      # alleles followed by is_insertion and is_deletion expectation
      (['A', 'C'], False, False),
      (['A', 'AT'], True, False),
      (['A', 'ATT'], True, False),
      (['AT', 'A'], False, True),
      (['ATT', 'A'], False, True),
      (['CAT', 'TCA'], False, False),

      # These are examples where ref is not simplified, such as could occur
      # a multi-allelic record, such as the following:
      # alleles = AT, A, ATT, CT (1 deletion, 1 insertion, 1 SNP)
      (['AT', 'A'], False, True),
      (['AT', 'ATT'], True, False),
      (['AT', 'CT'], False, False),
  )
  def test_is_insertion_deletion(self, alleles, is_insertion, is_deletion):
    self.assertEqual(variant_utils.is_insertion(*alleles), is_insertion)
    self.assertEqual(variant_utils.is_deletion(*alleles), is_deletion)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), False, False),
      (test_utils.make_variant(alleles=['A', 'C', 'T']), False, False),
      (test_utils.make_variant(alleles=['A', 'AT']), True, False),
      (test_utils.make_variant(alleles=['AT', 'A']), False, True),
      (test_utils.make_variant(alleles=['AT', 'A', 'ATT']), True, True),
      (test_utils.make_variant(alleles=['AT', 'A', 'CT']), False, True),
      (test_utils.make_variant(alleles=['A', 'C', 'AT']), True, False),
      (test_utils.make_variant(alleles=['A']), False, False),
      (test_utils.make_variant(alleles=['A', '.']), False, False),
  )
  def test_has_insertion_deletion(self, variant, has_insertion, has_deletion):
    self.assertEqual(variant_utils.has_insertion(variant), has_insertion)
    self.assertEqual(variant_utils.has_deletion(variant), has_deletion)

  @parameterized.parameters(
      (test_utils.make_variant(gt=None), False),
      (test_utils.make_variant(gt=[0, 0]), True),
      (test_utils.make_variant(gt=[0, 1]), True),
      (test_utils.make_variant(gt=[1, 1]), True),
      (test_utils.make_variant(gt=[-1, -1]), True),
      (variants_pb2.Variant(calls=[]), False),
      (variants_pb2.Variant(
          calls=[variants_pb2.VariantCall(call_set_name='no_geno')]), True),
      (variants_pb2.Variant(calls=[
          variants_pb2.VariantCall(call_set_name='no_geno'),
          variants_pb2.VariantCall(call_set_name='no_geno2'),
      ]), True),
  )
  def test_has_calls(self, variant, expected):
    self.assertEqual(variant_utils.has_calls(variant), expected)

  def test_has_calls_raises_with_bad_inputs(self):
    with self.assertRaises(Exception):
      variant_utils.has_calls(None)

  @parameterized.parameters(
      (test_utils.make_variant(gt=None), variant_utils.GenotypeType.no_call),
      (test_utils.make_variant(gt=[-1, -1]),
       variant_utils.GenotypeType.no_call),
      (test_utils.make_variant(gt=[0, 0]), variant_utils.GenotypeType.hom_ref),
      (test_utils.make_variant(gt=[0, 1]), variant_utils.GenotypeType.het),
      (test_utils.make_variant(gt=[1, 0]), variant_utils.GenotypeType.het),
      (test_utils.make_variant(gt=[0, 2]), variant_utils.GenotypeType.het),
      (test_utils.make_variant(gt=[2, 0]), variant_utils.GenotypeType.het),
      (test_utils.make_variant(gt=[1, 1]), variant_utils.GenotypeType.hom_var),
      (test_utils.make_variant(gt=[1, 2]), variant_utils.GenotypeType.het),
  )
  def test_genotype_type(self, variant, expected):
    self.assertEqual(variant_utils.genotype_type(variant), expected)

  def test_genotype_type_raises_with_bad_args(self):
    with self.assertRaises(Exception):
      variant_utils.genotype_type(None)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C'], gt=[0, 0]), ['A', 'A']),
      (test_utils.make_variant(alleles=['A', 'C'], gt=[0, 1]), ['A', 'C']),
      (test_utils.make_variant(alleles=['A', 'C'], gt=[1, 0]), ['C', 'A']),
      (test_utils.make_variant(alleles=['A', 'C'], gt=[1, 1]), ['C', 'C']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[0, 0]), ['A', 'A']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[0, 1]), ['A', 'C']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[0, 2]), ['A', 'T']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[1, 2]), ['C', 'T']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[2, 1]), ['T', 'C']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[1, 1]), ['C', 'C']),
      (test_utils.make_variant(alleles=['A', 'C', 'T'], gt=[2, 2]), ['T', 'T']),
      (test_utils.make_variant(alleles=['A', 'C'], gt=[-1, -1]), ['.', '.']),
  )
  def test_genotype_as_alleles(self, variant, expected):
    self.assertEqual(variant_utils.genotype_as_alleles(variant), expected)

  def test_genotype_as_alleles_raises_with_bad_inputs(self):
    with self.assertRaises(Exception):
      variant_utils.genotype_as_alleles(None)
    with self.assertRaises(Exception):
      variant_utils.genotype_as_alleles(test_utils.make_variant(gt=None))
    with self.assertRaises(Exception):
      variant_utils.genotype_as_alleles(
          test_utils.make_variant(alleles=['A', 'C'], gt=[0, 0]), call_ix=1)
    with self.assertRaises(Exception):
      variant_utils.genotype_type(None)

  @parameterized.parameters(
      # Ref without an alt isn't gVCF.
      (test_utils.make_variant(alleles=['A']), False),
      # SNPs and indels aren't gVCF records.
      (test_utils.make_variant(alleles=['A', 'T']), False),
      (test_utils.make_variant(alleles=['A', 'AT']), False),
      (test_utils.make_variant(alleles=['AT', 'T']), False),
      # These are gVCF records.
      (test_utils.make_variant(alleles=['A', '<*>']), True),
      (test_utils.make_variant(alleles=['A', '<*>'], filters='PASS'), True),
      (test_utils.make_variant(alleles=['A', '<*>'], filters='FAIL'), True),
      # These are close but not exactly gVCFs.
      (test_utils.make_variant(alleles=['A', '<*>', 'C']), False),
      (test_utils.make_variant(alleles=['A', '<*F>']), False),
      (test_utils.make_variant(alleles=['A', '<CNV>']), False),
  )
  def test_is_gvcf(self, variant, expected):
    self.assertEqual(variant_utils.is_gvcf(variant), expected)

  @parameterized.parameters(
      # Variants with one ref and one alt allele.
      (test_utils.make_variant(alleles=['A', 'C']), [(0, 0, 'A', 'A'),
                                                     (0, 1, 'A', 'C'),
                                                     (1, 1, 'C', 'C')]),
      # Variants with one ref and two alt alleles.
      (test_utils.make_variant(alleles=['A', 'C', 'G']), [(0, 0, 'A', 'A'),
                                                          (0, 1, 'A', 'C'),
                                                          (1, 1, 'C', 'C'),
                                                          (0, 2, 'A', 'G'),
                                                          (1, 2, 'C', 'G'),
                                                          (2, 2, 'G', 'G')]),
      # Variants with one ref and three alt alleles.
      (test_utils.make_variant(alleles=['A', 'C', 'G', 'T']),
       [(0, 0, 'A', 'A'), (0, 1, 'A', 'C'), (1, 1, 'C', 'C'), (0, 2, 'A', 'G'),
        (1, 2, 'C', 'G'), (2, 2, 'G', 'G'), (0, 3, 'A', 'T'), (1, 3, 'C', 'T'),
        (2, 3, 'G', 'T'), (3, 3, 'T', 'T')]),
  )
  def test_genotype_ordering_in_likelihoods(self, variant, expected):
    self.assertEqual(
        list(variant_utils.genotype_ordering_in_likelihoods(variant)), expected)

  @parameterized.parameters(
      # Haploid.
      dict(gls=[0.], allele_indices=[0], expected=0.),
      dict(gls=[-1, -2], allele_indices=[1], expected=-2),
      dict(gls=[-1, -2, -3], allele_indices=[2], expected=-3),
      # Diploid.
      dict(gls=[0.], allele_indices=[0, 0], expected=0.),
      dict(gls=[-1, -2, -3], allele_indices=[0, 0], expected=-1),
      dict(gls=[-1, -2, -3], allele_indices=[0, 1], expected=-2),
      dict(gls=[-1, -2, -3], allele_indices=[1, 0], expected=-2),
      dict(gls=[-1, -2, -3], allele_indices=[1, 1], expected=-3),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[0, 0], expected=-1),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[0, 1], expected=-2),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[1, 0], expected=-2),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[1, 1], expected=-3),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[0, 2], expected=-4),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[2, 0], expected=-4),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[1, 2], expected=-5),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[2, 1], expected=-5),
      dict(gls=[-1, -2, -3, -4, -5, -6], allele_indices=[2, 2], expected=-6),
      dict(gls=range(10), allele_indices=[0, 3], expected=6),
      dict(gls=range(10), allele_indices=[1, 3], expected=7),
      dict(gls=range(10), allele_indices=[2, 3], expected=8),
      dict(gls=range(10), allele_indices=[3, 3], expected=9),
  )
  def test_genotype_likelihood(self, gls, allele_indices, expected):
    variantcall = variants_pb2.VariantCall(genotype_likelihood=gls)
    actual = variant_utils.genotype_likelihood(variantcall, allele_indices)
    self.assertEqual(actual, expected)

  def test_unsupported_genotype_likelihood(self):
    variantcall = variants_pb2.VariantCall(genotype_likelihood=[-1, -2, -3])
    with self.assertRaisesRegexp(NotImplementedError,
                                 'only supports haploid and diploid'):
      variant_utils.genotype_likelihood(variantcall, [0, 1, 1])

  def test_haploid_allele_indices_for_genotype_likelihood_index(self):
    for aix in xrange(20):
      allele_indices = (aix,)
      ix = variant_utils.genotype_likelihood_index(allele_indices)
      actual = variant_utils.allele_indices_for_genotype_likelihood_index(
          ix, ploidy=1)
      self.assertEqual(actual, aix)

  def test_diploid_allele_indices_for_genotype_likelihood_index(self):
    for aix in xrange(20):
      for bix in xrange(20):
        allele_indices = (aix, bix)
        expected = tuple(sorted(allele_indices))
        ix = variant_utils.genotype_likelihood_index(allele_indices)
        actual = variant_utils.allele_indices_for_genotype_likelihood_index(
            ix, ploidy=2)
        self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(ploidy=-1),
      dict(ploidy=0),
      dict(ploidy=3),
  )
  def test_unsupported_allele_indices_for_genotype_likelihood_index(
      self, ploidy):
    with self.assertRaisesRegexp(NotImplementedError,
                                 'only supported for haploid and diploid'):
      variant_utils.allele_indices_for_genotype_likelihood_index(0, ploidy)

  @parameterized.parameters(
      dict(alt_bases=[], num_alts=0, expected=[(0, 0)]),
      dict(alt_bases=['A'], num_alts=0, expected=[(0, 0)]),
      dict(alt_bases=['A'], num_alts=1, expected=[(0, 1)]),
      dict(alt_bases=['A'], num_alts=2, expected=[(1, 1)]),
      dict(alt_bases=['A', 'C'], num_alts=0, expected=[(0, 0)]),
      dict(alt_bases=['A', 'C'], num_alts=1, expected=[(0, 1), (0, 2)]),
      dict(alt_bases=['A', 'C'], num_alts=2, expected=[(1, 1), (1, 2), (2, 2)]),
  )
  def test_allele_indices_with_num_alts(self, alt_bases, num_alts, expected):
    variant = variants_pb2.Variant(alternate_bases=alt_bases)
    actual = variant_utils.allele_indices_with_num_alts(
        variant, num_alts, ploidy=2)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      dict(alt_bases=['A'], num_alts=0, ploidy=1),
      dict(alt_bases=['A'], num_alts=0, ploidy=3),
      dict(alt_bases=['A'], num_alts=-1, ploidy=2),
      dict(alt_bases=['A'], num_alts=3, ploidy=2),
  )
  def test_invalid_allele_indices_with_num_alts(self, alt_bases, num_alts,
                                                ploidy):
    variant = variants_pb2.Variant(alternate_bases=alt_bases)
    with self.assertRaises((NotImplementedError, ValueError)):
      variant_utils.allele_indices_with_num_alts(variant, num_alts, ploidy)

  def test_variants_overlap(self):
    v1 = test_utils.make_variant(chrom='1', alleles=['A', 'C'], start=10)
    v2 = test_utils.make_variant(chrom='1', alleles=['A', 'C'], start=20)
    with mock.patch.object(ranges, 'ranges_overlap') as mock_overlap:
      mock_overlap.return_value = 'SENTINEL'
      self.assertEqual(variant_utils.variants_overlap(v1, v2), 'SENTINEL')
      mock_overlap.assert_called_once_with(
          variant_utils.variant_range(v1), variant_utils.variant_range(v2))

  @parameterized.parameters(
      # Degenerate cases - no and one variant.
      dict(sorted_variants=[],),
      dict(sorted_variants=[
          test_utils.make_variant(chrom='1', start=10),
      ],),
      # Two variants on the same chromosome.
      dict(
          sorted_variants=[
              test_utils.make_variant(chrom='1', start=10),
              test_utils.make_variant(chrom='1', start=15),
          ],),
      # The first variant has start > the second, but it's on a later chrom.
      dict(
          sorted_variants=[
              test_utils.make_variant(chrom='1', start=15),
              test_utils.make_variant(chrom='2', start=10),
          ],),
      # Make sure the end is respected.
      dict(
          sorted_variants=[
              test_utils.make_variant(chrom='1', start=10),
              test_utils.make_variant(chrom='1', start=15),
              test_utils.make_variant(chrom='1', alleles=['AA', 'A'], start=15),
          ],),
      # Complex example with multiple chromosomes, ends, etc.
      dict(
          sorted_variants=[
              test_utils.make_variant(chrom='1', start=10),
              test_utils.make_variant(chrom='2', start=5),
              test_utils.make_variant(chrom='2', alleles=['AA', 'A'], start=5),
              test_utils.make_variant(chrom='2', start=6),
              test_utils.make_variant(chrom='2', start=10),
              test_utils.make_variant(chrom='3', start=2),
          ],),
  )
  def test_sorted_variants(self, sorted_variants):
    for permutation in itertools.permutations(
        sorted_variants, r=len(sorted_variants)):

      # Check that sorting the permutations produced sorted.
      self.assertEqual(
          variant_utils.sorted_variants(permutation), sorted_variants)

      # Check that variants_are_sorted() is correct, which we detect if
      # the range_tuples of permutation == the range_tuples of sorted_variants.
      def _range_tuples(variants):
        return [variant_utils.variant_range_tuple(v) for v in variants]

      self.assertEqual(
          variant_utils.variants_are_sorted(permutation),
          _range_tuples(permutation) == _range_tuples(sorted_variants))

  @parameterized.parameters(
      dict(
          variant=test_utils.make_variant(
              chrom='1', start=10, alleles=['A', 'C']),
          expected_key='1:11:A->C'),
      dict(
          variant=test_utils.make_variant(
              chrom='1', start=10, alleles=['A', 'G', 'C']),
          sort_alleles=True,
          expected_key='1:11:A->C/G'),
      dict(
          variant=test_utils.make_variant(
              chrom='1', start=10, alleles=['A', 'G', 'C']),
          sort_alleles=False,
          expected_key='1:11:A->G/C'),
  )
  def test_variant_key(self, variant, expected_key, sort_alleles=True):
    self.assertEqual(
        variant_utils.variant_key(variant, sort_alleles=sort_alleles),
        expected_key)

  @parameterized.parameters(
      dict(
          field_name='AD',
          value=[23, 25],
          reader=None,
          expected=[
              struct_pb2.Value(int_value=23),
              struct_pb2.Value(int_value=25)
          ],
      ),
      dict(
          field_name='AA',
          value='C',
          reader=True,
          expected=[struct_pb2.Value(string_value='C')],
      ),
  )
  def test_set_info(self, field_name, value, reader, expected):
    if reader is not None:
      reader = mock.Mock()
      reader.field_access_cache.info_field_set_fn.return_value = (
          struct_utils.set_string_field)
    variant = variants_pb2.Variant()
    variant_utils.set_info(variant, field_name, value, reader)
    actual = variant.info[field_name].values
    self.assertEqual(len(actual), len(expected))
    for actual_elem, expected_elem in zip(actual, expected):
      self.assertEqual(actual_elem, expected_elem)

  @parameterized.parameters(
      dict(field_name='AD', reader=None, expected=[23, 25]),
      dict(field_name='AA', reader=True, expected='C'),
      dict(field_name='1000G', reader=None, expected=True),
  )
  def test_get_info(self, field_name, reader, expected):
    if reader is not None:
      reader = mock.Mock()
      reader.field_access_cache.info_field_get_fn.return_value = (
          functools.partial(
              struct_utils.get_string_field, is_single_field=True))
    variant = variants_pb2.Variant()
    variant_utils.set_info(variant, 'AD', [23, 25])
    variant_utils.set_info(variant, 'AA', 'C')
    variant_utils.set_info(variant, '1000G', True)
    variant_utils.set_info(variant, 'DB', False)
    actual = variant_utils.get_info(variant, field_name, vcf_object=reader)
    self.assertEqual(actual, expected)


if __name__ == '__main__':
  absltest.main()
