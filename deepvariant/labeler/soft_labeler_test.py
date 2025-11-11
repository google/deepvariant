# Copyright 2025 Google LLC.
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

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.labeler import soft_labeler
from third_party.nucleus.protos import variants_pb2


def _test_variant(chrom='20', start=10, alleles=('A', 'C'), gt=None):
  variant = variants_pb2.Variant(
      reference_name=chrom,
      start=start,
      end=start + len(alleles[0]),
      reference_bases=alleles[0],
      alternate_bases=alleles[1:],
  )

  if gt:
    variant.calls.add(genotype=gt)
  return variant


class ModifyTruthVariantTest(parameterized.TestCase):

  def test_modify_truth_variant_snp(self):
    ref = soft_labeler.ReferenceRegion('GATTACA', start=0)
    variant = _test_variant(chrom='1', start=1, alleles=['A', 'C'], gt=[0, 1])
    # SNPs are skipped.
    self.assertEmpty(list(soft_labeler.modify_truth_variant(variant, ref)))

  def test_modify_insertion(self):
    ref = soft_labeler.ReferenceRegion('GATTACA', start=0)
    # Insertion G -> GTT
    variant = _test_variant(chrom='1', start=0, alleles=['G', 'GTT'], gt=[0, 1])
    results = list(soft_labeler.modify_truth_variant(variant, ref))
    # i=-2: GTT -> GTT[:-2] = G
    # i=-1: GTT -> GTT[:-1] = GT
    # i=1: GTT -> GTT + 'T'*1 = GTTT
    # i=2: GTT -> GTT + 'T'*2 = GTTTT
    expected = [
        ('G', ('G',), 0.2),
        ('G', ('GT',), 0.1),
        ('G', ('GTTT',), 0.1),
        ('G', ('GTTTT',), 0.2),
    ]
    results_converted = []
    for variant, penalty in results:
      variant_alleles = tuple(variant.alternate_bases)
      results_converted.append(
          (variant.reference_bases, variant_alleles, penalty)
      )

    expected_converted = []
    for ref_bases, alt_bases, penalty in expected:
      expected_converted.append((ref_bases, tuple(alt_bases), penalty))

    self.assertCountEqual(results_converted, expected_converted)

  def test_modify_deletion(self):
    ref_seq = 'GATTACA'
    ref = soft_labeler.ReferenceRegion(ref_seq, start=0)
    # Deletion GATT -> G
    variant = _test_variant(
        chrom='1', start=0, alleles=['GATT', 'G'], gt=[0, 1]
    )
    # start=0, ref=GATT, end=4
    # i=-2: ref.end=4, bases_to_add=2. ref.bases(4, 6) = 'AC'
    #       ref_bases = GATT + AC = GATTAC
    # i=-1: ref.end=4, bases_to_add=1. ref.bases(4, 5) = 'A'
    #       ref_bases = GATT + A = GATTA
    # i=1: ref_bases[:-1] = GAT
    # i=2: ref_bases[:-2] = GA
    results = list(soft_labeler.modify_truth_variant(variant, ref))
    expected = [
        ('GATTAC', ('G',), 0.2),
        ('GATTA', ('G',), 0.1),
        ('GAT', ('G',), 0.1),
        ('GA', ('G',), 0.2),
    ]
    results_converted = []
    for variant, penalty in results:
      variant_alleles = tuple(variant.alternate_bases)
      results_converted.append(
          (variant.reference_bases, variant_alleles, penalty)
      )

    expected_converted = []
    for ref_bases, alt_bases, penalty in expected:
      expected_converted.append((ref_bases, tuple(alt_bases), penalty))

    self.assertCountEqual(results_converted, expected_converted)


class FindBestMatchingHaplotypesTest(parameterized.TestCase):

  def test_find_best_matching_haplotypes_no_variant(self):
    ref = soft_labeler.ReferenceRegion('GATTACA', start=0)
    result = soft_labeler.find_best_matching_haplotypes([], [], ref)
    self.assertNotEmpty(result)
    self.assertEqual(result[0].haplotypes, ['GATTACA'])
    self.assertEqual(result[0].candidate_genotypes, ())
    self.assertEqual(result[0].truth_genotypes, ())

  def test_find_best_matching_haplotypes_simple_snp(self):
    ref = soft_labeler.ReferenceRegion('GATTACA', start=0)
    candidates = [_test_variant(alleles=['A', 'C'], start=1)]
    truths = [_test_variant(alleles=['A', 'C'], start=1, gt=[0, 1])]
    result = soft_labeler.find_best_matching_haplotypes(candidates, truths, ref)
    self.assertNotEmpty(result)
    self.assertEqual(result[0].haplotypes, ['GATTACA', 'GCTTACA'])
    self.assertEqual(result[0].candidate_genotypes, ((0, 1),))
    self.assertEqual(result[0].truth_genotypes, ((0, 1),))

  def test_find_best_matching_haplotypes_matching_indel(self):
    ref = soft_labeler.ReferenceRegion('GATTACA', start=0)
    # Insertion G -> GTT
    candidates = [_test_variant(alleles=['G', 'GTT'], start=0)]
    truths = [_test_variant(alleles=['G', 'GTT'], start=0, gt=[0, 1])]
    result = soft_labeler.find_best_matching_haplotypes(candidates, truths, ref)
    self.assertNotEmpty(result)
    self.assertEqual(result[0].haplotypes, ['GATTACA', 'GTTATTACA'])
    self.assertEqual(result[0].candidate_genotypes, ((0, 1),))
    self.assertEqual(result[0].truth_genotypes, ((0, 1),))

  def test_find_best_matching_haplotypes_soft_indel(self):
    ref = soft_labeler.ReferenceRegion('GATTACA', start=0)
    # Insertion G -> GTT
    candidates = [_test_variant(alleles=['G', 'GTTT'], start=0)]
    truths = [_test_variant(alleles=['G', 'GTT'], start=0, gt=[0, 1])]
    result = soft_labeler.find_best_matching_haplotypes(candidates, truths, ref)
    self.assertNotEmpty(result)
    self.assertEqual(result[0].haplotypes, ['GATTACA', 'GTTTATTACA'])
    self.assertEqual(result[0].candidate_genotypes, ((0, 1),))
    self.assertEqual(result[0].truth_genotypes, ((0, 1),))
    self.assertEqual(result[0].truth_mod_penalties, [0.1])


if __name__ == '__main__':
  absltest.main()
