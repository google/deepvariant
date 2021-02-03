# Copyright 2020 Google LLC.
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
"""Tests for deepvariant.allele_frequency."""

from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import fasta
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import variants_pb2
from deepvariant import allele_frequency
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2


def setUpModule():
  testdata.init()


class AlleleFrequencyTest(parameterized.TestCase):

  @parameterized.parameters(
      # A SNP.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60168,
              end=60169,
              reference_bases='C',
              alternate_bases=['T']),
          reference_haplotype='GCACCT',
          reference_offset=60165,
          expected_return=[{
              'haplotype':
                  'GCATCT',
              'alt':
                  'T',
              'variant':
                  variants_pb2.Variant(
                      reference_name='chr20',
                      start=60168,
                      end=60169,
                      reference_bases='C',
                      alternate_bases=['T'])
          }]),
      # A deletion.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60284,
              end=60291,
              reference_bases='ATTCCAG',
              alternate_bases=['AT']),
          reference_haplotype='TTTCCATTCCAGTCCAT',
          reference_offset=60279,
          expected_return=[{
              'haplotype':
                  'TTTCCATTCCAT',
              'alt':
                  'AT',
              'variant':
                  variants_pb2.Variant(
                      reference_name='chr20',
                      start=60284,
                      end=60291,
                      reference_bases='ATTCCAG',
                      alternate_bases=['AT'])
          }]),
      # An insertion.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60279,
              end=60285,
              reference_bases='TTTCCA',
              alternate_bases=['TTTCCATTCCA']),
          reference_haplotype='TTTCCATTCCAGTCCAT',
          reference_offset=60279,
          expected_return=[{
              'haplotype':
                  'TTTCCATTCCATTCCAGTCCAT',
              'alt':
                  'TTTCCATTCCA',
              'variant':
                  variants_pb2.Variant(
                      reference_name='chr20',
                      start=60279,
                      end=60285,
                      reference_bases='TTTCCA',
                      alternate_bases=['TTTCCATTCCA'])
          }]),
      # A deletion.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60284,
              end=60291,
              reference_bases='ATTCCAG',
              alternate_bases=['AT']),
          reference_haplotype='TTTCCATTCCAG',
          reference_offset=60279,
          expected_return=[{
              'haplotype':
                  'TTTCCAT',
              'alt':
                  'AT',
              'variant':
                  variants_pb2.Variant(
                      reference_name='chr20',
                      start=60284,
                      end=60291,
                      reference_bases='ATTCCAG',
                      alternate_bases=['AT'])
          }]),
      # An insertion.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60279,
              end=60285,
              reference_bases='TTTCCA',
              alternate_bases=['TTTCCATTCCA']),
          reference_haplotype='TTTCCATTCCAG',
          reference_offset=60279,
          expected_return=[{
              'haplotype':
                  'TTTCCATTCCATTCCAG',
              'alt':
                  'TTTCCATTCCA',
              'variant':
                  variants_pb2.Variant(
                      reference_name='chr20',
                      start=60279,
                      end=60285,
                      reference_bases='TTTCCA',
                      alternate_bases=['TTTCCATTCCA'])
          }]))
  def test_update_haplotype(self, variant, reference_haplotype,
                            reference_offset, expected_return):
    list_hap_obj = allele_frequency.update_haplotype(variant,
                                                     reference_haplotype,
                                                     reference_offset)
    self.assertListEqual(list_hap_obj, expected_return)

  @parameterized.parameters([
      dict(
          dv_variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60284,
              end=60291,
              reference_bases='ATTCCAG',
              alternate_bases=['AT']),
          cohort_variants=[
              variants_pb2.Variant(
                  reference_name='chr20',
                  start=60279,
                  end=60285,
                  reference_bases='TTTCCA',
                  alternate_bases=['T', 'TTTCCATTCCA']),
              variants_pb2.Variant(
                  reference_name='chr20',
                  start=60285,
                  end=60291,
                  reference_bases='TTTCCA',
                  alternate_bases=['T']),
          ],
          expected_ref_haplotype='TTTCCATTCCAG',
          expected_ref_offset=60279)
  ])
  def test_get_ref_haplotype_and_offset(self, dv_variant, cohort_variants,
                                        expected_ref_haplotype,
                                        expected_ref_offset):
    ref_reader = fasta.IndexedFastaReader(testdata.GRCH38_FASTA)
    ref_haplotype, ref_offset = allele_frequency.get_ref_haplotype_and_offset(
        dv_variant, cohort_variants, ref_reader)
    self.assertEqual(ref_haplotype, expected_ref_haplotype)
    self.assertEqual(ref_offset, expected_ref_offset)

  @parameterized.parameters(
      # A matched SNP.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60168,
              end=60169,
              reference_bases='C',
              alternate_bases=['T']),
          expected_return=dict(C=0.9998, T=0.0002),
          label='matched_snp_1'),
      # A matched deletion.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60285,
              end=60291,
              reference_bases='TTCCAG',
              alternate_bases=['T']),
          expected_return=dict(T=0.001198, TTCCAG=0.998802),
          label='matched_del_1'),
      # A unmatched deletion.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60284,
              end=60291,
              reference_bases='ATTCCAG',
              alternate_bases=['A']),
          expected_return=dict(A=0, ATTCCAG=1),
          label='unmatched_del_1'),
      # A matched deletion, where the candidate is formatted differently.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60284,
              end=60291,
              reference_bases='ATTCCAG',
              alternate_bases=['AT']),
          expected_return=dict(AT=0.001198, ATTCCAG=0.998802),
          label='matched_del_2: diff representation'),
      # An unmatched SNP.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60150,
              end=60151,
              reference_bases='C',
              alternate_bases=['T']),
          expected_return=dict(C=1, T=0),
          label='unmatched_snp_1'),
      # A matched SNP and an unmatched SNP.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60168,
              end=60169,
              reference_bases='C',
              alternate_bases=['T', 'A']),
          expected_return=dict(C=0.9998, T=0.0002, A=0),
          label='mixed_snp_1'),
      # An unmatched SNP, where the REF allele frequency is not 1.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60168,
              end=60169,
              reference_bases='C',
              alternate_bases=['A']),
          expected_return=dict(C=0.9998, A=0),
          label='unmatched_snp_2: non-1 ref allele'),
      # A multi-allelic candidate at a multi-allelic locus.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60279,
              end=60285,
              reference_bases='TTTCCA',
              alternate_bases=['T', 'TTTCCATTCCA']),
          expected_return=dict(TTTCCA=0.999401, T=0.000399, TTTCCATTCCA=0.0002),
          label='matched_mult_1'),
      # A multi-allelic candidate at a multi-allelic locus.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60279,
              end=60285,
              reference_bases='TTTCCA',
              alternate_bases=['T', 'TATCCATTCCA']),
          expected_return=dict(TTTCCA=0.999401, T=0.000399, TATCCATTCCA=0),
          label='unmatched_mult_1'),
      # [Different representation]
      # A deletion where the cohort variant is represented differently.
      # In this case, REF frequency is calculated by going over all cohort ALTs.
      # Thus, the sum of all dict values is not equal to 1.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60295,
              end=60301,
              reference_bases='TTCCAT',
              alternate_bases=['T']),
          expected_return=dict(T=0.000399, TTCCAT=0.923922),
          label='matched_del_3: diff representation'),
      # [Non-candidate allele]
      # One allele of a multi-allelic cohort variant is not in candidate.
      # The non-candidate allele should be ignored.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=60279,
              end=60285,
              reference_bases='TTTCCA',
              alternate_bases=['T']),
          expected_return=dict(TTTCCA=0.999401, T=0.000399),
          label='matched_del_4: multi-allelic cohort'),
      # A left-align example.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=9074790,
              end=9074794,
              reference_bases='CT',
              alternate_bases=['C', 'CTTT']),
          expected_return=dict(C=0.167732, CTTT=0.215256, CT=0.442092),
          label='matched_mult_2: left align'),
      # A left-align example.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=9074790,
              end=9074794,
              reference_bases='C',
              alternate_bases=['CTTT']),
          expected_return=dict(CTTT=0.145367, C=0.442092),
          label='matched_ins_1: left align'),
      # A left-align example.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=9074790,
              end=9074793,
              reference_bases='CTT',
              alternate_bases=['CTTA']),
          expected_return=dict(CTTA=0, CTT=0.442092),
          label='unmatched_ins_1: left align'),
      # A matched mnps.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=61065,
              end=61066,
              reference_bases='T',
              alternate_bases=['C']),
          expected_return=dict(C=0.079872, T=0.919729),
          label='matched_mnps_1'),
      # A matched SNP.
      dict(
          variant=variants_pb2.Variant(
              reference_name='chr20',
              start=62022,
              end=62023,
              reference_bases='G',
              alternate_bases=['C', 'T']),
          expected_return=dict(G=0.996206, C=0.003594, T=0),
          label='matched_snp_2'))
  def test_find_matching_allele_frequency(self, variant, expected_return,
                                          label):
    ref_reader = fasta.IndexedFastaReader(testdata.GRCH38_FASTA)
    vcf_reader = vcf.VcfReader(testdata.VCF_WITH_ALLELE_FREQUENCIES)
    allele_frequencies = allele_frequency.find_matching_allele_frequency(
        variant, vcf_reader, ref_reader)
    # Compare keys.
    self.assertSetEqual(
        set(allele_frequencies.keys()), set(expected_return.keys()), msg=label)
    # Compare values (almost equal).
    for key in allele_frequencies.keys():
      self.assertAlmostEqual(
          allele_frequencies[key], expected_return[key], msg=label)

  def test_make_population_vcf_readers_with_multiple_vcfs(self):
    filenames = [testdata.AF_VCF_CHR20, testdata.AF_VCF_CHR21]

    output = allele_frequency.make_population_vcf_readers(filenames)

    self.assertIsInstance(output['chr20'], vcf.VcfReader)
    self.assertIsInstance(output['chr21'], vcf.VcfReader)
    self.assertEqual(next(output['chr20']).reference_name, 'chr20')
    self.assertEqual(next(output['chr21']).reference_name, 'chr21')
    # Check that chr22 has no reader rather than outputting another reader for
    # a different chromosome.
    self.assertIsNone(output['chr22'])

  def test_make_population_vcf_readers_with_one_vcf(self):
    filenames = [testdata.AF_VCF_CHR20_AND_21]

    output = allele_frequency.make_population_vcf_readers(filenames)

    self.assertIsInstance(output['chr20'], vcf.VcfReader)
    self.assertIsInstance(output['chr21'], vcf.VcfReader)
    self.assertIsInstance(output['chr22'], vcf.VcfReader)
    # All reference names should map to the same VCF that starts with chr20.
    self.assertEqual(next(output['chr20']).reference_name, 'chr20')
    self.assertEqual(next(output['chr21']).reference_name, 'chr20')
    self.assertEqual(next(output['chr22']).reference_name, 'chr20')

  def test_make_population_vcf_readers_raises_on_shared_chromosomes(self):
    filenames = [
        testdata.AF_VCF_CHR20, testdata.AF_VCF_CHR21,
        testdata.AF_VCF_CHR20_AND_21
    ]

    with self.assertRaisesRegex(
        expected_exception=ValueError,
        expected_regex='Variants on chr20 are included in multiple VCFs'):
      allele_frequency.make_population_vcf_readers(filenames)

  @parameterized.parameters(
      dict(
          dv_calls=iter([
              deepvariant_pb2.DeepVariantCall(
                  variant=variants_pb2.Variant(
                      reference_name='chr20',
                      start=60168,
                      end=60169,
                      reference_bases='C',
                      alternate_bases=['T']),
                  allele_support=None)
          ]),
          expected_return=dict(C=0.9998, T=0.0002),
          testcase='valid'),
      dict(
          dv_calls=iter([
              deepvariant_pb2.DeepVariantCall(
                  variant=variants_pb2.Variant(
                      reference_name='chrM',
                      start=10000,
                      end=10001,
                      reference_bases='T',
                      alternate_bases=['G']),
                  allele_support=None)
          ]),
          expected_return=dict(T=1, G=0),
          testcase='no VCF'))
  def test_add_allele_frequencies_to_candidates(self, dv_calls, expected_return,
                                                testcase):
    if testcase == 'valid':
      pop_vcf_reader = vcf.VcfReader(testdata.VCF_WITH_ALLELE_FREQUENCIES)
      ref_reader = fasta.IndexedFastaReader(testdata.GRCH38_FASTA)
    elif testcase == 'no VCF':
      pop_vcf_reader = None
      ref_reader = None
    else:
      raise ValueError('Invalid testcase for parameterized test.')
    updated_dv_call = list(
        allele_frequency.add_allele_frequencies_to_candidates(
            dv_calls, pop_vcf_reader, ref_reader))
    actual_frequency = updated_dv_call[0].allele_frequency
    # Compare keys.
    self.assertSetEqual(
        set(actual_frequency.keys()), set(expected_return.keys()))
    # Compare values (almost equal).
    for key in actual_frequency.keys():
      self.assertAlmostEqual(actual_frequency[key], expected_return[key])


if __name__ == '__main__':
  absltest.main()
