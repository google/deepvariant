# Copyright 2019 Google LLC.
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
r"""Tests for deepvariant .vcf_stats."""

import collections
import json
import os
import tempfile

from absl.testing import absltest
from absl.testing import parameterized
import six
import tensorflow as tf

from deepvariant import testdata
from deepvariant import vcf_stats
from third_party.nucleus.io import vcf
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils


def setUpModule():
  testdata.init()


class VcfStatsTest(parameterized.TestCase):

  def setUp(self):
    super(VcfStatsTest, self).setUp()
    self.variant = test_utils.make_variant(
        chrom='chr1', start=10, alleles=['A', 'G'], gt=[0, 1], gq=59)
    variantcall_utils.set_format(
        variant_utils.only_call(self.variant), 'DP', 20)

  @parameterized.parameters(
      (test_utils.make_variant(alleles=['A', 'C']), vcf_stats.BIALLELIC_SNP),
      (test_utils.make_variant(alleles=['A', 'C', '<*>']),
       vcf_stats.BIALLELIC_SNP),
      (test_utils.make_variant(alleles=['A', 'AG']),
       vcf_stats.BIALLELIC_INSERTION),
      (test_utils.make_variant(alleles=['A', 'AG', '<*>']),
       vcf_stats.BIALLELIC_INSERTION),
      (test_utils.make_variant(alleles=['AG', 'A']),
       vcf_stats.BIALLELIC_DELETION),
      (test_utils.make_variant(alleles=['AG', 'A', '<*>']),
       vcf_stats.BIALLELIC_DELETION),
      (test_utils.make_variant(alleles=['A', 'C', 'G']),
       vcf_stats.MULTIALLELIC_SNP),
      (test_utils.make_variant(alleles=['A', 'C', 'G', '<*>']),
       vcf_stats.MULTIALLELIC_SNP),
      (test_utils.make_variant(alleles=['A', 'AC', 'AG']),
       vcf_stats.MULTIALLELIC_INSERTION),
      (test_utils.make_variant(alleles=['A', 'AC', 'AG', '<*>']),
       vcf_stats.MULTIALLELIC_INSERTION),
      (test_utils.make_variant(alleles=['AGC', 'AC', 'A', 'AG']),
       vcf_stats.MULTIALLELIC_DELETION),
      (test_utils.make_variant(alleles=['AGC', 'AC', 'A', 'AG', '<*>']),
       vcf_stats.MULTIALLELIC_DELETION),
      (test_utils.make_variant(alleles=['AG', 'AC', 'A']),
       vcf_stats.MULTIALLELIC_COMPLEX),
      (test_utils.make_variant(alleles=['AG', 'AC', 'A', '<*>']),
       vcf_stats.MULTIALLELIC_COMPLEX),
      (test_utils.make_variant(alleles=['A', 'G', 'AT']),
       vcf_stats.MULTIALLELIC_COMPLEX),
      (test_utils.make_variant(alleles=['A', 'G', 'AT', '<*>']),
       vcf_stats.MULTIALLELIC_COMPLEX),
      (test_utils.make_variant(alleles=['AG', 'TC']), vcf_stats.BIALLELIC_MNP),
      (test_utils.make_variant(alleles=['AG', 'TC', '<*>']),
       vcf_stats.BIALLELIC_MNP),
      (test_utils.make_variant(alleles=['A']), vcf_stats.REFCALL),
      (test_utils.make_variant(alleles=['A', '<*>']), vcf_stats.REFCALL),
      (test_utils.make_variant(alleles=['A', 'G'],
                               filters='FAIL'), vcf_stats.REFCALL),
      (test_utils.make_variant(alleles=['A', 'G', '<*>'],
                               filters='FAIL'), vcf_stats.REFCALL),
      (test_utils.make_variant(alleles=['A', 'G'],
                               filters=['FAIL']), vcf_stats.REFCALL),
      (test_utils.make_variant(alleles=['A', '<*>'],
                               filters='.'), vcf_stats.REFCALL),
  )
  def test_get_variant_type(self, variant, expected_type):
    self.assertEqual(vcf_stats._get_variant_type(variant), expected_type)

  @parameterized.parameters(
      dict(
          alleles=['A', 'G'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_SNP,
          expected_transition=True,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'AG'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_INSERTION,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['AT', 'A'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_DELETION,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['AT', 'GC'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_MNP,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'G', 'T'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_SNP,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'AG', 'AT'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_INSERTION,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['AT', 'A', 'T'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_DELETION,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'GT', 'G'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_COMPLEX,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'G'],
          gt=[0, 0],
          expected_variant_type=vcf_stats.REFCALL,
          expected_transition=False,
          expected_transversion=False,
          expected_is_variant=False),
  )
  def test_get_variant_stats(self, alleles, gt, expected_variant_type,
                             expected_transition, expected_transversion,
                             expected_is_variant):
    variant = test_utils.make_variant(
        chrom='chr1', start=10, alleles=alleles, gt=gt, gq=59)
    variant_stats = vcf_stats._get_variant_stats(variant)
    self.assertEqual(
        variant_stats,
        vcf_stats.VariantStats(
            reference_name='chr1',
            position=11,
            reference_bases=alleles[0],
            alternate_bases=alleles[1:],
            variant_type=expected_variant_type,
            is_transition=expected_transition,
            is_transversion=expected_transversion,
            is_variant=expected_is_variant,
            depth=[],
            genotype_quality=59,
            genotype=str(gt),
            vaf=None,
            qual=0.0))

  def test_compute_variant_stats_for_charts(self):
    expected_keys = [
        'vaf_histograms_by_genotype', 'indel_sizes', 'base_changes',
        'qual_histogram', 'gq_histogram', 'variant_type_counts',
        'depth_histogram', 'titv_counts'
    ]
    vis_data = vcf_stats._compute_variant_stats_for_charts([self.variant])
    six.assertCountEqual(
        self,
        vis_data.keys(),
        expected_keys,
        msg='vis_data does not have the right keys')

  def test_vaf_histograms_by_genotype(self):
    variant_stats_lite = collections.namedtuple('variant_stats_lite',
                                                ['genotype', 'vaf'])
    variant_stats = [
        variant_stats_lite(genotype='[0, 0]', vaf=0),
        variant_stats_lite(genotype='[1, 1]', vaf=1),
        variant_stats_lite(genotype='[0, 1]', vaf=0.5),
        variant_stats_lite(genotype='[0, 1]', vaf=0.5),
        variant_stats_lite(genotype='[0, 0]', vaf=0.08),
        variant_stats_lite(genotype='[0, 0]', vaf=0.19),
        variant_stats_lite(genotype='[0, 1]', vaf=0.45),
        variant_stats_lite(genotype='[0, 1]', vaf=0.65)
    ]
    # s = bin_start, e = bin_end, c = count
    truth_histograms = """
    {
      "[0, 1]": [{"c": 0, "e": 0.1, "s": 0.0}, {"c": 0, "e": 0.2, "s": 0.1}, {"c": 0, "e": 0.3, "s": 0.2}, {"c": 0, "e": 0.4, "s": 0.3}, {"c": 1, "e": 0.5, "s": 0.4}, {"c": 2, "e": 0.6, "s": 0.5}, {"c": 1, "e": 0.7, "s": 0.6}, {"c": 0, "e": 0.8, "s": 0.7}, {"c": 0, "e": 0.9, "s": 0.8}, {"c": 0, "e": 1.0, "s": 0.9}],
      "[1, 1]": [{"c": 0, "e": 0.1, "s": 0.0}, {"c": 0, "e": 0.2, "s": 0.1}, {"c": 0, "e": 0.3, "s": 0.2}, {"c": 0, "e": 0.4, "s": 0.3}, {"c": 0, "e": 0.5, "s": 0.4}, {"c": 0, "e": 0.6, "s": 0.5}, {"c": 0, "e": 0.7, "s": 0.6}, {"c": 0, "e": 0.8, "s": 0.7}, {"c": 0, "e": 0.9, "s": 0.8}, {"c": 1, "e": 1.0, "s": 0.9}],
      "[0, 0]": [{"c": 2, "e": 0.1, "s": 0.0}, {"c": 1, "e": 0.2, "s": 0.1}, {"c": 0, "e": 0.3, "s": 0.2}, {"c": 0, "e": 0.4, "s": 0.3}, {"c": 0, "e": 0.5, "s": 0.4}, {"c": 0, "e": 0.6, "s": 0.5}, {"c": 0, "e": 0.7, "s": 0.6}, {"c": 0, "e": 0.8, "s": 0.7}, {"c": 0, "e": 0.9, "s": 0.8}, {"c": 0, "e": 1.0, "s": 0.9}],
      "[-1, -1]": [{"c": 0, "e": 0.5, "s": 0.0}, {"c": 0, "e": 1.0, "s": 0.5}],
      "[1, 2]": [{"c": 0, "e": 0.5, "s": 0.0}, {"c": 0, "e": 1.0, "s": 0.5}]
      }
    """
    self.assertEqual(
        vcf_stats._vaf_histograms_by_genotype(variant_stats),
        json.loads(truth_histograms))

  def test_format_histogram_for_vega(self):
    # s = bin_start, e = bin_end, c = count
    self.assertEqual(
        vcf_stats._format_histogram_for_vega(counts=[2, 2], bins=[1, 2.5, 4]),
        [{
            's': 1,
            'e': 2.5,
            'c': 2
        }, {
            's': 2.5,
            'e': 4,
            'c': 2
        }])

  def test_count_titv(self):
    variant_stats_lite = collections.namedtuple(
        'variant_stats_lite', ['is_transition', 'is_transversion'])
    variant_stats = [
        variant_stats_lite(is_transition=True, is_transversion=False),
        variant_stats_lite(is_transition=True, is_transversion=False),
        variant_stats_lite(is_transition=True, is_transversion=False),
        variant_stats_lite(is_transition=False, is_transversion=True),
        variant_stats_lite(is_transition=False, is_transversion=True),
        variant_stats_lite(is_transition=False, is_transversion=False)
    ]
    truth_counts = {'Transition': 3, 'Transversion': 2}
    self.assertEqual(vcf_stats._count_titv(variant_stats), truth_counts)

  def test_count_variant_types(self):
    variant_stats_lite = collections.namedtuple('variant_stats_lite',
                                                ['variant_type'])
    variant_stats = [
        variant_stats_lite(variant_type='A'),
        variant_stats_lite(variant_type='B'),
        variant_stats_lite(variant_type='C'),
        variant_stats_lite(variant_type='A'),
        variant_stats_lite(variant_type='B')
    ]
    truth_counts = {'A': 2, 'B': 2, 'C': 1}
    self.assertEqual(
        vcf_stats._count_variant_types(variant_stats), truth_counts)

  def test_count_base_changes_and_indel_sizes(self):
    variant_stats_lite = collections.namedtuple(
        'variant_stats_lite',
        ['reference_bases', 'alternate_bases', 'is_variant', 'variant_type'])
    variant_stats = [
        variant_stats_lite(
            reference_bases='A',
            alternate_bases=['G'],
            is_variant=True,
            variant_type=vcf_stats.BIALLELIC_SNP),
        variant_stats_lite(
            reference_bases='A',
            alternate_bases=['AGGG'],
            is_variant=True,
            variant_type=vcf_stats.BIALLELIC_INSERTION),
        variant_stats_lite(
            reference_bases='A',
            alternate_bases=['G'],
            is_variant=False,
            variant_type=vcf_stats.REFCALL),
        variant_stats_lite(
            reference_bases='A',
            alternate_bases=['G', 'T'],
            is_variant=True,
            variant_type=vcf_stats.MULTIALLELIC_COMPLEX)
    ]
    truth_base_changes = [['A', 'G', 1]]
    truth_indel_sizes = [[3, 1]]
    base_changes, indel_sizes = vcf_stats._count_base_changes_and_indel_sizes(
        variant_stats)
    self.assertEqual(base_changes, truth_base_changes)
    self.assertEqual(indel_sizes, truth_indel_sizes)

  def test_compute_qual_histogram(self):
    variant_stats_lite = collections.namedtuple('variant_stats_lite', ['qual'])
    variant_stats = [variant_stats_lite(qual=100), variant_stats_lite(qual=49)]
    hist = vcf_stats._compute_qual_histogram(variant_stats)
    # s = bin_start, e = bin_end, c = count
    self.assertEqual(hist, [{
        'c': 1,
        's': 49.0,
        'e': 50.0
    }, {
        'c': 1,
        's': 100.0,
        'e': 101.0
    }])

  def test_get_integer_counts(self):
    self.assertEqual(
        vcf_stats._get_integer_counts([1, 2, 2, 4]), [[1, 1], [2, 2], [4, 1]])

  def test_compute_gq_histogram(self):
    variant_stats_lite = collections.namedtuple('variant_stats_lite',
                                                ['genotype_quality'])
    variant_stats = [
        variant_stats_lite(genotype_quality=100),
        variant_stats_lite(genotype_quality=100),
        variant_stats_lite(genotype_quality=49)
    ]
    hist = vcf_stats._compute_gq_histogram(variant_stats)
    self.assertEqual(hist, [[49, 1], [100, 2]])

  def test_compute_depth_histogram(self):
    variant_stats_lite = collections.namedtuple('variant_stats_lite', ['depth'])
    variant_stats = [
        variant_stats_lite(depth=100),
        variant_stats_lite(depth=30),
        variant_stats_lite(depth=30)
    ]
    hist = vcf_stats._compute_depth_histogram(variant_stats)
    self.assertEqual(hist, [[30, 2], [100, 1]])

  def test_create_vcf_report(self):
    base_dir = tempfile.mkdtemp()
    outfile_base = os.path.join(base_dir, 'stats_test')
    sample_name = 'test_sample_name'
    with vcf.VcfReader(testdata.GOLDEN_POSTPROCESS_OUTPUT) as reader:
      vcf_stats.create_vcf_report(
          variants=reader.iterate(),
          output_basename=outfile_base,
          sample_name=sample_name,
          vcf_reader=reader)
    self.assertTrue(tf.io.gfile.exists(outfile_base + '.visual_report.html'))


if __name__ == '__main__':
  absltest.main()
