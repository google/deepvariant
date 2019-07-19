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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import collections
import json

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import testdata
from deepvariant import vcf_stats
from third_party.nucleus.io import vcf
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

VariantStatsLite = collections.namedtuple('VariantStatsLite',
                                          ['genotype', 'vaf'])


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
      dict(
          alleles=['A', 'G'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_SNP,
          expected_transition=True,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'AG'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_INSERTION,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['AT', 'A'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_DELETION,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['AT', 'GC'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.BIALLELIC_MNP,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'G', 'T'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_SNP,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'AG', 'AT'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_INSERTION,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['AT', 'A', 'T'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_DELETION,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'GT', 'G'],
          gt=[0, 1],
          expected_variant_type=vcf_stats.MULTIALLELIC_COMPLEX,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=True),
      dict(
          alleles=['A', 'G'],
          gt=[0, 0],
          expected_variant_type=vcf_stats.REFCALL,
          expected_transition=False,
          expection_transversion=False,
          expected_is_variant=False),
  )
  def test_get_variant_stats(self, alleles, gt, expected_variant_type,
                             expected_transition, expection_transversion,
                             expected_is_variant):
    variant = test_utils.make_variant(
        chrom='chr1', start=10, alleles=alleles, gt=gt, gq=59)
    variant_stats = vcf_stats.get_variant_stats(variant)
    self.assertEqual(
        variant_stats,
        vcf_stats.VariantStats(
            reference_name='chr1',
            position=11,
            reference_bases=alleles[0],
            alternate_bases=alleles[1:],
            variant_type=expected_variant_type,
            is_transition=expected_transition,
            is_transversion=expection_transversion,
            is_variant=expected_is_variant,
            depth=[],
            genotype_quality=59,
            genotype=str(gt),
            vaf=None))

  def test_summary_stats(self):
    with vcf.VcfReader(testdata.GOLDEN_POSTPROCESS_OUTPUT) as reader:
      single_stats = vcf_stats.single_variant_stats(reader.iterate())
      summary_stats = vcf_stats.summary_stats(single_stats)
      sum_variant_count = sum([
          summary_stats.snv_count,
          summary_stats.insertion_count,
          summary_stats.deletion_count,
          summary_stats.complex_count,
          summary_stats.mnp_count,
      ])
      self.assertEqual(summary_stats.variant_count, 71)
      self.assertEqual(summary_stats.variant_count, sum_variant_count)
      self.assertEqual(summary_stats.snv_count, 59)
      self.assertEqual(summary_stats.insertion_count, 7)
      self.assertEqual(summary_stats.deletion_count, 5)
      self.assertEqual(summary_stats.complex_count, 0)
      self.assertEqual(summary_stats.mnp_count, 0)
      self.assertEqual(summary_stats.record_count, 76)
      self.assertAlmostEqual(summary_stats.depth_mean, 47.289473684210527)
      self.assertAlmostEqual(summary_stats.depth_stdev, 8.8953207531791154)
      self.assertAlmostEqual(summary_stats.gq_mean, 40.236842105263158)
      self.assertAlmostEqual(summary_stats.gq_stdev, 14.59710535567045)

  def test_variants_to_stats_json(self):
    truth_stats_json = """
      {"alternate_bases":[["G"]],"depth":[20],"genotype_quality":[59],
      "is_transition":[true],"is_transversion":[false],"position":[11],
      "reference_bases":["A"],"reference_name":["chr1"],"is_variant":[true],
      "variant_type":["Biallelic_SNP"],"genotype":["[0, 1]"],"vaf":[null]}
      """

    truth_summary_json = """
      {"depth_mean": 20,"depth_stdev": 0,"gq_mean": 59,"gq_stdev": 0,
      "record_count": 1,"snv_count": 1,"insertion_count": 0,"deletion_count": 0,
      "complex_count": 0, "mnp_count": 0, "variant_count": 1,
      "transition_count": 1, "transversion_count": 0}
      """

    # Without a VcfReader containing VAF, the histogram is empty
    truth_histograms = """
      {"[0, 1]": [
        {"bin_end":0.1,"count":0,"bin_start": 0.0},
        {"bin_end":0.2,"count":0,"bin_start": 0.1},
        {"bin_end":0.3,"count":0,"bin_start": 0.2},
        {"bin_end":0.4,"count":0,"bin_start": 0.3},
        {"bin_end":0.5,"count":0,"bin_start": 0.4},
        {"bin_end":0.6,"count":0,"bin_start": 0.5},
        {"bin_end":0.7,"count":0,"bin_start": 0.6},
        {"bin_end":0.8,"count":0,"bin_start": 0.7},
        {"bin_end":0.9,"count":0,"bin_start": 0.8},
        {"bin_end":1.0,"count":0,"bin_start": 0.9}
      ]}
    """

    stats_json, summary_json, histograms = vcf_stats.variants_to_stats_json(
        [self.variant])
    self.assertEqual(json.loads(stats_json), json.loads(truth_stats_json))
    self.assertEqual(json.loads(summary_json), json.loads(truth_summary_json))
    self.assertEqual(json.loads(histograms), json.loads(truth_histograms))

  def test_vaf_histograms_by_genotype(self):
    variants = [
        VariantStatsLite(genotype='[0, 0]', vaf=0),
        VariantStatsLite(genotype='[1, 1]', vaf=1),
        VariantStatsLite(genotype='[0, 1]', vaf=0.5),
        VariantStatsLite(genotype='[0, 1]', vaf=0.5),
        VariantStatsLite(genotype='[0, 0]', vaf=0.08),
        VariantStatsLite(genotype='[0, 0]', vaf=0.19),
        VariantStatsLite(genotype='[0, 1]', vaf=0.45),
        VariantStatsLite(genotype='[0, 1]', vaf=0.65)
    ]
    truth_histograms = {
        '[0, 0]': [{
            'bin_end': 0.1,
            'count': 2,
            'bin_start': 0.0
        }, {
            'bin_end': 0.2,
            'count': 1,
            'bin_start': 0.1
        }, {
            'bin_end': 0.3,
            'count': 0,
            'bin_start': 0.2
        }, {
            'bin_end': 0.4,
            'count': 0,
            'bin_start': 0.3
        }, {
            'bin_end': 0.5,
            'count': 0,
            'bin_start': 0.4
        }, {
            'bin_end': 0.6,
            'count': 0,
            'bin_start': 0.5
        }, {
            'bin_end': 0.7,
            'count': 0,
            'bin_start': 0.6
        }, {
            'bin_end': 0.8,
            'count': 0,
            'bin_start': 0.7
        }, {
            'bin_end': 0.9,
            'count': 0,
            'bin_start': 0.8
        }, {
            'bin_end': 1.0,
            'count': 0,
            'bin_start': 0.9
        }],
        '[0, 1]': [{
            'bin_end': 0.1,
            'count': 0,
            'bin_start': 0.0
        }, {
            'bin_end': 0.2,
            'count': 0,
            'bin_start': 0.1
        }, {
            'bin_end': 0.3,
            'count': 0,
            'bin_start': 0.2
        }, {
            'bin_end': 0.4,
            'count': 0,
            'bin_start': 0.3
        }, {
            'bin_end': 0.5,
            'count': 1,
            'bin_start': 0.4
        }, {
            'bin_end': 0.6,
            'count': 2,
            'bin_start': 0.5
        }, {
            'bin_end': 0.7,
            'count': 1,
            'bin_start': 0.6
        }, {
            'bin_end': 0.8,
            'count': 0,
            'bin_start': 0.7
        }, {
            'bin_end': 0.9,
            'count': 0,
            'bin_start': 0.8
        }, {
            'bin_end': 1.0,
            'count': 0,
            'bin_start': 0.9
        }],
        '[1, 1]': [{
            'bin_end': 0.1,
            'count': 0,
            'bin_start': 0.0
        }, {
            'bin_end': 0.2,
            'count': 0,
            'bin_start': 0.1
        }, {
            'bin_end': 0.3,
            'count': 0,
            'bin_start': 0.2
        }, {
            'bin_end': 0.4,
            'count': 0,
            'bin_start': 0.3
        }, {
            'bin_end': 0.5,
            'count': 0,
            'bin_start': 0.4
        }, {
            'bin_end': 0.6,
            'count': 0,
            'bin_start': 0.5
        }, {
            'bin_end': 0.7,
            'count': 0,
            'bin_start': 0.6
        }, {
            'bin_end': 0.8,
            'count': 0,
            'bin_start': 0.7
        }, {
            'bin_end': 0.9,
            'count': 0,
            'bin_start': 0.8
        }, {
            'bin_end': 1.0,
            'count': 1,
            'bin_start': 0.9
        }]
    }
    self.assertEqual(
        vcf_stats.vaf_histograms_by_genotype(variants), truth_histograms)


if __name__ == '__main__':
  absltest.main()
