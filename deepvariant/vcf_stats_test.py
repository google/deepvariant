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


import json

from absl.testing import absltest

from deepvariant import testdata
from deepvariant import vcf_stats
from third_party.nucleus.testing import test_utils


def setUpModule():
  testdata.init()


class VcfStatsTest(absltest.TestCase):

  def setUp(self):
    super(VcfStatsTest, self).setUp()
    self.variant = test_utils.make_variant(
        chrom='chr1', start=10, alleles=['A', 'G'], gt=[0, 1])

  def test_get_variant_stats(self):
    variant_stats = vcf_stats.get_variant_stats(self.variant)
    self.assertEqual(
        variant_stats,
        vcf_stats.VariantStats(
            reference_name='chr1',
            position=11,
            reference_bases='A',
            alternate_bases=['G'],
            variant_type='SNP',
            is_transition=True,
            is_transversion=False,
            genotype_quality=[]))

  def test_variants_to_stats_json(self):
    stats_json = vcf_stats.variants_to_stats_json([self.variant])
    truth_stats_json = """
    {"alternate_bases":[["G"]],"genotype_quality":[[]],"is_transition":[true],
    "is_transversion":[false],"position":[11],"reference_bases":["A"],
    "reference_name":["chr1"],"variant_type":["SNP"]}
    """
    self.assertEqual(json.loads(stats_json), json.loads(truth_stats_json))


if __name__ == '__main__':
  absltest.main()
