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
"""Tests for deepvariant .vcf_stats_vis."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']



import os
import tempfile
from absl.testing import absltest
import altair as alt
import pandas as pd
import tensorflow as tf

from deepvariant import vcf_stats_vis

ALTAIR_CHART = "<class 'altair.vegalite.v3.api.Chart'>"

STATS_DATA = {
    "insertion_count": 23,
    "depth_mean": 27.972999999999999,
    "record_count": 1000,
    "snv_count": 309,
    "deletion_count": 30,
    "depth_stdev": 25.42967304154735,
    "variant_count": 362,
    "transition_count": 205,
    "mnp_count": 0,
    "complex_count": 0,
    "transversion_count": 104,
    "gq_mean": 20.869,
    "gq_stdev": 14.447139474650337
}

VIS_DATA = {
    "base_changes": [["G", "A", 56], ["T", "A", 17], ["C", "T", 47],
                     ["G", "C", 19], ["T", "C", 48], ["C", "A", 14],
                     ["A", "T", 9], ["A", "C", 15], ["T", "G", 9],
                     ["G", "T", 15], ["A", "G", 60], ["C", "G", 11]],
    "gq_histogram": [{
        "bin_start": 0,
        "bin_end": 0.5,
        "count": 3
    }, {
        "bin_start": 0.5,
        "bin_end": 1,
        "count": 24
    }],
    "indel_sizes": [[1, 6], [2, 4], [4, 2], [5, 2], [7, 2], [8, 1], [12, 1],
                    [-2, 6], [-5, 1], [-4, 7], [-3, 4], [-1, 11]],
    "qual_histogram": [{
        "bin_start": 0,
        "bin_end": 50,
        "count": 10
    }, {
        "bin_start": 50,
        "bin_end": 99,
        "count": 10
    }],
    "vaf_histograms_by_genotype": {
        "[-1, -1]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }],
        "[0, 0]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }],
        "[0, 1]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }],
        "[0, 2]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }],
        "[1, 1]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }],
        "[1, 2]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }],
        "[1, 3]": [{
            "bin_end": 0.5,
            "bin_start": 0,
            "count": 10
        }, {
            "bin_end": 1,
            "bin_start": 0.5,
            "count": 10
        }]
    },
    "variant_type_counts": {
        "Biallelic_SNP": 10,
        "RefCall": 3,
        "Multiallelic_Insertion": 1
    }
}


class VcfStatsVisTest(absltest.TestCase):

  def test_dict_to_dataframe(self):
    self.assertEqual("K", "K")
    self.assertEqual(
        vcf_stats_vis._dict_to_dataframe({
            "A": "a"
        }).to_dict("records"), [{
            "label": "A",
            "value": "a"
        }])

  def test_prettify_genotype(self):
    self.assertEqual(
        vcf_stats_vis._prettify_genotype("[0, 0]"), (vcf_stats_vis.REF, "main"))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype("[-1, -1]"),
        (vcf_stats_vis.UNCALLED, "others"))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype("[3, 3]"), (vcf_stats_vis.HOM, "main"))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype("[0, 3]"), (vcf_stats_vis.HET, "main"))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype("[6, 3]"),
        (vcf_stats_vis.HET_BOTH, "others"))

  def test_build_type_chart(self):
    stats = pd.DataFrame([["Deletion", 10], ["Insertion", 20]],
                         columns=["label", "value"])
    chart = vcf_stats_vis._build_type_chart(stats)
    self.assertEqual(
        str(type(chart)), "<class 'altair.vegalite.v3.api.LayerChart'>")

  def test_build_tt_chart(self):
    stats = pd.DataFrame([["Transition", 10], ["Transversion", 20]],
                         columns=["label", "value"])
    chart = vcf_stats_vis._build_tt_chart(stats, 2.0)
    self.assertEqual(str(type(chart)), ALTAIR_CHART)

  def test_build_qual_histogram(self):
    chart = vcf_stats_vis._build_qual_histogram(VIS_DATA)
    self.assertEqual(str(type(chart)), ALTAIR_CHART)

  def test_build_gq_histogram(self):
    chart = vcf_stats_vis._build_gq_histogram(VIS_DATA)
    self.assertEqual(str(type(chart)), ALTAIR_CHART)

  def test_build_vaf_histograms(self):
    chart = vcf_stats_vis._build_vaf_histograms(VIS_DATA)
    self.assertEqual(
        str(type(chart[0])), "<class 'altair.vegalite.v3.api.FacetChart'>")
    self.assertEqual(str(type(chart[1])), ALTAIR_CHART)

  def test_build_base_change_chart(self):
    chart = vcf_stats_vis._build_base_change_chart(VIS_DATA)
    self.assertEqual(
        str(type(chart)), "<class 'altair.vegalite.v3.api.FacetChart'>")

  def test_build_indel_size_chart(self):
    chart = vcf_stats_vis._build_indel_size_chart(VIS_DATA)
    self.assertEqual(
        str(type(chart)), "<class 'altair.vegalite.v3.api.VConcatChart'>")

  def test_build_all_charts(self):
    chart = vcf_stats_vis._build_all_charts(STATS_DATA, VIS_DATA)
    self.assertEqual(
        str(type(chart)), "<class 'altair.vegalite.v3.api.VConcatChart'>")

  def test_altair_chart_to_html(self):
    df = pd.DataFrame({"x": ["A", "B"], "y": [28, 55]})
    c = alt.Chart(df).mark_bar().encode(x="x", y="y")
    html_string = vcf_stats_vis._altair_chart_to_html(
        altair_chart=c, download_filename="TEST_DOWNLOAD_FILENAME")
    import_base = 'src="https://storage.googleapis.com/deepvariant/lib/vega/'
    self.assertNotEqual(
        html_string.find(import_base + 'vega@%s"' %
                         (vcf_stats_vis.VEGA_VERSION)), -1)
    self.assertNotEqual(
        html_string.find(import_base + 'vega-lite@%s"' %
                         (vcf_stats_vis.VEGA_LITE_VERSION)), -1)
    self.assertNotEqual(
        html_string.find(import_base + 'vega-embed@%s"' %
                         (vcf_stats_vis.VEGA_EMBED_VERSION)), -1)
    self.assertEqual(html_string.find("jsdelivr.net"), -1)
    self.assertNotEqual(html_string.find("TEST_DOWNLOAD_FILENAME"), -1)

  def test_create_visual_report(self):
    base_dir = tempfile.mkdtemp()
    outfile_base = os.path.join(base_dir, "stats_test")
    sample_name = "test_sample_name"
    vcf_stats_vis.create_visual_report(
        outfile_base, STATS_DATA, VIS_DATA, sample_name=sample_name)
    self.assertTrue(tf.io.gfile.exists(outfile_base + ".visual_report.html"))

if __name__ == "__main__":
  absltest.main()
