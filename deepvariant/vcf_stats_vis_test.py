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

import os
import tempfile
from absl.testing import absltest
import altair as alt
import pandas as pd
import six
import tensorflow as tf

from deepvariant import vcf_stats_vis

# Note: histograms all have keys s, e, and c, shortened versions of
# bin_start, bin_end, and count to save space in output HTML
VIS_DATA = {
    'base_changes': [['G', 'A', 56], ['T', 'A', 17], ['C', 'T', 47],
                     ['G', 'C', 19], ['T', 'C', 48], ['C', 'A', 14],
                     ['A', 'T', 9], ['A', 'C', 15], ['T', 'G', 9],
                     ['G', 'T', 15], ['A', 'G', 60], ['C', 'G', 11]],
    'gq_histogram': [[1, 3], [2, 24]],
    'indel_sizes': [[1, 6], [2, 4], [4, 2], [5, 2], [7, 2], [8, 1], [12, 1],
                    [-2, 6], [-5, 1], [-4, 7], [-3, 4], [-1, 11]],
    'qual_histogram': [{
        's': 0,
        'e': 50,
        'c': 10
    }, {
        's': 50,
        'e': 99,
        'c': 10
    }],
    'depth_histogram': [[0, 10], [1, 20]],
    'vaf_histograms_by_genotype': {
        '[-1, -1]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }],
        '[0, 0]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }],
        '[0, 1]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }],
        '[0, 2]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }],
        '[1, 1]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }],
        '[1, 2]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }],
        '[1, 3]': [{
            'e': 0.5,
            's': 0,
            'c': 10
        }, {
            'e': 1,
            's': 0.5,
            'c': 10
        }]
    },
    'variant_type_counts': {
        'Biallelic_SNP': 10,
        'RefCall': 3,
        'Multiallelic_Insertion': 1
    },
    'titv_counts': {
        'Transition': 20,
        'Transversion': 10
    }
}


def is_an_altair_chart(chart):
  # Chart type strings look like: "<class 'altair.vegalite.v3.api.FacetChart'>"
  # Chart, FacetChart, LayerChart, and VConcatChart.
  string_type = str(type(chart))
  return 'altair' in string_type and 'Chart' in string_type


class VcfStatsVisTest(absltest.TestCase):

  def test_dict_to_dataframe(self):
    self.assertEqual('K', 'K')
    self.assertEqual(
        vcf_stats_vis._dict_to_dataframe({
            'A': 'a'
        }).to_dict('records'), [{
            'label': 'A',
            'value': 'a'
        }])

  def test_prettify_genotype(self):
    self.assertEqual(
        vcf_stats_vis._prettify_genotype('[0, 0]'), (vcf_stats_vis.REF, 'main'))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype('[-1, -1]'),
        (vcf_stats_vis.UNCALLED, 'others'))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype('[3, 3]'), (vcf_stats_vis.HOM, 'main'))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype('[0, 3]'), (vcf_stats_vis.HET, 'main'))
    self.assertEqual(
        vcf_stats_vis._prettify_genotype('[6, 3]'),
        (vcf_stats_vis.HET_BOTH, 'others'))

  def test_integer_counts_to_histogram(self):
    test_input = [[1, 1], [2, 2], [4, 1]]
    expected_output = pd.DataFrame(
        data={
            'c': [1, 2, 1],
            's': [0.5, 1.5, 3.5],
            'e': [1.5, 2.5, 4.5]
        },
        columns=['c', 's', 'e'])
    observed_output = vcf_stats_vis._integer_counts_to_histogram(test_input)
    six.assertCountEqual(
        self,
        list(observed_output.columns),
        list(expected_output.columns),
        msg='Wrong column names')
    self.assertEqual(
        list(observed_output['c']),
        list(expected_output['c']),
        msg='column c differs')
    self.assertEqual(
        list(observed_output['s']),
        list(expected_output['s']),
        msg='column s differs')
    self.assertEqual(
        list(observed_output['e']),
        list(expected_output['e']),
        msg='column e differs')
    self.assertTrue((observed_output == expected_output).all().all())

  def test_chart_type_negative_control(self):
    self.assertFalse(is_an_altair_chart('some string'))
    self.assertFalse(is_an_altair_chart(None))

  def test_build_type_chart(self):
    chart = vcf_stats_vis._build_type_chart(VIS_DATA['variant_type_counts'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_tt_chart(self):
    chart = vcf_stats_vis._build_tt_chart(VIS_DATA['titv_counts'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_qual_histogram(self):
    chart = vcf_stats_vis._build_qual_histogram(VIS_DATA['qual_histogram'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_depth_histogram(self):
    chart = vcf_stats_vis._build_depth_histogram(VIS_DATA['depth_histogram'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_gq_histogram(self):
    chart = vcf_stats_vis._build_gq_histogram(VIS_DATA['gq_histogram'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_vaf_histograms(self):
    chart = vcf_stats_vis._build_vaf_histograms(
        VIS_DATA['vaf_histograms_by_genotype'])
    self.assertTrue(is_an_altair_chart(chart[0]))
    self.assertTrue(is_an_altair_chart(chart[1]))

  def test_build_base_change_chart(self):
    chart = vcf_stats_vis._build_base_change_chart(VIS_DATA['base_changes'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_indel_size_chart(self):
    chart = vcf_stats_vis._build_indel_size_chart(VIS_DATA['indel_sizes'])
    self.assertTrue(is_an_altair_chart(chart))

  def test_build_all_charts(self):
    chart = vcf_stats_vis._build_all_charts(VIS_DATA)
    self.assertTrue(is_an_altair_chart(chart))

  def test_altair_chart_to_html(self):
    df = pd.DataFrame({'x': ['A', 'B'], 'y': [28, 55]})
    c = alt.Chart(df).mark_bar().encode(x='x', y='y')
    html_string = vcf_stats_vis._altair_chart_to_html(
        altair_chart=c, download_filename='TEST_DOWNLOAD_FILENAME')
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
    self.assertEqual(html_string.find('jsdelivr.net'), -1)
    self.assertNotEqual(html_string.find('TEST_DOWNLOAD_FILENAME'), -1)

  def test_create_visual_report(self):
    base_dir = tempfile.mkdtemp()
    outfile_base = os.path.join(base_dir, 'stats_test')
    sample_name = 'test_sample_name'
    vcf_stats_vis.create_visual_report(
        outfile_base, VIS_DATA, sample_name=sample_name)
    self.assertTrue(tf.io.gfile.exists(outfile_base + '.visual_report.html'))


if __name__ == '__main__':
  absltest.main()
