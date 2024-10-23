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
"""Tests for DeepVariant runtime_by_region_vis visual report script."""

import io

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import runtime_by_region_vis
from deepvariant import testdata


def setUpModule():
  testdata.init()


def is_an_altair_chart(chart):
  # Chart type strings look like: "<class 'altair.vegalite.v3.api.FacetChart'>"
  # Chart, FacetChart, LayerChart, and VConcatChart.
  string_type = str(type(chart))
  return 'altair' in string_type and 'Chart' in string_type


class RuntimeByRegionVisTest(parameterized.TestCase):

  @classmethod
  def setUpClass(cls):
    super().setUpClass()
    cls.JSON_DF, cls.JSON_BY_TASK_DF = (
        runtime_by_region_vis.read_data_and_make_dataframes(
            testdata.RUNTIME_BY_REGION
        )
    )

  @parameterized.parameters(
      dict(sharded=False, expected_regions=101),
      dict(sharded=True, expected_regions=101),
  )
  def test_e2e(self, sharded, expected_regions):
    if sharded:
      input_path = testdata.RUNTIME_BY_REGION_SHARDED
    else:
      input_path = testdata.RUNTIME_BY_REGION

    html_output = io.StringIO()
    runtime_by_region_vis.make_report(
        input_path=input_path, title='my fancy title', html_output=html_output
    )
    html = html_output.getvalue()
    self.assertIn(
        'my fancy title', html, msg='The title is missing from the HTML.'
    )
    self.assertIn(
        '{} regions'.format(expected_regions),
        html,
        msg='The subtitle contains the number of regions.',
    )
    self.assertIn(
        'regions account for',
        html,
        msg='The Pareto curve may be missing or it changed title',
    )
    self.assertIn('bar', html, msg='Vega specs may be missing from the HTML')
    self.assertNotIn('sdlfkjdkjf', html, msg='Negative control failed')

  @parameterized.parameters(
      dict(raw_seconds=5, expected='5s'),
      dict(raw_seconds=3600, expected='1h'),
      dict(raw_seconds=62, expected='1m2s'),
      dict(raw_seconds=7200, expected='2h'),
      dict(raw_seconds=3661, expected='1h1m1s'),
      dict(raw_seconds=0.0001, expected='0.0s'),
      dict(raw_seconds=0.001, expected='0.001s'),
      dict(raw_seconds=0.1, expected='0.1s'),
  )
  def test_format_runtime_string(self, raw_seconds, expected):
    self.assertEqual(
        expected, runtime_by_region_vis.format_runtime_string(raw_seconds)
    )

  def test_chart_type_negative_control(self):
    self.assertFalse(is_an_altair_chart('some string'))
    self.assertFalse(is_an_altair_chart(None))

  def test_totals_by_stage(self):
    chart = runtime_by_region_vis.totals_by_stage(self.JSON_BY_TASK_DF)
    self.assertTrue(is_an_altair_chart(chart))

  def test_pareto_and_runtimes_by_task(self):
    chart = runtime_by_region_vis.pareto_and_runtimes_by_task(self.JSON_DF)
    self.assertTrue(is_an_altair_chart(chart))

  def test_stage_histogram(self):
    parameters = [
        dict(
            dataframe_json=self.JSON_BY_TASK_DF,
            msg='Histogram of tasks',
        ),
        dict(
            dataframe_json=self.JSON_DF,
            msg='Histogram of regions',
        ),
    ]
    for param in parameters:
      df = param['dataframe_json']
      chart = runtime_by_region_vis.stage_histogram(df, title='chart title')
      self.assertTrue(is_an_altair_chart(chart), msg=param['msg'])
      chart_json = chart.to_json()
      self.assertIn('chart title', chart_json)

  def test_selected_longest_and_median_regions(self):
    chart = runtime_by_region_vis.selected_longest_and_median_regions(
        self.JSON_DF
    )
    self.assertTrue(is_an_altair_chart(chart))

  def test_top_regions_producing_zero_examples(self):
    chart = runtime_by_region_vis.top_regions_producing_zero_examples(
        self.JSON_DF
    )
    self.assertTrue(is_an_altair_chart(chart))

  def test_correlation_scatter_charts(self):
    chart = runtime_by_region_vis.correlation_scatter_charts(
        self.JSON_DF, title='chart title'
    )
    self.assertTrue(is_an_altair_chart(chart))
    chart_json = chart.to_json()
    self.assertIn('chart title', chart_json)

  def test_individual_region_bars(self):
    chart = runtime_by_region_vis.individual_region_bars(
        self.JSON_DF, title='chart title'
    )
    self.assertTrue(is_an_altair_chart(chart))
    chart_json = chart.to_json()
    self.assertIn('chart title', chart_json)


if __name__ == '__main__':
  absltest.main()
