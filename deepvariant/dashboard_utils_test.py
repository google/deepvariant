# Copyright 2021 Google LLC.
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
"""Tests for deepvariant .dashboard_utils."""

import io

from absl.testing import absltest
from absl.testing import parameterized
import altair as alt
import pandas as pd

from deepvariant import dashboard_utils


class DashboardUtilsTest(parameterized.TestCase):

  def test_create_html_report(self):
    df = pd.DataFrame({'A': [1, 2], 'B': [10, 20]})
    chart = alt.Chart(df).mark_point().encode(x='A', y='B')
    charts = [{'id': 'my_chart_name', 'chart': chart}]

    html_output = io.StringIO()
    dashboard_utils.create_html_report(
        charts,
        html_output=html_output,
        title='my fancy title',
        subtitle='my fancy subtitle')

    html = html_output.getvalue()
    self.assertIn(
        'my_chart_name', html, msg='Chart ID is missing from the HTML.')
    self.assertIn('point', html, msg='Vega specs may be missing from the HTML.')
    self.assertIn(
        'my fancy title', html, msg='The title is missing from the HTML.')
    self.assertIn(
        'my fancy subtitle', html, msg='The subtitle is missing from the HTML.')
    self.assertIn('_blank', html, msg='Embed options missing.')


if __name__ == '__main__':
  absltest.main()
