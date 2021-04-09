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
"""Library code for creating visual reports or dashboards from Altair charts.

This is used by different dashboards that have multiple Altair charts to put
together into a single report or dashboard.
"""

from typing import Any, Dict, List, Text

import altair as alt

VEGA_URL = 'https://storage.googleapis.com/deepvariant/lib/vega'

CSS_STYLES = """
<style>
    body {
      font-family: sans-serif;
    }
    .chart-container {
      padding: 30px;
    }
</style>
"""


def create_html_report(charts: List[Dict[Text, alt.Chart]],
                       html_output: Any,
                       title: str = '',
                       subtitle: str = '',
                       charts_on_separate_lines: bool = False) -> None:
  """Makes the html report with all the charts inserted.

  Args:
    charts: A list of altair chart objects.
    html_output: A writable file object.
    title: The title to show at the top of the report.
    subtitle: The subtitle to show just below the title on the report.
    charts_on_separate_lines: Put charts on separate lines. If false, charts
      will set next to each other as space allows and flow to the next line,
      similar to text wrapping.

  Returns:
      None. Writes into the html_output file object.
  """

  chart_div_style = 'style="display:block"' if charts_on_separate_lines else ''

  # Start the HTML document.
  html_string = (
      f'<!DOCTYPE html>\n<html>\n<head>\n'
      # Add dependencies vega and vega-lite, which render the altair charts.
      f'<script type="text/javascript" src="{VEGA_URL}/vega@5"></script>\n'
      f'<script type="text/javascript" src="{VEGA_URL}/vega-lite@4.8.1">'
      '</script>\n'
      f'<script type="text/javascript" src="{VEGA_URL}/vega-embed@6">'
      '</script>\n'
      # Add styles (CSS).
      f'{CSS_STYLES}'
      '</head>\n<body>'
      # Titles
      f'<h1>{title}</h1>\n'
      f'<h2>{subtitle}</h2>\n'
      # Make a div containing all the charts.
      '<div>')
  for chart in charts:
    chart_id = chart['id']
    html_string += (f'<div class="chart-container" {chart_div_style} '
                    f'id="vis_{chart_id}"></div>\n')
  # End the chart container and star the JavaScript section.
  html_string += ('</div>' '<script>\n')

  # Add JSON vega specs and hook them up to the divs with VegaEmbed.
  for chart in charts:
    chart_id = chart['id']
    chart_json = chart['chart'].to_json()
    download_filename = '{}_{}'.format(title.replace(' ', '_'), chart['id'])
    embed_options = {'mode': 'vega-lite', 'downloadFileName': download_filename}
    html_string += (
        f'var spec_{chart_id} = {chart_json};\n'
        f'vegaEmbed("#vis_{chart_id}", spec_{chart_id}, {embed_options})\n')
  html_string += (
      '</script>\n'
      # Close HTML document.
      '</body></html>')

  html_output.write(html_string)
