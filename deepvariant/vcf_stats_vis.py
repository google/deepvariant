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
"""Create a visual report from a VCF file."""

import io
import json
import os

import altair as alt
import numpy as np
import pandas as pd
import tensorflow as tf

# Altair uses a lot of method chaining, such as
# chart.mark_bar().encode(...).properties(...), so allowing backslash
# continuation to break this into separate lines makes the code more readable.
# pylint: disable=g-backslash-continuation

OLD_LIB_BASE_URL = 'https://cdn.jsdelivr.net/npm//'
NEW_LIB_BASE_URL = 'https://storage.googleapis.com/deepvariant/lib/vega/'

VEGA_VERSION = '5'
VEGA_LITE_VERSION = '3.4.0'
VEGA_EMBED_VERSION = '4'

# "pretty" genotype strings:
REF = 'Ref (0/0)'
HET = 'Het (0/x)'
HOM = 'Hom (x/x)'
UNCALLED = 'Uncalled (./.)'
HET_BOTH = 'Het - two variants (x/y)'

# Establish ordering of bases to keep it consistent
BASES = ['A', 'G', 'T', 'C']

BAR_COLOR_DEPTH = '#4a1486'
BAR_COLOR_QUAL = '#0c2c84'
BAR_COLOR_GQ = '#0c2c84'

BIALLELIC_SNP = 'Biallelic_SNP'
BIALLELIC_INSERTION = 'Biallelic_Insertion'
BIALLELIC_DELETION = 'Biallelic_Deletion'
BIALLELIC_MNP = 'Biallelic_MNP'
MULTIALLELIC_SNP = 'Multiallelic_SNP'
MULTIALLELIC_INSERTION = 'Multiallelic_Insertion'
MULTIALLELIC_DELETION = 'Multiallelic_Deletion'
MULTIALLELIC_COMPLEX = 'Multiallelic_Complex'
REFCALL = 'RefCall'

ordered_variant_type_labels = [
    BIALLELIC_INSERTION,
    BIALLELIC_DELETION,
    BIALLELIC_SNP,
    BIALLELIC_MNP,
    MULTIALLELIC_INSERTION,
    MULTIALLELIC_DELETION,
    MULTIALLELIC_SNP,
    MULTIALLELIC_COMPLEX,
    REFCALL,
]


def _dict_to_dataframe(dictionary):
  """Turn a dict object into a dataframe of with label and value columns."""
  df = pd.DataFrame(
      {'label': list(dictionary.keys()), 'value': list(dictionary.values())}
  )
  return df


def _prettify_genotype(genotype):
  """Get more human-readable display name and grouping for a given genotype."""
  pretty = genotype
  group = 'others'
  alleles = json.loads(genotype)
  if len(alleles) == 2:
    g1, g2 = sorted(alleles)
    if g1 == 0 and g2 == 0:
      pretty = REF
      group = 'main'
    elif g1 == -1 and g2 == -1:
      pretty = UNCALLED
    elif g1 == 0 and g2 > 0:
      pretty = HET
      group = 'main'
    elif g1 == g2:
      pretty = HOM
      group = 'main'
    else:
      pretty = HET_BOTH
  return pretty, group


def _build_type_chart(variant_type_counts):
  """Create a chart of the counts of each variant type."""
  width = 400
  height = 200
  title = 'Variant types'
  variant_type_data = _dict_to_dataframe(variant_type_counts)
  type_chart = _placeholder_for_empty_chart(
      'No entries in VCF', width=width, height=height, title=title
  )
  if not variant_type_data.empty:
    bars = (
        alt.Chart(variant_type_data)
        .mark_bar()
        .encode(
            x=alt.X(
                'label',
                title=None,
                sort=ordered_variant_type_labels,
                axis=alt.Axis(labelAngle=-45),
            ),
            y=alt.Y('value', axis=alt.Axis(title='Count', format='s')),
            tooltip=alt.Tooltip('value', format='.4s'),
            color=alt.Color(
                'label',
                legend=None,
                scale=alt.Scale(
                    scheme='set1', domain=ordered_variant_type_labels
                ),
            ),
        )
    )
    labels = bars.mark_text(dy=-5).encode(text=alt.Text('value', format='.4s'))
    type_chart = (bars + labels).properties(
        width=width, height=height, title=title
    )
  return type_chart


def _build_qual_histogram(data):
  """Create the Quality(QUAL) histogram."""
  width = 200
  height = 200
  title = 'Quality score'
  qual_data = pd.DataFrame(data)
  qual_histogram = _placeholder_for_empty_chart(
      'No entries in VCF', width=width, height=height, title=title
  )
  if not qual_data.empty:
    # s = bin_start, e = bin_end, c = count
    domain = [min(0, data[0]['s']), max(150, data[-1]['e'])]
    qual_histogram = (
        alt.Chart(qual_data)
        .mark_bar(color=BAR_COLOR_QUAL)
        .encode(
            x=alt.X('s', title='QUAL', scale=alt.Scale(domain=domain)),
            x2='e',
            y=alt.Y('c', title='Count', stack=True, axis=alt.Axis(format='s')),
        )
        .properties(width=width, height=height, title=title)
        .interactive(bind_y=False)
    )
  return qual_histogram


def _build_gq_histogram(data):
  """Create the Genotype quality (GQ) histogram."""
  # gq = genotype quality, found at :GQ: in FORMAT column of VCF
  width = 200
  height = 200
  title = 'Genotype quality'
  gq_data = _integer_counts_to_histogram(data)
  gq_histogram = _placeholder_for_empty_chart(
      'No entries in VCF with GQ', width=width, height=height, title=title
  )
  if not gq_data.empty:
    # standardize x-axis limits across reports
    domain = [min(0, data[0][0]), max(150, data[-1][0])]
    # s = bin_start, e = bin_end, c = count
    gq_histogram = (
        alt.Chart(gq_data)
        .mark_bar(color=BAR_COLOR_GQ)
        .encode(
            x=alt.X('s', title='GQ', scale=alt.Scale(domain=domain)),
            x2='e',
            y=alt.Y('c', title='Count', stack=True, axis=alt.Axis(format='s')),
        )
        .properties(width=width, height=height, title=title)
        .interactive(bind_y=False)
    )
  return gq_histogram


def _build_vaf_histograms(histogram_json):
  """Create VAF histograms split by genotype."""
  guides = {REF: 0, HET: 0.5, HOM: 1}
  dfs = []
  for key in histogram_json:
    g = pd.DataFrame(histogram_json[key])
    pretty, group = _prettify_genotype(key)
    g['GT'] = pretty  # pretty genotype name
    g['g'] = group  # main/other genotypes
    g['l'] = guides.get(pretty, None)  # vertical line as guide
    dfs.append(g)

  hist_data = pd.concat(dfs)
  main_hist_data = hist_data[hist_data['g'] == 'main']
  other_hist_data = hist_data[hist_data['g'] == 'others']

  # Main genotypes (ref, het, hom-alt)
  # Histogram bars themselves
  # s = bin_start, e = bin_end, c = count
  bars = (
      alt.Chart(main_hist_data)
      .mark_bar()
      .encode(
          x=alt.X('s', title='VAF'),
          x2='e',
          y=alt.Y('c', title='Count', stack=True, axis=alt.Axis(format='s')),
      )
  )
  # Vertical lines
  guides = alt.Chart(main_hist_data).mark_rule().encode(x='l')
  # Facet into 3 plots by genotype
  vaf_histograms = (
      (bars + guides)
      .properties(width=200, height=200)
      .facet(
          column=alt.Column('GT', title='Main genotypes', sort=[REF, HET, HOM])
      )
      .resolve_scale(y='independent')
  )

  # Other genotypes (uncalled, het with two alt alleles)
  # s = bin_start, e = bin_end, c = count
  other_vaf_histograms = (
      alt.Chart(other_hist_data)
      .mark_bar()
      .encode(
          x=alt.X('s', title='VAF'),
          x2='e',
          y=alt.Y('c', title='Count', stack=True, axis=alt.Axis(format='s')),
          column=alt.Column('GT', title='Other genotypes'),
      )
      .properties(width=150, height=150)
      .resolve_scale(y='independent')
  )
  return vaf_histograms, other_vaf_histograms


def _placeholder_for_empty_chart(
    text_to_display, width=100, height=100, title=''
):
  chart = (
      alt.Chart({'values': [{'placeholder': text_to_display}]})
      .mark_text(size=14)
      .encode(text='placeholder:N')
      .properties(width=width, height=height, title=title)
  )
  return chart


def _build_base_change_chart(data):
  """Create the base change chart."""
  width = 100
  height = 200
  placeholder_width = (4 * width) + 80  # 4 charts, plus constant spacing
  title = 'Biallelic base changes from reference'
  base_change_data = pd.DataFrame(data, columns=['ref', 'alt', 'count'])

  base_change_chart = _placeholder_for_empty_chart(
      'No biallelic SNPs', width=placeholder_width, height=height, title=title
  )
  if not base_change_data.empty:
    bars = (
        alt.Chart(base_change_data)
        .mark_bar()
        .encode(
            x=alt.X('alt', title='to alt'),
            y=alt.Y('count', title='Count', axis=alt.Axis(format='s')),
            color=alt.Color(
                'alt',
                legend=None,
                sort=BASES,
                scale=alt.Scale(scheme='category20', domain=BASES),
            ),
            tooltip=alt.Tooltip('count', format='.4s'),
        )
    )
    labels = bars.mark_text(dy=-5, fontWeight='bold').encode(text='alt')

    base_change_chart = (
        (bars + labels)
        .properties(width=100, height=200)
        .facet(column=alt.Column('ref', title=title, sort=BASES))
    )

  return base_change_chart


def _integer_counts_to_histogram(num_count_pairs):
  """Turn paired numbers and their counts into data for a histogram.

  This centers the bars on the exact integer for clarity. For example, the bar
  for 3 is centered on 3 instead of being between 3 and 4 as in numpy's default
  histogram.

  Args:
    num_count_pairs: list of [num, count] pairs

  Returns:
    a pandas dataframe with num, count (bin count), s (bin start), e (bin end)
  """
  histogram_data = pd.DataFrame(num_count_pairs, columns=['num', 'c'])
  # For a proper histogram, use s and e to force each bar to cover
  # exactly one integer position:
  histogram_data['s'] = histogram_data['num'] - 0.5
  histogram_data['e'] = histogram_data['num'] + 0.5
  histogram_data = histogram_data.drop(columns=['num'])
  return histogram_data


def _build_indel_size_chart(data):
  """Create the indel size chart."""
  width = 400
  height = 100
  placeholder_height = (2 * height) + 20  # 2 charts, plus spacing
  title = 'Biallelic indel size distribution'
  ordered_labels = ['Insertion', 'Deletion']
  indel_size_data = _integer_counts_to_histogram(data)
  indel_size_data['type'] = np.where(
      indel_size_data['s'] > 0, 'Insertion', 'Deletion'
  )

  indel_size_chart = _placeholder_for_empty_chart(
      'No biallelic indels', width=width, height=placeholder_height, title=title
  )

  if not indel_size_data.empty:
    indels_linear = (
        alt.Chart(indel_size_data)
        .mark_bar()
        .encode(
            x=alt.X('s', title='size'),
            x2='e',
            y=alt.Y('c', title='Count', axis=alt.Axis(format='s')),
            color=alt.Color(
                'type', sort=ordered_labels, scale=alt.Scale(scheme='set1')
            ),
        )
        .properties(width=400, height=100, title=title)
        .interactive(bind_y=False)
    )

    indel_log = (
        alt.Chart(indel_size_data)
        .mark_bar()
        .encode(
            x=alt.X('s', title='size'),
            x2='e',
            y=alt.Y(
                'c',
                title='Count',
                axis=alt.Axis(format='s'),
                scale=alt.Scale(type='log', base=10),
            ),
            color=alt.Color(
                'type', sort=ordered_labels, scale=alt.Scale(scheme='set1')
            ),
        )
        .properties(width=400, height=100)
        .interactive(bind_y=False)
    )

    indel_size_chart = alt.vconcat(indels_linear, indel_log).resolve_scale(
        color='shared'
    )
  return indel_size_chart


def _build_depth_histogram(data):
  """Build histogram with depth (DP)."""
  width = 200
  height = 200
  title = 'Depth'
  depth_data = _integer_counts_to_histogram(data)
  depth_histogram = _placeholder_for_empty_chart(
      'No entries in VCF with DP', width=width, height=height, title=title
  )
  if not depth_data.empty:
    # s = bin_start, e = bin_end, c = count
    depth_histogram = (
        alt.Chart(depth_data)
        .mark_bar(color=BAR_COLOR_DEPTH)
        .encode(
            x=alt.X('s', title='Depth'),
            x2='e',
            y=alt.Y('c', title='Count', stack=True, axis=alt.Axis(format='s')),
        )
        .properties(width=width, height=height, title=title)
        .interactive(bind_y=False)
    )
  return depth_histogram


def _build_tt_chart(titv_counts):
  """Built chart showing counts of transitions and transversions."""
  width = 150
  height = 200

  ti = titv_counts['Transition']
  tv = titv_counts['Transversion']
  # Show TiTv ratio with fallback to avoid division by 0
  titv_ratio = '%.2f' % (float(ti) / tv) if tv > 0 else '%d / 0' % (ti)
  title = 'Biallelic Ti/Tv ratio: %s' % (titv_ratio)

  tt_chart = _placeholder_for_empty_chart(
      'No biallelic SNPs', width=width, height=height, title=title
  )
  tt_labels = ['Transition', 'Transversion']
  if sum([titv_counts[k] for k in titv_counts]) > 0:
    tt_data = _dict_to_dataframe(titv_counts)
    bars = (
        alt.Chart(tt_data)
        .mark_bar()
        .encode(
            x=alt.X(
                'label', sort=tt_labels, axis=alt.Axis(title=None, labelAngle=0)
            ),
            y=alt.Y('value', axis=alt.Axis(title='Count', format='s')),
            tooltip=alt.Tooltip('value', format='.4s'),
            color=alt.Color(
                'label',
                legend=None,
                sort=tt_labels,
                scale=alt.Scale(scheme='teals', domain=tt_labels),
            ),
        )
    )
    labels = bars.mark_text(dy=-5).encode(text=alt.Text('value', format='.4s'))
    tt_chart = (bars + labels).properties(
        title=title, width=width, height=height
    )
  return tt_chart


def _build_all_charts(vis_data, title=''):
  """Build all charts and combine into a single interface."""

  # Row 1
  type_chart = _build_type_chart(vis_data['variant_type_counts'])
  depth_chart = _build_depth_histogram(vis_data['depth_histogram'])
  qual_histogram = _build_qual_histogram(vis_data['qual_histogram'])
  gq_histogram = _build_gq_histogram(vis_data['gq_histogram'])
  row1 = alt.hconcat(
      type_chart, depth_chart, qual_histogram, gq_histogram
  ).resolve_scale(color='independent')

  # Row 2
  vaf_histograms, other_vaf_histograms = _build_vaf_histograms(
      vis_data['vaf_histograms_by_genotype']
  )
  row2 = alt.hconcat(vaf_histograms, other_vaf_histograms)

  # Row 3
  base_change_chart = _build_base_change_chart(vis_data['base_changes'])
  indel_size_chart = _build_indel_size_chart(vis_data['indel_sizes'])
  tt_chart = _build_tt_chart(vis_data['titv_counts'])
  row3 = alt.hconcat(
      base_change_chart, tt_chart, indel_size_chart
  ).resolve_scale(color='independent')

  # Putting it all together
  all_charts = alt.vconcat(row1, row2, row3)

  all_charts = (
      all_charts.properties(title=title, spacing=70)
      .configure_header(labelFontSize=16, titleFontSize=20)
      .configure_title(fontSize=20)
  )
  return all_charts


def _altair_chart_to_html(altair_chart, download_filename):
  """Write to a temporary string stand-in for the file to replace import URLs.

  Args:
    altair_chart: a chart object made by Altair.
    download_filename: string filename base for when users export images.

  Returns:
    HTML in string format.
  """
  temp_writer = io.StringIO()
  altair_chart.save(
      temp_writer,
      format='html',
      embed_options={'downloadFileName': download_filename},
      vegalite_version=VEGA_LITE_VERSION,
      vega_version=VEGA_VERSION,
      vegaembed_version=VEGA_EMBED_VERSION,
  )
  temp_html_string = temp_writer.getvalue()
  html_with_new_cdn = temp_html_string.replace(
      OLD_LIB_BASE_URL, NEW_LIB_BASE_URL
  )
  return html_with_new_cdn


def _save_html(basename, all_charts):
  """Save Altair chart as an HTML file."""
  output_path = basename + '.visual_report.html'
  image_download_filename = os.path.basename(basename) + '.visual_report'
  html_string = _altair_chart_to_html(
      altair_chart=all_charts, download_filename=image_download_filename
  )

  with tf.io.gfile.GFile(output_path, 'w') as writer:
    writer.write(html_string)


def create_visual_report(basename, vis_data, title=''):
  """Build visual report with several charts."""
  all_charts = _build_all_charts(vis_data, title=title)
  _save_html(basename, all_charts)
