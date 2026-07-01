#!/usr/bin/env python3
# Copyright 2026 Google LLC.
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
"""Aggregates runtimes, md5sums, and happy results from deepvariant runs."""

import glob

import pandas as pd

HAPPY_COLUMNS = [
    'uid',
    'sample',
    'dataset_name',
    'Type',
    'TRUTH.TOTAL',
    'TRUTH.TP',
    'TRUTH.FN',
    'QUERY.TOTAL',
    'QUERY.FP',
    'METRIC.Recall',
    'METRIC.Precision',
    'METRIC.F1_Score',
]


def df_to_markdown(dataframe: pd.DataFrame) -> str:
  """Convert a DataFrame to a markdown table string with aligned columns."""
  cols = list(dataframe.columns)
  # Convert all values to string. NaNs should have been handled by fillna
  # earlier if we wanted them to be empty.
  str_df = dataframe.astype(str)
  widths = {}
  for col in cols:
    max_val_width = str_df[col].str.len().max()
    # Handle empty dataframe case
    if pd.isna(max_val_width):
      max_val_width = 0
    max_col_width = max(len(str(col)), max_val_width)
    widths[col] = max_col_width

  header = '| ' + ' | '.join(str(c).ljust(widths[c]) for c in cols) + ' |'
  separator = (
      '| ' + ' | '.join('---'.ljust(widths[c], '-') for c in cols) + ' |'
  )
  rows = []
  for _, row in str_df.iterrows():
    rows.append('| ' + ' | '.join(row[c].ljust(widths[c]) for c in cols) + ' |')
  return '\n'.join([header, separator] + rows)


def move_columns_to_left(
    input_df: pd.DataFrame, cols: 'str | list[str]'
) -> pd.DataFrame:
  """Moves columns to the left of a DataFrame, preserving order."""
  cols = [cols] if isinstance(cols, str) else cols
  missing = [c for c in cols if c not in input_df.columns]
  if missing:
    raise ValueError(f'Columns not found: {missing}')
  remaining = [c for c in input_df.columns if c not in cols]
  return input_df[cols + remaining]


def secs_to_hruntime(secs: float) -> str:
  """Convert seconds to a human-readable duration string."""
  hours = int(secs // 3600)
  minutes = int((secs % 3600) // 60)
  seconds = int(secs % 60)
  parts = []
  if hours:
    parts.append(f'{hours}h')
  if minutes:
    parts.append(f'{minutes}m')
  parts.append(f'{seconds}s')
  return ' '.join(parts)


# ====================#
# Aggregate Runtimes #
# ====================#
SORT_ORDER = {
    'make_examples': 0,
    'call_variants': 1,
    'postprocess_variants': 2,
    'vcf_stats': 3,
    'labeled_examples_to_vcf': 4,
}

if glob.glob('*.runtimes.tsv'):
  # Read in runtime tsv
  dfs = []
  for f in sorted(glob.glob('*.runtimes.tsv')):
    df_part = pd.read_csv(
        f,
        sep='\t',
        header=None,
        names=['group', 'stage', 'hruntime', 'runtime'],
    )
    dfs.append(df_part)
  df_raw = pd.concat(dfs, ignore_index=True)

  split_group = df_raw['group'].str.rsplit('_', n=1)
  df_raw['uid'] = split_group.str[0].str.replace(
      r'-trial[0-9]+$', '', regex=True
  )
  df_raw['sample'] = split_group.str[1]

  # Aggregate by stage
  df_stage = (
      df_raw.groupby(['uid', 'sample', 'stage'])
      .agg(
          mean_runtime=('runtime', 'mean'),
          std_runtime=('runtime', 'std'),
          n_trials=('uid', 'count'),
      )
      .reset_index()
  )
  df_stage['mean_runtime'] = df_stage['mean_runtime'].round(2)
  df_stage['std_runtime'] = df_stage['std_runtime'].round(3)
  df_stage['mean_hruntime'] = df_stage['mean_runtime'].apply(secs_to_hruntime)
  df_stage['sort_order'] = df_stage['stage'].map(SORT_ORDER)

  # Aggregate total runtime per uid, sample
  df_no_vcf = df_raw[df_raw['stage'] != 'vcf_stats']
  df_group_totals = (
      df_no_vcf.groupby(['group', 'uid', 'sample'])['runtime']
      .sum()
      .round(2)
      .reset_index(name='runtime_sum')
  )
  df_total = (
      df_group_totals.groupby(['uid', 'sample'])
      .agg(
          mean_runtime=('runtime_sum', 'mean'),
          std_runtime=('runtime_sum', 'std'),
          n_trials=('uid', 'count'),
      )
      .reset_index()
  )
  df_total['mean_runtime'] = df_total['mean_runtime'].round(2)
  df_total['std_runtime'] = df_total['std_runtime'].round(3)
  df_total['mean_hruntime'] = df_total['mean_runtime'].apply(secs_to_hruntime)
  df_total['stage'] = 'total'
  df_total['sort_order'] = 999

  # Order by group, stage
  df_concat = pd.concat([df_stage, df_total], ignore_index=True)
  df = (
      df_concat.sort_values(['uid', 'sample', 'sort_order'])
      .drop(columns=['sort_order'])
      .reset_index(drop=True)
  )

  df.to_csv('runtimes.tsv', sep='\t', index=False)
  with open('runtimes.md', 'w') as f:
    f.write(df_to_markdown(df))

# ============================#
# Check uniqueness of hashes #
# ============================#
if glob.glob('*.md5sum.txt'):
  md5_frames = []
  for fpath in sorted(glob.glob('*.md5sum.txt')):
    df_part = pd.read_csv(
        fpath,
        sep=r'\s+',
        header=None,
        names=['md5sum', 'fname'],
        engine='python',
    )
    df_part['source_file'] = fpath
    md5_frames.append(df_part)

  md5 = pd.concat(md5_frames, ignore_index=True)
  md5['group'] = md5['source_file'].str.replace(
      r'\.md5sum\.txt$', '', regex=True
  )
  split_group = md5['group'].str.rsplit('_', n=1)
  md5['uid'] = split_group.str[0].str.replace(r'-[0-9]+', '', regex=True)
  md5['sample'] = split_group.str[1]
  md5 = md5.drop(columns=['group', 'source_file'])
  md5 = md5.drop_duplicates()
  md5 = md5.sort_values(['uid', 'sample']).reset_index(drop=True)
  md5 = md5.drop_duplicates(subset=['md5sum'], keep='last')

  md5_reordered = move_columns_to_left(md5, ['uid', 'sample'])

  md5_reordered.to_csv('md5sums.tsv', sep='\t', index=False)
  with open('md5sums.md', 'w') as f:
    f.write(df_to_markdown(md5_reordered))


# =======================#
# Combine Happy Results #
# =======================#

if glob.glob('*.happy.summary.csv') and glob.glob('*.info.tsv'):
  # Read in happy dataset info
  info_frames = []
  for fpath in sorted(glob.glob('*.info.tsv')):
    df_part = pd.read_csv(
        fpath,
        sep='\t',
        header=None,
        names=['uid', 'sample', 'dataset_name', 'source_file'],
    )
    info_frames.append(df_part)
  happy_info = pd.concat(info_frames, ignore_index=True)
  happy_info['uid'] = happy_info['uid'].astype(str)
  happy_info['dataset_name'] = happy_info['dataset_name'].fillna('')

  happy_frames = []
  for fpath in sorted(glob.glob('*.happy.summary.csv')):
    df_part = pd.read_csv(fpath)
    df_part['source_file'] = fpath
    happy_frames.append(df_part)
  happy = pd.concat(happy_frames, ignore_index=True)

  happy = happy.merge(happy_info, on='source_file')
  happy = happy.sort_values('uid')
  happy['group'] = happy['uid'].str.replace(r'-trial[0-9]+$', '', regex=True)
  happy = happy.drop_duplicates(
      subset=['group', 'sample', 'dataset_name', 'Filter', 'Type'],
      keep='first',
  )
  happy = happy.sort_values(['Filter', 'uid', 'dataset_name', 'Type'])
  happy = happy.drop(columns=['source_file', 'group'])

  happy_reordered = move_columns_to_left(
      happy, ['uid', 'sample', 'dataset_name']
  )
  happy_reordered.to_csv('happy.summary.tsv', sep='\t', index=False)

  happy_pass = happy_reordered[happy_reordered['Filter'] == 'PASS']
  happy_md = df_to_markdown(
      happy_pass[HAPPY_COLUMNS].rename(
          columns={
              'METRIC.Recall': 'Recall',
              'METRIC.Precision': 'Precision',
              'METRIC.F1_Score': 'F1_Score',
          }
      )
  )

  with open('happy.summary.md', 'w') as f:
    f.write(happy_md)

  # =====================#
  # Combine Happy Errors #
  # =====================#
  errors_pass = happy_reordered[happy_reordered['Filter'] == 'PASS'].copy()
  errors_pass['errors'] = errors_pass['TRUTH.FN'] + errors_pass['QUERY.FP']
  errors_df = errors_pass.pivot_table(
      index=['uid', 'sample', 'dataset_name'],
      columns='Type',
      values='errors',
      aggfunc='sum',
  ).reset_index()
  errors_df.columns.name = None

  # For very small tests, the INDEL or SNP column may not be present.
  snp_col = errors_df.get('SNP', 0)
  indel_col = errors_df.get('INDEL', 0)
  errors_df['total_errors'] = snp_col + indel_col
  errors_df = errors_df.rename(
      columns={'SNP': 'SNP_errors', 'INDEL': 'INDEL_errors'}
  )
  with open('happy.error_counts.md', 'w') as f:
    f.write(df_to_markdown(errors_df))

  errors_df.to_csv('happy.error_counts.tsv', sep='\t', index=False)

# ===============================#
# Combine Happy Extended Results #
# ===============================#

if glob.glob('*.happy.extended.csv'):
  ext_frames = []
  for fpath in sorted(glob.glob('*.happy.extended.csv')):
    df_part = pd.read_csv(fpath)
    df_part['source_file'] = fpath
    ext_frames.append(df_part)
  happy_ext = pd.concat(ext_frames, ignore_index=True)

  happy_ext['uid'] = (
      happy_ext['source_file']
      .str.replace(r'\.happy\.extended\.csv$', '', regex=True)
      .str.split('_')
      .str[0]
  )
  happy_ext['sample'] = (
      happy_ext['source_file']
      .str.replace(r'\.happy\.extended\.csv$', '', regex=True)
      .str.split('_')
      .str[1]
  )
  happy_ext = happy_ext.sort_values('uid')
  happy_ext['group'] = happy_ext['uid'].str.replace(
      r'-trial[0-9]+$', '', regex=True
  )
  happy_ext = happy_ext.drop_duplicates(
      subset=['group', 'sample', 'Filter', 'Type', 'Subtype', 'Subset'],
      keep='first',
  )
  happy_ext = happy_ext.sort_values(['Filter', 'uid', 'Type'])
  happy_ext = happy_ext.drop(columns=['source_file', 'group'])

  happy_ext = move_columns_to_left(happy_ext, ['uid', 'sample'])
  happy_ext.to_csv('happy.extended.tsv', sep='\t', index=False)
