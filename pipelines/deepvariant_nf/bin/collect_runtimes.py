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
"""This script is used to collect runtimes from DeepVariant/DeepSomatic logs."""

import enum
import glob
import os
import re
import sys


class Pipeline(enum.Enum):
  """Pipeline type."""

  DEEPVARIANT = 'deepvariant'
  ORACLE_INFERENCE = 'oracle_inference'


uid = sys.argv[1]
pipeline = Pipeline.DEEPVARIANT

if len(sys.argv) > 2:
  try:
    pipeline = Pipeline(sys.argv[2])
  except ValueError:
    print(f'Unknown pipeline type: {sys.argv[2]}', file=sys.stderr)
    sys.exit(1)


def ms_to_secs(ms_str):
  """Convert XmYs to seconds."""
  m = re.match('(.+)m(.+)s', ms_str)
  if m is None:
    raise ValueError(
        f"Invalid time format: {ms_str}. Expected format like 'XmYs'."
    )
  return int(m.groups()[0]) * 60 + float(m.groups()[1])


def secs_to_ms(secs):
  """Convert secs to XmYs."""
  m = int(secs // 60)
  s = secs - 60 * m
  return f'{m}m{s:.2f}s'

LOG_FILES_BY_PIPELINE = {
    Pipeline.DEEPVARIANT: {
        'make_examples': (
            glob.glob('make_examples*log')[0]
            if glob.glob('make_examples*log')
            else 'make_examples.log'
        ),
        'call_variants': 'call_variants.log',
        'postprocess_variants': 'postprocess_variants.log',
        'vcf_stats': 'vcf_stats_report.log',
    },
    Pipeline.ORACLE_INFERENCE: {
        'make_examples': 'make_examples.log',
        'labeled_examples_to_vcf': 'labeled_examples_to_vcf.log',
    },
}

times = {}
log_files = LOG_FILES_BY_PIPELINE[pipeline]
for log_name, log_file in log_files.items():
  if log_file == 'vcf_stats_report.log' and not os.path.exists(log_file):
    # vcf_stats_report.log is not always generated.
    print(f'Warning: Log file not found: {log_file}', file=sys.stderr)
    continue
  with open(log_file, 'r') as f:
    real_time = [
        ms_to_secs(x.split()[1])
        for x in f.readlines()
        if x.startswith('real\t')
    ][0]
    times[log_name] = real_time

with open(f'{uid}.runtimes.tsv', 'w') as f:
  for log_name, time in times.items():
    f.write(f'{uid}\t{log_name}\t{secs_to_ms(time)}\t{time}\n')
