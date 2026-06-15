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
"""This script is used to collect runtimes from DeepTrio logs."""

import re
import sys


LOG_FILE = 'deeptrio.log'
TIMES_IN_ORDER = [
    'make_examples',
    'call_variants_child',
    'call_variants_parent1',
    'call_variants_parent2',
    'postprocess_variants_child',
    'postprocess_variants_parent1',
    'postprocess_variants_parent2',
    'vcf_stats_report_child',
    'vcf_stats_report_parent1',
    'vcf_stats_report_parent2',
    'total',
]


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


uid = sys.argv[1]
times = {}
with open(LOG_FILE, 'r') as f:
  real_times = [
      ms_to_secs(x.split()[1]) for x in f.readlines() if x.startswith('real\t')
  ]
  for step, real_time in zip(TIMES_IN_ORDER, real_times):
    times[step] = real_time

with open(f'{uid}.runtimes.tsv', 'w') as f:
  for log_name, time in times.items():
    f.write(f'{uid}\t{log_name}\t{secs_to_ms(time)}\t{time}\n')
