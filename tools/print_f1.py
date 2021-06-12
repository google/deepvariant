# Copyright 2018 Google LLC.
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
r"""Parse and extract metrics from *.metrics files."""

import argparse
import json
import logging
from os import listdir
from os.path import isfile
from os.path import join
import re


def parse_cmdline(argv):
  """Parse the commandline."""
  parser = argparse.ArgumentParser()

  parser.add_argument(
      '--metrics_dir', help='Path to the directory with metrics files.')

  known_args, _ = parser.parse_known_args(argv)

  return known_args


def extract_checkpoint_number_from_metrics_filename(filename):
  match = re.search(r'ckpt-([\d]*)\.metrics', filename)
  if match:
    return int(match.group(1))


def read_metrics_file(path):
  """Reads metrics f and outputs metrics in a dict."""
  with open(path) as f:
    metrics = {
        key.replace('/', '_'): float(value)
        for key, value in json.loads(f.read()).items()
    }
  metrics['checkpoint'] = extract_checkpoint_number_from_metrics_filename(path)
  metrics['F1_All'] = 2 * metrics['TPs_All'] / (
      2 * metrics['TPs_All'] + metrics['FNs_All'] + metrics['FPs_All'])
  metrics['TPs+FNs_All'] = metrics['TPs_All'] + metrics['FNs_All']
  return metrics


def main(argv=None):
  """Main entry point."""
  known_args = parse_cmdline(argv)
  metrics_dir = known_args.metrics_dir
  metrics_files = [
      join(metrics_dir, f)
      for f in listdir(metrics_dir)
      if isfile(join(metrics_dir, f))
  ]
  metrics = [read_metrics_file(f) for f in metrics_files]
  for m in metrics:
    print('%s\t%s\t%s' % (m['checkpoint'], m['TPs+FNs_All'], m['F1_All']))


if __name__ == '__main__':
  logging.getLogger().setLevel(logging.INFO)
  main()
