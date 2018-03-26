# Copyright 2018 Google Inc.
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

"""Print an ASCII art pileup image."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

from absl import app

from deepvariant.util.io import sam
from deepvariant.util import ranges

ANSI_BOLD = '\033[1m'
ANSI_RED = '\033[91m'
ANSI_OFF = '\033[0m'

# redacted
# top if a reference fasta file is supplied.


def print_read(left_pos, start, highlight_position, seq):
  s = ' ' * (start - left_pos)
  i = highlight_position - start
  j = i + 1
  s += seq[:i] + ANSI_BOLD + ANSI_RED + seq[i:j] + ANSI_OFF + seq[j:]
  print(s)


def main(argv):
  if len(argv) != 3:
    print('Usage: {} <input_sam> <chromosome>:<position>'.format(argv[0]))
    sys.exit(-1)
  in_sam = argv[1]
  r = ranges.parse_literal(argv[2])
  position = r.start

  with sam.SamReader(in_sam) as sam_reader:
    reads = sam_reader.query(r)
    pos_seq_pairs = sorted(
        (read.alignment.position.position, read.aligned_sequence)
        for read in reads)
    if not pos_seq_pairs:
      print('No overlapping reads found for', argv[2])
      sys.exit(0)

    left_position = pos_seq_pairs[0][0]
    for start, seq in pos_seq_pairs:
      print_read(left_position, start, position, seq)


if __name__ == '__main__':
  app.run(main)
