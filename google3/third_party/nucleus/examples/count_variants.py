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
"""Counts variants in a VCF, both per length and per chromosome."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import sys

from absl import app
from deepvariant.util.io import vcf
from deepvariant.util import variant_utils


def main(argv):
  if len(argv) != 2:
    print('Usage: {} <input_vcf>'.format(argv[0]))
    sys.exit(-1)
  in_vcf = argv[1]

  total = 0
  by_type = collections.defaultdict(int)
  by_ref = collections.defaultdict(int)

  with vcf.VcfReader(in_vcf, use_index=False) as reader:
    for variant in reader:
      total += 1
      by_type[variant_utils.variant_type(variant)] += 1
      by_ref[variant.reference_name] += 1

  print('# variants: {}'.format(total))
  print('# ref variants: {}'.format(by_type[variant_utils.VariantType.ref]))
  print('# SNP variants: {}'.format(by_type[variant_utils.VariantType.snp]))
  print('# indel variants: {}'.format(by_type[variant_utils.VariantType.indel]))
  for k, v in sorted(by_ref.items()):
    print('# variants in {}: {}'.format(k, v))


if __name__ == '__main__':
  app.run(main)
