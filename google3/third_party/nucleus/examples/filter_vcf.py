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

"""Writes all the variants in a VCF file with a quality greater than 3.01."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

from absl import app

from deepvariant.util.io import vcf


def main(argv):
  if len(argv) != 3:
    print('Usage: {} <input_vcf> <output_vcf>'.format(argv[0]))
    sys.exit(-1)
  in_vcf = argv[1]
  out_vcf = argv[2]

  # Please try to keep the following part in sync with the documenation in
  # g3doc/overview.md.
  with vcf.VcfReader(in_vcf, use_index=False) as reader:
    print('Sample names in VCF: ', ' '.join(reader.header.sample_names))
    with vcf.VcfWriter(out_vcf, header=reader.header) as writer:
      for variant in reader:
        if variant.quality > 3.01:
          writer.write(variant)


if __name__ == '__main__':
  app.run(main)
