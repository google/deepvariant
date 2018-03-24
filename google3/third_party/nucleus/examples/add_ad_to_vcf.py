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

"""This example program adds the AD info field to a VCF file.

It assumes that the AD field of the individual variant calls is already
populated.

Sample usage:
  $ add_ad_to_vcf input.vcf.gz output.vcf.gz
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

from absl import app

from deepvariant.util.io import vcf
from deepvariant.util import variant_utils
from deepvariant.util import variantcall_utils
from deepvariant.util import vcf_constants


def get_variant_ad(variant):
  """Returns the allele depth for the Variant, calculated across its calls."""
  num_alleles = len(variant.alternate_bases) + 1
  call_ads = [variantcall_utils.get_format(vc, 'AD') for vc in
              variant.calls]
  assert(len(call_ad) == num_alleles for call_ad in call_ads)
  return [sum(call_ad[i] for call_ad in call_ads) for i in xrange(num_alleles)]


def main(argv):
  if len(argv) != 3:
    print('Usage: %s <input_vcf> <output_vcf>' % argv[0])
    sys.exit(-1)
  in_vcf = argv[1]
  out_vcf = argv[2]

  with vcf.VcfReader(in_vcf, use_index=False) as reader:
    if 'AD' in [info.id for info in reader.header.infos]:
      print('%s already contains AD field.' % in_vcf)
      sys.exit(-1)
    out_header = reader.header
    out_header.infos.extend([vcf_constants.reserved_info_field('AD')])

    with vcf.VcfWriter(out_vcf, header=out_header) as writer:
      for variant in reader:
        variant_utils.set_info(variant, 'AD', get_variant_ad(variant),
                               writer)
        writer.write(variant)


if __name__ == '__main__':
  app.run(main)
