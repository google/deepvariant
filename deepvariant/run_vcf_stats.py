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
r"""Creates a JSON summary of variants from a VCF file."""

from __future__ import absolute_import
from __future__ import division
from __future__ import google_type_annotations
from __future__ import print_function

import itertools

from absl import app
from absl import flags

from deepvariant import vcf_stats
from third_party.nucleus.io import vcf

FLAGS = flags.FLAGS

flags.DEFINE_string('input_vcf', None, 'Path to the input VCF.')
flags.DEFINE_string(
    'vcf_stats_outfile', None, 'Path to the output JSON file. '
    'The default value is input_vcf.json if this is not set.')
flags.DEFINE_integer(
    'num_records', -1, 'Maximum number of records to emit. If '
    'negative, emit all records.')


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  with vcf.VcfReader(FLAGS.input_vcf) as reader:
    if FLAGS.num_records == -1:
      variants = reader.iterate()
    else:
      variants = itertools.islice(reader.iterate(), FLAGS.num_records)
    stats = vcf_stats.variants_to_stats_json(variants)
    vcf_stats.write_json(stats, FLAGS.vcf_stats_outfile)


if __name__ == '__main__':
  flags.mark_flags_as_required(['input_vcf', 'vcf_stats_outfile'])
  app.run(main)
