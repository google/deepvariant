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
r"""Creates a visual HTML report about the variants from a VCF file."""

import itertools
from absl import flags
import tensorflow as tf

from deepvariant import vcf_stats
from third_party.nucleus.io import vcf
from third_party.nucleus.util import errors

FLAGS = flags.FLAGS

_INPUT_VCF = flags.DEFINE_string('input_vcf', None, 'Path to the input VCF.')
_OUTFILE_BASE = flags.DEFINE_string(
    'outfile_base',
    None,
    (
        'Base path for output. The HTML report will be created at '
        '<outfile_base>.visual_report.html.'
    ),
)
_TITLE = flags.DEFINE_string(
    'title',
    None,
    (
        'Title to show at the top of the HTML report. '
        'Default: Use the sample name from the VCF.'
    ),
)
_NUM_RECORDS = flags.DEFINE_integer(
    'num_records',
    -1,
    'Maximum number of VCF lines to read. If negative, read the whole VCF.',
)


def main(argv):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: vcf_stats_report does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv[1:])),
          errors.CommandLineError,
      )

  with vcf.VcfReader(_INPUT_VCF.value) as reader:
    sample_names = reader.header.sample_names
    if len(sample_names) != 1:
      raise ValueError(
          'There must be exactly one sample in VCF: {}'.format(_INPUT_VCF.value)
      )
    sample_name = sample_names[0]

    # Check for missing GT in VCF to avoid a confusing error downstream.
    vcf_columns = [col.id for col in reader.header.formats]
    if 'GT' not in vcf_columns:
      errors.log_and_raise('ERROR: No GT sub-column in VCF.')

    if _NUM_RECORDS.value < 0:
      variants = reader.iterate()
    else:
      variants = itertools.islice(reader.iterate(), FLAGS.num_records)

    vcf_stats.create_vcf_report(
        variants,
        output_basename=_OUTFILE_BASE.value,
        title=_TITLE.value or sample_name,
        vcf_reader=reader,
    )


if __name__ == '__main__':
  flags.mark_flags_as_required(['input_vcf', 'outfile_base'])
  tf.compat.v1.app.run()
