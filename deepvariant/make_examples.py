# Copyright 2017 Google LLC.
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
"""Step one of DeepVariant: creates tf.Example protos for training/calling."""

import os
from absl import app
from absl import flags

from deepvariant import dv_constants

from deepvariant import logging_level
from deepvariant import make_examples_core
from deepvariant import make_examples_options
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io.python import hts_verbose

from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils

# Sentinel command line flag value indicating no downsampling should occur.
NO_DOWNSAMPLING = 0.0

MAIN_SAMPLE_INDEX = 0  # 0 is the only sample.

FLAGS = flags.FLAGS

# Adopt more general flags from make_examples_options.
flags.adopt_module_key_flags(make_examples_options)

# Flags related to the sample in DeepVariant:
flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref. '
    'Can provide multiple BAMs (comma-separated).')
flags.DEFINE_string(
    'sample_name', '', 'Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads.')
flags.DEFINE_integer(
    'pileup_image_height', 0,
    'Height for the pileup image. If 0, uses the default height')
flags.DEFINE_float(
    'downsample_fraction', NO_DOWNSAMPLING,
    f'If not {NO_DOWNSAMPLING} must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input BAM. This argument makes it easy to create examples as '
    'though the input BAM had less coverage.')
flags.DEFINE_string(
    'proposed_variants', '',
    '(Only used when --variant_caller=vcf_candidate_importer.) '
    'Tabix-indexed VCF file containing the proposed positions and alts for '
    '`vcf_candidate_importer`. The GTs will be ignored.')


def one_sample_from_flags(add_flags=True, flags_obj=None):
  """Collect sample-related options into a list of samples."""
  # Sample-specific options.
  sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=flags_obj.sample_name, reads_filenames=flags_obj.reads)
  sample_options = deepvariant_pb2.SampleOptions(
      role='main_sample',
      name=sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=sample_name, flags_obj=flags_obj),
      order=[0],
      pileup_height=dv_constants.PILEUP_DEFAULT_HEIGHT)

  if add_flags:
    if flags_obj.reads:
      sample_options.reads_filenames.extend(flags_obj.reads.split(','))
    if flags_obj.downsample_fraction != NO_DOWNSAMPLING:
      sample_options.downsample_fraction = flags_obj.downsample_fraction
    if flags_obj.pileup_image_height:
      sample_options.pileup_height = flags_obj.pileup_image_height
    if flags_obj.proposed_variants:
      sample_options.proposed_variants_filename = flags_obj.proposed_variants
  samples_in_order = [sample_options]
  sample_role_to_train = sample_options.role
  return samples_in_order, sample_role_to_train


def default_options(add_flags=True, flags_obj=None):
  """Creates a MakeExamplesOptions proto populated with reasonable defaults.

  Args:
    add_flags: bool. defaults to True. If True, we will push the value of
      certain FLAGS into our options. If False, those option fields are left
      uninitialized.
    flags_obj: object.  If not None, use as the source of flags, else use global
      FLAGS.

  Returns:
    deepvariant_pb2.MakeExamplesOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  if not flags_obj:
    flags_obj = FLAGS

  samples_in_order, sample_role_to_train = one_sample_from_flags(
      add_flags=add_flags, flags_obj=flags_obj)

  options = make_examples_options.shared_flags_to_options(
      add_flags=add_flags,
      flags_obj=flags_obj,
      samples_in_order=samples_in_order,
      sample_role_to_train=sample_role_to_train,
      main_sample_index=MAIN_SAMPLE_INDEX)

  if add_flags:
    options.bam_fname = os.path.basename(flags_obj.reads)
  return options


def check_options_are_valid(options):
  """Checks that all the options chosen make sense together."""

  # Check for general flags (shared for DeepVariant and DeepTrio).
  make_examples_options.check_options_are_valid(
      options, main_sample_index=MAIN_SAMPLE_INDEX)


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: make_examples does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)), errors.CommandLineError)
    del argv  # Unused.

    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()
    hts_verbose.set(hts_verbose.htsLogLevel[FLAGS.hts_logging_level])

    # Set up options; may do I/O.
    options = default_options(add_flags=True, flags_obj=FLAGS)
    check_options_are_valid(options)

    # Run!
    make_examples_core.make_examples_runner(options)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'mode',
      'reads',
      'ref',
  ])
  app.run(main)
