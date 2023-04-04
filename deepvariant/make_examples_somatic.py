# Copyright 2021 Google LLC.
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
"""A prototype to create tumor-normal images (tf.Example protos)."""

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

MAIN_SAMPLE_INDEX = 1  # 1 is the tumor, 0 is the normal match.

FLAGS = flags.FLAGS

# Adopt more general flags from make_examples_options.
flags.adopt_module_key_flags(make_examples_options)

# Flags related to samples in DeepSomatic:
_READS_TUMOR = flags.DEFINE_string(
    'reads_tumor',
    None,
    (
        'Required. Reads from the tumor sample. '
        'Aligned, sorted, indexed BAM file. '
        'Should be aligned to a reference genome compatible with --ref. '
        'Can provide multiple BAMs (comma-separated).'
    ),
)
_READS_NORMAL = flags.DEFINE_string(
    'reads_normal',
    None,
    (
        'Required. Reads from the normal matched sample. '
        'Aligned, sorted, indexed BAM file. '
        'Should be aligned to a reference genome compatible with --ref. '
        'Can provide multiple BAMs (comma-separated).'
    ),
)
_SAMPLE_NAME_TUMOR = flags.DEFINE_string(
    'sample_name_tumor',
    '',
    (
        'Sample name for tumor to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_tumor.'
    ),
)
_SAMPLE_NAME_NORMAL = flags.DEFINE_string(
    'sample_name_normal',
    '',
    (
        'Sample name for normal match to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_normal.'
    ),
)
_DOWNSAMPLE_FRACTION_TUMOR = flags.DEFINE_float(
    'downsample_fraction_tumor',
    NO_DOWNSAMPLING,
    'If not '
    + str(NO_DOWNSAMPLING)
    + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input tumor sample BAM. This argument makes it easy to create '
    'examples as though the input BAM had less coverage.',
)
_DOWNSAMPLE_FRACTION_NORMAL = flags.DEFINE_float(
    'downsample_fraction_normal',
    NO_DOWNSAMPLING,
    'If not '
    + str(NO_DOWNSAMPLING)
    + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input normal matched BAMs. This argument makes it easy to create '
    'examples as though the input BAMs had less coverage.',
)
_PILEUP_IMAGE_HEIGHT_TUMOR = flags.DEFINE_integer(
    'pileup_image_height_tumor',
    0,
    (
        'Height for the part of the pileup image showing reads from the tumor. '
        'If 0, uses the default height'
    ),
)
_PILEUP_IMAGE_HEIGHT_NORMAL = flags.DEFINE_integer(
    'pileup_image_height_normal',
    0,
    (
        'Height for the part of the pileup image showing reads from the matched'
        ' normal. If 0, uses the default height'
    ),
)

# Change any flag defaults that differ for DeepSomatic.
# I'm setting this to float('inf') because we don't want to include any
# candidates from the non-target (i.e., normal) sample.
FLAGS.set_default('vsc_min_fraction_multiplier', float('inf'))


def tumor_normal_samples_from_flags(add_flags=True, flags_obj=None):
  """Collects sample-related options into a list of samples."""
  # Sample-specific options.
  tumor_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=flags_obj.sample_name_tumor,
      reads_filenames=flags_obj.reads_tumor,
  )

  normal_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=flags_obj.sample_name_normal,
      reads_filenames=flags_obj.reads_normal,
  )

  tumor_sample_options = deepvariant_pb2.SampleOptions(
      role='tumor',
      name=tumor_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=tumor_sample_name, flags_obj=flags_obj
      ),
      order=[0, 1],
      pileup_height=dv_constants.PILEUP_DEFAULT_HEIGHT,
  )
  normal_sample_options = deepvariant_pb2.SampleOptions(
      role='normal',
      name=normal_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=normal_sample_name, flags_obj=flags_obj
      ),
      order=[0, 1],
      pileup_height=dv_constants.PILEUP_DEFAULT_HEIGHT,
  )

  if add_flags:
    if flags_obj.reads_tumor:
      tumor_sample_options.reads_filenames.extend(
          flags_obj.reads_tumor.split(',')
      )
    if flags_obj.reads_normal:
      normal_sample_options.reads_filenames.extend(
          flags_obj.reads_normal.split(',')
      )

    if flags_obj.downsample_fraction_tumor != NO_DOWNSAMPLING:
      tumor_sample_options.downsample_fraction = (
          flags_obj.downsample_fraction_tumor
      )
    if flags_obj.downsample_fraction_normal != NO_DOWNSAMPLING:
      normal_sample_options.downsample_fraction = (
          flags_obj.downsample_fraction_normal
      )

    if flags_obj.pileup_image_height_tumor:
      tumor_sample_options.pileup_height = flags_obj.pileup_image_height_tumor
    if flags_obj.pileup_image_height_normal:
      normal_sample_options.pileup_height = flags_obj.pileup_image_height_normal

  # Ordering here determines the default order of samples, and when a sample
  # above has a custom .order, then this is the list those indices refer to.
  samples_in_order = [normal_sample_options, tumor_sample_options]
  sample_role_to_train = 'tumor'
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

  samples_in_order, sample_role_to_train = tumor_normal_samples_from_flags(
      add_flags=add_flags, flags_obj=flags_obj
  )

  options = make_examples_options.shared_flags_to_options(
      add_flags=add_flags,
      flags_obj=flags_obj,
      samples_in_order=samples_in_order,
      sample_role_to_train=sample_role_to_train,
      main_sample_index=MAIN_SAMPLE_INDEX,
  )

  if add_flags:
    options.bam_fname = f'{os.path.basename(flags_obj.reads_tumor)}|{os.path.basename(flags_obj.reads_normal)}'

  return options


def check_options_are_valid(options):
  """Checks that all the options chosen make sense together."""

  # Check for general flags (shared for DeepVariant and DeepTrio).
  make_examples_options.check_options_are_valid(
      options, main_sample_index=MAIN_SAMPLE_INDEX
  )

  normal = options.sample_options[MAIN_SAMPLE_INDEX]

  if normal.variant_caller_options.sample_name == _SAMPLE_NAME_TUMOR.value:
    errors.log_and_raise(
        (
            'Sample names of tumor and normal samples cannot be the same. Use '
            '--sample_name_tumor and --sample_name_normal with different names '
        ),
        errors.CommandLineError,
    )


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: make_examples does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)),
          errors.CommandLineError,
      )
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
      'reads_tumor',
      'reads_normal',
      'ref',
  ])
  app.run(main)
