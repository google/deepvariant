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

import logging
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

# 1 is the tumor, 0 is the normal match.
# Tumor sample is the "main" sample because the goal here is somatic calling.
NORMAL_SAMPLE_INDEX = (
    0  # If the normal sample is present, it will appear first.
)

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
    dv_constants.PILEUP_DEFAULT_HEIGHT,
    (
        'Height for the part of the pileup image showing reads from the tumor. '
        f'Uses {dv_constants.PILEUP_DEFAULT_HEIGHT} by default.'
    ),
)
_PILEUP_IMAGE_HEIGHT_NORMAL = flags.DEFINE_integer(
    'pileup_image_height_normal',
    dv_constants.PILEUP_DEFAULT_HEIGHT,
    (
        'Height for the part of the pileup image showing reads from the matched'
        f' normal. Uses {dv_constants.PILEUP_DEFAULT_HEIGHT} by default.'
    ),
)
_CANDIDATE_POSITIONS = flags.DEFINE_string(
    'candidate_positions',
    '',
    'Path to the binary file containing candidate positions.',
)

# Change any flag defaults that differ for DeepSomatic.
# I'm setting this to float('inf') because we don't want to include any
# candidates from the non-target (i.e., normal) sample.
FLAGS.set_default('vsc_min_fraction_multiplier', float('inf'))


def tumor_normal_samples_from_flags(flags_obj):
  """Collects sample-related options into a list of samples."""
  samples_in_order = []

  def setup_sample(role):
    pileup_height = (
        flags_obj.pileup_image_height_tumor
        if role == 'tumor'
        else flags_obj.pileup_image_height_normal
    )
    read_filenames = (
        flags_obj.reads_tumor if role == 'tumor' else flags_obj.reads_normal
    )
    sample_name_flag = (
        flags_obj.sample_name_tumor
        if role == 'tumor'
        else flags_obj.sample_name_normal
    )
    reads_filenames_split = read_filenames.split(',')
    sample_name = make_examples_core.assign_sample_name(
        sample_name_flag=sample_name_flag,
        reads_filenames=read_filenames,
    )
    skip_output_generation = True if role == 'normal' else False

    sample_options = deepvariant_pb2.SampleOptions(
        role=role,
        name=sample_name,
        variant_caller_options=make_examples_core.make_vc_options(
            sample_name=sample_name, flags_obj=flags_obj
        ),
        skip_output_generation=skip_output_generation,
        pileup_height=pileup_height,
        reads_filenames=reads_filenames_split,
    )
    if role == 'tumor':
      if flags_obj.reads_normal and flags_obj.reads_tumor:
        sample_options.order.extend([0, 1])
      else:
        sample_options.order.extend([0])

    downsample_fraction = (
        flags_obj.downsample_fraction_tumor
        if role == 'tumor'
        else flags_obj.downsample_fraction_normal
    )
    if downsample_fraction != NO_DOWNSAMPLING:
      sample_options.downsample_fraction = downsample_fraction
    samples_in_order.append(sample_options)

  if flags_obj.reads_normal:
    setup_sample('normal')

  setup_sample('tumor')

  return samples_in_order, 'tumor'


def default_options(main_sample_index, add_flags=True, flags_obj=None):
  """Creates a MakeExamplesOptions proto populated with reasonable defaults.

  Args:
    main_sample_index: int. Indicates the position of the tumor sample.
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
      flags_obj=flags_obj
  )
  samples_in_order[main_sample_index].candidate_positions = (
      flags_obj.candidate_positions
  )

  options = make_examples_options.shared_flags_to_options(
      add_flags=add_flags,
      flags_obj=flags_obj,
      samples_in_order=samples_in_order,
      sample_role_to_train=sample_role_to_train,
      main_sample_index=main_sample_index,
  )

  if _READS_NORMAL.value:
    options.bam_fname = f'{os.path.basename(flags_obj.reads_tumor)}|{os.path.basename(flags_obj.reads_normal)}'
  else:
    options.bam_fname = os.path.basename(flags_obj.reads_tumor)

  return options


def check_options_are_valid(options, main_sample_index):
  """Checks that all the options chosen make sense together."""

  # Check for general flags (shared for DeepVariant and DeepTrio).
  make_examples_options.check_options_are_valid(
      options, main_sample_index=main_sample_index
  )

  tumor = options.sample_options[main_sample_index]
  if _READS_NORMAL.value:
    normal = options.sample_options[NORMAL_SAMPLE_INDEX]

    if (
        tumor.variant_caller_options.sample_name
        == normal.variant_caller_options.sample_name
    ):
      errors.log_and_raise(
          (
              'Sample names of tumor and normal samples cannot be the same. Use'
              ' --sample_name_tumor and --sample_name_normal with different'
              ' names '
          ),
          errors.CommandLineError,
      )

  if options.sample_options[main_sample_index].candidate_positions:
    if options.max_reads_per_partition:
      logging.warning(
          'Since candidate_positions is set, we use '
          'max_reads_for_dynamic_bases_per_region instead of '
          'max_reads_per_partition. This is due to the dynamic nature of the '
          'partition size when candidate_positions is enabled, making a fixed '
          'max_reads_per_partition unsuitable.'
      )
      options.max_reads_for_dynamic_bases_per_region = (
          options.max_reads_per_partition
      )
      options.max_reads_per_partition = 0


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
    is_tumor_only = _READS_TUMOR.value and not _READS_NORMAL.value
    main_sample_index = 0 if is_tumor_only else 1
    options = default_options(main_sample_index, flags_obj=FLAGS)
    check_options_are_valid(options, main_sample_index)

    # Run!
    make_examples_core.make_examples_runner(options)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'mode',
      'reads_tumor',
      'ref',
  ])
  app.run(main)
