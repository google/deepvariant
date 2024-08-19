# Copyright 2020 Google LLC.
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
"""Step one of DeepTrio: creates tf.Example protos for training/calling."""

import os

from absl import app
from absl import flags

from deeptrio import dt_constants
from deepvariant import logging_level
from deepvariant import make_examples_core
from deepvariant import make_examples_options
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils

# Sentinel command line flag value indicating no downsampling should occur.
NO_DOWNSAMPLING = 0.0

MAIN_SAMPLE_INDEX = 1  # 1 is the child of the trio.

FLAGS = flags.FLAGS

# Adopt more general flags from make_examples_options.
flags.adopt_module_key_flags(make_examples_options)

# Flags related to samples in DeepTrio:
SAMPLE_NAME_TO_TRAIN_ = flags.DEFINE_string(
    'sample_name_to_train',
    None,
    (
        'Optional - if not set, default to the value in --sample_name, i.e. the'
        ' child. The default is set to be backward compatible. If set, it has'
        ' to match one of --sample_name, --sample_name_parent1, or'
        ' --sample_name_parent2. Only used for training. When run in calling'
        ' mode, this is unused because examples are generated for all 3 samples'
        ' together.'
    ),
)
READS_ = flags.DEFINE_string(
    'reads',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing reads from the '
        'child of the trio. '
        'Should be aligned to a reference genome compatible with --ref. '
        'Can provide multiple BAMs (comma-separated).'
    ),
)
READS_PARENT1_ = flags.DEFINE_string(
    'reads_parent1',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing reads from'
        ' parent 1 of the trio. Should be aligned to a reference genome'
        ' compatible with --ref. Can provide multiple BAMs (comma-separated).'
    ),
)
READS_PARENT2_ = flags.DEFINE_string(
    'reads_parent2',
    None,
    (
        'Aligned, sorted, indexed BAM file containing reads from parent 2 of'
        ' the trio. Should be aligned to a reference genome compatible with'
        ' --ref. Can provide multiple BAMs (comma-separated).'
    ),
)
DOWNSAMPLE_FRACTION_CHILD_ = flags.DEFINE_float(
    'downsample_fraction_child',
    NO_DOWNSAMPLING,
    'If not '
    + str(NO_DOWNSAMPLING)
    + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input child BAM. This argument makes it easy to create examples '
    'as though the input BAM had less coverage.',
)
DOWNSAMPLE_FRACTION_PARENTS_ = flags.DEFINE_float(
    'downsample_fraction_parents',
    NO_DOWNSAMPLING,
    'If not '
    + str(NO_DOWNSAMPLING)
    + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input parent BAMs. This argument makes it easy to create examples'
    ' as though the input BAMs had less coverage.',
)
SAMPLE_NAME_ = flags.DEFINE_string(
    'sample_name',
    '',
    (
        'Child sample name to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads.'
    ),
)
SAMPLE_NAME_PARENT1_ = flags.DEFINE_string(
    'sample_name_parent1',
    '',
    (
        'Parent1 Sample name to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_parent1.'
    ),
)
SAMPLE_NAME_PARENT2_ = flags.DEFINE_string(
    'sample_name_parent2',
    '',
    (
        'Parent2 Sample name to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_parent2.'
    ),
)
PILEUP_IMAGE_HEIGHT_PARENT_ = flags.DEFINE_integer(
    'pileup_image_height_parent',
    0,
    'Height for the parent pileup image. If 0, uses the default height',
)
PILEUP_IMAGE_HEIGHT_CHILD_ = flags.DEFINE_integer(
    'pileup_image_height_child',
    0,
    'Height for the child pileup image. If 0, uses the default height',
)
PROPOSED_VARIANTS_CHILD_ = flags.DEFINE_string(
    'proposed_variants_child',
    None,
    (
        '(Only used when --variant_caller=vcf_candidate_importer.) '
        'Tabix-indexed VCF file containing the proposed positions and alts for '
        '`vcf_candidate_importer` for the child. The GTs will be ignored.'
    ),
)
PROPOSED_VARIANTS_PARENT1_ = flags.DEFINE_string(
    'proposed_variants_parent1',
    None,
    (
        '(Only used when --variant_caller=vcf_candidate_importer.) '
        'Tabix-indexed VCF file containing the proposed positions and alts for '
        '`vcf_candidate_importer` for the parent 1. The GTs will be ignored.'
    ),
)
PROPOSED_VARIANTS_PARENT2_ = flags.DEFINE_string(
    'proposed_variants_parent2',
    None,
    (
        '(Only used when --variant_caller=vcf_candidate_importer.) '
        'Tabix-indexed VCF file containing the proposed positions and alts for '
        '`vcf_candidate_importer` for the parent 2. The GTs will be ignored.'
    ),
)
# We are using this flag for determining intervals for both child and parent
# models. In the future, we can consider extending into 3 samples.
CANDIDATE_POSITIONS_ = flags.DEFINE_string(
    'candidate_positions',
    None,
    (
        'Path to the binary file containing candidate positions used for '
        'make_examples partitioning by candidates. Currently this '
        'is only the child positions.'
    ),
)
_SKIP_PARENT_CALLING = flags.DEFINE_bool(
    'skip_parent_calling',
    False,
    'If True, parents will not be called. Default is False.',
)

# Change any flag defaults that differ for DeepTrio.
FLAGS.set_default('vsc_min_fraction_multiplier', 0.67)


def trio_samples_from_flags(add_flags=True, flags_obj=None):
  """Collects sample-related options into a list of samples."""
  # Sample-specific options.
  child_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=SAMPLE_NAME_.value, reads_filenames=READS_.value
  )

  parent1_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=SAMPLE_NAME_PARENT1_.value,
      reads_filenames=READS_PARENT1_.value,
  )

  parent2_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=SAMPLE_NAME_PARENT2_.value,
      reads_filenames=READS_PARENT2_.value,
  )

  parent1_options = deepvariant_pb2.SampleOptions(
      role='parent1',
      name=parent1_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=parent1_sample_name, flags_obj=flags_obj
      ),
      order=[0, 1, 2],
      pileup_height=dt_constants.PILEUP_DEFAULT_HEIGHT_PARENT,
      skip_output_generation=_SKIP_PARENT_CALLING.value,
  )
  child_options = deepvariant_pb2.SampleOptions(
      role='child',
      name=child_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=child_sample_name, flags_obj=flags_obj
      ),
      order=[0, 1, 2],
      pileup_height=dt_constants.PILEUP_DEFAULT_HEIGHT_CHILD,
  )
  parent2_options = deepvariant_pb2.SampleOptions(
      role='parent2',
      name=parent2_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=parent2_sample_name, flags_obj=flags_obj
      ),
      # Swap the two parents when calling on parent2.
      order=[2, 1, 0],
      pileup_height=dt_constants.PILEUP_DEFAULT_HEIGHT_PARENT,
      skip_output_generation=_SKIP_PARENT_CALLING.value,
  )

  # If --sample_name_to_train is not set, train on the child.
  # This is for backward compatibility.
  sample_role_to_train = 'child'

  if add_flags:
    if READS_.value:
      child_options.reads_filenames.extend(READS_.value.split(','))
    if READS_PARENT1_.value:
      parent1_options.reads_filenames.extend(READS_PARENT1_.value.split(','))
    if READS_PARENT2_.value:
      parent2_options.reads_filenames.extend(READS_PARENT2_.value.split(','))

    if CANDIDATE_POSITIONS_.value:
      child_options.candidate_positions = CANDIDATE_POSITIONS_.value

    if PROPOSED_VARIANTS_CHILD_.value:
      child_options.proposed_variants_filename = PROPOSED_VARIANTS_CHILD_.value
    if PROPOSED_VARIANTS_PARENT1_.value:
      parent1_options.proposed_variants_filename = (
          PROPOSED_VARIANTS_PARENT1_.value
      )
    if PROPOSED_VARIANTS_PARENT2_.value:
      parent2_options.proposed_variants_filename = (
          PROPOSED_VARIANTS_PARENT2_.value
      )

    if DOWNSAMPLE_FRACTION_CHILD_.value != NO_DOWNSAMPLING:
      child_options.downsample_fraction = DOWNSAMPLE_FRACTION_CHILD_.value
    if DOWNSAMPLE_FRACTION_PARENTS_.value != NO_DOWNSAMPLING:
      parent1_options.downsample_fraction = DOWNSAMPLE_FRACTION_PARENTS_.value
      parent2_options.downsample_fraction = DOWNSAMPLE_FRACTION_PARENTS_.value

    if PILEUP_IMAGE_HEIGHT_CHILD_.value:
      child_options.pileup_height = PILEUP_IMAGE_HEIGHT_CHILD_.value
    if PILEUP_IMAGE_HEIGHT_PARENT_.value:
      parent1_options.pileup_height = (
          parent2_options.pileup_height
      ) = PILEUP_IMAGE_HEIGHT_PARENT_.value

    if SAMPLE_NAME_TO_TRAIN_.value:
      if SAMPLE_NAME_TO_TRAIN_.value == SAMPLE_NAME_.value:
        sample_role_to_train = child_options.role
      elif SAMPLE_NAME_TO_TRAIN_.value == SAMPLE_NAME_PARENT1_.value:
        sample_role_to_train = parent1_options.role
      else:
        errors.log_and_raise(
            (
                '--sample_name_to_train must match either --sample_name or '
                '--sample_name_parent1, or it can be unset to default to '
                '--sample_name.'
            ),
            errors.CommandLineError,
        )

  # Ordering here determines the default order of samples, and when a sample
  # above has a custom .order, then this is the list those indices refer to.
  samples_in_order = [parent1_options, child_options, parent2_options]
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

  samples_in_order, sample_role_to_train = trio_samples_from_flags(
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
    options.bam_fname = (
        os.path.basename(READS_.value)
        + '|'
        + (
            os.path.basename(READS_PARENT1_.value)
            if READS_PARENT1_.value
            else 'None'
        )
        + '|'
        + (
            os.path.basename(READS_PARENT2_.value)
            if READS_PARENT2_.value
            else 'None'
        )
    )
    options.pic_options.sequencing_type = (
        deepvariant_pb2.PileupImageOptions.TRIO
    )
    if not options.pic_options.height:
      options.pic_options.height = dt_constants.PILEUP_DEFAULT_HEIGHT
    if not options.pic_options.width:
      options.pic_options.width = dt_constants.PILEUP_DEFAULT_WIDTH

  return options


def check_options_are_valid(options):
  """Checks that all the options chosen make sense together."""

  # Check for general flags (shared for DeepVariant and DeepTrio).
  make_examples_options.check_options_are_valid(
      options, main_sample_index=MAIN_SAMPLE_INDEX
  )

  child = options.sample_options[MAIN_SAMPLE_INDEX]

  # Sanity check the sample_names (specific to trio).
  if (
      child.variant_caller_options.sample_name == FLAGS.sample_name_parent1
      or child.variant_caller_options.sample_name == FLAGS.sample_name_parent2
  ):
    errors.log_and_raise(
        'The sample_name of the child is the same as one of the parents.',
        errors.CommandLineError,
    )

  if options.pic_options.alt_aligned_pileup == 'rows':
    errors.log_and_raise(
        '--alt_aligned_pileup="rows" cannot be used with '
        'DeepTrio because the pileup images would become '
        'too tall for InceptionV3.'
    )

  if (
      options.mode == deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP
      and child.candidate_positions is None
  ):
    errors.log_and_raise(
        '--candidate_positions is required when --positions_sweep_mode is set.'
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
    hts_verbose.set(hts_verbose.htsLogLevel.HTS_LOG_WARNING)

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
