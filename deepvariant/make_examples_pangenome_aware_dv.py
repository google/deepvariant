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
"""A prototype to create pangenome-aware deepvariant images (tf.Example protos)."""

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

# 1 is for "reads", 0 is for "pangenome"
# "reads" sample is the "main" sample because the goal here is calling variant
# based on read alignments.
MAIN_SAMPLE_INDEX = 1
PANGENOME_SAMPLE_INDEX = 0

FLAGS = flags.FLAGS

# Adopt more general flags from make_examples_options.
flags.adopt_module_key_flags(make_examples_options)

# Flags related to samples in Pangenome-aware DeepVariant:
_READS = flags.DEFINE_string(
    'reads',
    None,
    (
        'Required. Reads from the sample. '
        'Aligned, sorted, indexed BAM file. '
        'Should be aligned to a reference genome compatible with --ref. '
        'Can provide multiple BAMs (comma-separated).'
    ),
)
_PANGENOME = flags.DEFINE_string(
    'pangenome',
    None,
    (
        'Required. Pangenome haplotypes to aid the variant calling process.'
        'Aligned, sorted, indexed BAM file. It can also be a gbz file. If'
        'a gbz file is used make sure to set --sample_name_pangenome and '
        '--ref_name_pangenome'
        'Should be aligned to a reference genome compatible with --ref. '
        'Can provide multiple BAMs (comma-separated).'
    ),
)
_SAMPLE_NAME_READS = flags.DEFINE_string(
    'sample_name_reads',
    '',
    (
        'Sample name for reads to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads.'
    ),
)
_SAMPLE_NAME_PANGENOME = flags.DEFINE_string(
    'sample_name_pangenome',
    'pangenome',
    (
        'Sample name for pangenome panel to use for our sample_name in the'
        'output Variant/DeepVariantCall protos. If not specified, will be'
        'inferred from the header information from --pangenome'
        '(if the pangenome is saved as a bam file).'
    ),
)
_DOWNSAMPLE_FRACTION_READS = flags.DEFINE_float(
    'downsample_fraction_reads',
    NO_DOWNSAMPLING,
    'If not '
    + str(NO_DOWNSAMPLING)
    + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input BAM containing reads. This argument makes it easy to create'
    ' examples as though the input BAM had less coverage.',
)

_PILEUP_IMAGE_HEIGHT_READS = flags.DEFINE_integer(
    'pileup_image_height_reads',
    0,
    (
        'Height for the part of the pileup image showing reads from the sample.'
        'If 0, uses the default height'
    ),
)
_PILEUP_IMAGE_HEIGHT_PANGENOME = flags.DEFINE_integer(
    'pileup_image_height_pangenome',
    0,
    (
        'Height for the part of the pileup image showing haplotypes from the '
        ' pangenome panel. If 0, uses the default height'
    ),
)
_KEEP_ONLY_WINDOW_SPANNING_HAPLOTYPES = flags.DEFINE_bool(
    'keep_only_window_spanning_haplotypes',
    True,
    (
        '[optional] If True, it will trim all haplotypes overlapping the'
        'candidate variant and keep only the haplotypes that fully span'
        'the whole example window. This is an option for the pangenome-aware '
        'model.'
    ),
)
_VARIANT_TYPES_TO_BLANK = flags.DEFINE_list(
    'variant_types_to_blank',
    [],
    (
        '[optional] A comma-separated list of variant types to blank out in the'
        ' pangenome sample pileup image.'
        ' The possible values are INDEL, and SNP.'
    ),
)

_REF_NAME_PANGENOME = flags.DEFINE_string(
    'ref_name_pangenome',
    'GRCh38',
    (
        'The name of the reference genome in the pangenome gbz file.'
        'This reference should match the reference used for the reads. This'
        'flag is added since the exact name assigned to the pangenome reference'
        'can be different from the name of the reference fasta used for '
        'the reads.'
    ),
)

_REF_CHROM_PREFIX = flags.DEFINE_string(
    'ref_chrom_prefix',
    '',
    (
        'The prefix to add to the chromosome name in the pangenome gbz file. It'
        'is empty by default. However sometimes we need to add a prefix '
        '(like "GRCh38.") to the chromosome name in the pangenome gbz file to '
        'match the chromosome name in the reads.'
    ),
)

_USE_LOADED_GBZ_SHARED_MEMORY = flags.DEFINE_bool(
    'use_loaded_gbz_shared_memory',
    False,
    (
        'If enabled, the sequences of the gbz file are already loaded into'
        'shared memory using load_gbz_into_shared_memory.py and the SamReader'
        'reads the sequences from the shared memory.'
    ),
)

_GBZ_SHARED_MEMORY_NAME = flags.DEFINE_string(
    'gbz_shared_memory_name',
    'GBZ_SHARED_MEMORY',
    (
        'Name of the shared memory segment that contains the sequences of the'
        'gbz format. If --use_loaded_gbz_shared_memory is enabled, this flag '
        'must be set based on the name of the shared memory created by'
        'load_gbz_into_shared_memory.py'
    ),
)

_TRAINED_SMALL_MODEL_PATH = flags.DEFINE_string(
    'trained_small_model_path',
    '',
    'Path to a small model checkpoint directory.',
)

# Change any flag defaults that differ for Pangenome-aware DeepVariant.
# I'm setting this to float('inf') because we don't want to include any
# candidates from the non-target (i.e., pangenome) sample.
FLAGS.set_default('vsc_min_fraction_multiplier', float('inf'))
# trim_reads_for_pileup is always needed for the Pangenome input.
FLAGS.set_default('trim_reads_for_pileup', True)


def reads_and_pangenome_samples_from_flags(add_flags=True, flags_obj=None):
  """Collects sample-related options into a list of samples."""
  # Sample-specific options.
  reads_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=flags_obj.sample_name_reads,
      reads_filenames=flags_obj.reads,
  )

  pangenome_sample_name = make_examples_core.assign_sample_name(
      sample_name_flag=flags_obj.sample_name_pangenome,
      reads_filenames=flags_obj.pangenome,
  )

  reads_sample_options = deepvariant_pb2.SampleOptions(
      role='reads',
      name=reads_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=reads_sample_name, flags_obj=flags_obj
      ),
      order=[0, 1],
      pileup_height=dv_constants.PILEUP_DEFAULT_HEIGHT,
      small_model_path=_TRAINED_SMALL_MODEL_PATH.value,
  )
  def variant_type_string_to_enum(variant_type_string):
    variant_type_string = variant_type_string.upper()
    if variant_type_string == 'INDEL':
      return deepvariant_pb2.SampleOptions.VARIANT_TYPE_INDEL
    elif variant_type_string == 'SNP':
      return deepvariant_pb2.SampleOptions.VARIANT_TYPE_SNP
    else:
      raise ValueError('Invalid variant type: %s' % variant_type_string)

  pangenome_sample_options = deepvariant_pb2.SampleOptions(
      role='pangenome',
      name=pangenome_sample_name,
      variant_caller_options=make_examples_core.make_vc_options(
          sample_name=pangenome_sample_name, flags_obj=flags_obj
      ),
      skip_output_generation=True,
      pileup_height=dv_constants.PILEUP_DEFAULT_HEIGHT,
      keep_only_window_spanning_reads=flags_obj.keep_only_window_spanning_haplotypes,
      skip_phasing=True,
      skip_normalization=True,
      channels_enum_to_blank=[
          deepvariant_pb2.CH_HAPLOTYPE_TAG,
          deepvariant_pb2.CH_DIFF_CHANNELS_ALTERNATE_ALLELE_1,
          deepvariant_pb2.CH_DIFF_CHANNELS_ALTERNATE_ALLELE_2,
          deepvariant_pb2.CH_BASE_QUALITY,
          deepvariant_pb2.CH_MAPPING_QUALITY,
      ],
      variant_types_to_blank=[
          variant_type_string_to_enum(variant_type_string)
          for variant_type_string in _VARIANT_TYPES_TO_BLANK.value
      ],
  )

  if add_flags:
    if flags_obj.reads:
      reads_sample_options.reads_filenames.extend(flags_obj.reads.split(','))
    if flags_obj.pangenome:
      pangenome_sample_options.reads_filenames.extend(
          flags_obj.pangenome.split(',')
      )

    if flags_obj.downsample_fraction_reads != NO_DOWNSAMPLING:
      reads_sample_options.downsample_fraction = (
          flags_obj.downsample_fraction_reads
      )
    if flags_obj.pileup_image_height_reads:
      reads_sample_options.pileup_height = flags_obj.pileup_image_height_reads
    if flags_obj.pileup_image_height_pangenome:
      pangenome_sample_options.pileup_height = (
          flags_obj.pileup_image_height_pangenome
      )

  # Ordering here determines the default order of samples, and when a sample
  # above has a custom .order, then this is the list those indices refer to.
  samples_in_order = [pangenome_sample_options, reads_sample_options]
  sample_role_to_train = 'reads'
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

  samples_in_order, sample_role_to_train = (
      reads_and_pangenome_samples_from_flags(
          add_flags=add_flags, flags_obj=flags_obj
      )
  )

  options = make_examples_options.shared_flags_to_options(
      add_flags=add_flags,
      flags_obj=flags_obj,
      samples_in_order=samples_in_order,
      sample_role_to_train=sample_role_to_train,
      main_sample_index=MAIN_SAMPLE_INDEX,
  )
  if flags_obj.ref_name_pangenome:
    options.ref_name_pangenome = flags_obj.ref_name_pangenome
  if flags_obj.ref_chrom_prefix:
    options.ref_chrom_prefix = flags_obj.ref_chrom_prefix
  if flags_obj.use_loaded_gbz_shared_memory:
    options.use_loaded_gbz_shared_memory = (
        flags_obj.use_loaded_gbz_shared_memory
    )
  if flags_obj.gbz_shared_memory_name:
    options.gbz_shared_memory_name = flags_obj.gbz_shared_memory_name

  if add_flags:
    options.bam_fname = f'{os.path.basename(flags_obj.reads)}|{os.path.basename(flags_obj.pangenome)}'

  return options


def check_options_are_valid(options):
  """Checks that all the options chosen make sense together."""

  # Check for general flags (shared for DeepVariant and DeepTrio).
  make_examples_options.check_options_are_valid(
      options, main_sample_index=MAIN_SAMPLE_INDEX
  )

  reads = options.sample_options[MAIN_SAMPLE_INDEX]
  pangenome = options.sample_options[PANGENOME_SAMPLE_INDEX]

  if (
      reads.variant_caller_options.sample_name
      == pangenome.variant_caller_options.sample_name
  ):
    errors.log_and_raise(
        (
            'Sample names for reads and pangenome cannot be the same.'
            'Use --sample_name_reads and --sample_name_pangenome with '
            'different names '
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
      'pangenome',
      'ref',
  ])
  app.run(main)
