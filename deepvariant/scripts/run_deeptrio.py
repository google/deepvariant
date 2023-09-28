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
"""Runs all 3 steps to go from input DNA reads_child to output VCF/gVCF files.

This script currently provides the most common use cases and standard models.
If you want to access more flags that are available in `make_examples`,
`call_variants`, and `postprocess_variants`, you can also call them separately
using the binaries in the Docker image.
"""

import enum
import os
import subprocess
import sys
import tempfile
from typing import List, Optional

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

FLAGS = flags.FLAGS

# Required flags.
_MODEL_TYPE = flags.DEFINE_enum(
    'model_type',
    None,
    ['WGS', 'WES', 'PACBIO', 'HYBRID_PACBIO_ILLUMINA'],
    (
        'Required. Type of model to use for variant calling. Each '
        'model_type has an associated default model, which can be '
        'overridden by the --customized_model_{parent,child} flags.'
    ),
)
_REF = flags.DEFINE_string(
    'ref',
    None,
    (
        'Required. Genome reference to use. Must have an associated FAI index'
        ' as well. Supports text or gzipped references. Should match the'
        ' reference used to align the BAM file provided to --reads_child.'
    ),
)
_READS_CHILD = flags.DEFINE_string(
    'reads_child',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing the reads we'
        ' want to call. Should be aligned to a reference genome compatible with'
        ' --ref.'
    ),
)
_READS_PARENT1 = flags.DEFINE_string(
    'reads_parent1',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing parent 1 reads'
        ' of the person we want to call. Should be aligned to a reference'
        ' genome compatible with --ref.'
    ),
)
_READS_PARENT2 = flags.DEFINE_string(
    'reads_parent2',
    None,
    (
        'Aligned, sorted, indexed BAM file containing parent 2 reads of '
        'the person we want to call. Should be aligned to a reference genome '
        'compatible with --ref.'
    ),
)
_OUTPUT_VCF_CHILD = flags.DEFINE_string(
    'output_vcf_child',
    None,
    'Required. Path where we should write VCF file for the child.',
)
_OUTPUT_VCF_PARENT1 = flags.DEFINE_string(
    'output_vcf_parent1',
    None,
    'Required. Path where we should write VCF file for parent1.',
)
_OUTPUT_VCF_PARENT2 = flags.DEFINE_string(
    'output_vcf_parent2',
    None,
    'Required. Path where we should write VCF file for parent2.',
)
# Optional flags.
_DRY_RUN = flags.DEFINE_boolean(
    'dry_run',
    False,
    'Optional. If True, only prints out commands without executing them.',
)
_INTERMEDIATE_RESULTS_DIR = flags.DEFINE_string(
    'intermediate_results_dir',
    None,
    (
        'Optional. If specified, this should be an existing '
        'directory that is visible insider docker, and will be '
        'used to to store intermediate outputs.'
    ),
)
_VERSION = flags.DEFINE_boolean(
    'version',
    None,
    'Optional. If true, print out version number and exit.',
    allow_hide_cpp=True,
)
# TODO: Change to True as default before release.
_USE_SLIM_MODEL = flags.DEFINE_boolean(
    'use_slim_model',
    False,
    'Default to False. If True, the model provided has to be a Slim model.',
)

_LOGGING_DIR = flags.DEFINE_string(
    'logging_dir', None, 'Required. Directory where we should write log files.'
)
_RUNTIME_REPORT = flags.DEFINE_boolean(
    'runtime_report',
    False,
    (
        'Output make_examples runtime metrics '
        'and create a visual runtime report using runtime_by_region_vis. '
        'Only works with --logging_dir.'
    ),
)
# Optional flags for call_variants.
_CUSTOMIZED_MODEL_CHILD = flags.DEFINE_string(
    'customized_model_child',
    None,
    (
        'Optional. A path to a child model checkpoint to load for the'
        ' `call_variants` step. If not set, the default for each --model_type'
        ' will be used'
    ),
)
_CUSTOMIZED_MODEL_PARENT = flags.DEFINE_string(
    'customized_model_parent',
    None,
    (
        'Optional. A path to a parent model checkpoint to load for the'
        ' `call_variants` step. If not set, the default for each --model_type'
        ' will be used'
    ),
)
# Optional flags for make_examples.
_NUM_SHARDS = flags.DEFINE_integer(
    'num_shards', 1, 'Optional. Number of shards for make_examples step.'
)
_REGIONS = flags.DEFINE_string(
    'regions',
    None,
    (
        'Optional. Space-separated list of regions we want to process. Elements'
        ' can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE'
        ' files.'
    ),
)
_SAMPLE_NAME_CHILD = flags.DEFINE_string(
    'sample_name_child',
    None,
    (
        'Sample name to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_child.'
    ),
)
_SAMPLE_NAME_PARENT1 = flags.DEFINE_string(
    'sample_name_parent1',
    None,
    (
        'Parent1 Sample name to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_parent1.'
    ),
)
_SAMPLE_NAME_PARENT2 = flags.DEFINE_string(
    'sample_name_parent2',
    None,
    (
        'Parent2 Sample name to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads_parent2.'
    ),
)
_MAKE_EXAMPLES_EXTRA_ARGS = flags.DEFINE_string(
    'make_examples_extra_args',
    None,
    (
        'A comma-separated list of flag_name=flag_value. "flag_name" has to be'
        ' valid flags for make_examples.py. If the flag_value is boolean, it'
        ' has to be flag_name=true or flag_name=false.'
    ),
)
_CALL_VARIANTS_EXTRA_ARGS = flags.DEFINE_string(
    'call_variants_extra_args',
    None,
    (
        'A comma-separated list of flag_name=flag_value. "flag_name" has to be'
        ' valid flags for call_variants.py. If the flag_value is boolean, it'
        ' has to be flag_name=true or flag_name=false.'
    ),
)
_POSTPROCESS_VARIANTS_EXTRA_ARGS = flags.DEFINE_string(
    'postprocess_variants_extra_args',
    None,
    (
        'A comma-separated list of flag_name=flag_value. "flag_name" has to be'
        ' valid flags for postprocess_variants.py. If the flag_value is'
        ' boolean, it has to be flag_name=true or flag_name=false.'
    ),
)
_USE_CANDIDATE_PARTITION = flags.DEFINE_boolean(
    'use_candidate_partition',
    False,
    (
        '[This flag is experimental for internal testing. '
        'Do not set it to true.] '
        'Optional. If set, make_examples is run over partitions that contain an'
        ' equal number of candidates. Default value is False.'
    ),
)


# Optional flags for postprocess_variants.
_OUTPUT_GVCF_CHILD = flags.DEFINE_string(
    'output_gvcf_child',
    None,
    'Optional. Path where we should write gVCF file for child sample.',
)
_OUTPUT_GVCF_PARENT1 = flags.DEFINE_string(
    'output_gvcf_parent1',
    None,
    'Optional. Path where we should write gVCF file for parent1 sample.',
)
_OUTPUT_GVCF_PARENT2 = flags.DEFINE_string(
    'output_gvcf_parent2',
    None,
    'Optional. Path where we should write gVCF file for parent2 sample.',
)

# Optional flags for vcf_stats_report.
_VCF_STATS_REPORT = flags.DEFINE_boolean(
    'vcf_stats_report',
    True,
    (
        'Optional. Output a visual report (HTML) of '
        'statistics about each output VCF.'
    ),
)

MODEL_TYPE_MAP = {
    'WGS_child': '/opt/models/deeptrio/wgs/child/model.ckpt',
    'WGS_parent': '/opt/models/deeptrio/wgs/parent/model.ckpt',
    'WES_child': '/opt/models/deeptrio/wes/child/model.ckpt',
    'WES_parent': '/opt/models/deeptrio/wes/parent/model.ckpt',
    'PACBIO_child': '/opt/models/deeptrio/pacbio/child/model.ckpt',
    'PACBIO_parent': '/opt/models/deeptrio/pacbio/parent/model.ckpt',
}

# Current release version of DeepTrio.
# Should be the same in dv_vcf_constants.py.
DEEP_TRIO_VERSION = '1.4.0'
GLNEXUS_VERSION = 'v1.2.7'

DEEP_TRIO_WGS_PILEUP_HEIGHT_CHILD = 60
DEEP_TRIO_WGS_PILEUP_HEIGHT_PARENT = 40
DEEP_TRIO_WES_PILEUP_HEIGHT_CHILD = 100
DEEP_TRIO_WES_PILEUP_HEIGHT_PARENT = 100
DEEP_TRIO_PACBIO_PILEUP_HEIGHT_CHILD = 60
DEEP_TRIO_PACBIO_PILEUP_HEIGHT_PARENT = 40

CHILD = 'child'
PARENT1 = 'parent1'
PARENT2 = 'parent2'

CALL_VARIANTS_OUTPUT_COMMON_SUFFIX = 'tfrecord.gz'
NO_VARIANT_TFRECORD_SUFFIX = 'tfrecord.gz'
EXAMPLES_NAME_PATTERN = '{}_{}.{}'
CALL_VARIANTS_OUTPUT_PATTERN = '{}_{}.{}'
NO_VARIANT_TFRECORD_PATTERN = '{}_{}.{}'


@enum.unique
class CandidatePartitionCommand(enum.Enum):
  """make_examples mode for candidate partition."""

  SWEEP = enum.auto()  # Candidate sweep
  CANDIDATE_PARTITION_INFERENCE = (
      enum.auto()
  )  # Inference with candidate partition


def call_variants_output_common_prefix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'call_variants_output')


def examples_common_suffix(num_shards):
  return 'tfrecord@{}.gz'.format(num_shards)


def _candidate_positions_common_suffix(num_shards):
  return '@{}'.format(num_shards)


def examples_common_prefix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'make_examples')


def _candidate_positions_common_prefix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'candidate_positions')


def nonvariant_site_tfrecord_common_suffix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'gvcf')


def examples_common_name(intermediate_results_dir, num_shards):
  return '{}.{}'.format(
      examples_common_prefix(intermediate_results_dir),
      examples_common_suffix(num_shards),
  )


def _candidate_positions_common_name(intermediate_results_dir, num_shards):
  return '{}{}'.format(
      _candidate_positions_common_prefix(intermediate_results_dir),
      _candidate_positions_common_suffix(num_shards),
  )


def _is_quoted(value):
  if value.startswith('"') and value.endswith('"'):
    return True
  if value.startswith("'") and value.endswith("'"):
    return True
  return False


def _add_quotes(value):
  if isinstance(value, str) and _is_quoted(value):
    return value
  return '"{}"'.format(value)


def trim_suffix(string: str, suffix: str) -> str:
  if string.endswith(suffix):
    return string[: -len(suffix)]
  else:
    return string


def _extra_args_to_dict(extra_args):
  """Parses comma-separated list of flag_name=flag_value to dict."""
  args_dict = {}
  if extra_args is None:
    return args_dict
  for extra_arg in extra_args.split(','):
    (flag_name, flag_value) = extra_arg.split('=')
    flag_name = flag_name.strip('-')
    # Check for boolean values.
    if flag_value.lower() == 'true':
      flag_value = True
    elif flag_value.lower() == 'false':
      flag_value = False
    args_dict[flag_name] = flag_value
  return args_dict


def _extend_command_by_args_dict(command, extra_args):
  """Adds `extra_args` to the command string."""
  for key in sorted(extra_args):
    value = extra_args[key]
    if value is None:
      continue
    if isinstance(value, bool):
      added_arg = '' if value else 'no'
      added_arg += key
      command.extend(['--' + added_arg])
    else:
      command.extend(['--' + key, _add_quotes(value)])
  return command


def _update_kwargs_with_warning(kwargs, extra_args):
  """Updates `kwargs` with `extra_args`; gives a warning if values changed."""
  for k, v in extra_args.items():
    if k in kwargs:
      if kwargs[k] != v:
        print(
            '\nWarning: --{} is previously set to {}, now to {}.'.format(
                k, kwargs[k], v
            )
        )
    kwargs[k] = v
  return kwargs


def _make_examples_command(
    ref,
    reads_child,
    reads_parent1,
    reads_parent2,
    examples,
    sample_name_child,
    sample_name_parent1,
    sample_name_parent2,
    runtime_by_region_path,
    candidate_positions_path,
    extra_args,
    candidate_partition_mode=None,
    **kwargs,
):
  """Returns a make_examples command for subprocess.check_call.

  Args:
    ref: Input FASTA file.
    reads_child: Input BAM file for child.
    reads_parent1: Input BAM file for parent1.
    reads_parent2: Input BAM file for parent2.
    examples: Output tfrecord files suffix.
    sample_name_child: Sample name to use for child.
    sample_name_parent1: Sample name for parent1.
    sample_name_parent2: Sample name for parent2.
    runtime_by_region_path: Path for runtime statistics output.
    candidate_positions_path: Path to candidate positions file.
    extra_args: Comma-separated list of flag_name=flag_value.
    candidate_partition_mode: If set adds extra parameters to allow candidate
      partition.
    **kwargs: Additional arguments to pass in for make_examples.

  Returns:
    (string) A command to run.
  """
  command = [
      'time',
      'seq 0 {} |'.format(_NUM_SHARDS.value - 1),
      'parallel -q --halt 2 --line-buffer',
      '/opt/deepvariant/bin/deeptrio/make_examples',
  ]
  if candidate_partition_mode == CandidatePartitionCommand.SWEEP:
    command.extend(['--mode', 'candidate_sweep'])
  elif (
      candidate_partition_mode
      == CandidatePartitionCommand.CANDIDATE_PARTITION_INFERENCE
  ):
    command.extend(['--mode', 'calling'])
  elif candidate_partition_mode is None:
    command.extend(['--mode', 'calling'])
  else:
    raise ValueError('Invalid value of candidate_partition_mode.')

  command.extend(['--ref', '"{}"'.format(ref)])
  if _READS_PARENT1.value is not None:
    command.extend(['--reads_parent1', '"{}"'.format(reads_parent1)])
  if _READS_PARENT2.value is not None:
    command.extend(['--reads_parent2', '"{}"'.format(reads_parent2)])
  command.extend(['--reads', '"{}"'.format(reads_child)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--sample_name', '"{}"'.format(sample_name_child)])
  if _SAMPLE_NAME_PARENT1.value is not None:
    command.extend(
        ['--sample_name_parent1', '"{}"'.format(sample_name_parent1)]
    )
  if _SAMPLE_NAME_PARENT2.value is not None:
    command.extend(
        ['--sample_name_parent2', '"{}"'.format(sample_name_parent2)]
    )
  special_args = {}
  special_args['pileup_image_height_child'] = DEEP_TRIO_WGS_PILEUP_HEIGHT_CHILD
  special_args['pileup_image_height_parent'] = (
      DEEP_TRIO_WGS_PILEUP_HEIGHT_PARENT
  )

  if runtime_by_region_path is not None:
    command.extend(
        ['--runtime_by_region', '"{}"'.format(runtime_by_region_path)]
    )

  if _MODEL_TYPE.value == 'PACBIO':
    special_args['pileup_image_width'] = 199
    special_args['realign_reads'] = False
    special_args['vsc_min_fraction_indels'] = 0.12
    special_args['alt_aligned_pileup'] = 'diff_channels'
    special_args['add_hp_channel'] = True
    special_args['min_mapping_quality'] = 1
    special_args['track_ref_reads'] = True
    special_args['pileup_image_height_child'] = (
        DEEP_TRIO_PACBIO_PILEUP_HEIGHT_CHILD
    )
    special_args['pileup_image_height_parent'] = (
        DEEP_TRIO_PACBIO_PILEUP_HEIGHT_PARENT
    )
    # TODO make phasing optional.
    special_args['phase_reads'] = True
    if candidate_partition_mode != CandidatePartitionCommand.SWEEP:
      special_args['partition_size'] = 25000
    special_args['sort_by_haplotypes'] = True
    special_args['parse_sam_aux_fields'] = True
    special_args['discard_non_dna_regions'] = True
    special_args['max_reads_for_dynamic_bases_per_region'] = 200
    kwargs = _update_kwargs_with_warning(kwargs, special_args)

  if _MODEL_TYPE.value == 'WES':
    special_args['pileup_image_height_child'] = (
        DEEP_TRIO_WES_PILEUP_HEIGHT_CHILD
    )
    special_args['pileup_image_height_parent'] = (
        DEEP_TRIO_WES_PILEUP_HEIGHT_PARENT
    )
    special_args['channels'] = 'insert_size'

  if _MODEL_TYPE.value == 'WGS':
    special_args['channels'] = 'insert_size'

  if candidate_partition_mode == CandidatePartitionCommand.SWEEP:
    special_args['partition_size'] = 10000  # Should be approximately read
    # length to avoid having high
    # coverage intervals in multiple shards at a time
    special_args['candidate_positions'] = candidate_positions_path

  if (
      candidate_partition_mode
      == CandidatePartitionCommand.CANDIDATE_PARTITION_INFERENCE
  ):
    special_args['candidate_positions'] = candidate_positions_path

  if special_args:
    kwargs = _update_kwargs_with_warning(kwargs, special_args)

  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])

  if _LOGGING_DIR.value:
    log_filename = 'make_examples.log'
    if candidate_partition_mode == CandidatePartitionCommand.SWEEP:
      log_filename = 'make_examples_sweep.log'
    command.extend(
        ['2>&1 | tee {}/{}'.format(_LOGGING_DIR.value, log_filename)]
    )

  return ' '.join(command)


def call_variants_command(
    outfile: str,
    examples: str,
    model_ckpt: str,
    sample: str,
    extra_args: str,
    use_slim_model: bool = False,
) -> str:
  """Returns a call_variants command for subprocess.check_call."""
  binary_name = 'call_variants'
  if use_slim_model:
    binary_name = 'call_variants_slim'
  command = ['time', f'/opt/deepvariant/bin/{binary_name}']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(
      command, _extra_args_to_dict(extra_args)
  )
  if _LOGGING_DIR.value:
    command.extend(
        [
            '2>&1 | tee {}/{}_{}.log'.format(
                _LOGGING_DIR.value, binary_name, sample
            )
        ]
    )

  return ' '.join(command)


def postprocess_variants_command(
    ref,
    infile,
    outfile,
    sample,
    extra_args,
    nonvariant_site_tfrecord_path=None,
    gvcf_outfile=None,
):
  """Returns a postprocess_variants command for subprocess.check_call."""
  command = ['time', '/opt/deepvariant/bin/postprocess_variants']
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--infile', '"{}"'.format(infile)])
  command.extend(['--outfile', '"{}"'.format(outfile)])
  if nonvariant_site_tfrecord_path is not None:
    command.extend([
        '--nonvariant_site_tfrecord_path',
        '"{}"'.format(nonvariant_site_tfrecord_path),
    ])
  if gvcf_outfile is not None:
    command.extend(['--gvcf_outfile', '"{}"'.format(gvcf_outfile)])
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(
      command, _extra_args_to_dict(extra_args)
  )
  if _LOGGING_DIR.value:
    command.extend(
        [
            '2>&1 | tee {}/postprocess_variants_{}.log'.format(
                _LOGGING_DIR.value, sample
            )
        ]
    )

  return ' '.join(command)


def vcf_stats_report_command(vcf_path: str) -> str:
  """Returns a vcf_stats_report command for subprocess.

  Args:
    vcf_path: Path to VCF, which will be passed to --input_vcf and
      suffix-trimmed for --outfile_base.

  Returns:
    command string for subprocess
  """
  command = ['time', '/opt/deepvariant/bin/vcf_stats_report']
  command.extend(['--input_vcf', '"{}"'.format(vcf_path)])
  outfile_base = trim_suffix(trim_suffix(vcf_path, '.gz'), '.vcf')
  command.extend(['--outfile_base', '"{}"'.format(outfile_base)])

  return ' '.join(command)


def runtime_by_region_vis_command(runtime_by_region_path: str) -> str:
  """Returns a runtime_by_region_vis command for subprocess."""
  runtime_report = os.path.join(
      _LOGGING_DIR.value, 'make_examples_runtime_by_region_report.html'
  )

  command = ['time', '/opt/deepvariant/bin/runtime_by_region_vis']
  command.extend(['--input', '"{}"'.format(runtime_by_region_path)])
  command.extend(['--title', '"{}"'.format('DeepTrio')])
  command.extend(['--output', '"{}"'.format(runtime_report)])

  return ' '.join(command)


def check_or_create_intermediate_results_dir(
    intermediate_results_dir: Optional[str],
) -> str:
  """Checks or creates the path to the directory for intermediate results."""
  if intermediate_results_dir is None:
    intermediate_results_dir = tempfile.mkdtemp()
  if not os.path.isdir(intermediate_results_dir):
    logging.info(
        'Creating a directory for intermediate results in %s',
        intermediate_results_dir,
    )
    os.makedirs(intermediate_results_dir)
  else:
    logging.info(
        'Re-using the directory for intermediate results in %s',
        intermediate_results_dir,
    )
  return intermediate_results_dir


def model_exists(model_prefix: str) -> bool:
  if not tf.io.gfile.exists(
      model_prefix + '.data-00000-of-00001'
  ) or not tf.io.gfile.exists(model_prefix + '.index'):
    return False
  # If it's a Slim model, we also expect a .meta file.
  if not _USE_SLIM_MODEL.value and not tf.io.gfile.exists(
      model_prefix + '.meta'
  ):
    return False
  return True


def check_flags():
  """Additional logic to make sure flags are set appropriately."""
  if _CUSTOMIZED_MODEL_PARENT.value is not None:
    if not model_exists(_CUSTOMIZED_MODEL_PARENT.value):
      raise RuntimeError(
          'The model files {}* do not exist. Potentially '
          'relevant issue: '
          'https://github.com/google/deepvariant/blob/r1.5/docs/'
          'FAQ.md#why-cant-it-find-one-of-the-input-files-eg-'
          'could-not-open'.format(_CUSTOMIZED_MODEL_PARENT.value)
      )
    logging.info(
        (
            'You set --customized_model_parent. Instead of using the default '
            'model for %s, `call_variants` step will load %s* '
            'instead.'
        ),
        _MODEL_TYPE.value,
        _CUSTOMIZED_MODEL_PARENT.value,
    )

  if _CUSTOMIZED_MODEL_CHILD.value is not None:
    if not model_exists(_CUSTOMIZED_MODEL_CHILD.value):
      raise RuntimeError(
          'The model files {}* do not exist. Potentially '
          'relevant issue: '
          'https://github.com/google/deepvariant/blob/r1.5/docs/'
          'FAQ.md#why-cant-it-find-one-of-the-input-files-eg-'
          'could-not-open'.format(_CUSTOMIZED_MODEL_CHILD.value)
      )
    logging.info(
        (
            'You set --customized_model_child. Instead of using the default '
            'model for %s, `call_variants` step will load %s* '
            'instead.'
        ),
        _MODEL_TYPE.value,
        _CUSTOMIZED_MODEL_CHILD.value,
    )


def get_model_ckpt(model_type, customized_model):
  """Return the path to the model checkpoint based on the input args."""
  if customized_model is not None:
    return customized_model
  else:
    return MODEL_TYPE_MAP[model_type]


def generate_call_variants_command(
    sample, model_ckpt, intermediate_results_dir
):
  """Helper function to generate call_variants command line."""
  return call_variants_command(
      CALL_VARIANTS_OUTPUT_PATTERN.format(
          call_variants_output_common_prefix(intermediate_results_dir),
          sample,
          CALL_VARIANTS_OUTPUT_COMMON_SUFFIX,
      ),
      EXAMPLES_NAME_PATTERN.format(
          examples_common_prefix(intermediate_results_dir),
          sample,
          examples_common_suffix(_NUM_SHARDS.value),
      ),
      model_ckpt,
      sample,
      _CALL_VARIANTS_EXTRA_ARGS.value,
      use_slim_model=_USE_SLIM_MODEL.value,
  )


def generate_postprocess_variants_command(
    sample, intermediate_results_dir, output_vcf, output_gvcf
):
  """Helper function to generate post_process command line."""
  return postprocess_variants_command(
      ref=_REF.value,
      infile=CALL_VARIANTS_OUTPUT_PATTERN.format(
          call_variants_output_common_prefix(intermediate_results_dir),
          sample,
          CALL_VARIANTS_OUTPUT_COMMON_SUFFIX,
      ),
      outfile=output_vcf,
      sample=sample,
      extra_args=_POSTPROCESS_VARIANTS_EXTRA_ARGS.value,
      nonvariant_site_tfrecord_path=NO_VARIANT_TFRECORD_PATTERN.format(
          nonvariant_site_tfrecord_common_suffix(intermediate_results_dir),
          sample,
          examples_common_suffix(_NUM_SHARDS.value),
      ),
      gvcf_outfile=output_gvcf,
  )


def create_all_commands(intermediate_results_dir):
  """Creates commands for all stages to be executed later."""
  check_flags()
  commands = []
  post_process_commands = []
  report_commands = []

  # make_examples
  nonvariant_site_tfrecord_path = None
  if _OUTPUT_GVCF_CHILD.value is not None:
    nonvariant_site_tfrecord_path = '{}.{}'.format(
        nonvariant_site_tfrecord_common_suffix(intermediate_results_dir),
        examples_common_suffix(_NUM_SHARDS.value),
    )

  if _LOGGING_DIR.value and _RUNTIME_REPORT.value:
    runtime_directory = os.path.join(
        _LOGGING_DIR.value, 'make_examples_runtime_by_region'
    )
    if not os.path.isdir(runtime_directory):
      logging.info(
          'Creating a make_examples runtime by region directory in %s',
          runtime_directory,
      )
      os.makedirs(runtime_directory)
    # The path to runtime metrics output is sharded just like the examples.
    runtime_by_region_path = os.path.join(
        runtime_directory,
        'make_examples_runtime@{}.tsv'.format(_NUM_SHARDS.value),
    )
  else:
    runtime_by_region_path = None

  # If _USE_CANDIDATE_PARTITION is set add a call to make_examples to generate
  # candidates.
  # If _USE_CANDIDATE_PARTITION is set we generate two make_examples commands.
  # The first one to generate candidate_positions. The second command is for
  # generating DeepVariant examples. _USE_CANDIDATE_PARTITION is an option that
  # helps to better distribute the work between shards.
  candidate_partition_modes = []
  if _USE_CANDIDATE_PARTITION.value:
    candidate_partition_modes = [
        CandidatePartitionCommand.SWEEP,
        CandidatePartitionCommand.CANDIDATE_PARTITION_INFERENCE,
    ]
  else:
    candidate_partition_modes = [None]

  for candidate_partition_mode in candidate_partition_modes:
    commands.append(
        _make_examples_command(
            _REF.value,
            _READS_CHILD.value,
            _READS_PARENT1.value,
            _READS_PARENT2.value,
            examples_common_name(intermediate_results_dir, _NUM_SHARDS.value),
            _SAMPLE_NAME_CHILD.value,
            _SAMPLE_NAME_PARENT1.value,
            _SAMPLE_NAME_PARENT2.value,
            runtime_by_region_path=runtime_by_region_path,
            extra_args=_MAKE_EXAMPLES_EXTRA_ARGS.value,
            candidate_positions_path=_candidate_positions_common_name(
                intermediate_results_dir, _NUM_SHARDS.value
            ),
            candidate_partition_mode=candidate_partition_mode
            if _USE_CANDIDATE_PARTITION.value
            else None,
            # kwargs:
            gvcf=nonvariant_site_tfrecord_path,
            regions=_REGIONS.value,
        )
    )

  # Calling variants for child sample
  model_ckpt = get_model_ckpt(
      _MODEL_TYPE.value + '_child', _CUSTOMIZED_MODEL_CHILD.value
  )
  commands.append(
      generate_call_variants_command(
          CHILD, model_ckpt, intermediate_results_dir
      )
  )

  # Calling variants for parent1 sample
  model_ckpt = get_model_ckpt(
      _MODEL_TYPE.value + '_parent', _CUSTOMIZED_MODEL_PARENT.value
  )
  if _READS_PARENT1.value is not None:
    commands.append(
        generate_call_variants_command(
            PARENT1, model_ckpt, intermediate_results_dir
        )
    )
  if _READS_PARENT2.value is not None:
    commands.append(
        generate_call_variants_command(
            PARENT2, model_ckpt, intermediate_results_dir
        )
    )

  # postprocess_variants for child
  post_process_commands.append(
      generate_postprocess_variants_command(
          CHILD,
          intermediate_results_dir,
          _OUTPUT_VCF_CHILD.value,
          _OUTPUT_GVCF_CHILD.value,
      )
  )
  if _VCF_STATS_REPORT.value:
    report_commands.append(
        vcf_stats_report_command(vcf_path=_OUTPUT_VCF_CHILD.value)
    )

  if _READS_PARENT1.value is not None:
    post_process_commands.append(
        generate_postprocess_variants_command(
            PARENT1,
            intermediate_results_dir,
            _OUTPUT_VCF_PARENT1.value,
            _OUTPUT_GVCF_PARENT1.value,
        )
    )
    if _VCF_STATS_REPORT.value:
      report_commands.append(
          vcf_stats_report_command(vcf_path=_OUTPUT_VCF_PARENT1.value)
      )

  if _READS_PARENT2.value is not None:
    post_process_commands.append(
        generate_postprocess_variants_command(
            PARENT2,
            intermediate_results_dir,
            _OUTPUT_VCF_PARENT2.value,
            _OUTPUT_GVCF_PARENT2.value,
        )
    )
    if _VCF_STATS_REPORT.value:
      report_commands.append(
          vcf_stats_report_command(vcf_path=_OUTPUT_VCF_PARENT2.value)
      )

  # runtime-by-region
  if _LOGGING_DIR.value and _RUNTIME_REPORT.value:
    report_commands.append(
        runtime_by_region_vis_command(runtime_by_region_path)
    )

  return commands, post_process_commands, report_commands


def run_commands(
    commands: List[str], sequential: bool = True, dry_run: bool = False
) -> None:
  """Run commands using subprocess, sequentially or in parallel.

  Whether sequential or not, this function will finish when all commands have
  stopped running.

  Args:
    commands: List of string commands to run, to be executed by /bin/bash.
    sequential: True to execute commands one at a time, exiting if any fail.
      False to run all in parallel, exiting when all have succeeded or failed.
    dry_run: False to execute commands, True to only print commands.

  Raises:
    Exception: When one or more of the given commands have failed.
  """

  if sequential:
    for command in commands:
      print('\n***** Running the command:*****\n{}\n'.format(command))
      if not dry_run:
        try:
          subprocess.check_call(command, shell=True, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
          logging.info(e.output)
          raise Exception('Command(s) failed. See details above.') from e
  else:
    for command in commands:
      print('\n***** Running the command:*****\n{}\n'.format(command))
    if not _DRY_RUN.value:
      tasks = [
          subprocess.Popen(command, shell=True, executable='/bin/bash')
          for command in commands
      ]
      failed_task_indices = [
          i for i, task in enumerate(tasks) if task.wait() != 0
      ]
      if failed_task_indices:
        for i in failed_task_indices:
          logging.info('Failed command: %s', commands[i])
        raise Exception('Command(s) failed. See details above.')


def main(_):
  if _VERSION.value:
    print('DeepTrio version {}'.format(DEEP_TRIO_VERSION))
    return

  for flag_key in [
      'model_type',
      'ref',
      'reads_child',
      'output_vcf_child',
      'sample_name_child',
  ]:
    if FLAGS.get_flag_value(flag_key, None) is None:
      sys.stderr.write('--{} is required.\n'.format(flag_key))
      sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
      sys.exit(1)

  # Check flags consistency.
  # --reads_parent?, --output_vcf_parent?, --sample_name_parent? flags should
  # either all be set or all be unset.
  parent1_flags = [
      _READS_PARENT1.value,
      _OUTPUT_VCF_PARENT1.value,
      _SAMPLE_NAME_PARENT1.value,
  ]
  if any(parent1_flags) and not all(parent1_flags):
    sys.stderr.write(
        '--reads_parent1, --output_vcf_parent1, --sample_name_parent1 must be'
        ' set altogether\n'
    )
    sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
    sys.exit(1)

  parent2_flags = [
      _READS_PARENT2.value,
      _OUTPUT_VCF_PARENT2.value,
      _SAMPLE_NAME_PARENT2.value,
  ]
  if any(parent2_flags) and not all(parent2_flags):
    sys.stderr.write(
        '--reads_parent2, --output_vcf_parent2, --sample_name_parent2 must be'
        ' set altogether\n'
    )
    sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
    sys.exit(1)

  intermediate_results_dir = check_or_create_intermediate_results_dir(
      _INTERMEDIATE_RESULTS_DIR.value
  )

  commands, post_process_commands, report_commands = create_all_commands(
      intermediate_results_dir
  )
  print(
      '\n***** Intermediate results will be written to {} '
      'in docker. ****\n'.format(intermediate_results_dir)
  )
  # make_examples and call_variants:
  run_commands(commands=commands, sequential=True, dry_run=_DRY_RUN.value)
  run_commands(
      commands=post_process_commands, sequential=False, dry_run=_DRY_RUN.value
  )
  run_commands(
      commands=report_commands, sequential=False, dry_run=_DRY_RUN.value
  )


if __name__ == '__main__':
  app.run(main)
