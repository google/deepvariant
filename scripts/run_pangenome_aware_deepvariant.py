# Copyright 2023 Google LLC.
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
"""Runs all 3 steps to go from input DNA reads to output VCF/gVCF files.

This script is used to run pangenome-aware DeepVariant, which is a variation of
DeepVariant.
If you want to access more flags that are available in
`make_examples_pangenome_aware_dv`, `call_variants`, and `postprocess_variants`,
you can also call them separately using the binaries in the Docker image.
"""

import dataclasses
import enum
import os
import subprocess
import sys
import tempfile
from typing import Any, Optional

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

FLAGS = flags.FLAGS


class ModelType(enum.Enum):
  WGS = 'WGS'
  WES = 'WES'


# Required flags.
_MODEL_TYPE = flags.DEFINE_enum(
    'model_type',
    None,
    [m.value for m in ModelType],
    (
        'Required. Type of model to use for variant calling. Set this flag to'
        ' use the default model associated with each type, and it will set'
        ' necessary flags corresponding to each model. If you want to use a'
        ' customized model, add --customized_model flag in addition to this'
        ' flag.'
    ),
)
_REF = flags.DEFINE_string(
    'ref',
    None,
    (
        'Required. Genome reference to use. Must have an associated FAI index'
        ' as well. Supports text or gzipped references. Should match the'
        ' reference used to align the BAM file provided to --reads.'
    ),
)
_READS = flags.DEFINE_string(
    'reads',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing the reads we'
        ' want to call. Should be aligned to a reference genome compatible with'
        ' --ref.'
    ),
)
_PANGENOME = flags.DEFINE_string(
    'pangenome',
    None,
    (
        'Required. Pangenome haplotypes to aid the variant calling process.'
        'Aligned, sorted, indexed BAM file. '
        'Should be aligned to a reference genome compatible with --ref. '
        'Can provide multiple BAMs (comma-separated).'
    ),
)
_OUTPUT_VCF = flags.DEFINE_string(
    'output_vcf', None, 'Required. Path where we should write VCF file.'
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
_LOGGING_DIR = flags.DEFINE_string(
    'logging_dir',
    None,
    (
        'Optional. Directory where we should write log files '
        'for each stage and optionally runtime reports.'
    ),
)
_RUNTIME_REPORT = flags.DEFINE_boolean(
    'runtime_report',
    False,
    (
        'Output make_examples_pangenome_aware_dv runtime metrics '
        'and create a visual runtime report using runtime_by_region_vis. '
        'Only works with --logging_dir.'
    ),
)
_VERSION = flags.DEFINE_boolean(
    'version',
    None,
    'Optional. If true, print out version number and exit.',
    allow_hide_cpp=True,
)

# Optional flags for call_variants.
_CUSTOMIZED_MODEL = flags.DEFINE_string(
    'customized_model',
    None,
    (
        'Optional. A path to a model checkpoint to load for the `call_variants`'
        ' step. If not set, the default for each --model_type will be used'
    ),
)
_DISABLE_SMALL_MODEL = flags.DEFINE_boolean(
    'disable_small_model',
    True,
    'Optional. Disable the use of the small model to call variants during '
    'the `make_examples` step.',
)
_CUSTOMIZED_SMALL_MODEL = flags.DEFINE_string(
    'customized_small_model',
    None,
    (
        'Optional. A path to a small model checkpoint to call variants during '
        'the `make_examples` step.'
    ),
)
# Optional flags for make_examples_pangenome_aware_dv.
_NUM_SHARDS = flags.DEFINE_integer(
    'num_shards',
    1,
    'Optional. Number of shards for make_examples_pangenome_aware_dv step.',
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
_SAMPLE_NAME_READS = flags.DEFINE_string(
    'sample_name_reads',
    None,
    (
        'Sample name for reads to use for our sample_name in the output'
        ' Variant/DeepVariantCall protos. If not specified, will be inferred'
        ' from the header information from --reads.'
    ),
)
_SAMPLE_NAME_PANGENOME = flags.DEFINE_string(
    'sample_name_pangenome',
    'hprc_v1.1',
    (
        'Sample name for pangenome panel to use for our sample_name in the'
        'output Variant/DeepVariantCall protos. If not specified, will be'
        'inferred from the header information from --pangenome.'
        'The default here is corresponding to the default for --pangenome.'
    ),
)
_MAKE_EXAMPLES_EXTRA_ARGS = flags.DEFINE_string(
    'make_examples_extra_args',
    None,
    (
        'A comma-separated list of flag_name=flag_value. "flag_name" has to be'
        ' valid flags for make_examples_pangenome_aware_dv.py. If the '
        'flag_value is boolean, it has to be flag_name=true or flag_name=false.'
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
# Optional flag for postprocess variants
_POSTPROCESS_CPUS = flags.DEFINE_integer(
    'postprocess_cpus',
    None,
    'Optional. Number of cpus to use during'
    ' postprocess_variants. Set to 0 to disable multiprocessing. Default is'
    ' None which sets to num_shards.',
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

# Optional flags for postprocess_variants.
_OUTPUT_GVCF = flags.DEFINE_string(
    'output_gvcf', None, 'Optional. Path where we should write gVCF file.'
)

# Optional flags for vcf_stats_report.
_VCF_STATS_REPORT = flags.DEFINE_boolean(
    'vcf_stats_report',
    False,
    (
        'Optional. Output a visual report (HTML) of '
        'statistics about the output VCF.'
    ),
)
_REPORT_TITLE = flags.DEFINE_string(
    'report_title',
    None,
    (
        'Optional. Title for the VCF stats report (HTML).'
        'If not provided, the title will be the sample name.'
    ),
)

# Optional flags that are used for loading GBZ into shared memory.
_GBZ_SHARED_MEMORY_NAME = flags.DEFINE_string(
    'gbz_shared_memory_name',
    None,
    (
        'Optional. Name of the shared memory region to create. This flag will'
        ' be used by load_gbz_into_shared_memory binary and'
        ' make_examples_pangenome_aware_dv binary.'
    ),
)

_GBZ_SHARED_MEMORY_SIZE_GB = flags.DEFINE_integer(
    'gbz_shared_memory_size_gb',
    12,
    'Optional. Size of the shared memory region to create.',
)

MODEL_TYPE_MAP = {
    'WGS': '/opt/models/pangenome_aware_deepvariant/wgs',
    'WES': '/opt/models/pangenome_aware_deepvariant/wes',
}


@dataclasses.dataclass
class SmallModelConfig:
  small_model_checkpoint: str
  snp_gq_threshold: int
  indel_gq_threshold: int
  vaf_context_window: int


SMALL_MODEL_CONFIG_BY_MODEL_TYPE = {
    ModelType.WGS: SmallModelConfig(
        small_model_checkpoint='/opt/smallmodels/wgs',
        snp_gq_threshold=25,
        indel_gq_threshold=30,
        vaf_context_window=51,
    ),
}

# Current release version of DeepVariant.
# Should be the same in dv_vcf_constants.py.
DEEP_VARIANT_VERSION = '1.8.0'


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


def _extra_args_to_dict(extra_args: str) -> dict[str, Any]:
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


def _use_small_model() -> bool:
  """Determines if the small model is enabled based on flags and model type."""
  if _DISABLE_SMALL_MODEL.value:
    return False
  if _CUSTOMIZED_SMALL_MODEL.value:
    return True
  return ModelType(_MODEL_TYPE.value) in SMALL_MODEL_CONFIG_BY_MODEL_TYPE


def _set_small_model_config(
    special_args: dict[str, Any],
    model_type: ModelType,
    customized_small_model: str | None,
) -> None:
  """Sets small model config parameters."""
  if not _use_small_model():
    return
  special_args['call_small_model_examples'] = True
  config = SMALL_MODEL_CONFIG_BY_MODEL_TYPE.get(model_type)
  if customized_small_model:
    special_args['trained_small_model_path'] = customized_small_model
  elif config:
    special_args['trained_small_model_path'] = config.small_model_checkpoint
  if config:
    special_args['small_model_snp_gq_threshold'] = config.snp_gq_threshold
    special_args['small_model_indel_gq_threshold'] = config.indel_gq_threshold
    special_args['small_model_vaf_context_window_size'] = (
        config.vaf_context_window
    )


def load_gbz_into_shared_memory_command(
    gbz: str,
    gbz_shared_memory_name: str | None,
    gbz_shared_memory_size_gb: int,
) -> tuple[str, Optional[str]]:
  """Returns a load_gbz_into_shared_memory (command, logfile) for subprocess.

  Args:
    gbz: Input pangenome GBZ file(s).
    gbz_shared_memory_name: Name of the shared memory region to create.
    gbz_shared_memory_size_gb: Size of the shared memory region to create.

  Returns:
    (string, string) A command to run, and a log file to output to.
  """
  command = ['time', '/opt/deepvariant/bin/load_gbz_into_shared_memory']
  command.extend(['--pangenome_gbz', '"{}"'.format(gbz)])
  if gbz_shared_memory_name is not None:
    command.extend(
        ['--shared_memory_name', '"{}"'.format(gbz_shared_memory_name)]
    )
  command.extend(['--shared_memory_size_gb', str(gbz_shared_memory_size_gb)])
  command.extend(['--num_shards', '"{}"'.format(_NUM_SHARDS.value)])

  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/load_gbz_into_shared_memory.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def make_examples_pangenome_aware_dv_command(
    ref: str,
    reads: str,
    pangenome: str,
    examples: str,
    model_ckpt: str,
    extra_args: Optional[str],
    runtime_by_region_path: Optional[str] = None,
    **kwargs,
) -> tuple[str, Optional[str]]:
  """Returns a make_examples_pangenome_aware_dv (command, logfile) for subprocess.

  Args:
    ref: Input FASTA file.
    reads: Input BAM file.
    pangenome: Input pangenome BAM file(s).
    examples: Output tfrecord file containing tensorflow.Example files.
    model_ckpt: Path to the TensorFlow model checkpoint.
    extra_args: Comma-separated list of flag_name=flag_value.
    runtime_by_region_path: Output path for runtime by region metrics.
    **kwargs: Additional arguments to pass in for
      make_examples_pangenome_aware_dv.

  Returns:
    (string, string) A command to run, and a log file to output to.
  """
  command = [
      'time',
      'seq 0 {} |'.format(_NUM_SHARDS.value - 1),
      'parallel -q --halt 2 --line-buffer',
      '/opt/deepvariant/bin/make_examples_pangenome_aware_dv',
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads', '"{}"'.format(reads)])
  command.extend(['--pangenome', '"{}"'.format(pangenome)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])

  if runtime_by_region_path is not None:
    command.extend(
        ['--runtime_by_region', '"{}"'.format(runtime_by_region_path)]
    )

  special_args = {}
  model_type = ModelType(_MODEL_TYPE.value)
  if model_type == ModelType.WGS:
    # Specific flags that are not default can be added here.
    special_args['keep_only_window_spanning_haplotypes'] = True
    special_args['keep_supplementary_alignments'] = True
    special_args['sort_by_haplotypes'] = True
    special_args['min_mapping_quality'] = 0
    special_args['keep_legacy_allele_counter_behavior'] = True
    special_args['normalize_reads'] = True
    special_args['trim_reads_for_pileup'] = True
  elif model_type == ModelType.WES:
    special_args['keep_only_window_spanning_haplotypes'] = True
    special_args['keep_supplementary_alignments'] = True
    special_args['sort_by_haplotypes'] = True
    special_args['trim_reads_for_pileup'] = True
  else:
    raise ValueError('Invalid model_type: %s' % _MODEL_TYPE.value)

  _set_small_model_config(
      special_args, model_type, _CUSTOMIZED_SMALL_MODEL.value
  )
  # By default, use pangenome GBZ that's already into shared memory
  # load_gbz_into_shared_memory should be run
  # before make_examples_pangenome_aware_dv to have the GBZ loaded into shared
  # memory.
  if pangenome is not None and pangenome.endswith('.gbz'):
    special_args['use_loaded_gbz_shared_memory'] = True
  kwargs = _update_kwargs_with_warning(kwargs, special_args)

  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/make_examples_pangenome_aware_dv.log'.format(
        _LOGGING_DIR.value
    )
  return (' '.join(command), logfile)


def call_variants_command(
    outfile: str,
    examples: str,
    model_ckpt: str,
    extra_args: str,
    allow_empty_examples: bool,
) -> tuple[str, Optional[str]]:
  """Returns a call_variants (command, logfile) for subprocess."""
  binary_name = 'call_variants'
  command = ['time', f'/opt/deepvariant/bin/{binary_name}']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])
  if allow_empty_examples:
    command.extend(['--allow_empty_examples'])
  if extra_args and 'use_openvino' in extra_args:
    raise RuntimeError(
        'OpenVINO is not installed by default in DeepVariant '
        'Docker images. Please rerun without use_openvino flag.'
    )
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(
      command, _extra_args_to_dict(extra_args)
  )
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/{}.log'.format(_LOGGING_DIR.value, binary_name)
  return (' '.join(command), logfile)


def postprocess_variants_command(
    ref: str,
    infile: str,
    outfile: str,
    small_model_cvo_records: str,
    extra_args: str,
    **kwargs,
) -> tuple[str, Optional[str]]:
  """Returns a postprocess_variants (command, logfile) for subprocess."""
  cpus = _POSTPROCESS_CPUS.value
  if cpus is None:
    cpus = _NUM_SHARDS.value
    # WES does not benefit from multiprocessing.
    if ModelType(_MODEL_TYPE.value) == ModelType.WES:
      cpus = 0
  command = ['time', '/opt/deepvariant/bin/postprocess_variants']
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--infile', '"{}"'.format(infile)])
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--cpus', '"{}"'.format(cpus)])
  if _use_small_model():
    command.extend(
        ['--small_model_cvo_records', '"{}"'.format(small_model_cvo_records)]
    )

  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/postprocess_variants.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def vcf_stats_report_command(
    vcf_path: str, title: Optional[str] = None
) -> tuple[str, Optional[str]]:
  """Returns a vcf_stats_report (command, logfile) for subprocess.

  Args:
    vcf_path: Path to VCF, which will be passed to --input_vcf and
      suffix-trimmed for --outfile_base.
    title: Passed straight to command unless it's None.

  Returns:
    [command string for subprocess, optional log directory path]
  """
  command = ['time', '/opt/deepvariant/bin/vcf_stats_report']
  command.extend(['--input_vcf', '"{}"'.format(vcf_path)])
  outfile_base = trim_suffix(trim_suffix(vcf_path, '.gz'), '.vcf')
  command.extend(['--outfile_base', '"{}"'.format(outfile_base)])
  if title is not None:
    command.extend(['--title', '"{}"'.format(title)])

  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/vcf_stats_report.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def runtime_by_region_vis_command(
    runtime_by_region_path: str, title: str = 'DeepVariant'
) -> tuple[str, None]:
  """Returns a runtime_by_region_vis (command, logfile=None) for subprocess."""
  runtime_report = os.path.join(
      _LOGGING_DIR.value, 'make_examples_runtime_by_region_report.html'
  )

  command = ['time', '/opt/deepvariant/bin/runtime_by_region_vis']
  command.extend(['--input', '"{}"'.format(runtime_by_region_path)])
  command.extend(['--title', '"{}"'.format(title)])
  command.extend(['--output', '"{}"'.format(runtime_report)])

  return (' '.join(command), None)


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


def check_flags():
  """Additional logic to make sure flags are set appropriately."""
  if _CUSTOMIZED_MODEL.value is not None:
    if not tf.io.gfile.exists(
        _CUSTOMIZED_MODEL.value + '.data-00000-of-00001'
    ) or not tf.io.gfile.exists(_CUSTOMIZED_MODEL.value + '.index'):
      raise RuntimeError(
          'The model files {}* do not exist. Potentially '
          'relevant issue: '
          'https://github.com/google/deepvariant/blob/r1.8/docs/'
          'FAQ.md#why-cant-it-find-one-of-the-input-files-eg-'
          'could-not-open'.format(_CUSTOMIZED_MODEL.value)
      )
    logging.info(
        (
            'You set --customized_model. Instead of using the default '
            'model for %s, `call_variants` step will load %s* '
            'instead.'
        ),
        _MODEL_TYPE.value,
        _CUSTOMIZED_MODEL.value,
    )
  if _CUSTOMIZED_SMALL_MODEL.value is not None:
    logging.info(
        (
            'You set --customized_small_model. Instead of using the default'
            ' small model for %s, `make_examples` will load %s* instead. Make'
            ' sure you set the GQ thresholds explicitly via'
            ' make_examples_extra_args.'
        ),
        _MODEL_TYPE.value,
        _CUSTOMIZED_SMALL_MODEL.value,
    )


def get_model_ckpt(model_type, customized_model):  # pylint: disable=unused-argument
  """Return the path to the model checkpoint based on the input args."""
  if customized_model is not None:
    return customized_model
  else:
    return MODEL_TYPE_MAP[model_type]


def create_all_commands_and_logfiles(
    intermediate_results_dir: str, used_in_test: bool = False
):
  """Creates 3 (command, logfile) to be executed later."""
  if not used_in_test:
    check_flags()
  commands = []
  # make_examples_pangenome_aware_dv
  nonvariant_site_tfrecord_path = None
  if _OUTPUT_GVCF.value is not None:
    nonvariant_site_tfrecord_path = os.path.join(
        intermediate_results_dir,
        'gvcf.tfrecord@{}.gz'.format(_NUM_SHARDS.value),
    )

  examples = os.path.join(
      intermediate_results_dir,
      'make_examples_pangenome_aware_dv.tfrecord@{}.gz'.format(
          _NUM_SHARDS.value
      ),
  )
  small_model_cvo_records = os.path.join(
      intermediate_results_dir,
      'make_examples_pangenome_aware_dv_call_variant_outputs.tfrecord@{}.gz'
      .format(_NUM_SHARDS.value),
  )

  if _LOGGING_DIR.value and _RUNTIME_REPORT.value:
    runtime_directory = os.path.join(
        _LOGGING_DIR.value, 'make_examples_runtime_by_region'
    )
    if not os.path.isdir(runtime_directory):
      logging.info(
          'Creating a make_examples_pangenome_aware_dv runtime by region'
          ' directory in %s',
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

  # Load pangenome GBZ into shared memory.
  if _PANGENOME.value is not None and _PANGENOME.value.endswith('.gbz'):
    commands.append(
        load_gbz_into_shared_memory_command(
            gbz=_PANGENOME.value,
            gbz_shared_memory_name=_GBZ_SHARED_MEMORY_NAME.value,
            gbz_shared_memory_size_gb=_GBZ_SHARED_MEMORY_SIZE_GB.value,
        )
    )
  model_ckpt = get_model_ckpt(_MODEL_TYPE.value, _CUSTOMIZED_MODEL.value)
  commands.append(
      make_examples_pangenome_aware_dv_command(
          ref=_REF.value,
          reads=_READS.value,
          pangenome=_PANGENOME.value,
          examples=examples,
          model_ckpt=model_ckpt,
          runtime_by_region_path=runtime_by_region_path,
          extra_args=_MAKE_EXAMPLES_EXTRA_ARGS.value,
          # kwargs:
          gvcf=nonvariant_site_tfrecord_path,
          regions=_REGIONS.value,
          sample_name_reads=_SAMPLE_NAME_READS.value,
          sample_name_pangenome=_SAMPLE_NAME_PANGENOME.value,
          gbz_shared_memory_name=_GBZ_SHARED_MEMORY_NAME.value,
      )
  )

  # call_variants
  call_variants_output = os.path.join(
      intermediate_results_dir, 'call_variants_output.tfrecord.gz'
  )
  commands.append(
      call_variants_command(
          outfile=call_variants_output,
          examples=examples,
          model_ckpt=model_ckpt,
          extra_args=_CALL_VARIANTS_EXTRA_ARGS.value,
          allow_empty_examples=_CUSTOMIZED_SMALL_MODEL.value is not None,
      )
  )

  # postprocess_variants
  commands.append(
      postprocess_variants_command(
          ref=_REF.value,
          infile=call_variants_output,
          outfile=_OUTPUT_VCF.value,
          small_model_cvo_records=small_model_cvo_records,
          extra_args=_POSTPROCESS_VARIANTS_EXTRA_ARGS.value,
          nonvariant_site_tfrecord_path=nonvariant_site_tfrecord_path,
          gvcf_outfile=_OUTPUT_GVCF.value,
      )
  )

  # vcf_stats_report
  if _VCF_STATS_REPORT.value:
    commands.append(
        vcf_stats_report_command(
            vcf_path=_OUTPUT_VCF.value, title=_REPORT_TITLE.value
        )
    )

  # runtime-by-region
  if _LOGGING_DIR.value and _RUNTIME_REPORT.value:
    commands.append(
        runtime_by_region_vis_command(
            runtime_by_region_path, title=_REPORT_TITLE.value
        )
    )

  return commands


def main(_):
  if _VERSION.value:
    print(
        'Pangenome-aware DeepVariant: DeepVariant version {}'.format(
            DEEP_VARIANT_VERSION
        )
    )
    return

  for flag_key in [
      'model_type',
      'ref',
      'reads',
      'output_vcf',
  ]:
    if FLAGS.get_flag_value(flag_key, None) is None:
      sys.stderr.write('--{} is required.\n'.format(flag_key))
      sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
      sys.exit(1)

  intermediate_results_dir = check_or_create_intermediate_results_dir(
      _INTERMEDIATE_RESULTS_DIR.value
  )

  if _LOGGING_DIR.value and not os.path.isdir(_LOGGING_DIR.value):
    logging.info('Creating a directory for logs in %s', _LOGGING_DIR.value)
    os.makedirs(_LOGGING_DIR.value)

  commands_logfiles = create_all_commands_and_logfiles(intermediate_results_dir)
  print(
      '\n***** Intermediate results will be written to {} '
      'in docker. ****\n'.format(intermediate_results_dir)
  )
  for command, logfile in commands_logfiles:
    print('\n***** Running the command:*****\n{}\n'.format(command))
    if not _DRY_RUN.value:
      fp = open(logfile, 'w') if logfile is not None else None
      with subprocess.Popen(
          command,
          stdout=subprocess.PIPE,
          stderr=subprocess.STDOUT,
          bufsize=1,
          shell=True,
          executable='/bin/bash',
          universal_newlines=True,
      ) as proc:
        for line in proc.stdout:
          print(line, end='')
          if fp is not None:
            print(line, end='', file=fp)
      if fp is not None:
        fp.close()
      if proc.returncode != 0:
        sys.exit(proc.returncode)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'model_type',
      'output_vcf',
      'pangenome',
      'reads',
      'ref',
  ])
  app.run(main)
