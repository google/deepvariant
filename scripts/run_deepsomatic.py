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

This script is used to run DeepSomatic, which is an extension of DeepVariant.
If you want to access more flags that are available in `make_examples_somatic`,
`call_variants`, and `postprocess_variants`, you can also call them separately
using the binaries in the Docker image.
"""

import os
import re
import subprocess
import sys
import tempfile
from typing import Any, Dict, Optional, Tuple

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

FLAGS = flags.FLAGS

# Required flags.
_MODEL_TYPE = flags.DEFINE_enum(
    'model_type',
    'WGS',
    ['WGS', 'ONT_R104', 'PACBIO'],
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
_READS_TUMOR = flags.DEFINE_string(
    'reads_tumor',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing the reads from '
        ' the tumor sample we want to call. Should be aligned to a reference '
        'genome compatible with --ref.'
    ),
)
_READS_NORMAL = flags.DEFINE_string(
    'reads_normal',
    None,
    (
        'Required. Aligned, sorted, indexed BAM file containing the reads from'
        ' the normal sample. Should be aligned to a reference genome '
        'compatible with --ref.'
        'NOTE: It is possible to not provide this flag. However, it is not '
        'officially supported yet.'
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
        'Output make_examples_somatic runtime metrics '
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
# Optional flags for make_examples_somatic.
_NUM_SHARDS = flags.DEFINE_integer(
    'num_shards',
    1,
    'Optional. Number of shards for make_examples_somatic step.',
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
_SAMPLE_NAME_TUMOR = flags.DEFINE_string(
    'sample_name_tumor',
    None,
    (
        'Sample name to use instead of the sample name from the input '
        'reads_tumor BAM (SM tag in the header). This flag is used for both '
        'make_examples_somatic and postprocess_variants.'
    ),
)
_SAMPLE_NAME_NORMAL = flags.DEFINE_string(
    'sample_name_normal',
    None,
    (
        'Sample name to use instead of the sample name from the input '
        'reads_normal BAM (SM tag in the header). This flag is used for both '
        'make_examples_somatic and postprocess_variants.'
    ),
)
_MAKE_EXAMPLES_EXTRA_ARGS = flags.DEFINE_string(
    'make_examples_extra_args',
    None,
    (
        'A comma-separated list of flag_name=flag_value. "flag_name" has to be'
        ' valid flags for make_examples_somatic.py. If the flag_value is '
        'boolean, it has to be flag_name=true or flag_name=false.'
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

# Optional flags for postprocess_variants.
_OUTPUT_GVCF = flags.DEFINE_string(
    'output_gvcf', None, 'Optional. Path where we should write gVCF file.'
)

# Optional flags for vcf_stats_report.
_VCF_STATS_REPORT = flags.DEFINE_boolean(
    'vcf_stats_report',
    True,
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
_PROCESS_SOMATIC = flags.DEFINE_boolean(
    'process_somatic',
    True,
    'Optional. Is specified the input is treated as somatic.',
)


MODEL_TYPE_MAP = {
    'WGS': '/opt/models/deepsomatic/wgs',
    'PACBIO': '/opt/models/deepsomatic/pacbio',
}

# Current release version of DeepVariant.
# Should be the same in dv_vcf_constants.py.
DEEP_VARIANT_VERSION = '1.6.1'


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


def split_extra_args(input_string: str) -> list[str]:
  """Splits into strs that do not contain commas or are enclosed in quotes."""
  pattern = r"[^,]+=[\"'][^\"']*[\"']|[^,]+"
  return re.findall(pattern, input_string)


def _extra_args_to_dict(extra_args: str) -> Dict[str, Any]:
  """Parses comma-separated list of flag_name=flag_value to dict."""
  args_dict = {}
  if extra_args is None:
    return args_dict
  for extra_arg in split_extra_args(extra_args):
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


def make_examples_somatic_command(
    ref,
    reads_tumor,
    reads_normal,
    examples,
    model_ckpt,
    extra_args,
    runtime_by_region_path=None,
    **kwargs,
) -> Tuple[str, Optional[str]]:
  """Returns a make_examples_somatic (command, logfile) for subprocess.

  Args:
    ref: Input FASTA file.
    reads_tumor: Input tumor BAM file.
    reads_normal: Input normal BAM file.
    examples: Output tfrecord file containing tensorflow.Example files.
    model_ckpt: Path to the TensorFlow model checkpoint.
    extra_args: Comma-separated list of flag_name=flag_value.
    runtime_by_region_path: Output path for runtime by region metrics.
    **kwargs: Additional arguments to pass in for make_examples_somatic.

  Returns:
    (string, string) A command to run, and a log file to output to.
  """
  command = [
      'time',
      'seq 0 {} |'.format(_NUM_SHARDS.value - 1),
      'parallel -q --halt 2 --line-buffer',
      '/opt/deepvariant/bin/make_examples_somatic',
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads_tumor', '"{}"'.format(reads_tumor)])
  if reads_normal:
    command.extend(['--reads_normal', '"{}"'.format(reads_normal)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])

  if runtime_by_region_path is not None:
    command.extend(
        ['--runtime_by_region', '"{}"'.format(runtime_by_region_path)]
    )

  if _MODEL_TYPE.value == 'WGS':
    # Specific flags that are not default can be added here.
    special_args = {}
    special_args['vsc_min_fraction_indels'] = 0.05
    special_args['vsc_min_fraction_snps'] = 0.029
    special_args['vsc_max_fraction_indels_for_non_target_sample'] = 0.5
    special_args['vsc_max_fraction_snps_for_non_target_sample'] = 0.5
    kwargs = _update_kwargs_with_warning(kwargs, special_args)
  elif _MODEL_TYPE.value == 'PACBIO':
    special_args = {}
    special_args['alt_aligned_pileup'] = 'diff_channels'
    special_args['min_mapping_quality'] = 5
    special_args['parse_sam_aux_fields'] = True
    special_args['partition_size'] = 25000
    special_args['phase_reads'] = True
    special_args['pileup_image_width'] = 199
    special_args['realign_reads'] = False
    special_args['sort_by_haplotypes'] = True
    special_args['track_ref_reads'] = True
    special_args['vsc_max_fraction_indels_for_non_target_sample'] = 0.5
    special_args['vsc_max_fraction_snps_for_non_target_sample'] = 0.5
    special_args['vsc_min_count_snps'] = 1
    special_args['vsc_min_fraction_indels'] = 0.08
    special_args['vsc_min_fraction_snps'] = 0.02
    special_args['trim_reads_for_pileup'] = True
    kwargs = _update_kwargs_with_warning(kwargs, special_args)
  elif _MODEL_TYPE.value == 'ONT_R104':
    special_args = {}
    special_args['alt_aligned_pileup'] = 'diff_channels'
    special_args['min_mapping_quality'] = 5
    special_args['parse_sam_aux_fields'] = True
    special_args['partition_size'] = 25000
    special_args['phase_reads'] = True
    special_args['pileup_image_width'] = 199
    special_args['realign_reads'] = False
    special_args['sort_by_haplotypes'] = True
    special_args['track_ref_reads'] = True
    special_args['vsc_min_fraction_snps'] = 0.05
    special_args['vsc_min_fraction_indels'] = 0.08
    special_args['trim_reads_for_pileup'] = True
    kwargs = _update_kwargs_with_warning(kwargs, special_args)
  else:
    raise ValueError('Invalid model_type: %s' % _MODEL_TYPE.value)

  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/make_examples_somatic.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def call_variants_command(
    outfile: str,
    examples: str,
    model_ckpt: str,
    extra_args: str,
) -> Tuple[str, Optional[str]]:
  """Returns a call_variants (command, logfile) for subprocess."""
  binary_name = 'call_variants'
  command = ['time', f'/opt/deepvariant/bin/{binary_name}']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])
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
    process_somatic: bool,
    extra_args: str,
    **kwargs,
) -> Tuple[str, Optional[str]]:
  """Returns a postprocess_variants (command, logfile) for subprocess."""
  command = ['time', '/opt/deepvariant/bin/postprocess_variants']
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--infile', '"{}"'.format(infile)])
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--cpus 0'])
  if process_somatic:
    command.extend(['--process_somatic=true'])
  else:
    command.extend(['--noprocess_somatic'])
  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/postprocess_variants.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def vcf_stats_report_command(
    vcf_path: str, title: Optional[str] = None
) -> Tuple[str, Optional[str]]:
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
) -> Tuple[str, None]:
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
    if (
        not tf.io.gfile.exists(_CUSTOMIZED_MODEL.value + '.data-00000-of-00001')
        or not tf.io.gfile.exists(_CUSTOMIZED_MODEL.value + '.index')
    ):
      raise RuntimeError(
          'The model files {}* do not exist. Potentially '
          'relevant issue: '
          'https://github.com/google/deepvariant/blob/r1.6.1/docs/'
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
  # make_examples_somatic
  nonvariant_site_tfrecord_path = None
  if _OUTPUT_GVCF.value is not None:
    nonvariant_site_tfrecord_path = os.path.join(
        intermediate_results_dir,
        'gvcf.tfrecord@{}.gz'.format(_NUM_SHARDS.value),
    )

  examples = os.path.join(
      intermediate_results_dir,
      'make_examples_somatic.tfrecord@{}.gz'.format(_NUM_SHARDS.value),
  )

  if _LOGGING_DIR.value and _RUNTIME_REPORT.value:
    runtime_directory = os.path.join(
        _LOGGING_DIR.value, 'make_examples_runtime_by_region'
    )
    if not os.path.isdir(runtime_directory):
      logging.info(
          'Creating a make_examples_somatic runtime by region directory in %s',
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

  model_ckpt = get_model_ckpt(_MODEL_TYPE.value, _CUSTOMIZED_MODEL.value)
  commands.append(
      make_examples_somatic_command(
          ref=_REF.value,
          reads_tumor=_READS_TUMOR.value,
          reads_normal=_READS_NORMAL.value,
          examples=examples,
          model_ckpt=model_ckpt,
          runtime_by_region_path=runtime_by_region_path,
          extra_args=_MAKE_EXAMPLES_EXTRA_ARGS.value,
          # kwargs:
          gvcf=nonvariant_site_tfrecord_path,
          regions=_REGIONS.value,
          sample_name_tumor=_SAMPLE_NAME_TUMOR.value,
          sample_name_normal=_SAMPLE_NAME_NORMAL.value,
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
      )
  )

  # postprocess_variants
  commands.append(
      postprocess_variants_command(
          ref=_REF.value,
          infile=call_variants_output,
          outfile=_OUTPUT_VCF.value,
          process_somatic=_PROCESS_SOMATIC.value,
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
    print('DeepSomatic: DeepVariant version {}'.format(DEEP_VARIANT_VERSION))
    return

  for flag_key in [
      'model_type',
      'ref',
      'reads_tumor',
      'output_vcf',
  ]:
    if FLAGS.get_flag_value(flag_key, None) is None:
      sys.stderr.write('--{} is required.\n'.format(flag_key))
      sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
      sys.exit(1)

  for flag_key in ['reads_normal']:
    if FLAGS.get_flag_value(flag_key, None) is None:
      sys.stderr.write('--{} is not set.\n'.format(flag_key))
      sys.stderr.write('WARNING: This is not officially supported.\n')

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
  env = os.environ.copy()
  logging.info('env = %s', env)
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
          env=env,
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
  app.run(main)
