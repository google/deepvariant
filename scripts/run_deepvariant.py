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
"""Runs all 3 steps to go from input DNA reads to output VCF/gVCF files.

This script currently provides the most common use cases and standard models.
If you want to access more flags that are available in `make_examples`,
`call_variants`, and `postprocess_variants`, you can also call them separately
using the binaries in the Docker image.

For more details, see:
https://github.com/google/deepvariant/blob/r1.7/docs/deepvariant-quick-start.md
"""

import os
import re
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

FLAGS = flags.FLAGS

# Required flags.
_MODEL_TYPE = flags.DEFINE_enum(
    'model_type',
    None,
    ['WGS', 'WES', 'PACBIO', 'ONT_R104', 'HYBRID_PACBIO_ILLUMINA'],
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
_OUTPUT_VCF = flags.DEFINE_string(
    'output_vcf', None, 'Required. Path where we should write VCF file.'
)
# Optional flags.
_HAPLOID_CONTIGS = flags.DEFINE_string(
    'haploid_contigs',
    None,
    (
        'Optional list of non autosomal chromosomes. For all listed'
        ' chromosomes, HET probabilities are not considered. For samples with'
        ' XY karyotype it is expected to set --haploid_contigs="chrX,chrY" for'
        ' GRCh38 and --haploid_contigs="X,Y" for GRCh37. For samples with'
        ' XX karyotype --haploid_contigs flag should not be used.'
    ),
)

_PAR_REGIONS = flags.DEFINE_string(
    'par_regions_bed',
    None,
    (
        'Optional BED file containing Human Pseudoautosomal Region (PAR)'
        ' regions. This should be specific to the reference used. For example'
        ' GRCh38 PAR bed file would be different from GRCh37 bed file. Regions'
        ' in this bed file are treated as diploid, effectively subtracting them'
        ' from the --haploid_contigs.'
    ),
)

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
        'Output make_examples runtime metrics '
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
# TODO: Update once we have default trained small models.
_CUSTOMIZED_SMALL_MODEL = flags.DEFINE_string(
    'customized_small_model',
    None,
    (
        'Optional. A path to a small model checkpoint to call variants during '
        'the `make_examples` step. If not set, the small model is not used.'
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
_SAMPLE_NAME = flags.DEFINE_string(
    'sample_name',
    None,
    (
        'Sample name to use instead of the sample name from the input reads BAM'
        ' (SM tag in the header). This flag is used for both make_examples and'
        ' postprocess_variants.'
    ),
)
_USE_HP_INFORMATION = flags.DEFINE_boolean(
    'use_hp_information',
    None,
    (
        '(Deprecated in v1.4.0) Optional. If True, corresponding flags will be'
        ' set to properly use the HP information present in the BAM input.'
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


MODEL_TYPE_MAP = {
    'WGS': '/opt/models/wgs',
    'WES': '/opt/models/wes',
    'PACBIO': '/opt/models/pacbio',
    'ONT_R104': '/opt/models/ont_r104',
    'HYBRID_PACBIO_ILLUMINA': '/opt/models/hybrid_pacbio_illumina',
}

# Current release version of DeepVariant.
# Should be the same in dv_vcf_constants.py.
DEEP_VARIANT_VERSION = '1.7.0'


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


def split_extra_args(input_string: str) -> List[str]:
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


def make_examples_command(
    ref,
    reads,
    examples,
    model_ckpt,
    extra_args,
    runtime_by_region_path=None,
    **kwargs,
):
  """Returns a make_examples (command, logfile) for subprocess.

  Args:
    ref: Input FASTA file.
    reads: Input BAM file.
    examples: Output tfrecord file containing tensorflow.Example files.
    model_ckpt: Path to the TensorFlow model checkpoint.
    extra_args: Comma-separated list of flag_name=flag_value.
    runtime_by_region_path: Output path for runtime by region metrics.
    **kwargs: Additional arguments to pass in for make_examples.

  Returns:
    (string, string) A command to run, and a log file to output to.
  """
  command = [
      'time',
      'seq 0 {} |'.format(_NUM_SHARDS.value - 1),
      'parallel -q --halt 2 --line-buffer',
      '/opt/deepvariant/bin/make_examples',
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads', '"{}"'.format(reads)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])

  if runtime_by_region_path is not None:
    command.extend(
        ['--runtime_by_region', '"{}"'.format(runtime_by_region_path)]
    )

  special_args = {}
  if _MODEL_TYPE.value == 'PACBIO':
    special_args['alt_aligned_pileup'] = 'diff_channels'
    special_args['max_reads_per_partition'] = 600
    special_args['min_mapping_quality'] = 1
    special_args['parse_sam_aux_fields'] = True
    special_args['partition_size'] = 25000
    special_args['phase_reads'] = True
    special_args['pileup_image_width'] = 199
    special_args['realign_reads'] = False
    special_args['sort_by_haplotypes'] = True
    special_args['track_ref_reads'] = True
    special_args['vsc_min_fraction_indels'] = 0.12
    special_args['trim_reads_for_pileup'] = True
  elif _MODEL_TYPE.value == 'ONT_R104':
    special_args['alt_aligned_pileup'] = 'diff_channels'
    special_args['max_reads_per_partition'] = 600
    special_args['min_mapping_quality'] = 5
    special_args['parse_sam_aux_fields'] = True
    special_args['partition_size'] = 25000
    special_args['phase_reads'] = True
    special_args['pileup_image_width'] = 199
    special_args['realign_reads'] = False
    special_args['sort_by_haplotypes'] = True
    special_args['track_ref_reads'] = True
    special_args['vsc_min_fraction_snps'] = 0.08
    special_args['vsc_min_fraction_indels'] = 0.12
    special_args['trim_reads_for_pileup'] = True
  elif _MODEL_TYPE.value == 'HYBRID_PACBIO_ILLUMINA':
    special_args['trim_reads_for_pileup'] = True

  if _CUSTOMIZED_SMALL_MODEL.value:
    special_args['trained_small_model_path'] = _CUSTOMIZED_SMALL_MODEL.value
    special_args['call_small_model_examples'] = True
    # Note that the default GQ threshold can be overwritten by extra_args.
    special_args['small_model_snp_gq_threshold'] = 30
    special_args['small_model_indel_gq_threshold'] = 30

  kwargs = _update_kwargs_with_warning(kwargs, special_args)
  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/make_examples.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def call_variants_command(
    outfile: str,
    examples: str,
    model_ckpt: str,
    extra_args: str,
    allow_empty_examples: bool,
) -> Tuple[str, Optional[str]]:
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
) -> Tuple[str, Optional[str]]:
  """Returns a postprocess_variants (command, logfile) for subprocess."""
  cpus = _POSTPROCESS_CPUS.value
  if cpus is None:
    cpus = _NUM_SHARDS.value
  command = ['time', '/opt/deepvariant/bin/postprocess_variants']
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--infile', '"{}"'.format(infile)])
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--cpus', '"{}"'.format(cpus)])
  if _CUSTOMIZED_SMALL_MODEL.value:
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
    use_saved_model = tf.io.gfile.exists(
        _CUSTOMIZED_MODEL.value
    ) and tf.io.gfile.exists(f'{_CUSTOMIZED_MODEL.value}/saved_model.pb')

    if use_saved_model:
      logging.info('Using saved model: %s', str(use_saved_model))
    elif (
        not tf.io.gfile.exists(_CUSTOMIZED_MODEL.value + '.data-00000-of-00001')
        or not tf.io.gfile.exists(_CUSTOMIZED_MODEL.value + '.index')
    ):
      raise RuntimeError(
          'The model files {}* do not exist. Potentially '
          'relevant issue: '
          'https://github.com/google/deepvariant/blob/r1.7/docs/'
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


def get_model_ckpt(model_type, customized_model):
  """Return the path to the model checkpoint based on the input args."""
  if customized_model is not None:
    return customized_model
  else:
    return MODEL_TYPE_MAP[model_type]


def create_all_commands_and_logfiles(intermediate_results_dir):
  """Creates 3 (command, logfile) to be executed later."""
  check_flags()
  commands = []
  # make_examples
  nonvariant_site_tfrecord_path = None
  if _OUTPUT_GVCF.value is not None:
    nonvariant_site_tfrecord_path = os.path.join(
        intermediate_results_dir,
        'gvcf.tfrecord@{}.gz'.format(_NUM_SHARDS.value),
    )

  examples = os.path.join(
      intermediate_results_dir,
      'make_examples.tfrecord@{}.gz'.format(_NUM_SHARDS.value),
  )
  small_model_cvo_records = os.path.join(
      intermediate_results_dir,
      'make_examples_call_variant_outputs.tfrecord@{}.gz'.format(
          _NUM_SHARDS.value
      ),
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

  model_ckpt = get_model_ckpt(_MODEL_TYPE.value, _CUSTOMIZED_MODEL.value)
  commands.append(
      make_examples_command(
          ref=_REF.value,
          reads=_READS.value,
          examples=examples,
          model_ckpt=model_ckpt,
          runtime_by_region_path=runtime_by_region_path,
          extra_args=_MAKE_EXAMPLES_EXTRA_ARGS.value,
          # kwargs:
          gvcf=nonvariant_site_tfrecord_path,
          regions=_REGIONS.value,
          sample_name=_SAMPLE_NAME.value,
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
          sample_name=_SAMPLE_NAME.value,
          haploid_contigs=_HAPLOID_CONTIGS.value,
          par_regions_bed=_PAR_REGIONS.value,
          regions=_REGIONS.value,
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
  if _USE_HP_INFORMATION.value:
    raise NotImplementedError(
        'The --use_hp_information flag has been '
        'deprecated. DeepVariant now phases internally '
        'for PacBio mode.'
    )
  if _VERSION.value:
    print('DeepVariant version {}'.format(DEEP_VARIANT_VERSION))
    return

  for flag_key in ['model_type', 'ref', 'reads', 'output_vcf']:
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
