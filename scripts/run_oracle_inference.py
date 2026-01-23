# Copyright 2025 Google LLC.
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
"""Runs make_examples in training mode and labeled_examples_to_vcf.

This script runs make_examples in training mode to generate tf.Examples and
then runs labeled_examples_to_vcf to generate a VCF file from the labels in
tf.Examples.
"""

import enum
import os
import re
import subprocess
import sys
import tempfile
from typing import Any, Optional

from absl import app
from absl import flags
from absl import logging

FLAGS = flags.FLAGS


class ModelType(enum.Enum):
  WGS = 'WGS'
  WES = 'WES'
  PACBIO = 'PACBIO'
  ONT_R104 = 'ONT_R104'
  HYBRID_PACBIO_ILLUMINA = 'HYBRID_PACBIO_ILLUMINA'
  MASSEQ = 'MASSEQ'


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
_OUTPUT_VCF = flags.DEFINE_string(
    'output_vcf', None, 'Required. Path where we should write VCF file.'
)
_TRUTH_VARIANTS = flags.DEFINE_string(
    'truth_variants',
    None,
    'Required. VCF file containing truth variants.',
)
_CONFIDENT_REGIONS = flags.DEFINE_string(
    'confident_regions',
    None,
    'Required. BED file containing confident regions.',
)
_LABELER_ALGORITHM = flags.DEFINE_string(
    'labeler_algorithm',
    None,
    'Required. Labeler algorithm to use for calling variants.',
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
_VERSION = flags.DEFINE_boolean(
    'version',
    None,
    'Optional. If true, print out version number and exit.',
    allow_hide_cpp=True,
)
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
_MAKE_EXAMPLES_EXTRA_ARGS = flags.DEFINE_string(
    'make_examples_extra_args',
    None,
    (
        'A comma-separated list of flag_name=flag_value. "flag_name" has to be'
        ' valid flags for make_examples.py. If the flag_value is boolean, it'
        ' has to be flag_name=true or flag_name=false.'
    ),
)


# Current release version of DeepVariant.
# Should be the same in dv_vcf_constants.py.
DEEP_VARIANT_VERSION = '1.10.0'


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


def _extra_args_to_dict(extra_args: str) -> dict[str, Any]:
  """Parses comma-separated list of flag_name=flag_value to dict."""
  args_dict = {}
  if extra_args is None:
    return args_dict
  for extra_arg in split_extra_args(extra_args):
    flag_name, flag_value = extra_arg.split('=')
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
    ref: str,
    reads: str,
    examples: str,
    labeler_algorithm: str,
    extra_args: str | None,
    **kwargs,
) -> tuple[str, str | None]:
  """Returns a make_examples (command, logfile) for subprocess.

  Args:
    ref: Input FASTA file.
    reads: Input BAM file.
    examples: Output tfrecord file containing tensorflow.Example files.
    labeler_algorithm: Labeler algorithm to use for calling variants.
    extra_args: Comma-separated list of flag_name=flag_value.
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
  command.extend(['--mode', 'training'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads', '"{}"'.format(reads)])
  command.extend(['--labeler_algorithm', '"{}"'.format(labeler_algorithm)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--channel_list', '"BASE_CHANNELS"'])

  command.extend(['--max_reads_per_partition', '1500'])
  partition_size = 1000
  if _MODEL_TYPE.value in (ModelType.PACBIO.value, ModelType.ONT_R104.value):
    partition_size = 25000
  command.extend(['--partition_size', '"{}"'.format(partition_size)])

  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])
  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/make_examples.log'.format(_LOGGING_DIR.value)
  return (' '.join(command), logfile)


def labeled_examples_to_vcf_command(
    ref: str,
    examples: str,
    outfile: str,
    sample_name: str,
) -> tuple[str, Optional[str]]:
  """Returns a labeled_examples_to_vcf (command, logfile) for subprocess.

  Args:
    ref: Input FASTA file.
    examples: Input tfrecord file containing tensorflow.Example files.
    outfile: Output VCF file.
    sample_name: Sample name to use for the VCF.

  Returns:
    (string, string) A command to run, and a log file to output to.
  """
  binary_name = 'labeled_examples_to_vcf'
  command = ['time', f'/opt/deepvariant/bin/{binary_name}']
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--output_vcf', '"{}"'.format(outfile)])
  if sample_name:
    command.extend(['--sample_name', '"{}"'.format(sample_name)])

  logfile = None
  if _LOGGING_DIR.value:
    logfile = '{}/{}.log'.format(_LOGGING_DIR.value, binary_name)
  return (' '.join(command), logfile)


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
  pass


def create_all_commands_and_logfiles(
    intermediate_results_dir: str,
) -> list[tuple[str, str | None]]:
  """Creates 3 (command, logfile) to be executed later.

  Args:
    intermediate_results_dir: Directory to store intermediate results.

  Returns:
    A list of (command, logfile) tuples.
  """
  check_flags()
  commands = []
  # make_examples
  examples = os.path.join(
      intermediate_results_dir,
      'make_examples.tfrecord@{}.gz'.format(_NUM_SHARDS.value),
  )
  commands.append(
      make_examples_command(
          ref=_REF.value,
          reads=_READS.value,
          examples=examples,
          labeler_algorithm=_LABELER_ALGORITHM.value,
          extra_args=_MAKE_EXAMPLES_EXTRA_ARGS.value,
          # kwargs:
          truth_variants=_TRUTH_VARIANTS.value,
          confident_regions=_CONFIDENT_REGIONS.value,
          regions=_REGIONS.value,
          sample_name=_SAMPLE_NAME.value,
          haploid_contigs=_HAPLOID_CONTIGS.value,
          par_regions_bed=_PAR_REGIONS.value,
      )
  )

  # labeled_examples_to_vcf
  commands.append(
      labeled_examples_to_vcf_command(
          ref=_REF.value,
          examples=examples,
          outfile=_OUTPUT_VCF.value,
          sample_name=_SAMPLE_NAME.value,
      )
  )

  return commands


def main(_):
  if _VERSION.value:
    print('DeepVariant version {}'.format(DEEP_VARIANT_VERSION))
    return

  for flag_key in [
      'model_type',
      'ref',
      'reads',
      'output_vcf',
      'truth_variants',
      'confident_regions',
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
