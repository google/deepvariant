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
https://github.com/google/deepvariant/blob/r1.4/docs/deepvariant-quick-start.md
"""

import os
import subprocess
import sys
import tempfile

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

FLAGS = flags.FLAGS

# Required flags.
flags.DEFINE_enum(
    'model_type', None, ['WGS', 'WES', 'PACBIO', 'HYBRID_PACBIO_ILLUMINA'],
    'Required. Type of model to use for variant calling. Set this flag '
    'to use the default model associated with each type, and it will set '
    'necessary flags corresponding to each model. If you want to use a '
    'customized model, add --customized_model flag in addition to this flag.')
flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads.')
flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref.')
flags.DEFINE_string('output_vcf', None,
                    'Required. Path where we should write VCF file.')
# Optional flags.
flags.DEFINE_boolean(
    'dry_run', False,
    'Optional. If True, only prints out commands without executing them.')
flags.DEFINE_string(
    'intermediate_results_dir', None,
    'Optional. If specified, this should be an existing '
    'directory that is visible insider docker, and will be '
    'used to to store intermediate outputs.')
flags.DEFINE_string(
    'logging_dir', None, 'Optional. Directory where we should write log files '
    'for each stage and optionally runtime reports.')
flags.DEFINE_boolean(
    'runtime_report', False, 'Output make_examples runtime metrics '
    'and create a visual runtime report using runtime_by_region_vis. '
    'Only works with --logging_dir.')
flags.DEFINE_boolean(
    'version',
    None,
    'Optional. If true, print out version number and exit.',
    allow_hide_cpp=True)
# Optional flags for call_variants.
flags.DEFINE_string(
    'customized_model', None,
    'Optional. A path to a model checkpoint to load for the `call_variants` '
    'step. If not set, the default for each --model_type will be used')
# Optional flags for make_examples.
flags.DEFINE_integer('num_shards', 1,
                     'Optional. Number of shards for make_examples step.')
flags.DEFINE_string(
    'regions', None,
    'Optional. Space-separated list of regions we want to process. Elements '
    'can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.')
flags.DEFINE_string(
    'sample_name', None,
    'Sample name to use instead of the sample name from the input reads BAM '
    '(SM tag in the header). This flag is used for both make_examples and '
    'postprocess_variants.')
flags.DEFINE_boolean(
    'use_hp_information', None, '(Deprecated in v1.4.0) '
    'Optional. If True, corresponding flags will be set to properly use the HP '
    'information present in the BAM input.')
flags.DEFINE_string(
    'make_examples_extra_args', None,
    'A comma-separated list of flag_name=flag_value. "flag_name" has to be '
    'valid flags for make_examples.py. If the flag_value is boolean, it has to '
    'be flag_name=true or flag_name=false.')
flags.DEFINE_string(
    'call_variants_extra_args', None,
    'A comma-separated list of flag_name=flag_value. "flag_name" has to be '
    'valid flags for call_variants.py. If the flag_value is boolean, it has to '
    'be flag_name=true or flag_name=false.')
flags.DEFINE_string(
    'postprocess_variants_extra_args', None,
    'A comma-separated list of flag_name=flag_value. "flag_name" has to be '
    'valid flags for postprocess_variants.py. If the flag_value is boolean, '
    'it has to be flag_name=true or flag_name=false.')

# Optional flags for postprocess_variants.
flags.DEFINE_string('output_gvcf', None,
                    'Optional. Path where we should write gVCF file.')
flags.DEFINE_boolean(
    'vcf_stats_report', True, 'Optional. Output a visual report (HTML) of '
    'statistics about the output VCF.')

MODEL_TYPE_MAP = {
    'WGS': '/opt/models/wgs/model.ckpt',
    'WES': '/opt/models/wes/model.ckpt',
    'PACBIO': '/opt/models/pacbio/model.ckpt',
    'HYBRID_PACBIO_ILLUMINA': '/opt/models/hybrid_pacbio_illumina/model.ckpt',
}

# Current release version of DeepVariant.
# Should be the same in dv_vcf_constants.py.
DEEP_VARIANT_VERSION = '1.4.0'


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
        print('\nWarning: --{} is previously set to {}, now to {}.'.format(
            k, kwargs[k], v))
    kwargs[k] = v
  return kwargs


def make_examples_command(ref,
                          reads,
                          examples,
                          extra_args,
                          runtime_by_region_path=None,
                          **kwargs):
  """Returns a make_examples (command, logfile) for subprocess.

  Args:
    ref: Input FASTA file.
    reads: Input BAM file.
    examples: Output tfrecord file containing tensorflow.Example files.
    extra_args: Comma-separated list of flag_name=flag_value.
    runtime_by_region_path: Output path for runtime by region metrics.
    **kwargs: Additional arguments to pass in for make_examples.

  Returns:
    (string, string) A command to run, and a log file to output to.
  """
  command = [
      'time', 'seq 0 {} |'.format(FLAGS.num_shards - 1),
      'parallel -q --halt 2 --line-buffer', '/opt/deepvariant/bin/make_examples'
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads', '"{}"'.format(reads)])
  command.extend(['--examples', '"{}"'.format(examples)])

  if runtime_by_region_path is not None:
    command.extend(
        ['--runtime_by_region', '"{}"'.format(runtime_by_region_path)])

  if FLAGS.model_type == 'WGS' or FLAGS.model_type == 'WES':
    special_args = {}
    special_args['channels'] = 'insert_size'
    kwargs = _update_kwargs_with_warning(kwargs, special_args)
  elif FLAGS.model_type == 'PACBIO':
    special_args = {}
    special_args['add_hp_channel'] = True
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
    kwargs = _update_kwargs_with_warning(kwargs, special_args)
  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args))
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])
  logfile = None
  if FLAGS.logging_dir:
    logfile = '{}/make_examples.log'.format(FLAGS.logging_dir)
  return (' '.join(command), logfile)


def call_variants_command(outfile, examples, model_ckpt,
                          intermediate_results_dir, extra_args):
  """Returns a call_variants (command, logfile) for subprocess."""
  command = ['time', '/opt/deepvariant/bin/call_variants']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])
  # --openvino_model_dir will only be used if use_openvino is set to true for
  # call_variants. But it won't hurt to set it anyway, so setting it to the
  # intermediate_results_dir.
  command.extend(
      ['--openvino_model_dir', '"{}"'.format(intermediate_results_dir)])
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(command,
                                         _extra_args_to_dict(extra_args))
  logfile = None
  if FLAGS.logging_dir:
    logfile = '{}/call_variants.log'.format(FLAGS.logging_dir)
  return (' '.join(command), logfile)


def postprocess_variants_command(ref,
                                 infile,
                                 outfile,
                                 extra_args,
                                 nonvariant_site_tfrecord_path=None,
                                 gvcf_outfile=None,
                                 vcf_stats_report=True,
                                 sample_name=None):
  """Returns a postprocess_variants (command, logfile) for subprocess."""
  command = ['time', '/opt/deepvariant/bin/postprocess_variants']
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--infile', '"{}"'.format(infile)])
  command.extend(['--outfile', '"{}"'.format(outfile)])
  if nonvariant_site_tfrecord_path is not None:
    command.extend([
        '--nonvariant_site_tfrecord_path',
        '"{}"'.format(nonvariant_site_tfrecord_path)
    ])
  if gvcf_outfile is not None:
    command.extend(['--gvcf_outfile', '"{}"'.format(gvcf_outfile)])
  if not vcf_stats_report:
    command.extend(['--novcf_stats_report'])
  if sample_name is not None:
    command.extend(['--sample_name', '"{}"'.format(sample_name)])
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(command,
                                         _extra_args_to_dict(extra_args))
  logfile = None
  if FLAGS.logging_dir:
    logfile = '{}/postprocess_variants.log'.format(FLAGS.logging_dir)
  return (' '.join(command), logfile)


def runtime_by_region_vis_command(runtime_by_region_path: str):
  """Returns a runtime_by_region_vis (command, logfile=None) for subprocess."""
  runtime_report = os.path.join(FLAGS.logging_dir,
                                'make_examples_runtime_by_region_report.html')

  command = ['time', '/opt/deepvariant/bin/runtime_by_region_vis']
  command.extend(['--input', '"{}"'.format(runtime_by_region_path)])
  command.extend(['--title', '"{}"'.format('DeepVariant')])
  command.extend(['--output', '"{}"'.format(runtime_report)])

  return (' '.join(command), None)


def check_or_create_intermediate_results_dir(intermediate_results_dir):
  """Checks or creates the path to the directory for intermediate results."""
  if intermediate_results_dir is None:
    intermediate_results_dir = tempfile.mkdtemp()
  if not os.path.isdir(intermediate_results_dir):
    logging.info('Creating a directory for intermediate results in %s',
                 intermediate_results_dir)
    os.makedirs(intermediate_results_dir)
  else:
    logging.info('Re-using the directory for intermediate results in %s',
                 intermediate_results_dir)
  return intermediate_results_dir


def check_flags():
  """Additional logic to make sure flags are set appropriately."""
  if FLAGS.customized_model is not None:
    if (not tf.compat.v1.gfile.Exists(FLAGS.customized_model +
                                      '.data-00000-of-00001') or
        not tf.compat.v1.gfile.Exists(FLAGS.customized_model + '.index') or
        not tf.compat.v1.gfile.Exists(FLAGS.customized_model + '.meta')):
      raise RuntimeError('The model files {}* do not exist. Potentially '
                         'relevant issue: '
                         'https://github.com/google/deepvariant/blob/r1.4/docs/'
                         'FAQ.md#why-cant-it-find-one-of-the-input-files-eg-'
                         'could-not-open'.format(FLAGS.customized_model))
    logging.info(
        'You set --customized_model. Instead of using the default '
        'model for %s, `call_variants` step will load %s* '
        'instead.', FLAGS.model_type, FLAGS.customized_model)


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
  if FLAGS.output_gvcf is not None:
    nonvariant_site_tfrecord_path = os.path.join(
        intermediate_results_dir,
        'gvcf.tfrecord@{}.gz'.format(FLAGS.num_shards))

  examples = os.path.join(
      intermediate_results_dir,
      'make_examples.tfrecord@{}.gz'.format(FLAGS.num_shards))

  if FLAGS.logging_dir and FLAGS.runtime_report:
    runtime_directory = os.path.join(FLAGS.logging_dir,
                                     'make_examples_runtime_by_region')
    if not os.path.isdir(runtime_directory):
      logging.info('Creating a make_examples runtime by region directory in %s',
                   runtime_directory)
      os.makedirs(runtime_directory)
    # The path to runtime metrics output is sharded just like the examples.
    runtime_by_region_path = os.path.join(
        runtime_directory,
        'make_examples_runtime@{}.tsv'.format(FLAGS.num_shards))
  else:
    runtime_by_region_path = None

  commands.append(
      make_examples_command(
          ref=FLAGS.ref,
          reads=FLAGS.reads,
          examples=examples,
          runtime_by_region_path=runtime_by_region_path,
          extra_args=FLAGS.make_examples_extra_args,
          # kwargs:
          gvcf=nonvariant_site_tfrecord_path,
          regions=FLAGS.regions,
          sample_name=FLAGS.sample_name))

  # call_variants
  call_variants_output = os.path.join(intermediate_results_dir,
                                      'call_variants_output.tfrecord.gz')
  model_ckpt = get_model_ckpt(FLAGS.model_type, FLAGS.customized_model)
  commands.append(
      call_variants_command(call_variants_output, examples, model_ckpt,
                            intermediate_results_dir,
                            FLAGS.call_variants_extra_args))

  # postprocess_variants
  commands.append(
      postprocess_variants_command(
          FLAGS.ref,
          call_variants_output,
          FLAGS.output_vcf,
          FLAGS.postprocess_variants_extra_args,
          nonvariant_site_tfrecord_path=nonvariant_site_tfrecord_path,
          gvcf_outfile=FLAGS.output_gvcf,
          vcf_stats_report=FLAGS.vcf_stats_report,
          sample_name=FLAGS.sample_name))

  # runtime-by-region
  if FLAGS.logging_dir and FLAGS.runtime_report:
    commands.append(runtime_by_region_vis_command(runtime_by_region_path))

  return commands


def main(_):
  if FLAGS.use_hp_information:
    raise NotImplementedError('The --use_hp_information flag has been '
                              'deprecated. DeepVariant now phases internally '
                              'for PacBio mode.')
  if FLAGS.version:
    print('DeepVariant version {}'.format(DEEP_VARIANT_VERSION))
    return

  for flag_key in ['model_type', 'ref', 'reads', 'output_vcf']:
    if FLAGS.get_flag_value(flag_key, None) is None:
      sys.stderr.write('--{} is required.\n'.format(flag_key))
      sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
      sys.exit(1)

  intermediate_results_dir = check_or_create_intermediate_results_dir(
      FLAGS.intermediate_results_dir)

  if FLAGS.logging_dir and not os.path.isdir(FLAGS.logging_dir):
    logging.info('Creating a directory for logs in %s', FLAGS.logging_dir)
    os.makedirs(FLAGS.logging_dir)

  commands_logfiles = create_all_commands_and_logfiles(intermediate_results_dir)
  print('\n***** Intermediate results will be written to {} '
        'in docker. ****\n'.format(intermediate_results_dir))
  for command, logfile in commands_logfiles:
    print('\n***** Running the command:*****\n{}\n'.format(command))
    if not FLAGS.dry_run:
      fp = open(logfile, 'w') if logfile is not None else None
      with subprocess.Popen(
          command,
          stdout=subprocess.PIPE,
          stderr=subprocess.STDOUT,
          bufsize=1,
          shell=True,
          executable='/bin/bash',
          universal_newlines=True) as proc:
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
