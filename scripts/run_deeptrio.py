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
    'Required. Type of model to use for variant calling. Each '
    'model_type has an associated default model, which can be '
    'overridden by the --customized_model_{parent,child} flags.')
flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads_child.')
flags.DEFINE_string(
    'reads_child', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref.')
flags.DEFINE_string(
    'reads_parent1', None,
    'Required. Aligned, sorted, indexed BAM file containing parent 1 reads of '
    'the person we want to call. Should be aligned to a reference genome '
    'compatible with --ref.')
flags.DEFINE_string(
    'reads_parent2', None,
    'Aligned, sorted, indexed BAM file containing parent 2 reads of '
    'the person we want to call. Should be aligned to a reference genome '
    'compatible with --ref.')
flags.DEFINE_string(
    'output_vcf_child', None,
    'Required. Path where we should write VCF file for the child.')
flags.DEFINE_string(
    'output_vcf_parent1', None,
    'Required. Path where we should write VCF file for parent1.')
flags.DEFINE_string(
    'output_vcf_parent2', None,
    'Required. Path where we should write VCF file for parent2.')
# Optional flags.
flags.DEFINE_boolean(
    'dry_run', False,
    'Optional. If True, only prints out commands without executing them.')
flags.DEFINE_string(
    'intermediate_results_dir', None,
    'Optional. If specified, this should be an existing '
    'directory that is visible insider docker, and will be '
    'used to to store intermediate outputs.')
flags.DEFINE_boolean(
    'version',
    None,
    'Optional. If true, print out version number and exit.',
    allow_hide_cpp=True)
flags.DEFINE_string('logging_dir', None,
                    'Required. Directory where we should write log files.')
flags.DEFINE_boolean(
    'runtime_report', False, 'Output make_examples runtime metrics '
    'and create a visual runtime report using runtime_by_region_vis. '
    'Only works with --logging_dir.')
# Optional flags for call_variants.
# TODO: Support different child/parent models.
flags.DEFINE_string(
    'customized_model_child', None,
    'Optional. A path to a child model checkpoint to load for the '
    '`call_variants` step. If not set, the default for each --model_type will '
    'be used')
flags.DEFINE_string(
    'customized_model_parent', None,
    'Optional. A path to a parent model checkpoint to load for the '
    '`call_variants` step. If not set, the default for each --model_type will '
    'be used')
# Optional flags for make_examples.
flags.DEFINE_integer('num_shards', 1,
                     'Optional. Number of shards for make_examples step.')
flags.DEFINE_string(
    'regions', None,
    'Optional. Space-separated list of regions we want to process. Elements '
    'can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.')
flags.DEFINE_string(
    'sample_name_child', None,
    'Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads_child.')
flags.DEFINE_string(
    'sample_name_parent1', None,
    'Parent1 Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads_parent1.')
flags.DEFINE_string(
    'sample_name_parent2', None,
    'Parent2 Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads_parent2.')
flags.DEFINE_boolean(
    'use_hp_information', None,
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
flags.DEFINE_string(
    'output_gvcf_child', None,
    'Optional. Path where we should write gVCF file for child sample.')
flags.DEFINE_string(
    'output_gvcf_parent1', None,
    'Optional. Path where we should write gVCF file for parent1 sample.')
flags.DEFINE_string(
    'output_gvcf_parent2', None,
    'Optional. Path where we should write gVCF file for parent2 sample.')
flags.DEFINE_string(
    'output_gvcf_merged', None,
    'Optional. Path where we should write gVCF with all samples merged.')
flags.DEFINE_boolean(
    'vcf_stats_report', True, 'Optional. Output a visual report (HTML) of '
    'statistics about the output VCF.')

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

CHILD = 'child'
PARENT1 = 'parent1'
PARENT2 = 'parent2'

CALL_VARIANTS_OUTPUT_COMMON_SUFFIX = 'tfrecord.gz'
NO_VARIANT_TFRECORD_SUFFIX = 'tfrecord.gz'
EXAMPLES_NAME_PATTERN = '{}_{}.{}'
CALL_VARIANTS_OUTPUT_PATTERN = '{}_{}.{}'
NO_VARIANT_TFRECORD_PATTERN = '{}_{}.{}'


def call_variants_output_common_prefix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'call_variants_output')


def examples_common_suffix(num_shards):
  return 'tfrecord@{}.gz'.format(num_shards)


def examples_common_prefix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'make_examples')


def nonvariant_site_tfrecord_common_suffix(intermediate_results_dir):
  return os.path.join(intermediate_results_dir, 'gvcf')


def examples_common_name(intermediate_results_dir, num_shards):
  return '{}.{}'.format(
      examples_common_prefix(intermediate_results_dir),
      examples_common_suffix(num_shards))


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


def _update_kwargs_with_warning(kwargs, extra_args, conflict_args=None):
  """Updates `kwargs` with `extra_args`; crashes if `conflict_args` changed."""
  for k, v in extra_args.items():
    if k in kwargs:
      if conflict_args is not None and k in conflict_args and kwargs[k] != v:
        raise ValueError(
            'The extra_args "{}" might have conflicts with other flags. '
            'See '
            'https://github.com/google/deepvariant/blob/r1.4/docs/'
            'deepvariant-pacbio-model-case-study.md#clarification-'
            'of-the---use_hp_information-flag '
            'for an explanation, or report this issue on '
            'https://github.com/google/deepvariant/issues.'.format(k))
      if kwargs[k] != v:
        print('\nWarning: --{} is previously set to {}, now to {}.'.format(
            k, kwargs[k], v))
    kwargs[k] = v
  return kwargs


def make_examples_command(ref, reads_child, reads_parent1, reads_parent2,
                          examples, sample_name_child, sample_name_parent1,
                          sample_name_parent2, runtime_by_region_path,
                          extra_args, **kwargs):
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
    extra_args: Comma-separated list of flag_name=flag_value.
    **kwargs: Additional arguments to pass in for make_examples.

  Returns:
    (string) A command to run.
  """
  command = [
      'time', 'seq 0 {} |'.format(FLAGS.num_shards - 1),
      'parallel -q --halt 2 --line-buffer',
      '/opt/deepvariant/bin/deeptrio/make_examples'
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  if FLAGS.reads_parent1 is not None:
    command.extend(['--reads_parent1', '"{}"'.format(reads_parent1)])
  if FLAGS.reads_parent2 is not None:
    command.extend(['--reads_parent2', '"{}"'.format(reads_parent2)])
  command.extend(['--reads', '"{}"'.format(reads_child)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--sample_name', '"{}"'.format(sample_name_child)])
  if FLAGS.sample_name_parent1 is not None:
    command.extend(
        ['--sample_name_parent1', '"{}"'.format(sample_name_parent1)])
  if FLAGS.sample_name_parent2 is not None:
    command.extend(
        ['--sample_name_parent2', '"{}"'.format(sample_name_parent2)])
  special_args = {}
  special_args['pileup_image_height_child'] = DEEP_TRIO_WGS_PILEUP_HEIGHT_CHILD
  special_args[
      'pileup_image_height_parent'] = DEEP_TRIO_WGS_PILEUP_HEIGHT_PARENT

  if runtime_by_region_path is not None:
    command.extend(
        ['--runtime_by_region', '"{}"'.format(runtime_by_region_path)])

  conflict_args = None
  if FLAGS.model_type == 'PACBIO':
    special_args['pileup_image_width'] = 199
    special_args['realign_reads'] = False
    special_args['vsc_min_fraction_indels'] = 0.12
    special_args['alt_aligned_pileup'] = 'diff_channels'
    special_args['add_hp_channel'] = True
    special_args['sort_by_haplotypes'] = special_args[
        'parse_sam_aux_fields'] = bool(FLAGS.use_hp_information)
    kwargs = _update_kwargs_with_warning(kwargs, special_args)
    conflict_args = ['sort_by_haplotypes', 'parse_sam_aux_fields']

  if FLAGS.model_type == 'WES':
    special_args[
        'pileup_image_height_child'] = DEEP_TRIO_WES_PILEUP_HEIGHT_CHILD
    special_args[
        'pileup_image_height_parent'] = DEEP_TRIO_WES_PILEUP_HEIGHT_PARENT
    special_args['channels'] = 'insert_size'

  if FLAGS.model_type == 'WGS':
    special_args['channels'] = 'insert_size'

  kwargs = _update_kwargs_with_warning(kwargs, special_args)

  # Extend the command with all items in kwargs and extra_args.
  kwargs = _update_kwargs_with_warning(kwargs, _extra_args_to_dict(extra_args),
                                       conflict_args)
  command = _extend_command_by_args_dict(command, kwargs)

  command.extend(['--task {}'])

  if FLAGS.logging_dir:
    command.extend(
        ['2>&1 | tee {}/make_examples.log'.format(FLAGS.logging_dir)])

  return ' '.join(command)


def call_variants_command(outfile, examples, model_ckpt, sample, extra_args):
  """Returns a call_variants command for subprocess.check_call."""
  command = ['time', '/opt/deepvariant/bin/call_variants']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(command,
                                         _extra_args_to_dict(extra_args))
  if FLAGS.logging_dir:
    command.extend([
        '2>&1 | tee {}/call_variants_{}.log'.format(FLAGS.logging_dir, sample)
    ])

  return ' '.join(command)


def postprocess_variants_command(ref,
                                 infile,
                                 outfile,
                                 sample,
                                 extra_args,
                                 nonvariant_site_tfrecord_path=None,
                                 gvcf_outfile=None,
                                 vcf_stats_report=True):
  """Returns a postprocess_variants command for subprocess.check_call."""
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
  # Extend the command with all items in extra_args.
  command = _extend_command_by_args_dict(command,
                                         _extra_args_to_dict(extra_args))
  if FLAGS.logging_dir:
    command.extend([
        '2>&1 | tee {}/postprocess_variants_{}.log'.format(
            FLAGS.logging_dir, sample)
    ])

  return ' '.join(command)


def runtime_by_region_vis_command(runtime_by_region_path: str):
  """Returns a runtime_by_region_vis command for subprocess."""
  runtime_report = os.path.join(FLAGS.logging_dir,
                                'make_examples_runtime_by_region_report.html')

  command = ['time', '/opt/deepvariant/bin/runtime_by_region_vis']
  command.extend(['--input', '"{}"'.format(runtime_by_region_path)])
  command.extend(['--title', '"{}"'.format('DeepTrio')])
  command.extend(['--output', '"{}"'.format(runtime_report)])

  return ' '.join(command)


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
  if FLAGS.customized_model_parent is not None:
    if (not tf.compat.v1.gfile.Exists(FLAGS.customized_model_parent +
                                      '.data-00000-of-00001') or
        not tf.compat.v1.gfile.Exists(FLAGS.customized_model_parent + '.index')
        or
        not tf.compat.v1.gfile.Exists(FLAGS.customized_model_parent + '.meta')):
      raise RuntimeError('The model files {}* do not exist. Potentially '
                         'relevant issue: '
                         'https://github.com/google/deepvariant/blob/r1.4/docs/'
                         'FAQ.md#why-cant-it-find-one-of-the-input-files-eg-'
                         'could-not-open'.format(FLAGS.customized_model_parent))
    logging.info(
        'You set --customized_model_parent. Instead of using the default '
        'model for %s, `call_variants` step will load %s* '
        'instead.', FLAGS.model_type, FLAGS.customized_model_parent)

  if FLAGS.customized_model_child is not None:
    if (not tf.compat.v1.gfile.Exists(FLAGS.customized_model_child +
                                      '.data-00000-of-00001') or
        not tf.compat.v1.gfile.Exists(FLAGS.customized_model_child + '.index')
        or
        not tf.compat.v1.gfile.Exists(FLAGS.customized_model_child + '.meta')):
      raise RuntimeError('The model files {}* do not exist. Potentially '
                         'relevant issue: '
                         'https://github.com/google/deepvariant/blob/r1.4/docs/'
                         'FAQ.md#why-cant-it-find-one-of-the-input-files-eg-'
                         'could-not-open'.format(FLAGS.customized_model_child))
    logging.info(
        'You set --customized_model_child. Instead of using the default '
        'model for %s, `call_variants` step will load %s* '
        'instead.', FLAGS.model_type, FLAGS.customized_model_child)

  if FLAGS.use_hp_information and FLAGS.model_type != 'PACBIO':
    raise ValueError('--use_hp_information can only be used with '
                     '--model_type="PACBIO"')


def get_model_ckpt(model_type, customized_model):
  """Return the path to the model checkpoint based on the input args."""
  if customized_model is not None:
    return customized_model
  else:
    return MODEL_TYPE_MAP[model_type]


def generate_call_variants_command(sample, model_ckpt,
                                   intermediate_results_dir):
  """Helper function to generate call_variants command line."""
  return call_variants_command(
      CALL_VARIANTS_OUTPUT_PATTERN.format(
          call_variants_output_common_prefix(intermediate_results_dir), sample,
          CALL_VARIANTS_OUTPUT_COMMON_SUFFIX),
      EXAMPLES_NAME_PATTERN.format(
          examples_common_prefix(intermediate_results_dir), sample,
          examples_common_suffix(FLAGS.num_shards)), model_ckpt, sample,
      FLAGS.call_variants_extra_args)


def generate_postprocess_variants_command(sample, intermediate_results_dir,
                                          output_vcf, output_gvcf):
  """Helper function to generate post_process command line."""
  return postprocess_variants_command(
      FLAGS.ref,
      CALL_VARIANTS_OUTPUT_PATTERN.format(
          call_variants_output_common_prefix(intermediate_results_dir), sample,
          CALL_VARIANTS_OUTPUT_COMMON_SUFFIX),
      output_vcf,
      sample,
      FLAGS.postprocess_variants_extra_args,
      NO_VARIANT_TFRECORD_PATTERN.format(
          nonvariant_site_tfrecord_common_suffix(intermediate_results_dir),
          sample, examples_common_suffix(FLAGS.num_shards)),
      gvcf_outfile=output_gvcf,
      vcf_stats_report=FLAGS.vcf_stats_report)


def create_all_commands(intermediate_results_dir):
  """Creates 3 commands to be executed later."""
  check_flags()
  commands = []
  post_process_commands = []

  # make_examples
  nonvariant_site_tfrecord_path = None
  if FLAGS.output_gvcf_child is not None:
    nonvariant_site_tfrecord_path = '{}.{}'.format(
        nonvariant_site_tfrecord_common_suffix(intermediate_results_dir),
        examples_common_suffix(FLAGS.num_shards))

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
          FLAGS.ref,
          FLAGS.reads_child,
          FLAGS.reads_parent1,
          FLAGS.reads_parent2,
          examples_common_name(intermediate_results_dir, FLAGS.num_shards),
          FLAGS.sample_name_child,
          FLAGS.sample_name_parent1,
          FLAGS.sample_name_parent2,
          runtime_by_region_path=runtime_by_region_path,
          extra_args=FLAGS.make_examples_extra_args,
          # kwargs:
          gvcf=nonvariant_site_tfrecord_path,
          regions=FLAGS.regions))

  # Calling variants for child sample
  model_ckpt = get_model_ckpt(FLAGS.model_type + '_child',
                              FLAGS.customized_model_child)
  commands.append(
      generate_call_variants_command(CHILD, model_ckpt,
                                     intermediate_results_dir))

  # Calling variants for parent1 sample
  model_ckpt = get_model_ckpt(FLAGS.model_type + '_parent',
                              FLAGS.customized_model_parent)
  if FLAGS.reads_parent1 is not None:
    commands.append(
        generate_call_variants_command(PARENT1, model_ckpt,
                                       intermediate_results_dir))
  if FLAGS.reads_parent2 is not None:
    commands.append(
        generate_call_variants_command(PARENT2, model_ckpt,
                                       intermediate_results_dir))

  # postprocess_variants for child
  post_process_commands.append(
      generate_postprocess_variants_command(CHILD, intermediate_results_dir,
                                            FLAGS.output_vcf_child,
                                            FLAGS.output_gvcf_child))

  if FLAGS.reads_parent1 is not None:
    post_process_commands.append(
        generate_postprocess_variants_command(PARENT1, intermediate_results_dir,
                                              FLAGS.output_vcf_parent1,
                                              FLAGS.output_gvcf_parent1))

  if FLAGS.reads_parent2 is not None:
    post_process_commands.append(
        generate_postprocess_variants_command(PARENT2, intermediate_results_dir,
                                              FLAGS.output_vcf_parent2,
                                              FLAGS.output_gvcf_parent2))

  # runtime-by-region
  if FLAGS.logging_dir and FLAGS.runtime_report:
    commands.append(runtime_by_region_vis_command(runtime_by_region_path))

  return commands, post_process_commands


def main(_):
  if FLAGS.version:
    print('DeepTrio version {}'.format(DEEP_TRIO_VERSION))
    return

  for flag_key in [
      'model_type', 'ref', 'reads_child', 'output_vcf_child',
      'sample_name_child'
  ]:
    if FLAGS.get_flag_value(flag_key, None) is None:
      sys.stderr.write('--{} is required.\n'.format(flag_key))
      sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
      sys.exit(1)

  # Check flags consistency.
  # --reads_parent?, --output_vcf_parent?, --sample_name_parent? flags should
  # either all be set or all be unset.
  parent1_flags = [
      FLAGS.reads_parent1, FLAGS.output_vcf_parent1, FLAGS.sample_name_parent1
  ]
  if any(parent1_flags) and not all(parent1_flags):
    sys.stderr.write(
        '--reads_parent1, --output_vcf_parent1, --sample_name_parent1 must be set altogether\n'
    )
    sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
    sys.exit(1)

  parent2_flags = [
      FLAGS.reads_parent2, FLAGS.output_vcf_parent2, FLAGS.sample_name_parent2
  ]
  if any(parent2_flags) and not all(parent2_flags):
    sys.stderr.write(
        '--reads_parent2, --output_vcf_parent2, --sample_name_parent2 must be set altogether\n'
    )
    sys.stderr.write('Pass --helpshort or --helpfull to see help on flags.\n')
    sys.exit(1)

  intermediate_results_dir = check_or_create_intermediate_results_dir(
      FLAGS.intermediate_results_dir)

  commands, post_process_commands = create_all_commands(
      intermediate_results_dir)
  print('\n***** Intermediate results will be written to {} '
        'in docker. ****\n'.format(intermediate_results_dir))
  for command in commands:
    print('\n***** Running the command:*****\n{}\n'.format(command))
    if not FLAGS.dry_run:
      try:
        subprocess.check_call(command, shell=True, executable='/bin/bash')
      except subprocess.CalledProcessError as e:
        logging.info(e.output)
        raise
  for command in post_process_commands:
    print('\n***** Running the command:*****\n{}\n'.format(command))
  if not FLAGS.dry_run:
    return_codes = [
        subprocess.Popen(command, shell=True, executable='/bin/bash')
        for command in post_process_commands
    ]
    if any([x.wait() > 0 for x in return_codes]):
      raise Exception('One or more postprocess_variants failed.')


if __name__ == '__main__':
  app.run(main)
