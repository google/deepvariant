# Copyright 2019 Google LLC.
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
https://github.com/google/deepvariant/blob/r0.8/docs/deepvariant-quick-start.md
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import subprocess

from absl import app
from absl import flags
from absl import logging

FLAGS = flags.FLAGS

# Required flags.
flags.DEFINE_enum(
    'model_type', None, ['WGS', 'WES', 'PACBIO'],
    'Required. Type of model to use for variant calling. Each '
    'model_type has an associated default model, which can be '
    'overridden by the --customized_model flag.')
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
flags.DEFINE_string(
    'intermediate_results_dir', '/tmp/deepvariant_tmp_output',
    'Optional. If specified, this should be an existing '
    'directory that is visible insider docker, and will be '
    'used to to store intermediate outputs.')
# Optional flags for call_variants.
flags.DEFINE_string(
    'customized_model', None,
    'Optional. A path to a model checkpoint to load for the `call_variants` '
    'step. If not set, the default for each --model_type will be used')
# Optional flags for make_examples.
flags.DEFINE_boolean(
    'use_ref_for_cram', None,
    'If True, use --ref argument as the reference file for CRAM file passed to '
    '--reads. Reference must be on a local POSIX filesystem. Leaving this flag '
    'unspecified is equivalent to setting it to true, which means DeepVariant '
    'will use the --ref file to decode your CRAM file.')
flags.DEFINE_boolean('keep_secondary_alignments', None,
                     'If True, keep reads marked as secondary alignments.')
flags.DEFINE_boolean('keep_supplementary_alignments', None,
                     'If True, keep reads marked as supplementary alignments.')
flags.DEFINE_integer('num_shards', 1,
                     'Optional. Number of shards for make_examples step.')
flags.DEFINE_string(
    'regions', None,
    'Optional. Space-separated list of regions we want to process. Elements '
    'can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.')
flags.DEFINE_string(
    'sample_name', None,
    'Sample name to use instead of the sample name from the input reads BAM '
    '(SM tag in the header).')

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
}


def make_examples_command(ref, reads, examples, **kwargs):
  """Returns a make_examples command for subprocess.check_call.

  Args:
    ref: Input FASTA file.
    reads: Input BAM file.
    examples: Output tfrecord file containing tensorflow.Example files.
    **kwargs: Additional arguments to pass in for make_examples.

  Returns:
    (string) A command to run.
  """
  command = [
      'time', 'seq 0 {} |'.format(FLAGS.num_shards - 1),
      'parallel -k --line-buffer', '/opt/deepvariant/bin/make_examples'
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads', '"{}"'.format(reads)])
  command.extend(['--examples', '"{}"'.format(examples)])
  # Extend the command with all items in kwargs.
  for key in sorted(kwargs):
    value = kwargs[key]
    if value is None:
      continue
    if isinstance(value, bool):
      added_arg = '' if value else 'no'
      added_arg += key
      command.extend(['--' + added_arg])
    else:
      command.extend(['--' + key, '"{}"'.format(value)])
  command.extend(['--task {}'])
  return ' '.join(command)


def call_variants_command(outfile, examples, model_ckpt):
  """Returns a call_variants command for subprocess.check_call."""
  command = ['time', '/opt/deepvariant/bin/call_variants']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(model_ckpt)])
  return ' '.join(command)


def postprocess_variants_command(ref,
                                 infile,
                                 outfile,
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
  return ' '.join(command)


def check_or_create_intermediate_results_dir(intermediate_results_dir):
  """Checks or creates the path to the directory for intermediate results."""
  if not os.path.isdir(intermediate_results_dir):
    os.makedirs(intermediate_results_dir)


def check_flags():
  """Additional logic to make sure flags are set appropriately."""
  if FLAGS.customized_model is not None:
    logging.info(
        'You set --customized_model. Instead of using the default '
        'model for %s, `call_variants` step will load %s '
        'instead.', FLAGS.model_type, FLAGS.customized_model)


def get_model_ckpt(model_type, customized_model):
  """Return the path to the model checkpoint based on the input args."""
  if customized_model is not None:
    return customized_model
  else:
    return MODEL_TYPE_MAP[model_type]


def create_all_commands():
  """Creates 3 commands to be executed later."""
  commands = []
  # make_examples
  nonvariant_site_tfrecord_path = None
  if FLAGS.output_gvcf is not None:
    nonvariant_site_tfrecord_path = os.path.join(
        FLAGS.intermediate_results_dir,
        'gvcf.tfrecord@{}.gz'.format(FLAGS.num_shards))

  examples = os.path.join(
      FLAGS.intermediate_results_dir,
      'make_examples.tfrecord@{}.gz'.format(FLAGS.num_shards))

  commands.append(
      make_examples_command(
          FLAGS.ref,
          FLAGS.reads,
          examples,
          gvcf=nonvariant_site_tfrecord_path,
          regions=FLAGS.regions,
          realign_reads=False if FLAGS.model_type == 'PACBIO' else None,
          use_ref_for_cram=FLAGS.use_ref_for_cram,
          keep_secondary_alignments=FLAGS.keep_secondary_alignments,
          keep_supplementary_alignments=FLAGS.keep_supplementary_alignments,
          sample_name=FLAGS.sample_name))

  # call_variants
  call_variants_output = os.path.join(FLAGS.intermediate_results_dir,
                                      'call_variants_output.tfrecord.gz')
  model_ckpt = get_model_ckpt(FLAGS.model_type, FLAGS.customized_model)
  commands.append(
      call_variants_command(call_variants_output, examples, model_ckpt))

  # postprocess_variants
  commands.append(
      postprocess_variants_command(
          FLAGS.ref,
          call_variants_output,
          FLAGS.output_vcf,
          nonvariant_site_tfrecord_path=nonvariant_site_tfrecord_path,
          gvcf_outfile=FLAGS.output_gvcf,
          vcf_stats_report=FLAGS.vcf_stats_report))

  return commands


def main(_):
  check_or_create_intermediate_results_dir(FLAGS.intermediate_results_dir)
  check_flags()

  commands = create_all_commands()
  for command in commands:
    print('\n***** Running the command:*****\n{}\n'.format(command))
    subprocess.check_call(command, shell=True, executable='/bin/bash')


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'model_type',
      'ref',
      'reads',
      'output_vcf',
  ])
  app.run(main)
