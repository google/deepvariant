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

FLAGS = flags.FLAGS

flags.DEFINE_enum('model_type', None, ['WGS', 'WES', 'PacBio'],
                  'Required. Type of model to use for variant calling.')
flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads.')
flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref.')
flags.DEFINE_string(
    'regions', None,
    'Optional. Space-separated list of regions we want to process. Elements '
    'can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.')
flags.DEFINE_string('output_vcf', None,
                    'Required. Path where we should write VCF file.')
flags.DEFINE_string('output_gvcf', None,
                    'Optional. Path where we should write gVCF file.')
flags.DEFINE_string(
    'intermediate_results_dir', '/opt/tmp_output',
    'Optional. If specified, this should be an existing '
    'directory that is visible insider docker, and will be '
    'used to to store intermediate outputs.')
flags.DEFINE_integer('num_shards', 1,
                     'Optional. Number of shards for make_examples step.')

flags.mark_flag_as_required('model_type')
flags.mark_flag_as_required('ref')
flags.mark_flag_as_required('reads')
flags.mark_flag_as_required('output_vcf')

MODEL_TYPE_MAP = {
    'WGS': '/opt/models/wgs/model.ckpt',
    'WES': '/opt/models/wes/model.ckpt',
}


def make_examples_command(ref, reads, examples, gvcf=None, regions=None):
  """Returns a make_examples command for subprocess.check_call.

  Args:
    ref: Input FASTA file.
    reads: Input BAM file.
    examples: Output tfrecord file containing tensorflow.Example files.
    gvcf: (Optional) Output tfrecord file containing tensorflow.Example files.
    regions: (Optional) Input BED file or chromosome ranges.

  Returns:
    (string) A command to run.
  """
  command = [
      'seq 0 {} |'.format(FLAGS.num_shards - 1), 'parallel -k --line-buffer',
      '/opt/deepvariant/bin/make_examples'
  ]
  command.extend(['--mode', 'calling'])
  command.extend(['--ref', '"{}"'.format(ref)])
  command.extend(['--reads', '"{}"'.format(reads)])
  command.extend(['--examples', '"{}"'.format(examples)])
  if regions is not None:
    command.extend(['--regions', '"{}"'.format(regions)])
  if gvcf is not None:
    command.extend(['--gvcf', '"{}"'.format(gvcf)])
  command.extend(['--task {}'])
  return ' '.join(command)


def call_variants_command(outfile, examples, model_type):
  """Returns a call_variants command for subprocess.check_call."""
  checkpoint = MODEL_TYPE_MAP[model_type]
  command = ['/opt/deepvariant/bin/call_variants']
  command.extend(['--outfile', '"{}"'.format(outfile)])
  command.extend(['--examples', '"{}"'.format(examples)])
  command.extend(['--checkpoint', '"{}"'.format(checkpoint)])
  return ' '.join(command)


def postprocess_variants_command(ref,
                                 infile,
                                 outfile,
                                 nonvariant_site_tfrecord_path=None,
                                 gvcf_outfile=None):
  """Returns a postprocess_variants command for subprocess.check_call."""
  command = ['/opt/deepvariant/bin/postprocess_variants']
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
  return ' '.join(command)


def check_or_create_intermediate_results_dir(intermediate_results_dir):
  """Checks or creates the path to the directory for intermediate results."""
  if not os.path.isdir(intermediate_results_dir):
    os.makedirs(intermediate_results_dir)


def main(_):
  check_or_create_intermediate_results_dir(FLAGS.intermediate_results_dir)

  nonvariant_site_tfrecord_path = None
  if FLAGS.output_gvcf is not None:
    nonvariant_site_tfrecord_path = os.path.join(
        FLAGS.intermediate_results_dir,
        'gvcf.tfrecord@{}.gz'.format(FLAGS.num_shards))

  examples = os.path.join(
      FLAGS.intermediate_results_dir,
      'make_examples.tfrecord@{}.gz'.format(FLAGS.num_shards))
  command = make_examples_command(
      FLAGS.ref,
      FLAGS.reads,
      examples,
      gvcf=nonvariant_site_tfrecord_path,
      regions=FLAGS.regions)
  print('\n***** Running the command:*****\n{}\n'.format(command))
  subprocess.check_call(command, shell=True)

  call_variants_output = os.path.join(FLAGS.intermediate_results_dir,
                                      'call_variants_output.tfrecord.gz')
  command = call_variants_command(call_variants_output, examples,
                                  FLAGS.model_type)
  print('\n***** Running the command:*****\n{}\n'.format(command))
  subprocess.check_call(command, shell=True)

  command = postprocess_variants_command(
      FLAGS.ref,
      call_variants_output,
      FLAGS.output_vcf,
      nonvariant_site_tfrecord_path=nonvariant_site_tfrecord_path,
      gvcf_outfile=FLAGS.output_gvcf)
  print('\n***** Running the command:*****\n{}\n'.format(command))
  subprocess.check_call(command, shell=True)


if __name__ == '__main__':
  app.run(main)
