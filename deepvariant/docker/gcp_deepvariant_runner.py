# Copyright 2017 Google Inc.
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
r"""Runs the DeepVariant pipeline using the Google Genomics Pipelines API.

To run this script, you also need the pipelines tool in your $PATH.  You can
install it using:

$ go get github.com/googlegenomics/pipelines-tools/...

Sample run command (please run 'python gcp_deepvariant_runner.py --help' for
details on all available options):
$ python gcp_deepvariant_runner.py \
   --project alphanumeric_project_id \
   --zones 'us-*' \
   --docker_image gcr.io/path_to_deepvariant_cpu_docker_image \
   --outfile gs://bucket/output.vcf \
   --staging gs://bucket/staging \
   --model gs://path_to_deepvariant_model_folder \
   --bam gs://path_to_bam_file.bam \
   --ref gs://path_to_fasta_file.fasta
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import logging
import multiprocessing
import os
import subprocess

from subprocess import PIPE

_MAKE_EXAMPLES_JOB_NAME = 'make_examples'
_CALL_VARIANTS_JOB_NAME = 'call_variants'
_POSTPROCESS_VARIANTS_JOB_NAME = 'postprocess_variants'
_DEFAULT_BOOT_DISK_SIZE_GB = '50'

_MAKE_EXAMPLES_COMMAND = r"""
seq "${{SHARD_START_INDEX}}" "${{SHARD_END_INDEX}}" | parallel --halt 2
  /opt/deepvariant/bin/make_examples
    --mode calling
    --examples "${{EXAMPLES}}"/examples_output.tfrecord@"${{SHARDS}}".gz
    --reads "${{INPUT_BAM}}"
    --ref "${{INPUT_REF}}"
    --task {{}}
    {EXTRA_ARGS}
"""

_CALL_VARIANTS_COMMAND = r"""
seq -f "%05g" "${{SHARD_START_INDEX}}" "${{SHARD_END_INDEX}}" |
parallel --jobs "${{CONCURRENT_JOBS}}" --halt 2
/opt/deepvariant/bin/call_variants
  --examples "${{EXAMPLES}}"/examples_output.tfrecord-{{}}-of-"$(printf "%05d" "${{SHARDS}}")".gz
  --outfile "${{CALLED_VARIANTS}}"/call_variants_output.tfrecord-{{}}-of-"$(printf "%05d" "${{SHARDS}}")".gz
  --checkpoint "${{MODEL}}"/model.ckpt
"""

_POSTPROCESS_VARIANTS_COMMAND = r"""
/opt/deepvariant/bin/postprocess_variants
    --ref "${{INPUT_REF}}"
    --infile "${{CALLED_VARIANTS}}"/call_variants_output.tfrecord@"${{SHARDS}}".gz
    --outfile "${{OUTFILE}}"
    {EXTRA_ARGS}
"""


def _get_staging_examples_folder(pipeline_args, worker_index):
  """Returns the folder to store examples from make_examples job."""
  path_parts = [pipeline_args.staging, 'examples']
  # If the number of workers for both make_examples and call_variants is the
  # same, then we can optimize localization between make_examples and
  # call_variants by creating a sub-folder for each worker.
  if (pipeline_args.make_examples_workers ==
      pipeline_args.call_variants_workers):
    path_parts.append(str(worker_index))
  return os.path.join(*path_parts)


def _get_staging_gvcf_folder(pipeline_args):
  """Returns the folder to store gVCF TF records from make_examples job."""
  return os.path.join(pipeline_args.staging, 'gvcf')


def _get_staging_called_variants_folder(pipeline_args):
  """Returns the folder to store called variants from call_variants job."""
  return os.path.join(pipeline_args.staging, 'called_variants')


def _get_base_job_args(pipeline_args):
  """Base arguments that are common among all jobs."""
  pvm_attempts = 0
  if pipeline_args.preemptible:
    pvm_attempts = pipeline_args.max_preemptible_tries

  return [
      'pipelines', '--project', pipeline_args.project, 'run', '--attempts',
      str(pipeline_args.max_non_preemptible_tries), '--pvm-attempts',
      str(pvm_attempts), '--boot-disk-size', _DEFAULT_BOOT_DISK_SIZE_GB,
      '--zones'
  ] + pipeline_args.zones


def _run_job(run_args):
  """Runs a job using the pipelines CLI tool.

  Args:
    run_args: A list of arguments (type string) to pass to the pipelines tool.
  Raises:
    RuntimeError: if there was an error running the pipeline.
  """

  process = subprocess.Popen(
      run_args,
      stdin=PIPE,
      stdout=PIPE,
      stderr=PIPE,
      env={'PATH': os.environ['PATH']})

  try:
    _, stderr = process.communicate()
    if process.returncode == 0:
      return
  except KeyboardInterrupt:
    raise RuntimeError('Job cancelled by user')

  logging.error('Job failed with error %s. Job args: %s', stderr, run_args)
  raise RuntimeError('Job failed with error %s' % stderr)


def _run_make_examples(pipeline_args):
  """Runs the make_examples job."""

  def get_extra_args():
    """Optional arguments that are specific to make_examples binary."""
    extra_args = []
    if pipeline_args.gvcf_outfile:
      extra_args.extend(
          ['--gvcf', '"${GVCF}"/gvcf_output.tfrecord@"${SHARDS}".gz'])
    if pipeline_args.gvcf_gq_binsize:
      extra_args.extend(
          ['--gvcf_gq_binsize',
           str(pipeline_args.gvcf_gq_binsize)])
    if pipeline_args.regions:
      extra_args.extend(['--regions', ' '.join(pipeline_args.regions)])
    if pipeline_args.sample_name:
      extra_args.extend(['--sample_name', pipeline_args.sample_name])
    if pipeline_args.hts_block_size:
      extra_args.extend(['--hts_block_size', pipeline_args.hts_block_size])
    return extra_args

  command = _MAKE_EXAMPLES_COMMAND.format(EXTRA_ARGS=' '.join(get_extra_args()))

  machine_type = 'custom-{0}-{1}'.format(
      pipeline_args.make_examples_cores_per_worker,
      pipeline_args.make_examples_ram_per_worker_gb * 1024)

  num_workers = min(pipeline_args.make_examples_workers, pipeline_args.shards)
  shards_per_worker = pipeline_args.shards / num_workers
  threads = multiprocessing.Pool(num_workers)
  results = []
  for i in range(num_workers):
    inputs = [
        'INPUT_BAM=' + pipeline_args.bam,
        'INPUT_BAI=' + pipeline_args.bai,
        'INPUT_REF=' + pipeline_args.ref,
        'INPUT_REF_FAI=' + pipeline_args.ref_fai,
    ]
    if pipeline_args.ref_gzi:
      inputs.extend([pipeline_args.ref_gzi])

    outputs = [
        'EXAMPLES=' + _get_staging_examples_folder(pipeline_args, i) + '/*'
    ]
    if pipeline_args.gvcf_outfile:
      outputs.extend(['GVCF=' + _get_staging_gvcf_folder(pipeline_args) + '/*'])
    run_args = _get_base_job_args(pipeline_args) + [
        '--name', pipeline_args.job_name_prefix + _MAKE_EXAMPLES_JOB_NAME,
        '--image', pipeline_args.docker_image, '--output',
        os.path.join(pipeline_args.logging, _MAKE_EXAMPLES_JOB_NAME,
                     str(i)), '--inputs', ','.join(inputs), '--outputs',
        ','.join(outputs), '--machine-type', machine_type, '--disk-size',
        str(pipeline_args.make_examples_disk_per_worker_gb), '--set',
        'SHARDS=' + str(pipeline_args.shards), '--set',
        'SHARD_START_INDEX=' + str(int(i * shards_per_worker)), '--set',
        'SHARD_END_INDEX=' + str(int((i + 1) * shards_per_worker - 1)),
        '--command', command
    ]
    results.append(threads.apply_async(_run_job, [run_args]))

  _wait_for_results(threads, results)


def _wait_for_results(threads, results):
  threads.close()
  try:
    threads.join()
  except KeyboardInterrupt:
    raise RuntimeError('Cancelled')

  for result in results:
    if result:
      result.get()


def _run_call_variants(pipeline_args):
  """Runs the call_variants job."""

  command = _CALL_VARIANTS_COMMAND.format()

  machine_type = 'custom-{0}-{1}'.format(
      pipeline_args.call_variants_cores_per_worker,
      pipeline_args.call_variants_ram_per_worker_gb * 1024)

  num_workers = min(pipeline_args.call_variants_workers, pipeline_args.shards)
  shards_per_worker = pipeline_args.shards / num_workers
  threads = multiprocessing.Pool(processes=num_workers)
  results = []
  for i in range(num_workers):
    inputs = [
        'EXAMPLES=' + _get_staging_examples_folder(pipeline_args, i) + '/*'
    ]
    outputs = [
        'CALLED_VARIANTS=' + _get_staging_called_variants_folder(pipeline_args)
        + '/*'
    ]
    run_args = _get_base_job_args(pipeline_args) + [
        '--name', pipeline_args.job_name_prefix + _CALL_VARIANTS_JOB_NAME,
        '--output',
        os.path.join(pipeline_args.logging, _CALL_VARIANTS_JOB_NAME,
                     str(i)), '--image',
        (pipeline_args.docker_image_gpu if pipeline_args.gpu else
         pipeline_args.docker_image), '--inputs', ','.join(inputs), '--outputs',
        ','.join(outputs), '--machine-type', machine_type, '--disk-size',
        str(pipeline_args.call_variants_disk_per_worker_gb), '--set', 'MODEL=' +
        pipeline_args.model, '--set', 'SHARDS=' + str(pipeline_args.shards),
        '--set', 'SHARD_START_INDEX=' + str(int(i * shards_per_worker)),
        '--set', 'SHARD_END_INDEX=' + str(int((i + 1) * shards_per_worker - 1)),
        '--set', 'CONCURRENT_JOBS=' + ('1' if pipeline_args.gpu else (str(
            int(pipeline_args.call_variants_cores_per_worker /
                pipeline_args.call_variants_cores_per_shard)))), '--command',
        command
    ]
    if pipeline_args.gpu:
      run_args.extend(
          ['--gpu-type', pipeline_args.accelerator_type, '--gpus', '1'])
    results.append(threads.apply_async(_run_job, [run_args]))

  _wait_for_results(threads, results)


def _run_postprocess_variants(pipeline_args):
  """Runs the postprocess_variants job."""

  def get_extra_args():
    """Optional arguments that are specific to postprocess_variants binary."""
    extra_args = []
    if pipeline_args.gvcf_outfile:
      extra_args.extend([
          '--nonvariant_site_tfrecord_path',
          '"${GVCF}"/gvcf_output.tfrecord@"${SHARDS}".gz',
          '--gvcf_outfile',
          '"${GVCF_OUTFILE}"',
      ])
    return extra_args

  machine_type = 'custom-{0}-{1}'.format(
      pipeline_args.postprocess_variants_cores,
      pipeline_args.postprocess_variants_ram_gb * 1024)

  inputs = [
      'CALLED_VARIANTS=' + _get_staging_called_variants_folder(pipeline_args) +
      '/*',
      'INPUT_REF=' + pipeline_args.ref,
      'INPUT_REF_FAI=' + pipeline_args.ref_fai,
  ]
  outputs = ['OUTFILE=' + pipeline_args.outfile]

  if pipeline_args.ref_gzi:
    inputs.extend([pipeline_args.ref_gzi])

  if pipeline_args.gvcf_outfile:
    inputs.extend(['GVCF=' + _get_staging_gvcf_folder(pipeline_args) + '/*'])
    outputs.extend(['GVCF_OUTFILE=' + pipeline_args.gvcf_outfile])

  run_args = _get_base_job_args(pipeline_args) + [
      '--name', pipeline_args.job_name_prefix + _POSTPROCESS_VARIANTS_JOB_NAME,
      '--output',
      os.path.join(pipeline_args.logging,
                   _POSTPROCESS_VARIANTS_JOB_NAME), '--image',
      pipeline_args.docker_image, '--inputs', ','.join(inputs), '--outputs',
      ','.join(outputs), '--machine-type', machine_type, '--disk-size',
      str(pipeline_args.postprocess_variants_disk_gb), '--set',
      'SHARDS=' + str(pipeline_args.shards), '--command',
      _POSTPROCESS_VARIANTS_COMMAND.format(
          EXTRA_ARGS=' '.join(get_extra_args()))
  ]
  _run_job(run_args)


def _validate_and_complete_args(pipeline_args):
  """Validates pipeline arguments and fills some missing args (if any)."""
  # Basic validation logic. More detailed validation is done by pipelines API.
  if pipeline_args.preemptible and pipeline_args.max_preemptible_tries <= 0:
    raise ValueError('--max_preemptible_tries must be greater than zero.')
  if pipeline_args.max_non_preemptible_tries <= 0:
    raise ValueError('--max_non_preemptible_tries must be greater than zero.')
  if pipeline_args.make_examples_workers <= 0:
    raise ValueError('--make_examples_workers must be greater than zero.')
  if pipeline_args.call_variants_workers <= 0:
    raise ValueError('--call_variants_workers must be greater than zero.')
  if pipeline_args.shards <= 0:
    raise ValueError('--shards must be greater than zero.')
  if pipeline_args.shards % pipeline_args.make_examples_workers != 0:
    raise ValueError('--shards must be divisible by --make_examples_workers')
  if pipeline_args.shards % pipeline_args.call_variants_workers != 0:
    raise ValueError('--shards must be divisible by --call_variants_workers')
  if pipeline_args.gpu and not pipeline_args.docker_image_gpu:
    raise ValueError('--docker_image_gpu must be provided with --gpu')
  if (pipeline_args.call_variants_cores_per_worker <
      pipeline_args.call_variants_cores_per_shard):
    raise ValueError('--call_variants_cores_per_worker must be at least '
                     'as large as --call_variants_cores_per_shard')
  if (pipeline_args.gvcf_gq_binsize is not None and
      not pipeline_args.gvcf_outfile):
    raise ValueError('--gvcf_outfile must be provided with --gvcf_gq_binsize')
  if (pipeline_args.gvcf_gq_binsize is not None and
      pipeline_args.gvcf_gq_binsize < 1):
    raise ValueError('--gvcf_gq_binsize must be greater or equal to 1')

  # Automatically generate default values for missing args (if any).
  if not pipeline_args.logging:
    pipeline_args.logging = os.path.join(pipeline_args.staging, 'logs')
  if not pipeline_args.ref_fai:
    pipeline_args.ref_fai = pipeline_args.ref + '.fai'
  if not pipeline_args.ref_gzi and pipeline_args.ref.endswith('.gz'):
    pipeline_args.ref_gzi = pipeline_args.ref + '.gzi'
  if not pipeline_args.bai:
    pipeline_args.bai = pipeline_args.bam + '.bai'


def run(argv=None):
  """Runs the DeepVariant pipeline."""
  parser = argparse.ArgumentParser()

  # Required args.
  parser.add_argument(
      '--project',
      required=True,
      help='Cloud project ID in which to run the pipeline.')
  parser.add_argument(
      '--docker_image', required=True, help='DeepVariant docker image.')
  parser.add_argument(
      '--zones',
      required=True,
      nargs='+',
      help=('List of Google Compute Engine zones. Wildcard suffixes are '
            'supported, such as "us-central1-*" or "us-*".'))
  parser.add_argument(
      '--outfile',
      required=True,
      help=('Destination path in Google Cloud Storage where the resulting '
            'VCF file will be stored.'))
  parser.add_argument(
      '--staging',
      required=True,
      help=('A folder in Google Cloud Storage to use for storing intermediate '
            'files from the pipeline.'))
  parser.add_argument(
      '--model',
      required=True,
      help=('A folder in Google Cloud Storage that stores the TensorFlow '
            'model to use to evaluate candidate variant calls. It expects '
            'the files to be prefixed with "model.ckpt".'))
  parser.add_argument(
      '--bam',
      required=True,
      help='Path in Google Cloud Storage that stores the BAM file.')
  parser.add_argument(
      '--ref',
      required=True,
      help='Path in Google Cloud Storage that stores the reference file.')

  # Additional input args. These are required for the pipeline run.
  # Reasonable defaults would be chosen if unspecified (the generated paths
  # must map to valid files).
  parser.add_argument(
      '--bai', help='BAM index file. Defaults to --bam + ".bai" suffix.')
  parser.add_argument(
      '--ref_fai', help='FAI index file. Defaults to --ref + ".fai" suffix.')
  parser.add_argument(
      '--ref_gzi',
      help=('GZI index file. Required if --ref is gz. Defaults to '
            '--ref + ".gzi" suffix.'))
  parser.add_argument(
      '--logging',
      help=('A folder in Google Cloud Storage to use for storing logs. '
            'Defaults to --staging + "/logs".'))

  # Optinal make_examples args.
  parser.add_argument(
      '--sample_name',
      help=('By default, make_examples extracts sample_name from input BAM '
            'file. However, for BAM file with missing sample_name, this has to '
            'be manually set.'))
  parser.add_argument(
      '--hts_block_size',
      help=('Sets the htslib block size (in bytes). Zero or negative uses '
            'default htslib setting. Currently only applies to SAM/BAM '
            'reading.'))

  # Optional gVCF args.
  parser.add_argument(
      '--gvcf_outfile',
      help=('Destination path in Google Cloud Storage where the resulting '
            'gVCF file will be stored. This is optional, and gVCF file will '
            'only be generated if this is specified.'))
  parser.add_argument(
      '--gvcf_gq_binsize',
      type=int,
      help=('Bin size in which make_examples job quantizes gVCF genotype '
            'qualities. Larger bin size reduces the number of gVCF records '
            'at a loss of quality granularity.'))

  # Additional optional pipeline parameters.
  parser.add_argument(
      '--regions',
      default=None,
      nargs='+',
      help=('Optional space-separated list of regions to process. Elements can '
            'be region literals (chr20:10-20) or Google Cloud Storage paths '
            'to BED/BEDPE files.'))
  parser.add_argument(
      '--max_non_preemptible_tries',
      type=int,
      default=2,
      help=('Maximum number of times to try running each worker (within a job) '
            'with regular (non-preemptible) VMs. Regular VMs may still crash '
            'unexpectedly, so it may be worth to retry on transient failures. '
            'Note that if max_preemptible_tries is also specified, then '
            'the pipeline would first be run with preemptible VMs, and then '
            'with regular VMs following the value provided here.'))

  # Optional GPU args.
  parser.add_argument(
      '--gpu',
      default=False,
      action='store_true',
      help='Use GPUs for the call_variants step.')
  parser.add_argument(
      '--docker_image_gpu',
      help='DeepVariant docker image for GPUs. Required if --gpu is set.')
  parser.add_argument(
      '--accelerator_type',
      default='nvidia-tesla-k80',
      help=('GPU type defined by Compute Engine. Please see '
            'https://cloud.google.com/compute/docs/gpus/ for supported GPU '
            'types.'))

  # Optional preemptible args.
  parser.add_argument(
      '--preemptible',
      default=False,
      action='store_true',
      help=('Use preemptible VMs for the pipeline.'))
  parser.add_argument(
      '--max_preemptible_tries',
      type=int,
      default=3,
      help=('Maximum number of times to try running each worker (within a job) '
            'with preemptible VMs. Regular VMs will be used (for the '
            'particular shards assigned to that worker) after this many '
            'preemptions.'))

  # Optional pipeline sharding and machine shapes.
  parser.add_argument(
      '--shards',
      type=int,
      default=8,
      help=('Number of shards to use for the entire pipeline. The number of '
            'shards assigned to each worker is set by dividing --shards by '
            'the number of workers for each job.'))
  parser.add_argument(
      '--make_examples_workers',
      type=int,
      default=1,
      help=('Number of workers (machines) to use for running the make_examples '
            'job.'))
  parser.add_argument(
      '--make_examples_cores_per_worker',
      type=int,
      default=8,
      help='Number of cores for each worker in make_examples.')
  parser.add_argument(
      '--make_examples_ram_per_worker_gb',
      default=30,
      type=int,
      help='RAM (in GB) to use for each worker in make_examples.')
  parser.add_argument(
      '--make_examples_disk_per_worker_gb',
      type=int,
      default=50,
      help='Disk (in GB) to use for each worker in make_examples.')
  parser.add_argument(
      '--call_variants_workers',
      type=int,
      default=1,
      help=('Number of workers (machines) to use for running the call_variants '
            'job.'))
  parser.add_argument(
      '--call_variants_cores_per_worker',
      type=int,
      default=8,
      help='Number of cores for each worker in call_variants.')
  parser.add_argument(
      '--call_variants_cores_per_shard',
      type=int,
      default=4,
      help=('Number of cores for each shard in call_variants. Multiple shards '
            'may be assigned to each worker, so this flag effectively limits '
            'the number of parallel processes to run on each worker for '
            'call_variants.'))
  parser.add_argument(
      '--call_variants_ram_per_worker_gb',
      type=int,
      default=30,
      help='RAM (in GB) to use for each worker in call_variants.')
  parser.add_argument(
      '--call_variants_disk_per_worker_gb',
      type=int,
      default=30,
      help='Disk (in GB) to use for each worker in call_variants.')
  parser.add_argument(
      '--postprocess_variants_cores',
      type=int,
      default=8,
      help='Number of cores to use for postprocess_variants.')
  parser.add_argument(
      '--postprocess_variants_ram_gb',
      type=int,
      default=30,
      help='RAM (in GB) to use for postprocess_variants.')
  parser.add_argument(
      '--postprocess_variants_disk_gb',
      type=int,
      default=30,
      help='Disk (in GB) to use for postprocess_variants.')

  # Optional misc args.
  parser.add_argument(
      '--job_name_prefix',
      default='',
      help=('Prefix to add to the name of the jobs. Useful for distinguishing '
            'particular pipeline runs from others (e.g. in billing reports).'))
  parser.add_argument(
      '--jobs_to_run',
      nargs='+',
      default=[
          _MAKE_EXAMPLES_JOB_NAME, _CALL_VARIANTS_JOB_NAME,
          _POSTPROCESS_VARIANTS_JOB_NAME
      ],
      choices=[
          _MAKE_EXAMPLES_JOB_NAME, _CALL_VARIANTS_JOB_NAME,
          _POSTPROCESS_VARIANTS_JOB_NAME
      ],
      help=('DeepVariant jobs to run. The DeepVariant pipeline consists of 3 '
            'jobs. By default, the pipeline runs all 3 jobs (make_examples, '
            'call_variants, postprocess_variants) in sequence. '
            'This option may be used to run parts of the pipeline.'))

  pipeline_args = parser.parse_args(argv)
  _validate_and_complete_args(pipeline_args)

  if _MAKE_EXAMPLES_JOB_NAME in pipeline_args.jobs_to_run:
    logging.info('Running make_examples...')
    _run_make_examples(pipeline_args)
    logging.info('make_examples is done!')
  if _CALL_VARIANTS_JOB_NAME in pipeline_args.jobs_to_run:
    logging.info('Running call_variants...')
    _run_call_variants(pipeline_args)
    logging.info('call_variants is done!')
  if _POSTPROCESS_VARIANTS_JOB_NAME in pipeline_args.jobs_to_run:
    logging.info('Running postprocess_variants...')
    _run_postprocess_variants(pipeline_args)
    logging.info('postprocess_variants is done!')


if __name__ == '__main__':
  logging.basicConfig(
      level=logging.INFO,
      format='[%(asctime)s %(levelname)s %(filename)s] %(message)s',
      datefmt='%m/%d/%Y %H:%M:%S')
  run()
