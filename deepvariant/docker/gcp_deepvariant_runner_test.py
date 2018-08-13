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
"""Tests for gcp_deepvariant_runner.

To run the tests, first activate virtualenv and install mock:
$ virtualenv venv
$ . venv/bin/activate
$ pip install mock

Then run:
$ python gcp_deepvariant_runner_test.py
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import multiprocessing
import unittest

import gcp_deepvariant_runner
import gke_cluster

import mock


class _HasAllOf(object):
  """Helper class used in mock.call to check that all arguments are in a set."""

  def __init__(self, *values):
    self._values = set(values)

  def _expand(self, items):
    expanded = []
    explode = False
    for item in items:
      if explode:
        explode = False
        expanded.extend(item.split(','))
      else:
        expanded.append(item)

      if item == '--inputs' or item == '--outputs':
        explode = True

    return expanded

  def __eq__(self, other):
    return self._values.issubset(set(self._expand(other)))

  def __ne__(self, other):
    return not self.__eq__(other)

  def __repr__(self):
    return '<_HasAllOf({})>'.format(', '.join(repr(v) for v in self._values))


class AnyStringWith(str):
  """Helper class used in mocking to check string arguments."""

  def __eq__(self, other):
    return self in other


class DeepvariantRunnerTest(unittest.TestCase):

  def setUp(self):
    self._argv = [
        '--project',
        'project',
        '--docker_image',
        'gcr.io/dockerimage',
        '--zones',
        'zone-a',
        'zone-b',
        '--outfile',
        'gs://bucket/output.vcf',
        '--staging',
        'gs://bucket/staging',
        '--model',
        'gs://bucket/model',
        '--bam',
        'gs://bucket/bam',
        '--ref',
        'gs://bucket/ref',
    ]

  @mock.patch('gcp_deepvariant_runner._run_job')
  @mock.patch.object(multiprocessing, 'Pool')
  def testRunPipeline(self, mock_pool, mock_run_job):
    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None

    self._argv.extend(
        ['--make_examples_workers', '1', '--call_variants_workers', '1'])
    gcp_deepvariant_runner.run(self._argv)

    mock_apply_async.assert_has_calls([
        mock.call(mock_run_job, [
            _HasAllOf('make_examples', 'gcr.io/dockerimage',
                      'INPUT_BAM=gs://bucket/bam',
                      'INPUT_BAI=gs://bucket/bam.bai',
                      'INPUT_REF=gs://bucket/ref',
                      'INPUT_REF_FAI=gs://bucket/ref.fai',
                      'EXAMPLES=gs://bucket/staging/examples/0/*')
        ]),
        mock.call(mock_run_job, [
            _HasAllOf('call_variants', 'gcr.io/dockerimage',
                      'MODEL=gs://bucket/model',
                      'EXAMPLES=gs://bucket/staging/examples/0/*',
                      'CALLED_VARIANTS=gs://bucket/staging/called_variants/*')
        ]),
    ])
    mock_run_job.assert_called_once_with(
        _HasAllOf('postprocess_variants', 'gcr.io/dockerimage',
                  'CALLED_VARIANTS=gs://bucket/staging/called_variants/*',
                  'OUTFILE=gs://bucket/output.vcf'))

  @mock.patch('gcp_deepvariant_runner._run_job')
  @mock.patch.object(multiprocessing, 'Pool')
  def testRunPipeline_WithGVCFOutFile(self, mock_pool, mock_run_job):

    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None
    self._argv.extend([
        '--make_examples_workers', '1', '--call_variants_workers', '1',
        '--gvcf_outfile', 'gs://bucket/gvcf_output.vcf', '--gvcf_gq_binsize',
        '5'
    ])
    gcp_deepvariant_runner.run(self._argv)

    mock_apply_async.assert_has_calls([
        mock.call(mock_run_job, [
            _HasAllOf('make_examples', 'gcr.io/dockerimage',
                      'INPUT_BAM=gs://bucket/bam',
                      'INPUT_BAI=gs://bucket/bam.bai',
                      'INPUT_REF=gs://bucket/ref',
                      'INPUT_REF_FAI=gs://bucket/ref.fai',
                      'EXAMPLES=gs://bucket/staging/examples/0/*',
                      'GVCF=gs://bucket/staging/gvcf/*')
        ]),
        mock.call(mock_run_job, [
            _HasAllOf('call_variants', 'gcr.io/dockerimage',
                      'MODEL=gs://bucket/model',
                      'EXAMPLES=gs://bucket/staging/examples/0/*',
                      'CALLED_VARIANTS=gs://bucket/staging/called_variants/*')
        ]),
    ],)
    mock_run_job.assert_called_once_with(
        _HasAllOf('postprocess_variants', 'gcr.io/dockerimage',
                  'CALLED_VARIANTS=gs://bucket/staging/called_variants/*',
                  'OUTFILE=gs://bucket/output.vcf',
                  'GVCF=gs://bucket/staging/gvcf/*',
                  'GVCF_OUTFILE=gs://bucket/gvcf_output.vcf'))

  @mock.patch.object(multiprocessing, 'Pool')
  def testRunMakeExamples(self, mock_pool):
    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None
    self._argv.extend([
        '--jobs_to_run',
        'make_examples',
        '--make_examples_workers',
        '3',
        '--shards',
        '15',
        '--gpu',  # GPU should not have any effect.
        '--docker_image_gpu',
        'image_gpu',
        '--job_name_prefix',
        'prefix_',
    ])
    gcp_deepvariant_runner.run(self._argv)

    mock_apply_async.assert_has_calls(
        [
            mock.call(mock.ANY, [
                _HasAllOf('prefix_make_examples', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=0', 'SHARD_END_INDEX=4',
                          'EXAMPLES=gs://bucket/staging/examples/*')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('prefix_make_examples', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=5', 'SHARD_END_INDEX=9')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('prefix_make_examples', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=10', 'SHARD_END_INDEX=14')
            ]),
        ],
        any_order=True,
    )

  @mock.patch.object(multiprocessing, 'Pool')
  def testRunMakeExamples_WithGcsfuse(self, mock_pool):
    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None
    self._argv.extend([
        '--jobs_to_run',
        'make_examples',
        '--make_examples_workers',
        '3',
        '--shards',
        '15',
        '--gpu',  # GPU should not have any effect.
        '--docker_image_gpu',
        'image_gpu',
        '--job_name_prefix',
        'prefix_',
        '--gcsfuse',
    ])
    gcp_deepvariant_runner.run(self._argv)
    mock_apply_async.assert_has_calls(
        [
            mock.call(mock.ANY, [
                _HasAllOf('prefix_make_examples', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=0', 'SHARD_END_INDEX=4',
                          'EXAMPLES=gs://bucket/staging/examples/*',
                          'GCS_BUCKET=bucket', 'BAM=bam')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('prefix_make_examples', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=5', 'SHARD_END_INDEX=9',
                          'GCS_BUCKET=bucket', 'BAM=bam')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('prefix_make_examples', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=10', 'SHARD_END_INDEX=14',
                          'GCS_BUCKET=bucket', 'BAM=bam')
            ]),
        ],
        any_order=True,
    )

  @mock.patch.object(multiprocessing, 'Pool')
  def testRunCallVariants(self, mock_pool):
    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None
    self._argv.extend([
        '--jobs_to_run',
        'call_variants',
        '--call_variants_workers',
        '3',
        '--call_variants_cores_per_worker',
        '5',
        '--call_variants_cores_per_shard',
        '2',
        '--shards',
        '15',
    ])
    gcp_deepvariant_runner.run(self._argv)

    mock_apply_async.assert_has_calls(
        [
            mock.call(mock.ANY, [
                _HasAllOf('call_variants', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=0', 'SHARD_END_INDEX=4',
                          'CONCURRENT_JOBS=2')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('call_variants', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=5', 'SHARD_END_INDEX=9',
                          'CONCURRENT_JOBS=2')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('call_variants', 'gcr.io/dockerimage',
                          'SHARD_START_INDEX=10', 'SHARD_END_INDEX=14',
                          'CONCURRENT_JOBS=2')
            ]),
        ],
        any_order=True,
    )

  @mock.patch.object(multiprocessing, 'Pool')
  def testRunCallVariants_GPU(self, mock_pool):
    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None
    self._argv.extend([
        '--jobs_to_run',
        'call_variants',
        '--call_variants_workers',
        '3',
        '--call_variants_cores_per_worker',
        '5',
        '--call_variants_cores_per_shard',
        '5',
        '--shards',
        '15',
        '--gpu',
        '--docker_image_gpu',
        'gcr.io/dockerimage_gpu',
    ])
    gcp_deepvariant_runner.run(self._argv)

    mock_apply_async.assert_has_calls(
        [
            mock.call(mock.ANY, [
                _HasAllOf('call_variants', 'gcr.io/dockerimage_gpu',
                          'nvidia-tesla-k80', 'SHARD_START_INDEX=0',
                          'SHARD_END_INDEX=4', 'CONCURRENT_JOBS=1')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('call_variants', 'gcr.io/dockerimage_gpu',
                          'nvidia-tesla-k80', 'SHARD_START_INDEX=5',
                          'SHARD_END_INDEX=9', 'CONCURRENT_JOBS=1')
            ]),
            mock.call(mock.ANY, [
                _HasAllOf('call_variants', 'gcr.io/dockerimage_gpu',
                          'nvidia-tesla-k80', 'SHARD_START_INDEX=10',
                          'SHARD_END_INDEX=14', 'CONCURRENT_JOBS=1')
            ]),
        ],
        any_order=True,
    )

  @mock.patch.object(gke_cluster.GkeCluster, '__init__', return_value=None)
  @mock.patch.object(gke_cluster.GkeCluster, 'deploy_pod')
  def testRunCallVariants_TPU(self, mock_deploy_pod, mock_init):
    self._argv.extend([
        '--jobs_to_run',
        'call_variants',
        '--call_variants_workers',
        '1',
        '--shards',
        '15',
        '--tpu',
        '--gke_cluster_name',
        'foo-cluster',
        '--gke_cluster_zone',
        'us-central1-c',
        '--docker_image',
        'gcr.io/dockerimage',
    ])
    gcp_deepvariant_runner.run(self._argv)
    mock_init.assert_called_once_with('foo-cluster', None, 'us-central1-c')
    mock_deploy_pod.assert_called_once_with(
        pod_config=mock.ANY,
        pod_name=AnyStringWith('deepvariant-'),
        retries=1,
        wait=True)

  @mock.patch('gcp_deepvariant_runner._run_job')
  def testRunPostProcessVariants(self, mock_run_job):
    self._argv.extend([
        '--jobs_to_run',
        'postprocess_variants',
        '--shards',
        '15',
        '--gpu',  # GPU should not have any effect.
        '--docker_image_gpu',
        'gcr.io/dockerimage_gpu',
    ])
    gcp_deepvariant_runner.run(self._argv)
    mock_run_job.assert_called_once_with(
        _HasAllOf('postprocess_variants', 'gcr.io/dockerimage',
                  'CALLED_VARIANTS=gs://bucket/staging/called_variants/*',
                  'INPUT_REF=gs://bucket/ref',
                  'INPUT_REF_FAI=gs://bucket/ref.fai',
                  'OUTFILE=gs://bucket/output.vcf'))


if __name__ == '__main__':
  unittest.main()
