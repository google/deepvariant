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

To run the tests, first activate virtualenv and install dsub:
$ virtualenv venv
$ . venv/bin/activate
$ pip install dsub

Then run:
$ python gcp_deepvariant_runner_test.py
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import multiprocessing
import unittest

from dsub.lib import dsub_errors
import gcp_deepvariant_runner

import mock
from mock import patch


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
        'gs://output.vcf',
        '--staging',
        'gs://staging',
        '--model',
        'gs://model',
        '--bam',
        'gs://bam',
        '--ref',
        'gs://ref',
    ]

  @patch('dsub.commands.dsub.call')
  @patch.object(multiprocessing, 'Pool')
  def testRunPipeline(self, mock_pool, mock_dsub_call):
    mock_apply_async = mock_pool.return_value.apply_async
    mock_apply_async.return_value = None
    self._argv.extend(
        ['--make_examples_workers', '1', '--call_variants_workers', '1'])
    gcp_deepvariant_runner.run(self._argv)

    mock_apply_async.assert_has_calls([
        mock.call(
            func=mock.ANY,
            args=(mock.HASALLOF(
                'make_examples', 'gcr.io/dockerimage', 'INPUT_BAM=gs://bam',
                'INPUT_BAI=gs://bam.bai', 'INPUT_REF=gs://ref',
                'INPUT_REF_FAI=gs://ref.fai',
                'EXAMPLES=gs://staging/examples/0'), mock.ANY, 0)),
        mock.call(
            func=mock.ANY,
            args=(mock.HASALLOF(
                'call_variants', 'gcr.io/dockerimage', 'MODEL=gs://model',
                'EXAMPLES=gs://staging/examples/0',
                'CALLED_VARIANTS=gs://staging/called_variants'), mock.ANY, 0)),
    ],)
    mock_dsub_call.assert_called_once_with(
        mock.HASALLOF('postprocess_variants', 'gcr.io/dockerimage',
                      'CALLED_VARIANTS=gs://staging/called_variants',
                      'OUTFILE=gs://output.vcf'))

  @patch.object(multiprocessing, 'Pool')
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
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF(
                    'prefix_make_examples', 'gcr.io/dockerimage',
                    'SHARD_START_INDEX=0', 'SHARD_END_INDEX=4',
                    'EXAMPLES=gs://staging/examples'), mock.ANY, 0)),
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('prefix_make_examples',
                                    'gcr.io/dockerimage', 'SHARD_START_INDEX=5',
                                    'SHARD_END_INDEX=9'), mock.ANY, 1)),
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF(
                    'prefix_make_examples', 'gcr.io/dockerimage',
                    'SHARD_START_INDEX=10', 'SHARD_END_INDEX=14'), mock.ANY,
                      2)),
        ],
        any_order=True,
    )

  @patch.object(multiprocessing, 'Pool')
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
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('call_variants', 'gcr.io/dockerimage',
                                    'SHARD_START_INDEX=0', 'SHARD_END_INDEX=4',
                                    'CONCURRENT_JOBS=2'), mock.ANY, 0)),
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('call_variants', 'gcr.io/dockerimage',
                                    'SHARD_START_INDEX=5', 'SHARD_END_INDEX=9',
                                    'CONCURRENT_JOBS=2'), mock.ANY, 1)),
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('call_variants', 'gcr.io/dockerimage',
                                    'SHARD_START_INDEX=10',
                                    'SHARD_END_INDEX=14', 'CONCURRENT_JOBS=2'),
                      mock.ANY, 2)),
        ],
        any_order=True,
    )

  @patch.object(multiprocessing, 'Pool')
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
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('call_variants', 'gcr.io/dockerimage_gpu',
                                    'nvidia-tesla-k80', 'SHARD_START_INDEX=0',
                                    'SHARD_END_INDEX=4', 'CONCURRENT_JOBS=1'),
                      mock.ANY, 0)),
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('call_variants', 'gcr.io/dockerimage_gpu',
                                    'nvidia-tesla-k80', 'SHARD_START_INDEX=5',
                                    'SHARD_END_INDEX=9', 'CONCURRENT_JOBS=1'),
                      mock.ANY, 1)),
            mock.call(
                func=mock.ANY,
                args=(mock.HASALLOF('call_variants', 'gcr.io/dockerimage_gpu',
                                    'nvidia-tesla-k80', 'SHARD_START_INDEX=10',
                                    'SHARD_END_INDEX=14', 'CONCURRENT_JOBS=1'),
                      mock.ANY, 2)),
        ],
        any_order=True,
    )

  @patch('dsub.commands.dsub.call')
  def testRunPostProcessVariants(self, mock_dsub_call):
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
    mock_dsub_call.assert_called_once_with(
        mock.HASALLOF('postprocess_variants', 'gcr.io/dockerimage',
                      'CALLED_VARIANTS=gs://staging/called_variants',
                      'INPUT_REF=gs://ref', 'INPUT_REF_FAI=gs://ref.fai',
                      'OUTFILE=gs://output.vcf'))

  @patch('dsub.commands.dsub.call')
  def testRunWithPreemptibles(self, mock_dsub_call):
    self._argv.extend([
        '--jobs_to_run',
        'postprocess_variants',
        '--preemptible',
        '--max_preemptible_tries',
        '2',
    ])
    preemptible_exception = dsub_errors.JobExecutionError(
        'error', ['Error in job - code 10: 14: VM stopped unexpectedly'])
    mock_dsub_call.side_effect = [
        preemptible_exception, preemptible_exception, None
    ]
    gcp_deepvariant_runner.run(self._argv)

    mock_dsub_call.assert_has_calls([
        mock.call(mock.HASALLOF('postprocess_variants', '--preemptible')),
        mock.call(mock.HASALLOF('postprocess_variants', '--preemptible')),
        mock.call(mock.HASALLOF('postprocess_variants'))
    ])

  @patch('dsub.commands.dsub.call')
  def testRunWithPreemptibles_NonPreemptibleFailure(self, mock_dsub_call):
    self._argv.extend([
        '--jobs_to_run',
        'postprocess_variants',
        '--preemptible',
        '--max_preemptible_tries',
        '2',
    ])
    preemptible_exception = dsub_errors.JobExecutionError(
        'error', ['Error in job - code 10: 13: VM stopped unexpectedly'])
    # Note that the error code changed to '130' from '13'.
    non_preemptible_exception = dsub_errors.JobExecutionError(
        'error', ['Error in job - code 10: 130: VM stopped unexpectedly'])
    mock_dsub_call.side_effect = [
        preemptible_exception, non_preemptible_exception
    ]

    try:
      gcp_deepvariant_runner.run(self._argv)
      self.fail('Non-preemptible failures should throw exception.')
    except RuntimeError:
      pass

    # Two preemptible tries should have still happened.
    mock_dsub_call.assert_has_calls([
        mock.call(mock.HASALLOF('postprocess_variants', '--preemptible')),
        mock.call(mock.HASALLOF('postprocess_variants', '--preemptible'))
    ])

  @patch('dsub.commands.dsub.call')
  def testRunUnexpectedStopFromRegularVm(self, mock_dsub_call):
    self._argv.extend([
        '--jobs_to_run',
        'postprocess_variants',
        '--max_non_preemptible_tries',
        '3',
    ])
    unexpected_stop_exception = dsub_errors.JobExecutionError(
        'error', ['Error in job - code 10: 13: VM stopped unexpectedly'])
    mock_dsub_call.side_effect = [
        unexpected_stop_exception, unexpected_stop_exception, None
    ]
    gcp_deepvariant_runner.run(self._argv)

    mock_dsub_call.assert_has_calls([
        mock.call(mock.HASALLOF('postprocess_variants')),
        mock.call(mock.HASALLOF('postprocess_variants')),
        mock.call(mock.HASALLOF('postprocess_variants'))
    ])

  @patch('dsub.commands.dsub.call')
  def testRunUnexpectedStopFromRegularVm_WithPreemptibles(self, mock_dsub_call):
    self._argv.extend([
        '--jobs_to_run',
        'postprocess_variants',
        '--preemptible',
        '--max_preemptible_tries',
        '1',
        '--max_non_preemptible_tries',
        '2',
    ])
    preemptible_exception = dsub_errors.JobExecutionError(
        'error', ['Error in job - code 10: 14: VM stopped unexpectedly'])
    unexpected_stop_exception = dsub_errors.JobExecutionError(
        'error', ['Error in job - code 10: 13: VM stopped unexpectedly'])
    mock_dsub_call.side_effect = [
        preemptible_exception, unexpected_stop_exception,
        unexpected_stop_exception
    ]
    try:
      gcp_deepvariant_runner.run(self._argv)
      self.fail(
          'Too many unexpected stops from regular VMs should throw exception.')
    except RuntimeError:
      pass

    mock_dsub_call.assert_has_calls([
        mock.call(mock.HASALLOF('postprocess_variants', '--preemptible')),
        mock.call(mock.HASALLOF('postprocess_variants')),
        mock.call(mock.HASALLOF('postprocess_variants'))
    ])


if __name__ == '__main__':
  unittest.main()
