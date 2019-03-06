# Copyright 2018 Google LLC.
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
"""Tests for gke_cluster.py.

To run the tests, first activate virtualenv and install required packages:
$ virtualenv venv
$ . venv/bin/activate
$ pip install -r requirements.txt
$ pip install mock


Then run:
$ python gke_cluster_test.py
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import unittest
import gke_cluster
import mock


class GkeClusterTest(unittest.TestCase):
  """Tests for GkeCluster class."""

  @mock.patch('process_util.run_command')
  def test_create_new_cluster(self, mock_call):
    gke_cluster.GkeCluster(
        'foo-cluster', cluster_zone='foo-zone', alpha_cluster=True)
    mock_call.assert_called_with(
        [
            'gcloud', 'alpha', 'container', 'clusters', 'create', 'foo-cluster',
            '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)

  @mock.patch(
      'process_util.run_command',
      side_effect=(KeyboardInterrupt, 'foo-cluster', None))
  @mock.patch(
      'gke_cluster.GkeCluster._cluster_exists', side_effect=(False, True))
  def test_create_new_cluster_with_keyboard_interrupt(
      self, unused_mock_cluster_exists, mock_call):
    gke_cluster.GkeCluster(
        'foo-cluster', cluster_zone='foo-zone', alpha_cluster=True)
    mock_call.assert_any_call(
        [
            'gcloud', 'alpha', 'container', 'clusters', 'create', 'foo-cluster',
            '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'describe', 'foo-cluster',
            '--format=value(status)', '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'describe', 'foo-cluster',
            '--format=value(status)', '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)

  @mock.patch(
      'process_util.run_command', side_effect=('foo-cluster', None, 'RUNNING'))
  def test_reuse_existing_cluster(self, mock_call):
    gke_cluster.GkeCluster('foo-cluster', cluster_zone='foo-zone')
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'list', '--format=value(name)',
            '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'get-credentials', 'foo-cluster',
            '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'describe', 'foo-cluster',
            '--format=value(status)', '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)

  @mock.patch(
      'process_util.run_command',
      side_effect=('foo-cluster', None, 'RUNNING', 'RUNNING'))
  def test_get_cluster_status(self, mock_call):
    self.assertEqual(
        gke_cluster.GkeCluster('foo-cluster',
                               cluster_zone='foo-zone')._get_cluster_status(),
        gke_cluster.ClusterStatus.RUNNING)

    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'list', '--format=value(name)',
            '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'get-credentials', 'foo-cluster',
            '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)
    mock_call.assert_any_call(
        [
            'gcloud', 'container', 'clusters', 'describe', 'foo-cluster',
            '--format=value(status)', '--zone', 'foo-zone'
        ],
        retry_delay_sec=1,
        retries=1)

  @mock.patch('process_util.run_command', return_value='foo')
  def test_get_cluster_unknown_status(self, unused_mock_call):
    self.assertEqual(
        gke_cluster.GkeCluster('foo-cluster',
                               cluster_zone='foo-zone')._get_cluster_status(),
        gke_cluster.ClusterStatus.UNKNOWN)

  @mock.patch('time.sleep')
  @mock.patch('process_util.run_command')
  def test_delete_cluster_in_provisioning_state(self, mock_call,
                                                unused_mock_sleep):
    mock_call.side_effect = [
        'foo-cluster', None, 'RUNNING', 'foo-cluster', 'PROVISIONING',
        'RUNNING', 'RUNNING'
    ]
    gke_cluster.GkeCluster(
        'foo-cluster', cluster_zone='foo-zone').delete_cluster(wait=False)
    self.assertEqual(mock_call.call_count, 7)

  @mock.patch('process_util.run_command', return_value='foo-cluster\n')
  @mock.patch('gke_cluster.GkeCluster._create_cluster')
  def test_cluster_exists(self, unused_mock_call, unused_mock_create_cluster):
    self.assertEqual(
        gke_cluster.GkeCluster('foo-cluster',
                               cluster_zone='foo-zone')._cluster_exists(), True)

  @mock.patch('process_util.run_command', return_value='bar-cluster\n')
  @mock.patch('gke_cluster.GkeCluster._create_cluster')
  def test_cluster_exists_for_non_existent_cluster(
      self, unused_mock_call, unused_unused_mock_create_cluster):
    self.assertEqual(
        gke_cluster.GkeCluster(
            'foo-cluster', cluster_zone='foo-zone')._cluster_exists(), False)

  @mock.patch('time.sleep')
  @mock.patch('process_util.run_command')
  @mock.patch(
      'gke_cluster.GkeCluster.get_pod_status',
      return_value=gke_cluster.PodStatus.SUCCEEDED)
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_deploy_pod(self, unused_mock_cluster_exists,
                      unused_mock_get_pod_status, mock_call, unused_mock_sleep):
    gke_cluster.GkeCluster(
        'foo-cluster', cluster_zone='foo-zone').deploy_pod(
            pod_config='foo-config', pod_name='foo-pod')
    mock_call.assert_any_call(
        ['kubectl', 'create', '-f', '-'],
        std_input='foo-config',
        retries=1,
        retry_delay_sec=1)

  @mock.patch('time.sleep')
  @mock.patch('process_util.run_command')
  @mock.patch(
      'gke_cluster.GkeCluster.get_pod_status',
      return_value=gke_cluster.PodStatus.FAILED)
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_deploy_pod_fails(self, unused_mock_cluster_exists,
                            unused_mock_get_pod_status, unused_mock_call,
                            unused_mock_sleep):
    with self.assertRaises(RuntimeError):
      gke_cluster.GkeCluster(
          'foo-cluster', cluster_zone='foo-zone').deploy_pod(
              pod_config='foo-config', pod_name='foo-pod')

  @mock.patch('time.sleep')
  @mock.patch(
      'time.time', side_effect=(0, gke_cluster._PENDING_STATE_TIMEOUT_SEC + 1))
  @mock.patch('process_util.run_command')
  @mock.patch(
      'gke_cluster.GkeCluster.get_pod_status',
      return_value=gke_cluster.PodStatus.PENDING)
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_deploy_pod_timeout_in_pending_state(
      self, unused_mock_cluster_exists, unused_mock_get_pod_status,
      unused_mock_call, unused_mock_time, unused_mock_sleep):
    with self.assertRaises(RuntimeError):
      gke_cluster.GkeCluster(
          'foo-cluster', cluster_zone='foo-zone').deploy_pod(
              pod_config='foo-config', pod_name='foo-pod')

  @mock.patch('time.sleep')
  @mock.patch('process_util.run_command')
  @mock.patch('gke_cluster.GkeCluster.get_pod_status')
  @mock.patch('gke_cluster.GkeCluster.delete_pod')
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_deploy_pod_retries(self, unused_mock_cluster_exists,
                              unused_mock_delete_pod, mock_get_pod_status,
                              mock_call, unused_mock_sleep):
    mock_get_pod_status.side_effect = [
        gke_cluster.PodStatus.FAILED, gke_cluster.PodStatus.SUCCEEDED
    ]
    gke_cluster.GkeCluster(
        'foo-cluster', cluster_zone='foo-zone').deploy_pod(
            pod_config='foo-config', pod_name='foo-pod', retries=1)
    mock_call.assert_any_call(
        ['kubectl', 'create', '-f', '-'],
        std_input='foo-config',
        retries=1,
        retry_delay_sec=1)
    mock_call.assert_any_call(
        ['kubectl', 'replace', '--force', '-f', '-'],
        std_input='foo-config',
        retries=1,
        retry_delay_sec=1)

  @mock.patch('process_util.run_command', return_value='Succeeded')
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_get_pod_status(self, unused_mock_cluster_exists, mock_call):
    self.assertEqual(
        gke_cluster.GkeCluster(
            'foo-cluster',
            cluster_zone='foo-zone').get_pod_status(pod_name='foo-pod'),
        gke_cluster.PodStatus.SUCCEEDED)
    mock_call.assert_any_call(
        ['kubectl', 'get', 'pods', 'foo-pod', '-o', 'jsonpath={.status.phase}'],
        std_input=None,
        retries=0,
        retry_delay_sec=1)

  # First two side effects are for creating cluster and getting its status, and
  # needed but irrelevant to this test.
  @mock.patch(
      'process_util.run_command',
      side_effect=('foo-cluster', 'RUNNING', 'Pending', 'ImagePullBackOff'))
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_get_pod_status_pulling_image_failed(self, unused_mock_cluster_exists,
                                               mock_call):
    self.assertEqual(
        gke_cluster.GkeCluster(
            'foo-cluster',
            cluster_zone='foo-zone').get_pod_status(pod_name='foo-pod'),
        gke_cluster.PodStatus.FAILED)
    mock_call.assert_any_call(
        ['kubectl', 'get', 'pods', 'foo-pod', '-o', 'jsonpath={.status.phase}'],
        std_input=None,
        retries=0,
        retry_delay_sec=1)
    mock_call.assert_any_call([
        'kubectl', 'get', 'pods', 'foo-pod', '-o',
        'go-template="{{range .status.containerStatuses}}{{if not .ready}}'
        '{{.state.waiting.reason}}{{end}}{{end}}"'
    ],
                              std_input=None,
                              retries=0,
                              retry_delay_sec=1)

  @mock.patch('process_util.run_command', return_value='foo-status')
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_get_pod_status_unknown(self, unused_mock_cluster_exists, mock_call):
    self.assertEqual(
        gke_cluster.GkeCluster(
            'foo-cluster',
            cluster_zone='foo-zone').get_pod_status(pod_name='foo-pod'),
        gke_cluster.PodStatus.UNKNOWN)
    mock_call.assert_any_call(
        ['kubectl', 'get', 'pods', 'foo-pod', '-o', 'jsonpath={.status.phase}'],
        std_input=None,
        retries=0,
        retry_delay_sec=1)

  @mock.patch('process_util.run_command')
  @mock.patch('gke_cluster.GkeCluster._pod_exists', return_value=True)
  @mock.patch('gke_cluster.GkeCluster._cluster_exists', return_value=True)
  def test_delete_pod(self, unused_mock_cluster_exists, unused_mock_pod_exists,
                      mock_call):
    gke_cluster.GkeCluster(
        'foo-cluster', cluster_zone='foo-zone').delete_pod(
            pod_name='foo-pod', wait=True)
    mock_call.assert_any_call(
        ['kubectl', 'delete', 'pod', 'foo-pod', '--wait'],
        std_input=None,
        retries=1,
        retry_delay_sec=0)


if __name__ == '__main__':
  unittest.main()
