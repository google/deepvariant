# Copyright 2018 Google Inc.
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
"""GKE class for creating/reusing a GKE cluster."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import logging
import time
import enum

import process_util

# Number of times a gcloud CLI should be retried before reporting failure.
_GCLOUD_RETRIES = 1

# Delay (in seconds) between gcloud CLI retries.
_GCLOUD_RETRY_DELAY_SEC = 1

# Number of times a gcloud CLI should be retried before reporting failure.
_KUBECTL_RETRIES = 1

# Delay (in seconds) between gcloud CLI retries.
_KUBECTL_RETRY_DELAY_SEC = 1


@enum.unique
class ClusterStatus(enum.Enum):
  """Enums for GKE cluster status."""
  UNKNOWN = 0
  STATUS_UNSPECIFIED = 1
  PROVISIONING = 2
  RUNNING = 3
  RECONCILING = 4
  STOPPING = 5
  ERROR = 6
  DEGRADED = 7


_CLUSTER_STATUS_MAP = {
    # GKE cluster status from https://goo.gl/Ejieso.
    'STATUS_UNSPECIFIED': ClusterStatus.STATUS_UNSPECIFIED,
    'PROVISIONING': ClusterStatus.PROVISIONING,
    'RUNNING': ClusterStatus.RUNNING,
    'RECONCILING': ClusterStatus.RECONCILING,
    'STOPPING': ClusterStatus.STOPPING,
    'ERROR': ClusterStatus.ERROR,
    'DEGRADED': ClusterStatus.DEGRADED,
}


@enum.unique
class PodStatus(enum.Enum):
  """Enums for kubernetes pod status."""
  UNKNOWN = 0
  FAILED = 1
  PENDING = 2
  RUNNING = 3
  SUCCEEDED = 4


_POD_STATUS_MAP = {
    # Pod status from https://goo.gl/WnSPbQ.
    'Pending': PodStatus.PENDING,
    'Running': PodStatus.RUNNING,
    'Succeeded': PodStatus.SUCCEEDED,
    'Failed': PodStatus.FAILED,
}


class GkeCluster(object):
  """Helper for provisioning GKE cluster."""

  def __init__(self,
               cluster_name,
               cluster_region=None,
               cluster_zone=None,
               alpha_cluster=False,
               extra_create_args=None):
    """Helper for provisioning GKE cluster.

    On initialization, creates a new cluster or reuses if already exists. Uses
    gcloud's config file implicitly (for determining project, account, and etc).

    Args:
      cluster_name: (str) GKE cluster name. Must be unique within a project.
      cluster_region: (str) GCP region.
      cluster_zone: (str) GCP zone.
      alpha_cluster: (bool) whether the cluster is an alpha cluster.
      extra_create_args: (list) list of additional args (type str) to be used
          when creating a new cluster. E.g. --num-nodes=1, --enable-tpu, etc.

    Raises:
      ValueError: if both or neither of cluster region and zone is provided.
          Also if zone or region is specified in extra_creare_args.
    """
    if not bool(cluster_region) ^ bool(cluster_zone):
      raise ValueError('Either zone or region must be specified for a cluster.')
    if extra_create_args and ('--region' in extra_create_args or
                              '--zone' in extra_create_args):
      raise ValueError('extra_create_args cannot have --zone or --region. '
                       'Specify them explicitly via cluster_region or '
                       'cluster_zone.')

    self._cluster_name = cluster_name
    self._cluster_region = cluster_region
    self._cluster_zone = cluster_zone
    self._alpha_cluster = alpha_cluster
    self._extra_create_args = extra_create_args

    if self._cluster_exists():
      self._reuse_cluster()
    else:
      self._create_cluster()

  def _create_cluster(self):
    """Creates a kubernetes cluster.

    Takes care of clean up in case of interruption.

    Raises:
      RuntimeError: if interrupted.
    """
    args = []
    if self._alpha_cluster:
      args = ['gcloud', 'alpha']
    else:
      args = ['gcloud']
    args.extend(['container', 'clusters', 'create', self._cluster_name])
    if self._extra_create_args:
      args.extend(self._extra_create_args)

    logging.info('Creating GKE cluster: %s ...', self._cluster_name)
    try:
      self._gcloud_call(args)
    except KeyboardInterrupt:
      logging.error(
          'GKE Cluster creation interrupted. Deallocating the cluster %s ...',
          self._cluster_name)
      self.delete_cluster(wait=False)

  def _reuse_cluster(self):
    """Configures kubernetes tools (kubectl) to use existing cluster.

    Raises:
      RuntimeError: if cluster does not exist or is not reachable.
    """
    self._store_cluster_credentials()
    cluster_status = self._get_cluster_status()
    # Wait for cluster to become ready.
    while cluster_status == ClusterStatus.PROVISIONING:
      time.sleep(_GCLOUD_RETRY_DELAY_SEC)
      cluster_status = self._get_cluster_status()
    if cluster_status in [
        ClusterStatus.STOPPING, ClusterStatus.ERROR, ClusterStatus.DEGRADED
    ]:
      raise RuntimeError('Provided GKE cluster is not reachable. Cluster: %s' %
                         self._cluster_name)

  def _store_cluster_credentials(self):
    """Refreshes credentials for the cluster and stores it in kubectl config.

       See http://goo.gl/cFyUDz for more details.
    """
    args = [
        'gcloud', 'container', 'clusters', 'get-credentials', self._cluster_name
    ]
    logging.info('Storing credential for GKE cluster: %s ...',
                 self._cluster_name)
    self._gcloud_call(args)

  def _get_cluster_status(self):
    """Returns cluster's status."""
    args = [
        'gcloud', 'container', 'clusters', 'describe', self._cluster_name,
        '--format=value(status)'
    ]
    status_str = self._gcloud_call(args).strip()
    if status_str in _CLUSTER_STATUS_MAP:
      return _CLUSTER_STATUS_MAP[status_str]
    return ClusterStatus.UNKNOWN

  def delete_cluster(self, wait=False):
    """Deletes GKE cluster.

    Args:
      wait: (bool) if set blocks on completion.

    Raises:
      ValueError: if cluster does not exist.
      RuntimeError: if fails to delete the cluster.
    """
    if not self._cluster_exists():
      raise ValueError(
          'Cannot delete a non-existent cluster: %s' % self._cluster_name)
    status = self._get_cluster_status()
    if status == ClusterStatus.STOPPING:
      logging.warning(
          'Cannot delete GKE cluster %s. Cluster is being already deleted.',
          self._cluster_name)
      return
    # Cannot delete in PROVISIONING state.
    while status == ClusterStatus.PROVISIONING:
      time.sleep(_GCLOUD_RETRY_DELAY_SEC)
      status = self._get_cluster_status()

    args = [
        'gcloud', 'container', 'clusters', 'delete', self._cluster_name,
        '--quiet'
    ]
    if not wait:
      args = args + ['--async']
    try:
      self._gcloud_call(args)
    except RuntimeError:
      # redacted
      raise RuntimeError('Failed to delete cluster: %s. Please delete it '
                         'manually on GCP console.' % self._cluster_name)

  def _cluster_exists(self):
    """Returns true iff the cluster exists (not deleted)."""
    args = ['gcloud', 'container', 'clusters', 'list', '--format=value(name)']
    return self._cluster_name in self._gcloud_call(args).splitlines()

  def _gcloud_call(self,
                   args,
                   retry_delay_sec=_GCLOUD_RETRY_DELAY_SEC,
                   retries=_GCLOUD_RETRIES):
    """Make a gcloud CLI call.

    Args:
      args: (list)  A list of arguments (type string) to pass to Popen  call.
      retry_delay_sec: (int) delay in retries.
      retries: (int) number of retries.

    Returns:
      stdout of process call.
    """
    if self._cluster_region:
      args.extend(['--region', self._cluster_region])
    if self._cluster_zone:
      args.extend(['--zone', self._cluster_zone])
    return process_util.run_command(
        args, retry_delay_sec=retry_delay_sec, retries=retries)

  def deploy_pod(self, pod_config, pod_name, retries=0, wait=True):
    """Deploy a pod into Kubernetes cluster.

    Args:
      pod_config: (str) pod config in json format.
      pod_name: (str) Pod's name. Must be the same as the one used in the file.
      retries: (int) retries this number of times if pod status is
        PodStatus.Failure.
      wait: (bool) Whether to wait on completion. If retries is positive, it
        waits on completion regardless.

    Raises:
      RuntimeError: if pod fails after all retries.
    """
    def get_args(is_first_try=False):
      """Returns kubectl CLI args.

      Args:
        is_first_try: (bool) whether it is the first try.
      """
      if is_first_try:
        return ['kubectl', 'create', '-f', '-']
      else:
        # With replace arg, the pod is first deleted and then re-deployed.
        return ['kubectl', 'replace', '--force', '-f', '-']

    logging.info('Deploying pod %s to GKE cluster %s.', pod_name,
                 self._cluster_name)
    status = None
    for i in range(retries + 1):
      is_first_try = (i == 0)
      self._kubectl_call(get_args(is_first_try), std_input=pod_config)
      if not wait and not retries:
        return
      status = self.get_pod_status(pod_name)
      while status != PodStatus.SUCCEEDED and status != PodStatus.FAILED:
        time.sleep(_KUBECTL_RETRY_DELAY_SEC)
        status = self.get_pod_status(pod_name)
      if status == PodStatus.SUCCEEDED:
        self.delete_pod(pod_name, wait=True)
        return
      if i < retries:
        logging.warning('Retrying deploying pod (attempt %d/%d): ', i + 1,
                        retries)

    self.delete_pod(pod_name, wait=True)
    raise RuntimeError(
        'Pod %s failed after %d attempts.' % (pod_name, retries + 1))

  def get_pod_status(self, pod_name):
    """Returns given pod's status.

    Args:
      pod_name: (str) name of pod.

    Returns:
    """
    args = [
        'kubectl', 'get', 'pods', pod_name, '-o', 'jsonpath={.status.phase}'
    ]
    status_str = self._kubectl_call(args).strip()
    return _POD_STATUS_MAP.get(status_str, PodStatus.UNKNOWN)

  def delete_pod(self, pod_name, wait=True):
    """Deletes the given pod.

    Args:
      pod_name: (str) pod's name.
      wait: (bool) whether to wait on completion.
    """
    if not self._pod_exists(pod_name):
      return
    args = ['kubectl', 'delete', 'pod', pod_name]
    if wait:
      args += ['--wait']
      # Retry right away if pod fails and we have already waited for it.
      self._kubectl_call(args, retry_delay_sec=0)
    else:
      self._kubectl_call(args)

  def _pod_exists(self, pod_name):
    """Returns true iff the pod exists (not deleted)."""
    args = [
        'kubectl', 'get', 'pods', '-o',
        'jsonpath={.items[*].spec.containers[*].name}'
    ]
    return pod_name in self._kubectl_call(args).split()

  def _kubectl_call(self,
                    args,
                    std_input=None,
                    retry_delay_sec=_KUBECTL_RETRY_DELAY_SEC,
                    retries=_KUBECTL_RETRIES):
    """Make a kubectl CLI call.

    Args:
      args: (list)  A list of arguments (type string) to pass to Popen call.
      std_input: (str) standard input to be passed.
      retry_delay_sec: (int) delay in retries.
      retries: (int) number of retries.

    Returns:
      stdout of process call.
    """
    return process_util.run_command(
        args,
        std_input=std_input,
        retry_delay_sec=retry_delay_sec,
        retries=retries)
