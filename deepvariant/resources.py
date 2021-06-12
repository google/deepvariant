# Copyright 2017 Google LLC.  All Rights Reserved.
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
"""Library to gather runtime performance metrics.

This module exposes the ResourceMonitor class, which client code can use to
gather resource usage metrics about their program. An example usage would look
something like:

  with ResourceMonitor() as monitor:
    ... do work ...
    metrics = monitor.metrics()
"""

import platform
import resource
import time

import psutil

from deepvariant.protos import resources_pb2


class ResourceMonitor(object):
  """Class for collecting resource usage info from this or child process."""

  def __init__(self):
    """Constructs a ResourceMonitor object."""
    self.wall_start = None
    self.metrics_pb = self._initial_metrics_protobuf()

  def _initial_metrics_protobuf(self):
    """Returns an initialized ResourceMetrics proto.

    This function also fills in the "constant" fields of the ResourceMetrics
    proto that don't depend on the actual running commands, such as host_name.

    Returns:
      learning.genomics.deepvariant.ResourceMetrics proto.
    """
    return resources_pb2.ResourceMetrics(
        host_name=_get_host_name(),
        cpu_frequency_mhz=_get_cpu_frequency(),
        physical_core_count=_get_cpu_count(),
        total_memory_mb=_get_total_memory())

  def __enter__(self):
    return self.start()

  def __exit__(self, unused_type, unused_value, unused_traceback):
    pass

  def start(self):
    """Starts timers associated with resource collection.

    This method must be called before metrics().

    Returns:
      self to enable the idiom `monitor = ResourceMonitor().start()`.
    """
    self.wall_start = time.time()
    return self

  def metrics(self):
    """Collects and return runtime metrics as a ResourceMetrics proto.

    This method can be called multiple times, but wall clock time is always
    reckoned from the time of the last start() call.

    Returns:
      A learning.genomics.deepvariant.ResourceMetrics proto message.

    Raises:
      RuntimeError: if start() was not called previously.
    """
    if self.wall_start is None:
      raise RuntimeError('start() must be called prior to metrics()')

    self.metrics_pb.wall_time_seconds = time.time() - self.wall_start

    # Consider using psutil.cpu_times() instead to get more detailed information
    # about the usage in self and all children.
    try:
      rusage = resource.getrusage(resource.RUSAGE_SELF)
      self.metrics_pb.cpu_user_time_seconds = rusage.ru_utime
      self.metrics_pb.cpu_system_time_seconds = rusage.ru_stime
      self.metrics_pb.memory_peak_rss_mb = int(rusage.ru_maxrss / 1024)
    except resource.error:
      # The OS call to get rusage failed, so just don't set the field values,
      # leaving them as the default values of 0.
      pass

    # Create a psutil.Process pointed at the current process.
    process = psutil.Process()
    io_counters = process.io_counters()
    self.metrics_pb.read_bytes = io_counters.read_bytes
    self.metrics_pb.write_bytes = io_counters.write_bytes

    return self.metrics_pb


# ------------------------------------------------------------------------------
# Simple functions for getting host_name, cpu count, etc. Isolated here to make
# them mockable.
# ------------------------------------------------------------------------------


def _get_host_name():
  """Gets the host name of this machine."""
  return platform.node()


def _get_cpu_count():
  """Gets the number of physical cores in this machine.

  Returns:
    int >= 1 if the call to get the cpu_count succeeded, or 0 if not.
  """
  return psutil.cpu_count(logical=False) or 0


def _get_cpu_frequency():
  """Gets the frequency in MHz of the cpus in this machine.

  Returns:
    float > 0 if the call to get the cpu_frequency succeeded. This information
    may not be available on all systems, in which case we return 0.0.
  """
  try:
    freq = psutil.cpu_freq()
    return freq.current if freq is not None else 0.0
  except NotImplementedError:
    return 0.0


def _get_total_memory():
  """Gets the total memory in megabytes in this machine."""
  return int(psutil.virtual_memory().total / (1024 * 1024))
