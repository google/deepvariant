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
"""Tests for learning.genomics.deepvariant.resources."""

import resource
from unittest import mock
from absl.testing import absltest
# https://stackoverflow.com/questions/34630393/python2-7-contextlib-exitstack-equivalent
import contextlib2

from deepvariant import resources


class ResourcesTest(absltest.TestCase):

  def test_monitor_real_function_calls(self):
    # We want to actually make all of the real function calls under test, but
    # we of course don't know their values and can only do sanity checks.
    with resources.ResourceMonitor() as monitor:
      metrics = monitor.metrics()
      self.assertNotEmpty(metrics.host_name)
      self.assertGreater(metrics.physical_core_count, 0)
      self.assertGreater(metrics.total_memory_mb, 0)
      self.assertGreater(metrics.cpu_user_time_seconds, 0)
      self.assertGreater(metrics.cpu_system_time_seconds, 0)
      self.assertGreater(metrics.memory_peak_rss_mb, 0)

      # We unfortunately cannot make sure that read_bytes and write_bytes is
      # greater than zero, so these tests are commented out.
      # self.assertGreater(metrics.read_bytes, 0)
      # self.assertGreater(metrics.write_bytes, 0)

      # CPU frequency may not be available on all systems, so the value is
      # either a real frequency (> 0) or the magic value of 0.0 indicating that
      # the value could not be determined.
      self.assertGreaterEqual(metrics.cpu_frequency_mhz, 0.0)

  def test_monitor_gets_expected_metric_values(self):
    patchers = [
        mock.patch.object(resources, '_get_host_name', return_value='hostname'),
        mock.patch.object(resources, '_get_cpu_count', return_value=3),
        mock.patch.object(resources, '_get_cpu_frequency', return_value=2312.4),
        mock.patch.object(resources, '_get_total_memory', return_value=1234),
        mock.patch.object(
            resources.resource,
            'getrusage',
            return_value=mock.Mock(
                ru_utime=1.0,
                ru_stime=2.0,
                ru_maxrss=3 * 1024  # 3 Megabytes
            )),
    ]

    process_mock = mock.Mock()
    process_mock.io_counters.return_value = mock.Mock(
        read_bytes=12, write_bytes=34)
    patchers.append(
        mock.patch.object(
            resources.psutil, 'Process', return_value=process_mock))

    with contextlib2.ExitStack() as stack:
      for ctx in patchers:
        stack.enter_context(ctx)
      with resources.ResourceMonitor() as monitor:
        metrics = monitor.metrics()

        # Environment metrics; all mocked out.
        self.assertEqual(metrics.host_name, 'hostname')
        self.assertEqual(metrics.physical_core_count, 3)
        self.assertEqual(metrics.cpu_frequency_mhz, 2312.4)
        self.assertEqual(metrics.total_memory_mb, 1234)

        # Runtime metrics; they are all mocked out.
        self.assertEqual(metrics.cpu_user_time_seconds, 1.0)
        self.assertEqual(metrics.cpu_system_time_seconds, 2.0)
        self.assertEqual(metrics.memory_peak_rss_mb, 3)
        self.assertEqual(metrics.read_bytes, 12)
        self.assertEqual(metrics.write_bytes, 34)

  def test_start_returns_self(self):
    monitor = resources.ResourceMonitor()
    self.assertIs(monitor.start(), monitor)

  def test_metrics_without_start_raises_exception(self):
    monitor = resources.ResourceMonitor()
    with self.assertRaises(RuntimeError):
      monitor.metrics()

  def test_metrics_is_ok_if_rusage_fails(self):
    # Some psutil functions, such as cpu_freq(), can return None depending on
    # the environment; make sure we don't crash when that occurs.
    with mock.patch.object(
        resources.resource, 'getrusage', side_effect=resource.error):
      with resources.ResourceMonitor() as monitor:
        self.assertEqual(monitor.metrics().cpu_user_time_seconds, 0)
        self.assertEqual(monitor.metrics().cpu_system_time_seconds, 0)
        self.assertEqual(monitor.metrics().memory_peak_rss_mb, 0)

  def test_metrics_is_ok_when_cpu_freq_returns_none(self):
    # Some psutil functions, such as cpu_freq(), can return None depending on
    # the environment; make sure we don't crash when that occurs.
    with mock.patch.object(resources.psutil, 'cpu_freq', return_value=None):
      with resources.ResourceMonitor() as monitor:
        self.assertEqual(monitor.metrics().cpu_frequency_mhz, 0)

  def test_metrics_is_ok_when_cpu_freq_throws(self):
    # psutil.cpu_freq() can throw NotImplementedError in certain environments;
    # make sure we don't crash when that occurs.
    with mock.patch.object(
        resources.psutil, 'cpu_freq', side_effect=NotImplementedError):
      with resources.ResourceMonitor() as monitor:
        self.assertEqual(monitor.metrics().cpu_frequency_mhz, 0)

  def test_metrics_is_ok_when_cpu_count_returns_none(self):
    # Some psutil functions, such as cpu_freq(), can return None depending on
    # the environment; make sure we don't crash when that occurs.
    with mock.patch.object(resources.psutil, 'cpu_count', return_value=None):
      with resources.ResourceMonitor() as monitor:
        self.assertEqual(monitor.metrics().physical_core_count, 0)


if __name__ == '__main__':
  absltest.main()
