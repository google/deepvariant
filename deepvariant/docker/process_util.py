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
"""Helpers for running a command."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import logging
import subprocess
import time


def run_command(args, std_input=None, retry_delay_sec=1, retries=0):
  """Runs a command with optional retry behaviour.

  Args:
    args: (list) A list of arguments (type string) to pass to Popen.
    std_input: (str) will be passed as stdin to the command.
    retry_delay_sec: (int) delay in retries.
    retries: (int) number of retries.

  Returns:
    stdout.

  Raises:
    ValueError: if number of retries is less than zero.
    RuntimeError: if process returns a non-zero exit code in all retries.
  """
  if retries < 0:
    raise ValueError('Number of retries cannot be negative.')

  logging.debug('Calling command: %s', ' '.join(args))
  for i in range(retries + 1):
    if std_input:
      process = subprocess.Popen(
          args,
          stdin=subprocess.PIPE,
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
      stdout, stderr = process.communicate(input=std_input)
    else:
      process = subprocess.Popen(
          args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      stdout, stderr = process.communicate()
    if process.returncode == 0:
      if stderr:
        logging.info('%s succeeded with stderr: %s', ' '.join(args), stderr)
      return stdout
    logging.warning('%s exited with code %d: %s', ' '.join(args),
                    process.returncode, (stdout + stderr).strip())
    if i < retries:
      logging.warning('Retrying command (attempt %d/%d): ', i + 1, retries)
      time.sleep(retry_delay_sec)
  raise RuntimeError(
      '%s failed after %d attempts.' % (' '.join(args), retries + 1))
