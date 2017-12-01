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
"""Tests for deepvariant .core.errors."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import errno
import sys



from absl.testing import absltest
from absl.testing import parameterized
import mock
from absl import logging
from deepvariant.core import errors


class ErrorsTest(parameterized.TestCase):

  @parameterized.parameters(
      ('empty flag', errors.CommandLineError),
      ('bad value', ValueError),
      ('base error', errors.Error),
  )
  def test_log_and_raise(self, msg, cls):
    with mock.patch.object(logging, 'error') as mock_logging:
      with self.assertRaisesRegexp(cls, msg):
        errors.log_and_raise(msg, cls)
      mock_logging.assert_called_once_with(msg)

  @parameterized.parameters(
      (ValueError, 'ValueError exception'),
      (IOError, 'IOError exception'),
  )
  def test_clean_commandline_error_exit_raise_non_allowed(self, exc_type, msg):
    with self.assertRaisesRegexp(exc_type, msg):
      with errors.clean_commandline_error_exit():
        raise exc_type(msg)

  @parameterized.parameters(
      (errors.CommandLineError, errno.ENOENT),
      (errors.Error, errno.EINVAL),
  )
  def test_clean_commandline_error_exit_clean_exit(self, exc_type, exit_value):
    with mock.patch.object(sys, 'exit') as mock_exit:
      with errors.clean_commandline_error_exit(exit_value=exit_value):
        raise exc_type()
    mock_exit.assert_called_once_with(exit_value)


if __name__ == '__main__':
  absltest.main()
