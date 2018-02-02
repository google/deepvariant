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
"""Library of application-specific errors.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import contextlib
import errno
import sys

from absl import logging


class Error(Exception):
  """Base class for core error types."""
  pass


class CommandLineError(Error):
  """Exception class related to invalid command-line flags."""
  pass


def log_and_raise(msg, exception_class=Error):
  """Logs the given message at ERROR level and raises exception.

  Args:
    msg: [`string`]. The message to log and use in the raised exception.
    exception_class: [`Exception`]. The class of exception to raise.

  Raises:
    Error: An exception of the type specified by the input exception_class.
  """
  logging.error(msg)
  raise exception_class(msg)


@contextlib.contextmanager
def clean_commandline_error_exit(
    allowed_exceptions=(Error, CommandLineError), exit_value=errno.ENOENT):
  """Wraps commands to capture certain exceptions and exit without stacktraces.

  This function is intended to wrap all code within main() of Python binaries
  to provide a mechanism for user errors to exit abnormally without causing
  exceptions to be thrown. Any exceptions that are subclasses of those listed
  in `allowed_exceptions` will be caught and the program will quietly exit with
  `exit_value`. Other exceptions are propagated normally. It should only be used
  as a context manager and its usage should be limited to main().

  Args:
    allowed_exceptions: [`tuple of Exception`]. A tuple of Exception classes
        that should not be raised, but instead quietly caused to exit the
        program.
    exit_value: [`int`]. The value to return upon program exit.

  Yields:
    The yield in this function is used to allow the block nested in the with
    statement to be executed.
  """
  try:
    yield
  except allowed_exceptions:
    sys.exit(exit_value)
