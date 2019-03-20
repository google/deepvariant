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
"""A Python interface for files."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import six

from third_party.nucleus.io.python import gfile


def Exists(filename):
  return gfile.Exists(filename)


def Glob(pattern):
  return gfile.Glob(pattern)


class ReadableFile(six.Iterator):
  """Wraps gfile.ReadableFile to add iteration, enter/exit and readlines."""

  def __init__(self, filename):
    self._file = gfile.ReadableFile.New(filename)

  def __enter__(self):
    return self

  def __exit__(self, type_, value, traceback):
    self._file.__exit__()

  def __iter__(self):
    return self

  def __next__(self):
    ok, line = self._file.Readline()
    if ok:
      return line
    else:
      raise StopIteration

  def readlines(self):
    lines = []
    while True:
      ok, line = self._file.Readline()
      if ok:
        lines.append(line)
      else:
        break
    return lines


def Open(filename, mode="r"):
  if mode.startswith("r"):
    return ReadableFile(filename)
  elif mode.startswith("w"):
    return gfile.WritableFile.New(filename)
  else:
    raise ValueError("Unsupported mode '{}' for Open".format(mode))
