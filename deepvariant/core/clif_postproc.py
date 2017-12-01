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
"""CLIF postprocessors."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


def ValueErrorOnFalse(ok, *args):
  """Returns None / arg / (args,...) if ok."""
  if not isinstance(ok, bool):
    raise TypeError('Use ValueErrorOnFalse only on bool return value')
  if not ok:
    raise ValueError('CLIF wrapped call returned False')
  # Plain return args will turn 1 into (1,)  and None into () which is unwanted.
  if args:
    return args if len(args) > 1 else args[0]
  return None


class WrappedCppIterable(object):
  """This class gives Python iteration semantics on top of a C++ 'Iterable'."""

  def __init__(self, cc_iterable):
    self._cc_iterable = cc_iterable

  def __enter__(self):
    self._cc_iterable.__enter__()
    return self

  def __exit__(self, type_, value, traceback):
    self._cc_iterable.__exit__(type_, value, traceback)

  def __iter__(self):
    while True:
      not_done, record = self._cc_iterable.Next()
      if not_done:
        yield record
      else:
        return
