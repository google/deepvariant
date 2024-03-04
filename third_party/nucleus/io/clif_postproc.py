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
"""CLIF postprocessors."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import abc
from typing import Any

import six

from third_party.nucleus.protos import bed_pb2
from third_party.nucleus.protos import bedgraph_pb2
from third_party.nucleus.protos import fastq_pb2
from third_party.nucleus.protos import gff_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import variants_pb2


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


def ValueErrorOnInaccurate(accuracy: int, alignment: Any):
  """Returns Alignment if accuracy is good (equals 0).

  See: third_party/ssw/src/ssw_cpp.h;l=137
  """
  if accuracy == 0:
    return alignment
  else:
    raise ValueError(f"Alignment is not accurate, returned '{accuracy}.'")


class WrappedCppIterable(six.Iterator):
  """This class gives Python iteration semantics on top of a C++ 'Iterable'."""

  __metaclass__ = abc.ABCMeta

  def __init__(self, cc_iterable):
    self._cc_iterable = cc_iterable

  def __enter__(self):
    self._cc_iterable.__enter__()
    return self

  def __exit__(self, type_, value, traceback):
    self._cc_iterable.__exit__(type_, value, traceback)

  def __iter__(self):
    return self

  @abc.abstractmethod
  def _raw_next(self):
    """Sub-classes should implement __next__ in this method."""

  def __next__(self):
    try:
      record, not_done = self._raw_next()
    except AttributeError:
      if self._cc_iterable is None:
        raise ValueError('No underlying iterable. This may occur if trying to '
                         'create multiple concurrent iterators from the same '
                         'reader. Try wrapping your call to the iterator in a '
                         '`with` block or materializing the entire iterable '
                         'explicitly.')
      else:
        raise
    if not_done:
      return record
    else:
      raise StopIteration


class WrappedBedIterable(WrappedCppIterable):

  def _raw_next(self):
    record = bed_pb2.BedRecord()
    not_done = self._cc_iterable.PythonNext(record)
    return record, not_done


class WrappedBedGraphIterable(WrappedCppIterable):

  def _raw_next(self):
    record = bedgraph_pb2.BedGraphRecord()
    not_done = self._cc_iterable.PythonNext(record)
    return record, not_done


class WrappedFastqIterable(WrappedCppIterable):

  def _raw_next(self):
    record = fastq_pb2.FastqRecord()
    not_done = self._cc_iterable.PythonNext(record)
    return record, not_done


class WrappedGffIterable(WrappedCppIterable):

  def _raw_next(self):
    record = gff_pb2.GffRecord()
    not_done = self._cc_iterable.PythonNext(record)
    return record, not_done


class WrappedReferenceIterable(WrappedCppIterable):

  def _raw_next(self):
    not_done, record = self._cc_iterable.Next()
    return record, not_done


class WrappedSamIterable(WrappedCppIterable):

  def _raw_next(self):
    record = reads_pb2.Read()
    not_done = self._cc_iterable.PythonNext(record)
    return record, not_done


class WrappedVariantIterable(WrappedCppIterable):

  def _raw_next(self):
    record = variants_pb2.Variant()
    not_done = self._cc_iterable.PythonNext(record)
    return record, not_done
