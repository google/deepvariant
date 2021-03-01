# Copyright 2020 Google LLC.
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
"""Shareable functionality for make_examples."""

from typing import Sequence

from third_party.nucleus.io import sam
from deepvariant.protos import deepvariant_pb2


# redacted
class Sample(object):
  """Sample organizes sample-level properties and sam readers."""
  name: str = None
  nickname: str = None
  reads_filenames: Sequence[str] = None
  sam_readers: Sequence[sam.SamReader] = None
  in_memory_sam_reader: sam.InMemorySamReader = None
  variant_caller_options: deepvariant_pb2.VariantCallerOptions = None
  pileup_height: int = None

  def __init__(self,
               sam_readers=None,
               in_memory_sam_reader=None,
               name=None,
               nickname='default',
               pileup_height=None,
               reads_filenames=None,
               variant_caller_options=None):
    self.name = name
    self.nickname = nickname
    self.reads_filenames = reads_filenames
    self.sam_readers = sam_readers
    self.in_memory_sam_reader = in_memory_sam_reader
    self.variant_caller_options = variant_caller_options

    if pileup_height is not None:
      # Downstream, None defaults to the PileupImageOptions height.
      if isinstance(pileup_height, int):
        self.pileup_height = pileup_height
      else:
        raise ValueError('pileup_height must an integer.')

  def __repr__(self):
    return '<Sample {}>'.format(str(self.__dict__))
