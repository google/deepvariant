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

from typing import Optional, Sequence

from third_party.nucleus.io import sam
from deepvariant.protos import deepvariant_pb2


# redacted
class Sample(object):
  """Sample organizes sample-level properties and sam readers.

  name: This is the same name as set by --sample_name in make_examples.
  role: A string to identify the role of this sample in the analysis, e.g. in
      trios 'child', 'parent1', or 'parent2'. Importantly, `role` strings should
      not be checked inside make_examples_core.py. For example, instead of
      checking whether sample.role == "child" to set pileup height, instead add
      the pileup height to the sample, adding new properties to this class as
      needed. This keeps make_examples_core.py functioning for multiple samples
      without having to reason about sample roles that belong to each
      application.
  reads_filenames: Filenames of the bam files corresponding to this sample. A
      list of just one is most common, but multiple are supported.
  sam_readers: SamReader objects with handles on the `reads_filenames`.
  in_memory_sam_reader: InMemorySamReader for this sample, which stores the
      alignments for this sample that have been read into memory from the
      sam_readers.
  variant_caller_options: Options for variant-calling on this sample.
  pileup_height: Integer height of the pileup image desired for this sample.
  order: List of integers indicating the order in which samples should be shown
      in the pileup image when calling on this sample. The indices refer to the
      list of samples in the regionprocessor.
  downsample_fraction: Fraction by which to downsample the reads.
  """
  name: Optional[str] = None
  role: Optional[str] = None
  reads_filenames: Optional[Sequence[str]] = None
  sam_readers: Optional[Sequence[sam.SamReader]] = None
  in_memory_sam_reader: Optional[sam.InMemorySamReader] = None
  variant_caller_options: Optional[deepvariant_pb2.VariantCallerOptions] = None
  pileup_height: Optional[int] = None
  order: Optional[Sequence[int]] = None
  downsample_fraction: Optional[float] = None

  def __init__(self,
               name=None,
               role='default_role',
               reads_filenames=None,
               sam_readers=None,
               in_memory_sam_reader=None,
               variant_caller_options=None,
               pileup_height=None,
               order=None,
               downsample_fraction=None):
    self.name = name
    self.role = role
    self.reads_filenames = reads_filenames
    self.sam_readers = sam_readers
    self.in_memory_sam_reader = in_memory_sam_reader
    self.variant_caller_options = variant_caller_options
    self.order = order
    self.downsample_fraction = downsample_fraction

    if pileup_height is not None:
      # Downstream, None defaults to the PileupImageOptions height.
      if isinstance(pileup_height, int):
        self.pileup_height = pileup_height
      else:
        raise ValueError('pileup_height must an integer.')

  def __repr__(self):
    return '<Sample {}>'.format(str(self.__dict__))
