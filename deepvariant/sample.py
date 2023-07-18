# Copyright 2023 Google LLC.
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
"""Sample dataclass."""

import dataclasses
from typing import List, Optional, Sequence

from deepvariant import variant_caller as vc_base
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from third_party.nucleus.io import sam
from third_party.nucleus.protos import reads_pb2


# ---------------------------------------------------------------------------
# Working with samples
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class Sample:
  """Organizes sample-level properties.

  options: A SampleOptions proto containing instructions for how to treat the
      sample, most of which will be set from flags.
  sam_readers: SamReader objects with handles on the `reads_filenames` from the
      options.
  in_memory_sam_reader: InMemorySamReader for this sample, which stores the
      alignments for this sample that have been read into memory from the
      sam_readers.
  reads: A list of reads queried from the sam readers.
  allele_counter: An allele counter object for the sample.
  variant_caller: A variant caller for the sample, should be instantiated using
      the options.variant_caller_options.
  """

  options: deepvariant_pb2.SampleOptions
  sam_readers: Optional[Sequence[sam.SamReader]] = None
  in_memory_sam_reader: Optional[sam.InMemorySamReader] = None
  reads: Optional[List[reads_pb2.Read]] = None
  allele_counter: Optional[allelecounter.AlleleCounter] = None
  variant_caller: Optional[vc_base.VariantCaller] = None

  def __repr__(self):
    return '<Sample {}>'.format(str(self.__dict__))
