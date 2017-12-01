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
"""Tests for VariantCalling CLIF python wrappers."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest

from deepvariant import test_utils

from deepvariant.core import genomics_io
from deepvariant.core import ranges
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter as _allelecounter
from deepvariant.python import variant_calling


def setUpModule():
  test_utils.init()


class WrapVariantCallingTest(absltest.TestCase):

  def test_call_from_allele_counter(self):
    ref = genomics_io.make_ref_reader(test_utils.CHR20_FASTA)
    sam_reader = genomics_io.make_sam_reader(test_utils.CHR20_BAM)
    size = 1000
    region = ranges.make_range('chr20', 10000000, 10000000 + size)
    allele_counter = _allelecounter.AlleleCounter(
        ref, region, deepvariant_pb2.AlleleCounterOptions(partition_size=size))
    caller = variant_calling.VariantCaller(
        deepvariant_pb2.VariantCallerOptions(
            min_count_snps=2,
            min_count_indels=2,
            min_fraction_snps=0.12,
            min_fraction_indels=0.12,
            sample_name='sample_name',
            p_error=0.001,
            max_gq=50,
            gq_resolution=1,
            ploidy=2))

    # Grab all of the reads in our region and add them to the allele_counter.
    reads = list(sam_reader.query(region))
    self.assertNotEmpty(reads)
    for read in reads:
      allele_counter.add(read)

    # Get the candidates records for this whole region.
    candidates = caller.calls_from_allele_counter(allele_counter)

    # We should have at least some candidates and some gvcf records.
    self.assertNotEmpty(candidates)

    # Each candidate should be a DeepVariantCall.
    for candidate in candidates:
      self.assertIsInstance(candidate, deepvariant_pb2.DeepVariantCall)


if __name__ == '__main__':
  absltest.main()
