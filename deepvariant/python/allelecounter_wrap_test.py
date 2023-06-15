# Copyright 2017 Google LLC.
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
"""Tests for AlleleCounter CLIF python wrappers."""



from absl.testing import absltest

from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.util import ranges
from deepvariant import testdata

from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter as _allelecounter


def setUpModule():
  testdata.init()


class WrapAlleleCounterTest(absltest.TestCase):

  def test_wrap(self):
    ref = fasta.IndexedFastaReader(testdata.CHR20_FASTA)
    sam_reader = sam.SamReader(testdata.CHR20_BAM)
    size = 100
    region = ranges.make_range('chr20', 10000000, 10000000 + size)
    options = deepvariant_pb2.AlleleCounterOptions(partition_size=size)
    allele_counter = _allelecounter.AlleleCounter(
        ref.c_reader, region, [], options
    )
    reads = list(sam_reader.query(region))
    self.assertGreater(len(reads), 0)
    for read in reads:
      allele_counter.add(read, 'sample_id')
    counts = allele_counter.counts()
    self.assertLen(counts, size)


if __name__ == '__main__':
  absltest.main()
