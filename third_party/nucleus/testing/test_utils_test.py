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
"""Tests for nucleus's testing.test_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest

from third_party.nucleus.protos import cigar_pb2
from third_party.nucleus.testing import test_utils


class TestUtilsTests(absltest.TestCase):

  def test_make_read(self):
    bases = 'ACG'
    quals = [30, 40, 50]
    cigar = '3M'
    mapq = 42
    chrom = 'chr10'
    start = 123
    name = 'myname'
    read = test_utils.make_read(
        bases,
        quals=quals,
        cigar=cigar,
        mapq=mapq,
        chrom=chrom,
        start=start,
        name=name)

    self.assertEqual(read.aligned_sequence, bases)
    self.assertEqual(read.aligned_quality, quals)
    self.assertEqual(list(read.alignment.cigar), [
        cigar_pb2.CigarUnit(
            operation_length=3, operation=cigar_pb2.CigarUnit.ALIGNMENT_MATCH)
    ])
    self.assertEqual(read.alignment.mapping_quality, mapq)
    self.assertEqual(read.alignment.position.reference_name, chrom)
    self.assertEqual(read.alignment.position.position, start)
    self.assertEqual(read.fragment_name, name)

  def test_make_read_produces_unique_read_names(self):
    start = 0
    read1 = test_utils.make_read('A', start=start)
    read2 = test_utils.make_read('A', start=start)
    self.assertGreater(len(read1.fragment_name), 0)
    self.assertGreater(len(read2.fragment_name), 0)
    self.assertNotEqual(read1.fragment_name, read2.fragment_name)


if __name__ == '__main__':
  absltest.main()
