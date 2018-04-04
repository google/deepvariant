# Copyright 2018 Google Inc.
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
"""Tests for bed_reader CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import clif_postproc
from third_party.nucleus.io.python import bed_reader
from third_party.nucleus.protos import bed_pb2
from third_party.nucleus.testing import test_utils


class BedReaderTest(parameterized.TestCase):

  def setUp(self):
    self.bed = test_utils.genomics_core_testdata('test_regions.bed')
    self.zipped_bed = test_utils.genomics_core_testdata('test_regions.bed.gz')
    self.options = bed_pb2.BedReaderOptions()
    self.first = bed_pb2.BedRecord(
        reference_name='chr1',
        start=10,
        end=20,
        name='first',
        score=100,
        strand=bed_pb2.BedRecord.FORWARD_STRAND,
        thick_start=12,
        thick_end=18,
        item_rgb='255,124,1',
        block_count=3,
        block_sizes='2,6,2',
        block_starts='10,12,18')

  def test_bed_iterate(self):
    with bed_reader.BedReader.from_file(self.bed, self.options) as reader:
      self.assertEqual(reader.header.num_fields, 12)
      iterable = reader.iterate()
      self.assertIsInstance(iterable, clif_postproc.WrappedCppIterable)
      actual = list(iterable)
      self.assertLen(actual, 2)
      self.assertEqual(actual[0], self.first)

    zreader = bed_reader.BedReader.from_file(
        self.zipped_bed,
        bed_pb2.BedReaderOptions(
            compression_type=bed_pb2.BedReaderOptions.GZIP))
    self.assertEqual(zreader.header.num_fields, 12)
    with zreader:
      ziterable = zreader.iterate()
      self.assertIsInstance(ziterable, clif_postproc.WrappedCppIterable)
      zactual = list(ziterable)
      self.assertLen(zactual, 2)
      self.assertEqual(zactual[0], self.first)

  def test_from_file_raises_with_missing_bed(self):
    with self.assertRaisesRegexp(ValueError,
                                 'Not found: Could not open missing.bed'):
      bed_reader.BedReader.from_file('missing.bed', self.options)

  def test_ops_on_closed_reader_raise(self):
    reader = bed_reader.BedReader.from_file(self.bed, self.options)
    with reader:
      pass
    # At this point the reader is closed.
    with self.assertRaisesRegexp(ValueError, 'Cannot Iterate a closed'):
      reader.iterate()

  @parameterized.parameters('malformed.bed', 'malformed2.bed')
  def test_bed_iterate_raises_on_malformed_record(self, filename):
    malformed = test_utils.genomics_core_testdata(filename)
    reader = bed_reader.BedReader.from_file(malformed, self.options)
    iterable = iter(reader.iterate())
    self.assertIsNotNone(next(iterable))
    with self.assertRaises(ValueError):
      list(iterable)


if __name__ == '__main__':
  absltest.main()
