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
"""Tests for third_party.nucleus.io.beam_readers."""

from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import beam_readers
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.testing import test_utils

# pylint:disable=g-bad-import-order
# Beam needs to be imported after nucleus for this to work in open source.
from apache_beam.testing.test_pipeline import TestPipeline
from apache_beam.testing import util as beam_testing_util
# pylint:enable=g-bad-import-order


class TestReadBed(parameterized.TestCase):
  """Tests for reading BED files using Beam."""

  def _records_match(self, expected):
    def _equal(records):
      self.assertLen(records, len(expected))
      self.assertEqual([(r.reference_name, r.start, r.end) for r in records],
                       expected)
    return _equal

  def _len_match(self, expected):
    def _equal(records):
      self.assertLen(records, expected)
    return _equal

  @parameterized.parameters('test_regions.bed', 'test_regions.bed.gz',
                            'test_regions.bed.tfrecord',
                            'test_regions.bed.tfrecord.gz')
  def test_iterate_bed_reader(self, bed_filename):
    bed_path = test_utils.genomics_core_testdata(bed_filename)
    expected = [('chr1', 10, 20), ('chr1', 100, 200)]
    bed_path = test_utils.genomics_core_testdata(bed_filename)
    with TestPipeline() as p:
      records = (p | beam_readers.ReadBed(bed_path))
      beam_testing_util.assert_that(records, self._records_match(expected))

  def test_process_single_file_with_num_fields(self):
    # This BED file has 8 columns, but we will only read in three.
    bed_path = test_utils.genomics_core_testdata('test_regions.bed')
    with TestPipeline() as p:
      records = (p | beam_readers.ReadBed(bed_path, num_fields=3))
      beam_testing_util.assert_that(records, self._len_match(2))

  @parameterized.parameters(1, 2, 7, 10, 11, 13)
  def test_invalid_num_fields(self, invalid_num_fields):
    bed_path = test_utils.genomics_core_testdata('test_regions.bed')
    with self.assertRaisesRegexp(ValueError, 'Invalid requested number of fie'):
      with TestPipeline() as p:
        _ = (p | beam_readers.ReadBed(bed_path, num_fields=invalid_num_fields))

  def test_read_multiple_files(self):
    bed_path = test_utils.genomics_core_testdata('test.bed*')
    # test.bed and test.bed.gz each have four records.
    with TestPipeline() as p:
      records = (p | beam_readers.ReadBed(bed_path))
      beam_testing_util.assert_that(
          records, self._len_match(8))


class TestReadSam(parameterized.TestCase):
  """Tests for reading SAM/BAM files using Beam."""

  def _len_match(self, expected):
    def _equal(records):
      self.assertLen(records, expected)
    return _equal

  @parameterized.parameters('test.sam', 'test.sam.golden.tfrecord')
  def test_sam_iterate(self, sam_filename):
    with TestPipeline() as p:
      records = (
          p
          | beam_readers.ReadSam(
              test_utils.genomics_core_testdata(sam_filename)))
      beam_testing_util.assert_that(records, self._len_match(6))

  def test_bam_iterate(self):
    with TestPipeline() as p:
      records = (
          p
          | beam_readers.ReadSam(test_utils.genomics_core_testdata('test.bam')))
      beam_testing_util.assert_that(records, self._len_match(106))

  @parameterized.parameters(True, False)
  def test_read_multiple_files(self, keep_unaligned):
    file_pattern = test_utils.genomics_core_testdata('test.*am')
    read_requirements = reads_pb2.ReadRequirements(
        keep_unaligned=keep_unaligned)
    # test.sam contains 6 records and test.bam has 106 records. Both files
    # contain 1 unmapped record.
    expected_count = 112 if keep_unaligned else 110
    with TestPipeline() as p:
      result = (
          p
          | beam_readers.ReadSam(
              file_pattern, read_requirements=read_requirements))
      beam_testing_util.assert_that(result, self._len_match(expected_count))


if __name__ == '__main__':
  absltest.main()
