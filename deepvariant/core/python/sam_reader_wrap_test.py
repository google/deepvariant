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
"""Tests for sam_reader CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest

from deepvariant.core import clif_postproc
from deepvariant.core import ranges
from deepvariant.core import test_utils
from deepvariant.core.protos import core_pb2
from deepvariant.core.python import sam_reader


class SamReaderTest(absltest.TestCase):

  def setUp(self):
    self.bam = test_utils.genomics_core_testdata('test.bam')
    self.options = core_pb2.SamReaderOptions()
    self.indexed_options = core_pb2.SamReaderOptions(
        index_mode=core_pb2.INDEX_BASED_ON_FILENAME)

  def test_bam_iterate(self):
    reader = sam_reader.SamReader.from_file(self.bam, self.options)
    with reader:
      iterable = reader.iterate()
      self.assertIsInstance(iterable, clif_postproc.WrappedCppIterable)
      self.assertEqual(test_utils.iterable_len(iterable), 106)

  def test_bam_query(self):
    reader = sam_reader.SamReader.from_file(self.bam, self.indexed_options)
    expected = [(ranges.parse_literal('chr20:10,000,000-10,000,100'), 106),
                (ranges.parse_literal('chr20:10,000,000-10,000,000'), 45)]
    with reader:
      for interval, n_expected in expected:
        with reader.query(interval) as iterable:
          self.assertIsInstance(iterable, clif_postproc.WrappedCppIterable)
          self.assertEqual(test_utils.iterable_len(iterable), n_expected)

  def test_bam_samples(self):
    reader = sam_reader.SamReader.from_file(self.bam, self.options)
    with reader:
      self.assertIsInstance(reader.samples, set)
      self.assertEqual(reader.samples, {'NA12878'})

  def test_sam_contigs(self):
    reader = sam_reader.SamReader.from_file(self.bam, self.options)
    with reader:
      self.assertEqual([
          core_pb2.ContigInfo(name='chrM', pos_in_fasta=0, n_bases=16571),
          core_pb2.ContigInfo(name='chr1', pos_in_fasta=1, n_bases=249250621),
          core_pb2.ContigInfo(name='chr2', pos_in_fasta=2, n_bases=243199373),
          core_pb2.ContigInfo(name='chr3', pos_in_fasta=3, n_bases=198022430),
          core_pb2.ContigInfo(name='chr4', pos_in_fasta=4, n_bases=191154276),
          core_pb2.ContigInfo(name='chr5', pos_in_fasta=5, n_bases=180915260),
          core_pb2.ContigInfo(name='chr6', pos_in_fasta=6, n_bases=171115067),
          core_pb2.ContigInfo(name='chr7', pos_in_fasta=7, n_bases=159138663),
          core_pb2.ContigInfo(name='chr8', pos_in_fasta=8, n_bases=146364022),
          core_pb2.ContigInfo(name='chr9', pos_in_fasta=9, n_bases=141213431),
          core_pb2.ContigInfo(name='chr10', pos_in_fasta=10, n_bases=135534747),
          core_pb2.ContigInfo(name='chr11', pos_in_fasta=11, n_bases=135006516),
          core_pb2.ContigInfo(name='chr12', pos_in_fasta=12, n_bases=133851895),
          core_pb2.ContigInfo(name='chr13', pos_in_fasta=13, n_bases=115169878),
          core_pb2.ContigInfo(name='chr14', pos_in_fasta=14, n_bases=107349540),
          core_pb2.ContigInfo(name='chr15', pos_in_fasta=15, n_bases=102531392),
          core_pb2.ContigInfo(name='chr16', pos_in_fasta=16, n_bases=90354753),
          core_pb2.ContigInfo(name='chr17', pos_in_fasta=17, n_bases=81195210),
          core_pb2.ContigInfo(name='chr18', pos_in_fasta=18, n_bases=78077248),
          core_pb2.ContigInfo(name='chr19', pos_in_fasta=19, n_bases=59128983),
          core_pb2.ContigInfo(name='chr20', pos_in_fasta=20, n_bases=63025520),
          core_pb2.ContigInfo(name='chr21', pos_in_fasta=21, n_bases=48129895),
          core_pb2.ContigInfo(name='chr22', pos_in_fasta=22, n_bases=51304566),
          core_pb2.ContigInfo(name='chrX', pos_in_fasta=23, n_bases=155270560),
          core_pb2.ContigInfo(name='chrY', pos_in_fasta=24, n_bases=59373566),
      ], reader.contigs)

  def test_context_manager(self):
    """Test that we can use context manager to do two queries in sequence."""
    reader = sam_reader.SamReader.from_file(self.bam, self.indexed_options)
    region = ranges.parse_literal('chr20:10,000,000-10,000,100')
    with reader:
      with reader.query(region) as query_iterable1:
        self.assertIsNotNone(query_iterable1)
        self.assertIsInstance(query_iterable1, clif_postproc.WrappedCppIterable)
      with reader.query(region) as query_iterable2:
        self.assertIsNotNone(query_iterable2)
        self.assertIsInstance(query_iterable2, clif_postproc.WrappedCppIterable)

  def test_from_file_raises_with_missing_bam(self):
    with self.assertRaisesRegexp(ValueError,
                                 'Not found: Could not open missing.bam'):
      sam_reader.SamReader.from_file('missing.bam', self.options)

  def test_from_file_raises_with_missing_index(self):
    with self.assertRaisesRegexp(ValueError, 'Not found: No index found for'):
      sam_reader.SamReader.from_file(
          test_utils.genomics_core_testdata('unindexed.bam'),
          self.indexed_options)

  def test_ops_on_closed_reader_raise(self):
    reader = sam_reader.SamReader.from_file(self.bam, self.indexed_options)
    with reader:
      pass
    # At this point the reader is closed.
    with self.assertRaisesRegexp(ValueError, 'Cannot Iterate a closed'):
      reader.iterate()
    with self.assertRaisesRegexp(ValueError, 'Cannot Query a closed'):
      reader.query(ranges.parse_literal('chr20:10,000,000-10,000,100'))

  def test_query_on_unindexed_reader_raises(self):
    with sam_reader.SamReader.from_file(self.bam, self.options) as reader:
      with self.assertRaisesRegexp(ValueError, 'Cannot query without an index'):
        reader.query(ranges.parse_literal('chr20:10,000,000-10,000,100'))

  def test_query_raises_with_bad_range(self):
    with sam_reader.SamReader.from_file(self.bam,
                                        self.indexed_options) as reader:
      with self.assertRaisesRegexp(ValueError, 'Unknown reference_name'):
        reader.query(ranges.parse_literal('XXX:1-10'))
      with self.assertRaisesRegexp(ValueError, 'unknown reference interval'):
        reader.query(ranges.parse_literal('chr20:10-5'))

  def test_sam_iterate_raises_on_malformed_record(self):
    malformed = test_utils.genomics_core_testdata('malformed.sam')
    reader = sam_reader.SamReader.from_file(malformed, self.options)
    iterable = iter(reader.iterate())
    self.assertIsNotNone(next(iterable))
    with self.assertRaises(ValueError):
      list(iterable)


if __name__ == '__main__':
  absltest.main()
