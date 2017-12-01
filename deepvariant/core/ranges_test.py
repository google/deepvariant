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
"""Tests for ranges.py."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools



from absl.testing import absltest
from absl.testing import parameterized
import mock

from deepvariant.core import ranges
from deepvariant.core import test_utils
from deepvariant.core.genomics import position_pb2
from deepvariant.core.genomics import variants_pb2
from deepvariant.core.protos import core_pb2

_TEST_BED_REGIONS = [
    ranges.make_range('chr1', 1, 10),
    ranges.make_range('chr2', 20, 30),
    ranges.make_range('chr2', 40, 60),
    ranges.make_range('chr3', 80, 90),
]


class RangesTests(parameterized.TestCase):

  def test_ranges_overlaps(self):

    def check_overlaps(chr1, start1, end1, chr2, start2, end2, expected):
      i1 = ranges.make_range(chr1, start1, end1)
      i2 = ranges.make_range(chr2, start2, end2)
      self.assertEquals(ranges.ranges_overlap(i1, i2), expected)
      self.assertEquals(ranges.ranges_overlap(i2, i1), expected)

    check_overlaps('chr1', 0, 3, 'chr1', 4, 10, False)
    check_overlaps('chr1', 0, 3, 'chr1', 3, 10, False)
    check_overlaps('chr1', 0, 3, 'chr1', 2, 10, True)
    check_overlaps('chr1', 0, 3, 'chr1', 1, 10, True)
    check_overlaps('chr1', 0, 3, 'chr1', 0, 10, True)
    check_overlaps('chr1', 0, 3, 'chr1', 0, 1, True)
    check_overlaps('chr1', 0, 3, 'chr1', 0, 2, True)
    check_overlaps('chr1', 0, 3, 'chr1', 0, 3, True)
    check_overlaps('chr1', 0, 3, 'chr1', 1, 2, True)
    check_overlaps('chr1', 0, 3, 'chr1', 1, 3, True)
    check_overlaps('chr1', 0, 3, 'chr1', 2, 3, True)
    check_overlaps('chr1', 0, 3, 'chr1', 3, 3, False)
    check_overlaps('chr1', 1, 3, 'chr1', 0, 4, True)
    check_overlaps('chr1', 1, 3, 'chr1', 1, 4, True)

  def test_detector_no_ranges(self):
    range_set = ranges.RangeSet()
    # don't have any ranges by default
    self.assertEqual(bool(range_set), False)
    # make sure we can call overlaps without any ranges
    self.assertFalse(range_set.overlaps('chr1', 10))

  def test_from_regions_not_empty(self):
    literals = ['chr1', 'chr2:10-20']
    contig_map = {
        'chr1': core_pb2.ContigInfo(name='chr1', n_bases=10),
        'chr2': core_pb2.ContigInfo(name='chr2', n_bases=100),
    }
    self.assertItemsEqual(
        [ranges.make_range('chr1', 0, 10),
         ranges.make_range('chr2', 9, 20)],
        ranges.RangeSet.from_regions(literals, contig_map))

  def test_from_regions_empty_literals(self):
    range_set = ranges.RangeSet.from_regions([], contig_map=None)
    # The set is empty.
    self.assertItemsEqual([], range_set)
    self.assertFalse(range_set)

  @parameterized.parameters(
      # Overlapping intervals get merged.
      (['1:1-5', '1:3-8'], ['1:1-8']),
      (['1:1-5', '1:3-8', '1:6-9'], ['1:1-9']),
      # Adjacent intervals are merged.
      (['1:1-5', '1:5-8'], ['1:1-8']),
      (['1:1-5', '1:5-8', '1:8-10'], ['1:1-10']),
      # Sanity check that non-overlapping aren't merged.
      (['1:1-5', '1:6-8'], ['1:1-5', '1:6-8']),
  )
  def test_overlapping_and_adjacent_ranges_are_merged(self, regions, expected):
    self.assertCountEqual(
        ranges.RangeSet.from_regions(expected),
        ranges.RangeSet.from_regions(regions))

  def test_detector_ranges(self):
    test_ranges = [
        ranges.make_range('chr1', 0, 5),
        ranges.make_range('chr1', 8, 10),
        ranges.make_range('chr1', 12, 13),
        ranges.make_range('chr2', 2, 5),
    ]
    range_set = ranges.RangeSet(test_ranges)
    self.assertEqual(bool(range_set), True)
    self.assertEqual(len(range_set), 4)

    self.assertEqual(range_set.overlaps('chr1', 0), True)
    self.assertEqual(range_set.overlaps('chr1', 1), True)
    self.assertEqual(range_set.overlaps('chr1', 2), True)
    self.assertEqual(range_set.overlaps('chr1', 3), True)
    self.assertEqual(range_set.overlaps('chr1', 4), True)
    self.assertEqual(range_set.overlaps('chr1', 5), False)
    self.assertEqual(range_set.overlaps('chr1', 6), False)
    self.assertEqual(range_set.overlaps('chr1', 7), False)
    self.assertEqual(range_set.overlaps('chr1', 8), True)
    self.assertEqual(range_set.overlaps('chr1', 9), True)
    self.assertEqual(range_set.overlaps('chr1', 10), False)
    self.assertEqual(range_set.overlaps('chr1', 11), False)
    self.assertEqual(range_set.overlaps('chr1', 12), True)
    self.assertEqual(range_set.overlaps('chr1', 13), False)
    self.assertEqual(range_set.overlaps('chr1', 100), False)
    self.assertEqual(range_set.overlaps('chr1', 1000), False)
    self.assertEqual(range_set.overlaps('chr2', 0), False)
    self.assertEqual(range_set.overlaps('chr2', 1), False)
    self.assertEqual(range_set.overlaps('chr2', 2), True)
    self.assertEqual(range_set.overlaps('chr2', 3), True)
    self.assertEqual(range_set.overlaps('chr2', 4), True)
    self.assertEqual(range_set.overlaps('chr2', 5), False)
    self.assertEqual(range_set.overlaps('chr2', 6), False)
    self.assertEqual(range_set.overlaps('chr3', 3), False)

  def test_overlaps_variant_with_ranges(self):
    variant = variants_pb2.Variant(reference_name='chr2', start=10, end=11)
    range_set = ranges.RangeSet([ranges.make_range('chr1', 0, 5)])
    with mock.patch.object(range_set, 'overlaps') as mock_overlaps:
      mock_overlaps.return_value = True
      self.assertEqual(range_set.variant_overlaps(variant), True)
      mock_overlaps.assert_called_once_with('chr2', 10)

  def test_overlaps_variant_empty_range(self):
    variant = variants_pb2.Variant(reference_name='chr2', start=10, end=11)
    empty_set = ranges.RangeSet()
    self.assertEqual(
        empty_set.variant_overlaps(variant, empty_set_return_value='foo'),
        'foo')

  @parameterized.parameters(
      (ranges.make_range('1', 10, 50), '1', 9, False),
      (ranges.make_range('1', 10, 50), '1', 10, True),
      (ranges.make_range('1', 10, 50), '2', 10, False),
      (ranges.make_range('1', 10, 50), '1', 30, True),
      (ranges.make_range('1', 10, 50), '2', 30, False),
      (ranges.make_range('1', 10, 50), '1', 49, True),
      (ranges.make_range('1', 10, 50), '1', 50, False),
      (ranges.make_range('1', 10, 50), '1', 51, False),
  )
  def test_position_overlaps(self, interval, chrom, pos, expected):
    self.assertEqual(ranges.position_overlaps(chrom, pos, interval), expected)

  def test_make_position(self):
    self.assertEqual(
        ranges.make_position('chr1', 10),
        position_pb2.Position(
            reference_name='chr1', position=10, reverse_strand=False))
    self.assertEqual(
        ranges.make_position('chr2', 100, reverse_strand=True),
        position_pb2.Position(
            reference_name='chr2', position=100, reverse_strand=True))

  def test_make_range(self):
    interval = ranges.make_range('chr1', 1, 10)
    self.assertEqual(interval.reference_name, 'chr1')
    self.assertEqual(interval.start, 1)
    self.assertEqual(interval.end, 10)

  def test_to_literal(self):
    self.assertEqual(
        ranges.to_literal(ranges.make_range('chr1', 0, 20)), 'chr1:1-20')

  @parameterized.parameters(['chr1', '1', 'MT', 'chrM', 'chrX', 'X', 'Y'])
  def test_parse_literal_chromosomes(self, chrom):
    self.assertEqual(
        ranges.parse_literal(chrom + ':1-20'), ranges.make_range(chrom, 0, 20))

  @parameterized.parameters(
      ('chr1:{}-{}'.format(start_str, end_str), start_val, end_val)
      for start_str, start_val in [('12', 11), ('1,234', 1233)]
      for end_str, end_val in [('56789', 56789), ('56,789', 56789)])
  def test_parse_literal_numerics(self, literal, start_val, end_val):
    self.assertEqual(
        ranges.parse_literal(literal),
        ranges.make_range('chr1', start_val, end_val))

  def test_parse_literal_one_bp(self):
    self.assertEqual(
        ranges.parse_literal('1:10'), ranges.make_range('1', 9, 10))
    self.assertEqual(
        ranges.parse_literal('1:100'), ranges.make_range('1', 99, 100))
    self.assertEqual(
        ranges.parse_literal('1:1,000'), ranges.make_range('1', 999, 1000))

  @parameterized.parameters(['x', 'chr1', 'chr1:', 'chr1:10-', 'chr1:-1-10'])
  def test_parse_literal_bad(self, bad_literal):
    with self.assertRaises(ValueError):
      ranges.parse_literal(bad_literal)

  def test_from_bed(self):
    source = test_utils.genomics_core_testdata('test.bed')
    self.assertCountEqual([
        ranges.make_range('chr1', 1, 10),
        ranges.make_range('chr2', 20, 30),
        ranges.make_range('chr2', 40, 60),
        ranges.make_range('chr3', 80, 90),
    ], ranges.RangeSet.from_bed(source))

  @parameterized.parameters(
      dict(regions=[], expected=[]),
      dict(regions=['chr1:10-20'], expected=[ranges.make_range('chr1', 9, 20)]),
      dict(regions=['test.bed'], expected=_TEST_BED_REGIONS),
      dict(
          regions=['test.bed', 'test.bed'],
          expected=_TEST_BED_REGIONS + _TEST_BED_REGIONS),
      dict(
          regions=['chr1:10-20', 'test.bed'],
          expected=[ranges.make_range('chr1', 9, 20)] + _TEST_BED_REGIONS),
      dict(
          regions=['test.bed', 'chr1:10-20'],
          expected=_TEST_BED_REGIONS + [ranges.make_range('chr1', 9, 20)]),
      dict(
          regions=['chr1:9-19', 'test.bed', 'chr1:10-20'],
          expected=([ranges.make_range('chr1', 8, 19)] + _TEST_BED_REGIONS +
                    [ranges.make_range('chr1', 9, 20)])),
  )
  def test_from_regions(self, regions, expected):
    # For convenience we allow 'test.bed' in our regions but the actual file
    # path is in our testdata directory.
    for i in range(len(regions)):
      if regions[i] == 'test.bed':
        regions[i] = test_utils.genomics_core_testdata('test.bed')

    self.assertEqual(list(ranges.from_regions(regions)), expected)

  @parameterized.parameters(
      # Intersection with 1, 2, 3 identical RangeSets produces the original set.
      ([['1:1-10']], ['1:1-10']),
      ([['1:1-10'], ['1:1-10']], ['1:1-10']),
      ([['1:1-10'], ['1:1-10'], ['1:1-10']], ['1:1-10']),
      # Test some simple overlap configurations.
      ([['1:1-10'], ['1:11-15']], []),
      ([['1:1-10'], ['1:10-15']], ['1:10']),
      ([['1:1-10'], ['1:9-15']], ['1:9-10']),
      ([['1:5-10'], ['1:1-15']], ['1:5-10']),
      ([['1:5-10'], ['1:1-4']], []),
      ([['1:5-10'], ['1:1-5']], ['1:5']),
      # Check cutting a single interval into multiple pieces.
      ([['1:5-15'], ['1:6-8', '1:10-12']], ['1:6-8', '1:10-12']),
      ([['1:5-15'], ['1:3-8', '1:10-12']], ['1:5-8', '1:10-12']),
      ([['1:5-15'], ['1:3-8', '1:10-20']], ['1:5-8', '1:10-15']),
      # We have multiple overlapping intervals; make sure we merge intervals.
      ([['1:5-15'], ['1:3-8', '1:6-10']], ['1:5-10']),
      ([['1:5-15'], ['1:3-8', '1:6-10', '1:13']], ['1:5-10', '1:13']),
      # Check that multiple intervals work.
      ([['1:5-15', '1:20-25'], ['1:3-8', '1:16-23']], ['1:5-8', '1:20-23']),
      ([['1:5-15', '1:20-25'], ['1:3-8', '1:50-60']], ['1:5-8']),
      ([['1:5-15', '1:20-25'], ['1:3-4', '1:16-23']], ['1:20-23']),
      # Check that multiple sets can be intersected.
      ([['1:10-20'], ['1:5-15']], ['1:10-15']),
      ([['1:10-20'], ['1:5-15'], ['1:13-30']], ['1:13-15']),
      ([['1:10-20'], ['1:5-15'], ['1:25-30']], []),
      # Check that different chromosomes are kept separate.
      ([['1:10-20'], ['2:10-20']], []),
      ([['1:10-20', '2:11-14'], ['1:11-14']], ['1:11-14']),
      ([['1:10-20', '2:11-14'], ['2:10-20']], ['2:11-14']),
  )
  def test_intersection(self, regions, expected):
    regions_list = [ranges.RangeSet.from_regions(r) for r in regions]
    copies = [ranges.RangeSet(rs) for rs in regions_list]

    # Check that the intersection is as expected.
    self.assertCountEqual(
        ranges.RangeSet.from_regions(expected),
        regions_list[0].intersection(*regions_list[1:]))

    # Check that no one was modified.
    for pre, post in zip(copies, regions_list):
      self.assertCountEqual(pre, post)

  @parameterized.parameters(('chr1', ranges.make_range('chr1', 0, 10)),
                            ('chr2', ranges.make_range('chr2', 0, 5)))
  def test_parse_literal_with_contig_map(self, contig_name, expected):
    contig_map = {
        'chr1': core_pb2.ContigInfo(name='chr1', n_bases=10),
        'chr2': core_pb2.ContigInfo(name='chr2', n_bases=5),
    }
    self.assertEqual(
        ranges.parse_literal(contig_name, contig_map=contig_map), expected)

  @parameterized.parameters(['x', 'chr1:', 'chr1:10-', 'chr1:-1-10'])
  def test_parse_literal_with_contig_map_and_bad_input_raises_exception(
      self, bad_literal):
    with self.assertRaises(ValueError):
      ranges.parse_literal(
          bad_literal,
          contig_map={
              'chr1': core_pb2.ContigInfo(name='chr1', n_bases=10)
          })

  def test_from_contigs(self):
    contigs = [
        core_pb2.ContigInfo(name='chr1', n_bases=10),
        core_pb2.ContigInfo(name='chr2', n_bases=5),
    ]
    self.assertCountEqual([
        ranges.make_range('chr1', 0, 10),
        ranges.make_range('chr2', 0, 5),
    ], ranges.RangeSet.from_contigs(contigs))

  @parameterized.parameters(
      # Chop our contigs into 50 bp pieces.
      (50, [('chrM', 0, 50), ('chrM', 50, 100), ('chr1', 0, 50),
            ('chr1', 50, 76), ('chr2', 0, 50), ('chr2', 50, 100),
            ('chr2', 100, 121)]),
      # Chop our contigs in 120 bp pieces, leaving a 1 bp fragment in chr2.
      (120, [('chrM', 0, 100), ('chr1', 0, 76), ('chr2', 0, 120),
             ('chr2', 120, 121)]),
      # A 500 max size spans each of our contigs fully.
      (500, [('chrM', 0, 100), ('chr1', 0, 76), ('chr2', 0, 121)]),
  )
  def test_partitions(self, interval_size, expected):
    rangeset = ranges.RangeSet([
        ranges.make_range('chrM', 0, 100),
        ranges.make_range('chr1', 0, 76),
        ranges.make_range('chr2', 0, 121),
    ])
    self.assertCountEqual([ranges.make_range(*args) for args in expected],
                          rangeset.partition(interval_size))

  def test_partitions_bad_interval_size_raises(self):
    # list() is necessary to force the generator to execute.
    with self.assertRaisesRegexp(ValueError, 'max_size'):
      list(ranges.RangeSet([ranges.make_range('chrM', 0, 100)]).partition(-10))
    with self.assertRaisesRegexp(ValueError, 'max_size'):
      list(ranges.RangeSet([ranges.make_range('chrM', 0, 100)]).partition(0))

  @parameterized.parameters(
      (10, [('1', 0, 10), ('1', 20, 30), ('1', 30, 40), ('1', 45, 50)]),
      (7, [('1', 0, 7), ('1', 7, 10), ('1', 20, 27), ('1', 27, 34),
           ('1', 34, 40), ('1', 45, 50)]),
      (50, [('1', 0, 10), ('1', 20, 40), ('1', 45, 50)]),
  )
  def test_partition_of_multiple_intervals(self, interval_size, expected):
    rangeset = ranges.RangeSet([
        ranges.make_range('1', 0, 10),
        ranges.make_range('1', 20, 40),
        ranges.make_range('1', 45, 50),
    ])
    self.assertCountEqual([ranges.make_range(*args) for args in expected],
                          rangeset.partition(interval_size))

  def test_unknown_filetype(self):
    with self.assertRaises(ValueError):
      ranges.parse_lines([], file_format='png')

  def test_bed_parser(self):
    data = [
        'chr20\t61724611\t61725646',
        'chr20\t61304163\t61305182',
        'chr20\t61286467\t61286789',
    ]
    self.assertEqual(
        list(ranges.parse_lines(data, 'bed')), [
            ranges.make_range('chr20', 61724611, 61725646),
            ranges.make_range('chr20', 61304163, 61305182),
            ranges.make_range('chr20', 61286467, 61286789),
        ])

  def test_bedpe_parser(self):
    # pylint: disable=line-too-long
    data = [
        'chr20\t25763416\t25765517\tchr20\t25825181\t25826882\tP2_PM_20_1549\t63266\t+\tTYPE:DELETION',
        'chr20\t25972820\t25972991\tchr20\t26045347\t26045538\tP2_PM_20_696\t72548\t+\tTYPE:DELETION',
        'chr20\t23719873\t23721974\tchr20\t23794822\t23796523\tP2_PM_20_1548\t76450\t+\tTYPE:DELETION',
    ]
    self.assertEqual(
        list(ranges.parse_lines(data, 'bedpe')), [
            ranges.make_range('chr20', 25763416, 25826882),
            ranges.make_range('chr20', 25972820, 26045538),
            ranges.make_range('chr20', 23719873, 23796523),
        ])

  def test_bedpe_parser_skips_cross_chr_events(self):
    # pylint: disable=line-too-long
    data = [
        'chr20\t25763416\t25765517\tchr21\t25825181\t25826882\tP2_PM_20_1549\t63266\t+\tTYPE:DELETION',
        'chr20\t25972820\t25972991\tchr20\t26045347\t26045538\tP2_PM_20_696\t72548\t+\tTYPE:DELETION',
        'chr20\t23719873\t23721974\tchr20\t23794822\t23796523\tP2_PM_20_1548\t76450\t+\tTYPE:DELETION',
    ]
    self.assertEqual(
        list(ranges.parse_lines(data, 'bedpe')), [
            ranges.make_range('chr20', 25972820, 26045538),
            ranges.make_range('chr20', 23719873, 23796523),
        ])

  def test_contigs_n_bases(self):
    c1 = core_pb2.ContigInfo(name='c', n_bases=100, pos_in_fasta=0)
    c2 = core_pb2.ContigInfo(name='a', n_bases=50, pos_in_fasta=1)
    c3 = core_pb2.ContigInfo(name='b', n_bases=25, pos_in_fasta=2)
    self.assertEqual(100, ranges.contigs_n_bases([c1]))
    self.assertEqual(50, ranges.contigs_n_bases([c2]))
    self.assertEqual(25, ranges.contigs_n_bases([c3]))
    self.assertEqual(150, ranges.contigs_n_bases([c1, c2]))
    self.assertEqual(125, ranges.contigs_n_bases([c1, c3]))
    self.assertEqual(175, ranges.contigs_n_bases([c1, c2, c3]))

  def test_sort_ranges(self):
    contigs = [
        core_pb2.ContigInfo(name='c', n_bases=100, pos_in_fasta=0),
        core_pb2.ContigInfo(name='a', n_bases=76, pos_in_fasta=1),
        core_pb2.ContigInfo(name='b', n_bases=121, pos_in_fasta=2),
    ]
    unsorted = ranges.parse_literals(
        ['a:10', 'c:20', 'b:30', 'b:10-15', 'b:10', 'a:5'])

    # Without contigs we sort the contigs by name lexicographically.
    self.assertEqual(
        ranges.parse_literals(
            ['a:5', 'a:10', 'b:10', 'b:10-15', 'b:30', 'c:20']),
        ranges.sorted_ranges(unsorted))

    # With contigs we sort by the position of the contigs themselves.
    self.assertEqual(
        ranges.parse_literals(
            ['c:20', 'a:5', 'a:10', 'b:10', 'b:10-15', 'b:30']),
        ranges.sorted_ranges(unsorted, contigs))

  @parameterized.parameters(
      (ranges.make_range('1', 0, 10), ranges.make_range('2', 0, 10), 0),
      (ranges.make_range('1', 0, 10), ranges.make_range('1', 10, 20), 0),
      (ranges.make_range('1', 0, 10), ranges.make_range('1', 100, 200), 0),
      (ranges.make_range('1', 10, 10), ranges.make_range('1', 0, 20), 0),
      (ranges.make_range('1', 0, 100), ranges.make_range('1', 50, 99), 49),
      # Check that the overlap handles a few key edge cases.
      (ranges.make_range('1', 0, 10), ranges.make_range('1', 0, 1), 1),
      (ranges.make_range('1', 0, 10), ranges.make_range('1', 0, 2), 2),
      (ranges.make_range('1', 1, 10), ranges.make_range('1', 0, 1), 0),
  )
  def test_overlap_len(self, region_1, region_2, expected_overlap):
    """Test ReadAssigner.overlap_len()."""
    self.assertEqual(expected_overlap, ranges.overlap_len(region_1, region_2))
    self.assertEqual(expected_overlap, ranges.overlap_len(region_2, region_1))

  @parameterized.parameters(
      # No search_regions produces None.
      dict(
          query_range=ranges.make_range('1', 20, 30),
          search_ranges=[],
          expected=None),

      # Read overlaps with none of the ranges returns None.
      dict(
          query_range=ranges.make_range('1', 20, 30),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('1', 5, 10)
          ],
          expected=None),

      # Read has longer overlap with the first range.
      dict(
          query_range=ranges.make_range('1', 4, 10),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('1', 5, 10)
          ],
          expected=0),

      # Read has longer overlap with the second range.
      dict(
          query_range=ranges.make_range('1', 9, 20),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('1', 5, 15)
          ],
          expected=1),

      # Read has the maximum overlap with the third range.
      dict(
          query_range=ranges.make_range('1', 9, 20),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('1', 0, 15),
              ranges.make_range('1', 5, 20)
          ],
          expected=2),

      # Read has the maximum overlap with the middle range.
      dict(
          query_range=ranges.make_range('1', 5, 13),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('1', 0, 15),
              ranges.make_range('1', 10, 20)
          ],
          expected=1),

      # Read has a different reference_name with other ranges.
      dict(
          query_range=ranges.make_range('2', 0, 10),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('2', 5, 15),
              ranges.make_range('3', 0, 10)
          ],
          expected=1),

      # Read has equal overlap in two ranges.
      dict(
          query_range=ranges.make_range('1', 5, 15),
          search_ranges=[
              ranges.make_range('1', 0, 10),
              ranges.make_range('1', 10, 20),
              ranges.make_range('1', 12, 20)
          ],
          expected=0),
  )
  def test_find_max_overlapping(self, query_range, search_ranges, expected):
    actual = ranges.find_max_overlapping(query_range, search_ranges)
    self.assertEqual(expected, actual)

  def test_find_max_overlapping_allows_unordered_search_ranges(self):
    query_range = ranges.make_range('1', 4, 12)
    search_ranges = [
        ranges.make_range('1', 0, 10),
        ranges.make_range('1', 10, 20),
        ranges.make_range('1', 12, 20)
    ]
    max_overlapping_range = search_ranges[0]

    for permutated_ranges in itertools.permutations(search_ranges):
      self.assertEqual(
          permutated_ranges.index(max_overlapping_range),
          ranges.find_max_overlapping(query_range, permutated_ranges))

  def test_find_max_overlapping_returns_least_index(self):
    query_range = ranges.make_range('1', 0, 10)
    search_ranges = [
        ranges.make_range('1', 0, 5),
        ranges.make_range('1', 5, 10)
    ]

    for to_search in [search_ranges, list(reversed(search_ranges))]:
      self.assertEqual(0, ranges.find_max_overlapping(query_range, to_search))


if __name__ == '__main__':
  absltest.main()
