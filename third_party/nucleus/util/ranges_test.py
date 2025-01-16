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
"""Tests for ranges.py."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools

from absl.testing import absltest
from absl.testing import parameterized
import mock

from third_party.nucleus.protos import position_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges

_TEST_BED_REGIONS = [
    ranges.make_range('chr1', 1, 10),
    ranges.make_range('chr2', 20, 30),
    ranges.make_range('chr2', 40, 60),
    ranges.make_range('chr3', 80, 90),
]

_TEST_CONTIGS = [
    reference_pb2.ContigInfo(name='chr1', n_bases=10, pos_in_fasta=0),
    reference_pb2.ContigInfo(name='chr2', n_bases=100, pos_in_fasta=1),
    reference_pb2.ContigInfo(name='chr3', n_bases=500, pos_in_fasta=2),
]


class RangesTests(parameterized.TestCase):

  def test_ranges_overlaps(self):

    def check_overlaps(chr1, start1, end1, chr2, start2, end2, expected):
      i1 = ranges.make_range(chr1, start1, end1)
      i2 = ranges.make_range(chr2, start2, end2)
      self.assertEqual(ranges.ranges_overlap(i1, i2), expected)
      self.assertEqual(ranges.ranges_overlap(i2, i1), expected)

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
    self.assertItemsEqual(
        [ranges.make_range('chr1', 0, 10),
         ranges.make_range('chr2', 9, 20)],
        ranges.RangeSet.from_regions(
            literals, ranges.contigs_dict(_TEST_CONTIGS)))

  def test_from_regions_empty_literals(self):
    range_set = ranges.RangeSet.from_regions([])
    # The set is empty.
    self.assertItemsEqual([], range_set)
    self.assertFalse(range_set)

  def test_unrecognized_contig_triggers_exception(self):
    with self.assertRaises(ValueError):
      _ = ranges.RangeSet([ranges.make_range('bogus_chromosome', 1, 10)],
                          _TEST_CONTIGS)

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

  def test_envelops(self):
    start_ix = 5
    end_ix = 10
    start_ix2 = end_ix + 1
    end_ix2 = end_ix + 5
    range_set = ranges.RangeSet([
        ranges.make_range('chr1', start_ix, end_ix),
        ranges.make_range('chr1', start_ix2, end_ix2)
    ])

    # No start position before the first start range is enveloped.
    for i in range(start_ix):
      self.assertFalse(range_set.envelops('chr1', i, start_ix + 1))

    # All regions within a single record are enveloped.
    for six in range(start_ix, end_ix):
      for eix in range(six, end_ix + 1):
        self.assertTrue(
            range_set.envelops('chr1', six, eix),
            'chr1 {} {} not enveloped'.format(six, eix))

    # Bridging across two ranges is not enveloped.
    for six in range(start_ix, end_ix):
      for eix in range(start_ix2, end_ix2 + 1):
        self.assertFalse(range_set.envelops('chr1', six, eix))

    # Other chromosome is not spanned.
    self.assertFalse(range_set.envelops('chr2', start_ix, start_ix + 1))

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
    with self.assertRaisesRegex(ValueError, bad_literal):
      ranges.parse_literal(bad_literal)

  @parameterized.parameters('test.bed', 'test.bed.gz')
  def test_from_bed(self, bed_filename):
    source = test_utils.genomics_core_testdata(bed_filename)
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

    # Check that the intersection is as expected even if we do it in a different
    # direction.
    self.assertCountEqual(
        ranges.RangeSet.from_regions(expected),
        regions_list[-1].intersection(*regions_list[:-1]))

    # Check that no one was modified.
    for pre, post in zip(copies, regions_list):
      self.assertCountEqual(pre, post)

  @parameterized.parameters(
      dict(lhs=['1:1-100'], rhs=['1:10-20'], expected=['1:1-9', '1:21-100']),
      dict(lhs=['1:1-100'], rhs=[], expected=['1:1-100']),
      dict(lhs=['1:1-100', '2:1-10'], rhs=['2:1-100'], expected=['1:1-100']),
      dict(
          lhs=['1:1-100'],
          rhs=['1:10-20', '1:15-30'],
          expected=['1:1-9', '1:31-100']),
      dict(
          lhs=['1:1-100'],
          rhs=['1:10-20', '1:30-40'],
          expected=['1:1-9', '1:21-29', '1:41-100']),
      # Excluding regions not in lhs has no impact.
      dict(lhs=['1:1-100'], rhs=['2:1-100'], expected=['1:1-100']),
      # Check that excluding the whole region results in an empty RangeSet.
      dict(lhs=['1:1-100'], rhs=['1:1-100'], expected=[]),
      # An empty tree remains empty.
      dict(lhs=[], rhs=['1:1-100'], expected=[]),
  )
  def test_exclude_regions(self, lhs, rhs, expected):
    lhs = ranges.RangeSet.from_regions(lhs)
    rhs = ranges.RangeSet.from_regions(rhs)
    # Mutating operation returns None.
    self.assertIsNone(lhs.exclude_regions(rhs))
    self.assertCountEqual(ranges.RangeSet.from_regions(expected), lhs)

  @parameterized.parameters(('chr1', ranges.make_range('chr1', 0, 10)),
                            ('chr2', ranges.make_range('chr2', 0, 5)))
  def test_parse_literal_with_contig_map(self, contig_name, expected):
    contig_map = {
        'chr1': reference_pb2.ContigInfo(name='chr1', n_bases=10),
        'chr2': reference_pb2.ContigInfo(name='chr2', n_bases=5),
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
              'chr1': reference_pb2.ContigInfo(name='chr1', n_bases=10)
          })

  def test_from_contigs(self):
    contigs = [
        reference_pb2.ContigInfo(name='chr1', n_bases=10),
        reference_pb2.ContigInfo(name='chr2', n_bases=5),
    ]
    self.assertCountEqual([
        ranges.make_range('chr1', 0, 10),
        ranges.make_range('chr2', 0, 5),
    ], ranges.RangeSet.from_contigs(contigs))

  @parameterized.parameters(
      # Chop our contigs into 50 bp pieces.
      (50, [('chr1', 0, 50), ('chr1', 50, 76), ('chr2', 0, 50),
            ('chr2', 50, 100), ('chr2', 100, 121), ('chrM', 0, 50),
            ('chrM', 50, 100)]),
      # Chop our contigs in 120 bp pieces, leaving a 1 bp fragment in chr2.
      (120, [('chr1', 0, 76), ('chr2', 0, 120), ('chr2', 120, 121),
             ('chrM', 0, 100)]),
      # A 500 max size spans each of our contigs fully.
      (500, [('chr1', 0, 76), ('chr2', 0, 121), ('chrM', 0, 100)]),
  )
  def test_partitions(self, interval_size, expected):
    rangeset = ranges.RangeSet([
        ranges.make_range('chrM', 0, 100),
        ranges.make_range('chr1', 0, 76),
        ranges.make_range('chr2', 0, 121),
    ])
    self.assertEqual([ranges.make_range(*args) for args in expected],
                     list(rangeset.partition(interval_size)))

  def test_partitions_bad_interval_size_raises(self):
    # list() is necessary to force the generator to execute.
    with self.assertRaisesRegex(ValueError, 'max_size'):
      list(ranges.RangeSet([ranges.make_range('chrM', 0, 100)]).partition(-10))
    with self.assertRaisesRegex(ValueError, 'max_size'):
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

  def test_bed_parser(self):
    test_bed_path = test_utils.test_tmpfile(
        'test_bed_parser.bed', '\n'.join([
            'chr20\t61724611\t61725646', 'chr20\t61304163\t61305182',
            'chr20\t61286467\t61286789'
        ]))
    self.assertEqual(
        list(ranges.bed_parser(test_bed_path)), [
            ranges.make_range('chr20', 61724611, 61725646),
            ranges.make_range('chr20', 61304163, 61305182),
            ranges.make_range('chr20', 61286467, 61286789),
        ])

  def test_bedpe_parser(self):
    # pylint: disable=line-too-long
    data = '\n'.join([
        'chr20\t25763416\t25765517\tchr20\t25825181\t25826882\tP2_PM_20_1549\t63266\t+\tTYPE:DELETION',
        'chr20\t25972820\t25972991\tchr20\t26045347\t26045538\tP2_PM_20_696\t72548\t+\tTYPE:DELETION',
        'chr20\t23719873\t23721974\tchr20\t23794822\t23796523\tP2_PM_20_1548\t76450\t+\tTYPE:DELETION',
    ])
    test_bedpe_path = test_utils.test_tmpfile('test_bedpe_parser.bedpe', data)
    self.assertEqual(
        list(ranges.bedpe_parser(test_bedpe_path)), [
            ranges.make_range('chr20', 25763416, 25826882),
            ranges.make_range('chr20', 25972820, 26045538),
            ranges.make_range('chr20', 23719873, 23796523),
        ])

  def test_bedpe_parser_skips_cross_chr_events(self):
    # pylint: disable=line-too-long
    data = '\n'.join([
        'chr20\t25763416\t25765517\tchr21\t25825181\t25826882\tP2_PM_20_1549\t63266\t+\tTYPE:DELETION',
        'chr20\t25972820\t25972991\tchr20\t26045347\t26045538\tP2_PM_20_696\t72548\t+\tTYPE:DELETION',
        'chr20\t23719873\t23721974\tchr20\t23794822\t23796523\tP2_PM_20_1548\t76450\t+\tTYPE:DELETION',
    ])
    test_bedpe_path = test_utils.test_tmpfile('test_bedpe_parser2.bedpe', data)
    self.assertEqual(
        list(ranges.bedpe_parser(test_bedpe_path)), [
            ranges.make_range('chr20', 25972820, 26045538),
            ranges.make_range('chr20', 23719873, 23796523),
        ])

  def test_contigs_n_bases(self):
    c1 = reference_pb2.ContigInfo(name='c', n_bases=100, pos_in_fasta=0)
    c2 = reference_pb2.ContigInfo(name='a', n_bases=50, pos_in_fasta=1)
    c3 = reference_pb2.ContigInfo(name='b', n_bases=25, pos_in_fasta=2)
    self.assertEqual(100, ranges.contigs_n_bases([c1]))
    self.assertEqual(50, ranges.contigs_n_bases([c2]))
    self.assertEqual(25, ranges.contigs_n_bases([c3]))
    self.assertEqual(150, ranges.contigs_n_bases([c1, c2]))
    self.assertEqual(125, ranges.contigs_n_bases([c1, c3]))
    self.assertEqual(175, ranges.contigs_n_bases([c1, c2, c3]))

  def test_rangeset_iteration_order(self):
    contigs = [
        reference_pb2.ContigInfo(name='c', n_bases=100, pos_in_fasta=0),
        reference_pb2.ContigInfo(name='b', n_bases=121, pos_in_fasta=2),
        reference_pb2.ContigInfo(name='a', n_bases=76, pos_in_fasta=1),
    ]
    unsorted = ranges.parse_literals(
        ['a:10', 'c:20', 'b:30', 'b:10-15', 'a:5'])

    # Iteration order over a RangeSet instantiated with a contigs list is
    # determined by pos_in_fasta, start, end.
    range_set_with_contigs = ranges.RangeSet(unsorted, contigs)
    self.assertEqual(
        ranges.parse_literals(
            ['c:20', 'a:5', 'a:10', 'b:10-15', 'b:30']),
        [range_ for range_ in range_set_with_contigs])

    # For a RangeSet instantiated *without* a contig map, the iteration order
    # is determined by reference_name, start, end.
    range_set_no_contigs = ranges.RangeSet(unsorted)
    self.assertEqual(
        ranges.parse_literals(
            ['a:5', 'a:10', 'b:10-15', 'b:30', 'c:20']),
        [range_ for range_ in range_set_no_contigs])

  def test_sort_ranges(self):
    contigs = [
        reference_pb2.ContigInfo(name='c', n_bases=100, pos_in_fasta=0),
        reference_pb2.ContigInfo(name='a', n_bases=76, pos_in_fasta=1),
        reference_pb2.ContigInfo(name='b', n_bases=121, pos_in_fasta=2),
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

  @parameterized.parameters(
      dict(
          regions=[
              ranges.make_range('1', 1, 10),
          ],
          expected_span=ranges.make_range('1', 1, 10),
      ),
      dict(
          regions=[
              ranges.make_range('1', 1, 10),
              ranges.make_range('1', 10, 100),
          ],
          expected_span=ranges.make_range('1', 1, 100),
      ),
      dict(
          regions=[
              ranges.make_range('1', 1, 10),
              ranges.make_range('1', 10, 100),
              ranges.make_range('1', 2, 20),
          ],
          expected_span=ranges.make_range('1', 1, 100),
      ),
      # potential edge cases:
      # same start, different ends.
      dict(
          regions=[
              ranges.make_range('1', 1, 10),
              ranges.make_range('1', 1, 100),
          ],
          expected_span=ranges.make_range('1', 1, 100),
      ),
      # same end, different starts.
      dict(
          regions=[
              ranges.make_range('1', 1, 10),
              ranges.make_range('1', 2, 10),
          ],
          expected_span=ranges.make_range('1', 1, 10),
      ),
  )
  def test_span_computes_span_correctly(self, regions, expected_span):
    for permutation in itertools.permutations(regions, len(regions)):
      self.assertEqual(expected_span, ranges.span(permutation))

  @parameterized.parameters(
      dict(regions=[], regexp='empty'),
      dict(
          regions=[
              ranges.make_range('1', 0, 2),
              ranges.make_range('2', 0, 2),
          ],
          regexp='regions must be all on the same contig'),
      dict(
          regions=[
              ranges.make_range('1', 0, 2),
              ranges.make_range('1', 0, 3),
              ranges.make_range('2', 0, 2),
          ],
          regexp='regions must be all on the same contig'),
  )
  def test_span_raises_on_bad_input(self, regions, regexp):
    with self.assertRaisesRegex(ValueError, regexp):
      ranges.span(regions)

  @parameterized.parameters(
      dict(
          region=ranges.make_range('1', 10, 20),
          n_bp=n_bp,
          contig_map=None,
          expected=ranges.make_range('1', 10 - n_bp, 20 + n_bp),
      ) for n_bp in range(10))
  def test_expand_is_correct(self, region, n_bp, contig_map, expected):
    self.assertEqual(expected, ranges.expand(region, n_bp, contig_map))

  @parameterized.parameters(
      # Check that we don't create Ranges with negative starts.
      dict(
          region=ranges.make_range('1', 10, 20),
          n_bp=20,
          contig_map=None,
          expected=ranges.make_range('1', 0, 40),
      ),
      # Check that we respect n_bp if contig_map is provided.
      dict(
          region=ranges.make_range('1', 10, 20),
          n_bp=40,
          contig_map={
              '1': reference_pb2.ContigInfo(name='1', n_bases=50),
          },
          expected=ranges.make_range('1', 0, 50),
      ),
  )
  def test_expand_handles_boundaries(self, region, n_bp, contig_map, expected):
    self.assertEqual(expected, ranges.expand(region, n_bp, contig_map))

  def test_expand_raises_on_negative_n_bp(self):
    with self.assertRaisesRegex(ValueError, 'n_bp must be >= 0 but got -10'):
      ranges.expand(ranges.make_range('1', 10, 20), -10)

  def test_expand_raises_with_missing_contig_in_map(self):
    # Empty contig_map should raise.
    with self.assertRaises(KeyError):
      ranges.expand(ranges.make_range('1', 10, 20), 1, contig_map={})

    # Missing '1' from the contig map should raise.
    with self.assertRaises(KeyError):
      ranges.expand(
          ranges.make_range('1', 10, 20),
          1,
          contig_map={
              '2': reference_pb2.ContigInfo(name='2', n_bases=50),
          })

  @parameterized.parameters(
      dict(
          region=ranges.make_range(chrom, start, start + length),
          expected_length=length,
      )
      for length in range(10)
      for start in [10, 20, 1000]
      for chrom in ['1', '20']
  )
  def test_length_is_correct(self, region, expected_length):
    self.assertEqual(expected_length, ranges.length(region))


if __name__ == '__main__':
  absltest.main()
