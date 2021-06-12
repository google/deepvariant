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
"""Tests for deepvariant.realigner.window_selector."""



from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import fasta
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import window_selector


class AlleleCountLinearWindowSelectorTest(parameterized.TestCase):

  def setUp(self):
    window_selector_model = realigner_pb2.WindowSelectorModel(
        model_type=realigner_pb2.WindowSelectorModel.ALLELE_COUNT_LINEAR,
        allele_count_linear_model=realigner_pb2.WindowSelectorModel
        .AlleleCountLinearModel(
            bias=0,
            coeff_soft_clip=0,
            coeff_substitution=-0.5,
            coeff_insertion=1,
            coeff_deletion=1,
            coeff_reference=-0.5,
            decision_boundary=0))
    self.config = realigner_pb2.WindowSelectorOptions(
        min_mapq=20,
        min_base_quality=20,
        min_windows_distance=4,
        region_expansion_in_bp=20,
        window_selector_model=window_selector_model)

  def assertCandidatesFromReadsEquals(self,
                                      reads,
                                      expected,
                                      start=None,
                                      end=None,
                                      ref=None):
    chrom = reads[0].alignment.position.reference_name
    start = 0 if start is None else start
    end = 20 if end is None else end
    region = ranges.make_range(chrom, start, end)

    if ref is None:
      ref = 'A' * (ranges.length(region) + 512)

    ref_reader = fasta.InMemoryFastaReader([(chrom, 0, ref)])
    if isinstance(expected, type) and issubclass(expected, Exception):
      with self.assertRaises(expected):
        window_selector._candidates_from_reads(self.config, ref_reader, reads,
                                               region)
    else:
      actual = window_selector._candidates_from_reads(self.config, ref_reader,
                                                      reads, region)
      self.assertEqual(actual, expected)

  @parameterized.parameters(
      # ------------------------------------------------------------------------
      # These reads are all simple and just test the basic position calculation.
      # ------------------------------------------------------------------------
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='4M', quals=[64] * 4),
          expected=[]),
      dict(
          read=test_utils.make_read(
              'AAGTA', start=10, cigar='2M2I1M', quals=[64] * 5),
          expected=[10, 11, 12, 13]),
      dict(
          read=test_utils.make_read(
              'AAA', start=10, cigar='2M2D1M', quals=[64] * 3),
          expected=[12, 13]),
      dict(
          read=test_utils.make_read(
              'TGATAC', start=10, cigar='2S3M1S', quals=[64] * 6),
          expected=[]),
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='2M1X1M', quals=[64] * 4),
          expected=[]),
  )
  def test_candidates_from_one_read(self, read, expected):
    """Test WindowSelector.process_read() with reads of low quality."""
    self.assertCandidatesFromReadsEquals(reads=[read], expected=expected)

  @parameterized.parameters(
      # --------------------------------------------------
      # Systematic combination of simple CIGAR operations.
      # --------------------------------------------------
      dict(
          reads=[
              test_utils.make_read(
                  'AAGA', start=10, cigar='4M', quals=[64] * 4),
              test_utils.make_read(
                  'AAAA', start=10, cigar='4M', quals=[64] * 4),
          ],
          expected=[]),
      dict(
          reads=[
              test_utils.make_read(
                  'AAAA', start=10, cigar='4M', quals=[64] * 4),
              test_utils.make_read(
                  'AAA', start=10, cigar='3M1D', quals=[64] * 3),
          ],
          expected=[13]),
      dict(
          reads=[
              test_utils.make_read(
                  'AAGA', start=10, cigar='4M', quals=[64] * 4),
              test_utils.make_read(
                  'AAA', start=10, cigar='3M1D', quals=[64] * 3),
          ],
          expected=[13]),
      dict(
          reads=[
              test_utils.make_read(
                  'AAAA', start=10, cigar='4M', quals=[64] * 4),
              test_utils.make_read(
                  'AAAAT', start=10, cigar='4M1I', quals=[64] * 5),
          ],
          expected=[13, 14]),
      dict(
          reads=[
              test_utils.make_read(
                  'AAAT', start=10, cigar='3M1S', quals=[64] * 4),
              test_utils.make_read(
                  'AAAAT', start=10, cigar='4M1I', quals=[64] * 5),
          ],
          expected=[13, 14]),
  )
  def test_candidates_from_reads(self, reads, expected):
    """Test WindowSelector.process_read() with reads of low quality."""
    self.assertCandidatesFromReadsEquals(reads=reads, expected=expected)


class WindowSelectorTest(parameterized.TestCase):

  def setUp(self):
    window_selector_model = realigner_pb2.WindowSelectorModel(
        model_type=realigner_pb2.WindowSelectorModel.VARIANT_READS,
        variant_reads_model=realigner_pb2.WindowSelectorModel
        .VariantReadsThresholdModel(
            min_num_supporting_reads=1, max_num_supporting_reads=10))
    self.config = realigner_pb2.WindowSelectorOptions(
        min_mapq=20,
        min_base_quality=20,
        min_windows_distance=4,
        region_expansion_in_bp=20,
        window_selector_model=window_selector_model)

  def assertCandidatesFromReadsEquals(self,
                                      reads,
                                      expected,
                                      start=None,
                                      end=None,
                                      ref=None):
    chrom = reads[0].alignment.position.reference_name
    start = 0 if start is None else start
    end = 20 if end is None else end
    region = ranges.make_range(chrom, start, end)

    if ref is None:
      ref = 'A' * (ranges.length(region) + 512)

    ref_reader = fasta.InMemoryFastaReader([(chrom, 0, ref)])
    if isinstance(expected, type) and issubclass(expected, Exception):
      with self.assertRaises(expected):
        window_selector._candidates_from_reads(self.config, ref_reader, reads,
                                               region)
    else:
      actual = window_selector._candidates_from_reads(self.config, ref_reader,
                                                      reads, region)
      self.assertEqual(actual, expected)

  @parameterized.parameters(
      # ------------------------------------------------------------------------
      # These reads are all simple and just test the basic position calculation.
      # ------------------------------------------------------------------------
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='4M', quals=[64] * 4),
          expected=[12]),
      dict(
          read=test_utils.make_read(
              'AAGTA', start=10, cigar='2M2I1M', quals=[64] * 5),
          expected=[10, 11, 12, 13]),
      dict(
          read=test_utils.make_read(
              'AAA', start=10, cigar='2M2D1M', quals=[64] * 3),
          expected=[12, 13]),
      dict(
          read=test_utils.make_read(
              'TGATAC', start=10, cigar='2S3M1S', quals=[64] * 6),
          expected=[8, 9, 10, 11, 12, 13]),
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='2M1X1M', quals=[64] * 4),
          expected=[12]),
      # ------------------------------------------------------------------------
      # These reads test that we correctly ignore bases with low qualities.
      # ------------------------------------------------------------------------
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='4M', quals=[64, 64, 10, 30]),
          expected=[]),
      # Only insertions/soft clips where all bases have above our minimum base
      # quality are included.
      dict(
          read=test_utils.make_read(
              'AAGTA', start=10, cigar='2M2I1M', quals=[64, 64, 10, 21, 64]),
          expected=[]),
      # The left 2S operator is ignored because one base has a 10 quality.
      dict(
          read=test_utils.make_read(
              'TGATAC',
              start=10,
              cigar='2S3M1S',
              quals=[21, 10, 64, 64, 64, 64]),
          expected=[11, 12, 13]),
      # The right 1S operator is ignored because one base has a 10 quality.
      dict(
          read=test_utils.make_read(
              'TGATAC',
              start=10,
              cigar='2S3M1S',
              quals=[64, 64, 64, 64, 64, 10]),
          expected=[8, 9, 10, 11]),
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='2M1X1M', quals=[64, 64, 30, 10]),
          expected=[12]),
  )
  def test_candidates_from_one_read(self, read, expected):
    """Test WindowSelector.process_read() with reads of low quality."""
    self.assertCandidatesFromReadsEquals(reads=[read], expected=expected)

  # Systematically test all combinations of cigar operations and positions in a
  # read.
  @parameterized.parameters(
      # Check that the M operator works. We have to look at the bases on the
      # genome to decide if it generates a position at 10.
      dict(bases='A', cigar='1M', expected=[]),
      dict(bases='C', cigar='1M', expected=[10]),
      # The mismatch operator X indicates that a position mismatches the
      # reference regardless of whether it actually matches the genome or not.
      # The window selector inspects the actual reference genome bases, though,
      # and generates candidate positions only if we genuinely mismatch.
      dict(bases='A', cigar='1X', expected=[]),
      dict(bases='C', cigar='1X', expected=[10]),
      # The match operator = indicates that a position matches the reference
      # even if that base actually mismatches the reference genome. The window
      # selector inspects the actual reference genome bases, though, and
      # generates candidate positions only if we genuinely mismatch.
      dict(bases='A', cigar='1=', expected=[]),
      dict(bases='C', cigar='1=', expected=[10]),
      # The deletion operator generates positions at start for operator length
      # in the 5' direction starting at the base after the deletion.
      dict(bases='A', cigar='1M1D', expected=[11]),
      dict(bases='A', cigar='1M2D', expected=[11, 12]),
      dict(bases='A', cigar='1M3D', expected=[11, 12, 13]),
      dict(bases='A', cigar='1M4D', expected=[11, 12, 13, 14]),
      # The insertion operator generates positions at start for + length
      # basepairs in the 5' direction and length - 1 in the 3' direction.
      dict(bases='AA', cigar='1M1I', expected=[10, 11]),
      dict(bases='AAA', cigar='1M2I', expected=[9, 10, 11, 12]),
      dict(bases='AAAA', cigar='1M3I', expected=[8, 9, 10, 11, 12, 13]),
      # The soft-clip operator generates positions at the start for operator
      # length bases.
      dict(bases='AA', cigar='1M1S', expected=[10, 11]),
      dict(bases='AAA', cigar='1M2S', expected=[9, 10, 11, 12]),
      dict(bases='AAAA', cigar='1M3S', expected=[8, 9, 10, 11, 12, 13]),
      dict(bases='AA', cigar='1S1M', expected=[9, 10]),
      dict(bases='AAA', cigar='2S1M', expected=[8, 9, 10, 11]),
      dict(bases='AAAA', cigar='3S1M', expected=[7, 8, 9, 10, 11, 12]),
      # The skip (N) and hard clip (H) operators are both ignored.
      dict(bases='AA', cigar='1M1N1M', expected=[]),
      dict(bases='AA', cigar='1M2N1M', expected=[]),
      dict(bases='A', cigar='1M1H', expected=[]),
      dict(bases='A', cigar='1M1H', expected=[]),
      dict(bases='A', cigar='1H1M', expected=[]),
      dict(bases='A', cigar='1H1M', expected=[]),
      # The python version raises an exception when seeing a PAD, which is ok
      # but isn't strictly necessary. The C++ implementation handles PADs when
      # counting alleles, so we've commented out this test.
      # C++ version:
      # dict(bases='AA', cigar='1M1P1M', expected=[]),
      # dict(bases='AA', cigar='1M2P1M', expected=[]),
      # Python version:
      # dict(bases='AA', cigar='1M1P1M', expected=ValueError),
      # dict(bases='AA', cigar='1M2P1M', expected=ValueError),
  )
  def test_candidates_from_reads_all_cigars(self, bases, cigar, expected):
    """Test WindowSelector.process_read() with reads of low quality."""
    read = test_utils.make_read(
        bases, start=10, cigar=cigar, quals=[64] * len(bases))
    self.assertCandidatesFromReadsEquals(reads=[read], expected=expected)

  @parameterized.parameters(
      dict(
          read=test_utils.make_read(
              'AGA', start=read_start, cigar='3M', quals=[64] * 3),
          region_start=region_start,
          region_end=region_start + 100,
          expected=[read_start + 1],
      ) for region_start in range(10) for read_start in range(region_start, 10))
  def test_candidates_from_reads_position_invariance(self, read, region_start,
                                                     region_end, expected):
    # Tests that a read with a mismatch at position read_start + 1 produces a
    # single candidate position at read_start + 1 regardless of where it occurs
    # within a single region spanning region_start - region_end.
    self.assertCandidatesFromReadsEquals(
        reads=[read], expected=expected, start=region_start, end=region_end)

  # Our region is 5-8 and we are testing that the read's mismatch is only
  # included when it's within the region and not when it's outside.
  # Expected region boundaries are extended according to region_expansion_in_bp
  # flag. region_expansion_in_bp is set to 20 by default,
  # so 5 to 8  becomes 5 - 20 to 8 + 20 <=> 0 to 28
  @parameterized.parameters(
      dict(
          read=test_utils.make_read('G', start=start, cigar='1M', quals=[64]),
          expected=[start] if 0 <= start < 28 else [],
      ) for start in range(10))
  def test_candidates_from_reads_respects_region(self, read, expected):
    self.assertCandidatesFromReadsEquals(
        reads=[read], expected=expected, start=5, end=8)

  # Our region is 5-8 and we have a 4 basepair deletion in our read. We expect
  # a mismatch count of one for each position in the deletion that overlaps the
  # interval.
  # Expected region boundaries are extended according to region_expansion_in_bp
  # flag. region_expansion_in_bp is set to 20 by default,
  # so 5 to 8  becomes 5 - 20 to 8 + 20 <=> 0 to 28
  @parameterized.parameters(
      dict(
          read=test_utils.make_read(
              'AA', start=start, cigar='1M4D1M', quals=[64, 64]),
          expected=[
              pos for pos in range(start + 1, start + 5) if 0 <= pos < 28
          ],
      ) for start in range(10))
  def test_candidates_from_reads_respects_region_deletion(self, read, expected):
    self.assertCandidatesFromReadsEquals(
        reads=[read], expected=expected, start=5, end=8, ref='A' * 100)

  def test_candidates_from_reads_counts_overlapping_events(self):
    # This read has a mismatch at position 2 and a 2 bp insertion at position 4,
    # so we need to double count the candidate positions from the mismatch and
    # insertion at position 2.
    read = test_utils.make_read(
        'AAGACCAAA', start=0, cigar='4M2I3M', quals=[64] * 9)
    expected = [2, 3, 4, 5]
    self.assertCandidatesFromReadsEquals(reads=[read], expected=expected)

  @parameterized.parameters(
      dict(
          read_mapq=read_mapq,
          min_mapq=min_mapq,
          expect_read_to_be_included=read_mapq >= min_mapq)
      for read_mapq in range(10, 15)
      for min_mapq in range(8, 17))
  def test_candidates_from_reads_respects_mapq(self, read_mapq, min_mapq,
                                               expect_read_to_be_included):
    read = test_utils.make_read(
        'AGA', start=10, cigar='3M', quals=[64] * 3, mapq=read_mapq)
    self.config.min_mapq = min_mapq
    self.assertCandidatesFromReadsEquals(
        reads=[read], expected=[11] if expect_read_to_be_included else [])

  @parameterized.parameters(
      dict(
          candidates=[100, 200],
          expected_ranges=[
              ranges.make_range('ref', 96, 104),
              ranges.make_range('ref', 196, 204),
          ]),
      # Check that this works with 3 isolated regions.
      dict(
          candidates=[100, 200, 300],
          expected_ranges=[
              ranges.make_range('ref', 96, 104),
              ranges.make_range('ref', 196, 204),
              ranges.make_range('ref', 296, 304),
          ]),
      # Check a simple example where we have two candidates from the same
      # region:
      dict(
          candidates=[2, 8],
          expected_ranges=[
              ranges.make_range('ref', -2, 12),
          ]),
      # Check a simple example where we have candidates from two regions:
      dict(
          candidates=[2, 14],
          expected_ranges=[
              ranges.make_range('ref', -2, 6),
              ranges.make_range('ref', 10, 18),
          ]),
      # Check boundary conditions for merging windows: should merge.
      dict(
          candidates=[2, 10],
          expected_ranges=[
              ranges.make_range('ref', -2, 14),
          ]),
      # Check boundary conditions for merging windows: should not merge.
      dict(
          candidates=[2, 11],
          expected_ranges=[
              ranges.make_range('ref', -2, 6),
              ranges.make_range('ref', 7, 15),
          ]),
  )
  def test_candidates_to_windows(self, candidates, expected_ranges):
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        expected_ranges)

  @parameterized.parameters(range(1, 20))
  def test_candidates_to_windows_window_size(self, size):
    # We have a single candidate at position 100 with a 5 count.
    candidates = [100]
    # We expect the created window to be +/- size from 100.
    expected = ranges.make_range('ref', 100 - size, 100 + size)
    self.config.min_windows_distance = size
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        [expected])

  @parameterized.parameters(range(1, 20))
  def test_candidates_to_windows_min_window_distance(self, distance):
    candidates = [
        # We one candidate at position 100 with a 5 count.
        100,
        # We have another candidate at outside of our distance with a 5 count,
        # so it should produce a candidate but not be joined with our our
        # candidate at 100.
        100 - 2 * distance - 1,
        # Finally, we have another variant that is exactly distance away from
        # 100. It should be joined with the candidate at 100 to produce a single
        # larger window.
        100 + distance
    ]
    expected = [
        # Our first window is for the 100 - 2 * distance one.
        ranges.make_range('ref', 100 - 3 * distance - 1, 100 - distance - 1),
        # Our second window starts at 100 (- distance for the window size) and
        # ends at 100 + distance + distance (again for window size).
        ranges.make_range('ref', 100 - distance, 100 + 2 * distance),
    ]
    self.config.min_windows_distance = distance
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        expected)

  @parameterized.parameters(range(1, 20))
  def test_candidates_to_windows_merged_close_candidates(self, distance):
    # Create five candidates separated by exactly distance from each other:
    # 100, 101, 102, 103, 104 for distance == 1
    # 100, 102, 104, 106, 108 for distance == 2
    candidates = [100 + i * distance for i in range(5)]
    # Which should all be merged together into one giant window.
    expected = [
        ranges.make_range('ref', 100 - distance,
                          max(candidates) + distance),
    ]
    self.config.min_windows_distance = distance
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        expected)

  def test_select_windows(self):
    # Simple end-to-end test of the high-level select_windows function. We give
    # it a few reads with a single candidate at 100 and we expect a window back
    # centered at 100.
    reads = [
        test_utils.make_read('AGA', start=99, cigar='3M', quals=[64] * 3),
        test_utils.make_read('AGA', start=99, cigar='3M', quals=[63] * 3),
        test_utils.make_read('AGA', start=99, cigar='3M', quals=[62] * 3),
    ]
    chrom = reads[0].alignment.position.reference_name
    ref_reader = fasta.InMemoryFastaReader([(chrom, 0, 'A' * 300)])
    region = ranges.make_range(chrom, 0, 200)

    self.assertEqual(
        window_selector.select_windows(self.config, ref_reader, reads, region),
        [ranges.make_range(chrom, 96, 104)])

  def test_select_windows_returns_empty_list_when_no_reads(self):
    self.assertEqual([],
                     window_selector.select_windows(
                         self.config,
                         ref_reader=fasta.InMemoryFastaReader([('chr1', 0,
                                                                'A' * 500)]),
                         reads=[],
                         region=ranges.make_range('chr1', 1, 100)))


if __name__ == '__main__':
  absltest.main()
