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
"""Tests for deepvariant.realigner.window_selector."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import window_selector


class WindowSelectorTest(parameterized.TestCase):

  def setUp(self):
    self.config = realigner_pb2.RealignerOptions.WindowSelectorOptions(
        min_num_supporting_reads=2,
        max_num_supporting_reads=10,
        min_mapq=20,
        min_base_quality=20,
        min_windows_distance=4)

  @parameterized.parameters(
      # ------------------------------------------------------------------------
      # These reads are all simple and just test the basic position calculation.
      # ------------------------------------------------------------------------
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='4M', quals=[64] * 4),
          expected_candidate_positions=[12]),
      dict(
          read=test_utils.make_read(
              'AAGTA', start=10, cigar='2M2I1M', quals=[64] * 5),
          expected_candidate_positions=[10, 11, 12, 13]),
      dict(
          read=test_utils.make_read(
              'AAA', start=10, cigar='2M2D1M', quals=[64] * 3),
          expected_candidate_positions=[12, 13]),
      dict(
          read=test_utils.make_read(
              'TGATAC', start=10, cigar='2S3M1S', quals=[64] * 6),
          expected_candidate_positions=[8, 9, 11, 13]),
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='2M1X1M', quals=[64] * 4),
          expected_candidate_positions=[12]),
      # ------------------------------------------------------------------------
      # These reads test that we correctly ignore bases with low qualities.
      # ------------------------------------------------------------------------
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='4M', quals=[64, 64, 10, 30]),
          expected_candidate_positions=[]),
      dict(
          read=test_utils.make_read(
              'AAGTA', start=10, cigar='2M2I1M', quals=[64, 64, 10, 30, 64]),
          expected_candidate_positions=[11, 13]),
      dict(
          read=test_utils.make_read(
              'TGATAC',
              start=10,
              cigar='2S3M1S',
              quals=[64, 10, 64, 64, 64, 64]),
          expected_candidate_positions=[8, 11, 13]),
      dict(
          read=test_utils.make_read(
              'AAGA', start=10, cigar='2M1X1M', quals=[64, 64, 30, 10]),
          expected_candidate_positions=[12]),
  )
  def test_candidates_from_one_read(self, read, expected_candidate_positions):
    """Test WindowSelector.process_read() with reads of low quality."""
    region = ranges.make_range(read.alignment.position.reference_name, 0, 100)
    ref = 'A' * ranges.length(region)
    self.assertEqual(
        window_selector._candidates_from_reads(self.config, ref, [read],
                                               region),
        {pos: 1 for pos in expected_candidate_positions})

  @parameterized.parameters(
      # If we have no candidates, we have no windows to assemble.
      dict(candidates={}, expected_ranges=[]),
      # Our min count for candidates is 2, so we don't get any ranges back when
      # we have only a single candidate with a count of 1.
      dict(candidates={4: 1}, expected_ranges=[]),
      # Our max count for candidates is 10, so we don't get any ranges back when
      # we have only a single candidate with a count of 11.
      dict(candidates={4: 11}, expected_ranges=[]),
      # Check that this works with 2 isolated regions.
      dict(
          candidates={
              100: 5,
              200: 5,
          },
          expected_ranges=[
              ranges.make_range('ref', 96, 104),
              ranges.make_range('ref', 196, 204),
          ]),
      # Check that this works with 3 isolated regions.
      dict(
          candidates={
              100: 5,
              200: 5,
              300: 5,
          },
          expected_ranges=[
              ranges.make_range('ref', 96, 104),
              ranges.make_range('ref', 196, 204),
              ranges.make_range('ref', 296, 304),
          ]),
      # redacted
      # redacted
      # Check a simple example where we have candidates from two regions:
      dict(
          candidates={
              0: 2,
              2: 4,
              3: 11,
              8: 3
          },
          expected_ranges=[
              ranges.make_range('ref', -4, 6),
              ranges.make_range('ref', 4, 12),
          ]),
  )
  def test_candidates_to_windows(self, candidates, expected_ranges):
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        expected_ranges)

  @parameterized.parameters(range(1, 20))
  def test_candidates_to_windows_window_size(self, size):
    # We have a single candidate at position 100 with a 5 count.
    candidates = {100: 5}
    # We expect the created window to be +/- size from 100.
    expected = ranges.make_range('ref', 100 - size, 100 + size)
    self.config.min_windows_distance = size
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        [expected])

  @parameterized.parameters(range(1, 20))
  def test_candidates_to_windows_min_window_distance(self, distance):
    candidates = {
        # We one candidate at position 100 with a 5 count.
        100:
            5,
        # We have another candidate at outside of our distance with a 5 count,
        # so it should produce a candidate but not be joined with our our
        # candidate at 100.
        100 - 2 * distance:
            5,
        # Finally, we have another variant that is exactly distance away from
        # 100. It should be joined with the candidate at 100 to produce a single
        # larger window.
        100 + distance:
            5,
    }
    expected = [
        # Our first window is for the 100 - 2 * distance one.
        ranges.make_range('ref', 100 - 3 * distance, 100 - distance),
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
    candidates = {100 + i * distance: 5 for i in range(5)}
    # Which should all be merged together into one giant window.
    expected = [
        ranges.make_range('ref', 100 - distance,
                          max(candidates) + distance),
    ]
    self.config.min_windows_distance = distance
    self.assertEqual(
        window_selector._candidates_to_windows(self.config, candidates, 'ref'),
        expected)


if __name__ == '__main__':
  absltest.main()
