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
"""Tests for calling_regions_utils.py."""

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import calling_regions_utils
from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.util import ranges


class CallingRegionsUtilsTest(parameterized.TestCase):
  # Total bps: 2100
  CONTIGS = [
      reference_pb2.ContigInfo(
          name='chr1',
          n_bases=1000,
          pos_in_fasta=0,
      ),
      reference_pb2.ContigInfo(
          name='chr2',
          n_bases=500,
          pos_in_fasta=1,
      ),
      reference_pb2.ContigInfo(
          name='chr3',
          n_bases=300,
          pos_in_fasta=2,
      ),
      reference_pb2.ContigInfo(
          name='chr4',
          n_bases=200,
          pos_in_fasta=3,
      ),
      reference_pb2.ContigInfo(
          name='chr5',
          n_bases=100,
          pos_in_fasta=4,
      ),
  ]

  @parameterized.parameters(
      dict(includes=[], excludes=[], expected=['1:1-100', '2:1-200']),
      dict(includes=['1'], excludes=[], expected=['1:1-100']),
      # Check that excludes work as expected.
      dict(includes=[], excludes=['1'], expected=['2:1-200']),
      dict(includes=[], excludes=['2'], expected=['1:1-100']),
      dict(includes=[], excludes=['1', '2'], expected=[]),
      # Check that excluding pieces works. The main checks on taking the
      # difference between two RangeSets live in ranges.py so here we are just
      # making sure some basic logic works.
      dict(includes=['1'], excludes=['1:1-10'], expected=['1:11-100']),
      # Check that includes and excludes work together.
      dict(
          includes=['1', '2'],
          excludes=['1:5-10', '1:20-50', '2:10-20'],
          expected=['1:1-4', '1:11-19', '1:51-100', '2:1-9', '2:21-200'],
      ),
      dict(
          includes=['1'],
          excludes=['1:5-10', '1:20-50', '2:10-20'],
          expected=['1:1-4', '1:11-19', '1:51-100'],
      ),
      dict(
          includes=['2'],
          excludes=['1:5-10', '1:20-50', '2:10-20'],
          expected=['2:1-9', '2:21-200'],
      ),
      # A complex example of including and excluding.
      dict(
          includes=['1:10-20', '2:50-60', '2:70-80'],
          excludes=['1:1-13', '1:19-50', '2:10-65'],
          expected=['1:14-18', '2:70-80'],
      ),
  )
  def test_build_calling_regions(self, includes, excludes, expected):
    contigs = [
        reference_pb2.ContigInfo(name='1', n_bases=100, pos_in_fasta=0),
        reference_pb2.ContigInfo(name='2', n_bases=200, pos_in_fasta=1),
    ]
    actual = calling_regions_utils.build_calling_regions(
        contigs, includes, excludes, None
    )
    self.assertCountEqual(actual, ranges.parse_literals(expected))

  @parameterized.parameters(
      dict(
          num_partitions=1,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=1000),
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=2,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=1000),
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=3,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=700),
                  range_pb2.Range(reference_name='chr1', start=700, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=4,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=525),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=525, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=10,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=210),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=210, end=420),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=420, end=630),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=630, end=840),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=840, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=210),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=210, end=420),
                  range_pb2.Range(reference_name='chr2', start=420, end=500),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=0, end=210),
                  range_pb2.Range(reference_name='chr3', start=210, end=300),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
              ],
              [
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=11,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=190, end=380),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=380, end=570),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=570, end=760),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=760, end=950),
                  range_pb2.Range(reference_name='chr1', start=950, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=190, end=380),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=380, end=500),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=190, end=300),
                  range_pb2.Range(reference_name='chr4', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=190, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
  )
  def test_partition_calling_regions(
      self, num_partitions, expected_partition_groups
  ):
    partition_groups = calling_regions_utils.partition_calling_regions(
        calling_regions=ranges.RangeSet.from_contigs(self.CONTIGS),
        num_partitions=num_partitions,
    )
    self.assertLen(partition_groups, num_partitions)
    self.assertEqual(partition_groups, expected_partition_groups)

  @parameterized.parameters(
      dict(
          regions_to_include=['chr1', 'chr2', 'chr3'],
          regions_to_exclude=[],
          num_partitions=4,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=450),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=450, end=900),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=900, end=1000),
                  range_pb2.Range(reference_name='chr2', start=0, end=450),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=450, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
              ],
          ],
      ),
      dict(
          regions_to_include=['chr1', 'chr2'],
          regions_to_exclude=['chr1:0-200'],
          num_partitions=3,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=200, end=633),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=633, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=433),
                  range_pb2.Range(reference_name='chr2', start=433, end=500),
              ],
          ],
      ),
  )
  def test_build_and_partition_calling_regions(
      self,
      regions_to_include,
      regions_to_exclude,
      num_partitions,
      expected_partition_groups,
  ):
    """Tests multiple functions to mimic use in postprocess_variants."""
    calling_regions = calling_regions_utils.build_calling_regions(
        contigs=self.CONTIGS,
        regions_to_include=regions_to_include,
        regions_to_exclude=regions_to_exclude,
        ref_n_regions=[],
    )
    partition_groups = calling_regions_utils.partition_calling_regions(
        calling_regions=calling_regions,
        num_partitions=num_partitions,
    )
    self.assertLen(partition_groups, num_partitions)
    self.assertEqual(partition_groups, expected_partition_groups)


if __name__ == '__main__':
  absltest.main()
