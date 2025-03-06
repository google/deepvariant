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
"""Tests for deepvariant.make_examples_core."""

from unittest import mock

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
import numpy as np

from deepvariant import dv_constants
from deepvariant import make_examples
from deepvariant import make_examples_core
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.small_model import make_small_model_examples
from third_party.nucleus.io import fasta
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import errors
from third_party.nucleus.util import ranges

FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


def _make_contigs(specs):
  """Makes ContigInfo protos from specs.

  Args:
    specs: A list of 2- or 3-tuples. All tuples should be of the same length. If
      2-element, these should be the name and length in basepairs of each
      contig, and their pos_in_fasta will be set to their index in the list. If
      the 3-element, the tuple should contain name, length, and pos_in_fasta.

  Returns:
    A list of ContigInfo protos, one for each spec in specs.
  """
  if specs and len(specs[0]) == 3:
    return [
        reference_pb2.ContigInfo(name=name, n_bases=length, pos_in_fasta=i)
        for name, length, i in specs
    ]
  else:
    return [
        reference_pb2.ContigInfo(name=name, n_bases=length, pos_in_fasta=i)
        for i, (name, length) in enumerate(specs)
    ]


def _from_literals_list(literals, contig_map=None):
  """Makes a list of Range objects from literals."""
  return ranges.parse_literals(literals, contig_map)


def _from_literals(literals, contig_map=None):
  """Makes a RangeSet of intervals from literals."""
  return ranges.RangeSet.from_regions(literals, contig_map)


class MakeExamplesCoreUnitTest(parameterized.TestCase):

  def test_read_write_run_info(self):
    def _read_lines(path):
      with open(path) as fin:
        return list(fin.readlines())

    golden_actual = make_examples_core.read_make_examples_run_info(
        testdata.GOLDEN_MAKE_EXAMPLES_RUN_INFO
    )
    # We don't really want to inject too much knowledge about the golden right
    # here, so we only use a minimal test that (a) the run_info_filename is
    # a non-empty string and (b) the number of candidates sites in the labeling
    # metrics field is greater than 0. Any reasonable golden output will have at
    # least one candidate variant, and the reader should have filled in the
    # value.
    self.assertNotEmpty(golden_actual.options.run_info_filename)
    self.assertEqual(
        golden_actual.labeling_metrics.n_candidate_variant_sites,
        testdata.N_GOLDEN_TRAINING_EXAMPLES,
    )

    # Check that reading + writing the data produces the same lines:
    tmp_output = test_utils.test_tmpfile('written_run_info.pbtxt')
    make_examples_core.write_make_examples_run_info(golden_actual, tmp_output)
    self.assertEqual(
        _read_lines(testdata.GOLDEN_MAKE_EXAMPLES_RUN_INFO),
        _read_lines(tmp_output),
    )

  @parameterized.parameters(
      dict(
          flag_value='CALLING',
          expected=deepvariant_pb2.MakeExamplesOptions.CALLING,
      ),
      dict(
          flag_value='TRAINING',
          expected=deepvariant_pb2.MakeExamplesOptions.TRAINING,
      ),
      dict(
          flag_value='CANDIDATE_SWEEP',
          expected=deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP,
      ),
  )
  def test_parse_proto_enum_flag(self, flag_value, expected):
    enum_pb2 = deepvariant_pb2.MakeExamplesOptions.Mode
    self.assertEqual(
        make_examples_core.parse_proto_enum_flag(enum_pb2, flag_value), expected
    )

  def test_parse_proto_enum_flag_error_handling(self):
    with self.assertRaisesRegex(
        ValueError,
        (
            'Unknown enum option "foo". Allowed options are'
            ' CALLING,CANDIDATE_SWEEP,TRAINING'
        ),
    ):
      make_examples_core.parse_proto_enum_flag(
          deepvariant_pb2.MakeExamplesOptions.Mode, 'foo'
      )

  def test_extract_sample_name_from_reads_single_sample(self):
    mock_sample_reader = mock.Mock()
    mock_sample_reader.header = reads_pb2.SamHeader(
        read_groups=[reads_pb2.ReadGroup(sample_id='sample_name')]
    )
    self.assertEqual(
        make_examples_core.extract_sample_name_from_sam_reader(
            mock_sample_reader
        ),
        'sample_name',
    )

  @parameterized.parameters(
      # No samples could be found in the reads.
      dict(samples=[], expected_sample_name=dv_constants.DEFAULT_SAMPLE_NAME),
      # Check that we detect an empty sample name and use default instead.
      dict(samples=[''], expected_sample_name=dv_constants.DEFAULT_SAMPLE_NAME),
      # We have more than one sample in the reads.
      dict(samples=['sample1', 'sample2'], expected_sample_name='sample1'),
  )
  def test_extract_sample_name_from_reads_uses_default_when_necessary(
      self, samples, expected_sample_name
  ):
    mock_sample_reader = mock.Mock()
    mock_sample_reader.header = reads_pb2.SamHeader(
        read_groups=[
            reads_pb2.ReadGroup(sample_id=sample) for sample in samples
        ]
    )
    self.assertEqual(
        expected_sample_name,
        make_examples_core.extract_sample_name_from_sam_reader(
            mock_sample_reader
        ),
    )

  @flagsaver.flagsaver
  def test_confident_regions(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)

    options = make_examples.default_options(add_flags=True)
    confident_regions = make_examples_core.read_confident_regions(options)

    # Our expected intervals, inlined from CONFIDENT_REGIONS_BED.
    expected = _from_literals_list([
        'chr20:10000847-10002407',
        'chr20:10002521-10004171',
        'chr20:10004274-10004964',
        'chr20:10004995-10006386',
        'chr20:10006410-10007800',
        'chr20:10007825-10008018',
        'chr20:10008044-10008079',
        'chr20:10008101-10008707',
        'chr20:10008809-10008897',
        'chr20:10009003-10009791',
        'chr20:10009934-10010531',
    ])
    # Our confident regions should be exactly those found in the BED file.
    self.assertCountEqual(expected, list(confident_regions))

  @flagsaver.flagsaver
  def test_invalid_examples_filename_extension(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = 'out.invalid_extension'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)

    options = make_examples.default_options(add_flags=True)
    with self.assertRaisesRegex(
        ValueError, 'Unsupported file extension: out.invalid_extension'
    ):
      make_examples_core.OutputsWriter(options)

  @flagsaver.flagsaver
  def test_gvcf_output_enabled_is_false_without_gvcf_flag(self):
    FLAGS.mode = 'training'
    FLAGS.gvcf = ''
    FLAGS.reads = ''
    FLAGS.ref = ''
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    options = make_examples.default_options(add_flags=True)
    self.assertFalse(make_examples_core.gvcf_output_enabled(options))

  @flagsaver.flagsaver
  def test_gvcf_output_enabled_is_true_with_gvcf_flag(self):
    FLAGS.mode = 'training'
    FLAGS.gvcf = '/tmp/foo.vcf'
    FLAGS.reads = ''
    FLAGS.ref = ''
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    options = make_examples.default_options(add_flags=True)
    self.assertTrue(make_examples_core.gvcf_output_enabled(options))

  @flagsaver.flagsaver
  def test_phase_reads_without_track_ref_reads_error(self):
    FLAGS.mode = 'training'
    FLAGS.gvcf = '/tmp/foo.vcf'
    FLAGS.reads = ''
    FLAGS.ref = ''
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.track_ref_reads = False
    FLAGS.phase_reads = True
    with self.assertRaisesRegex(
        errors.CommandLineError,
        '--track_ref_reads must be set to True when --phase_reads is set.',
    ):
      make_examples.default_options(add_flags=True)

  def test_validate_ref_contig_coverage(self):
    ref_contigs = _make_contigs([('1', 100), ('2', 100)])

    # Fully covered reference contigs don't trigger an error.
    for threshold in [0.5, 0.9, 1.0]:
      self.assertIsNone(
          make_examples_core.validate_reference_contig_coverage(
              ref_contigs, ref_contigs, threshold
          )
      )

    # No common contigs always blows up.
    for threshold in [0.0, 0.1, 0.5, 0.9, 1.0]:
      with self.assertRaisesRegex(ValueError, 'span 200'):
        make_examples_core.validate_reference_contig_coverage(
            ref_contigs, [], threshold
        )

    # Dropping either contig brings up below our 0.9 threshold.
    with self.assertRaisesRegex(ValueError, 'span 200'):
      make_examples_core.validate_reference_contig_coverage(
          ref_contigs, _make_contigs([('1', 100)]), 0.9
      )

    with self.assertRaisesRegex(ValueError, 'span 200'):
      make_examples_core.validate_reference_contig_coverage(
          ref_contigs, _make_contigs([('2', 100)]), 0.9
      )

    # Our actual overlap is 50%, so check that we raise when appropriate.
    with self.assertRaisesRegex(ValueError, 'span 200'):
      make_examples_core.validate_reference_contig_coverage(
          ref_contigs, _make_contigs([('2', 100)]), 0.6
      )
    self.assertIsNone(
        make_examples_core.validate_reference_contig_coverage(
            ref_contigs, _make_contigs([('2', 100)]), 0.4
        )
    )

  @parameterized.parameters(
      # all intervals are shared.
      ([[('chrM', 10)], [('chrM', 10)]], [('chrM', 10)]),
      # No common intervals.
      ([[('chrM', 10)], [('chr1', 10)]], []),
      # The names are the same but sizes are different, so not common.
      ([[('chrM', 10)], [('chrM', 20)]], []),
      # One common interval and one not.
      (
          [[('chrM', 10), ('chr1', 20)], [('chrM', 10), ('chr2', 30)]],
          [('chrM', 10)],
      ),
      # Check that the order doesn't matter.
      (
          [[('chr1', 20), ('chrM', 10)], [('chrM', 10), ('chr2', 30)]],
          [('chrM', 10, 1)],
      ),
      # Three-way merges.
      (
          [
              [('chr1', 20), ('chrM', 10)],
              [('chrM', 10), ('chr2', 30)],
              [('chr2', 30), ('chr3', 30)],
          ],
          [],
      ),
  )
  def test_common_contigs(self, contigs_list, expected):
    self.assertEqual(
        _make_contigs(expected),
        make_examples_core.common_contigs(
            [_make_contigs(contigs) for contigs in contigs_list]
        ),
    )

  @parameterized.named_parameters(
      # calling_regions = contigs & calling_regions.
      # Note that _from_literals() decreases start coordinate of an interval
      # by 1. But, candidate positions are given in real coordinates.
      dict(
          testcase_name='one_interval',
          regions=['1:1-10'],
          candidate_positions=[2, 4, 5, make_examples_core.END_OF_REGION],
          max_size=2,
          expected=['1:1-5', '1:6-10'],
      ),  # One interval 2 partitions.
      dict(
          testcase_name='two_intervals',
          regions=['1:1-10', '1:15-20'],
          candidate_positions=[
              2,
              4,
              5,
              make_examples_core.END_OF_REGION,
              16,
              19,
              make_examples_core.END_OF_REGION,
          ],
          max_size=2,
          expected=['1:1-5', '1:6-10', '1:15-20'],
      ),  # 2 intervals.
      dict(
          testcase_name='two_intervals_different_contigs',
          regions=['1:1-10', '2:1-20'],
          candidate_positions=[
              2,
              4,
              5,
              make_examples_core.END_OF_REGION,
              2,
              3,
              12,
              make_examples_core.END_OF_REGION,
          ],
          max_size=2,
          expected=['1:1-5', '1:6-10', '2:1-4', '2:5-20'],
      ),  # 2 intervals, different contigs.
      dict(
          testcase_name='candidate_far_apart',
          regions=['1:1-1000200'],
          candidate_positions=[3, 5, 7, make_examples_core.END_OF_REGION],
          max_size=2,
          expected=['1:1-6', '1:7-1000006', '1:1000007-1000200'],
      ),
  )  # Distance between candidates > max partition
  def test_partition_by_candidates(
      self, regions, candidate_positions, max_size, expected
  ):
    partitioned = make_examples_core.partition_by_candidates(
        regions=_from_literals(regions),
        candidate_positions=candidate_positions,
        max_size=max_size,
    )
    expected = _from_literals_list(expected)
    print(partitioned)
    self.assertCountEqual(expected, partitioned)

  @parameterized.named_parameters(
      dict(
          testcase_name='simple',
          position_arrays=[
              np.array([
                  1,
                  2,
                  3,
                  make_examples_core.END_OF_PARTITION,
                  7,
                  8,
                  9,
                  make_examples_core.END_OF_PARTITION,
                  make_examples_core.END_OF_REGION,
              ]),
              np.array([4, 5, 6, make_examples_core.END_OF_PARTITION]),
          ],
          expected=[
              1,
              2,
              3,
              4,
              5,
              6,
              7,
              8,
              9,
              make_examples_core.END_OF_REGION,
          ],
          check_failure=False,
      ),
      dict(
          testcase_name='one_partition_in_each_shard',
          position_arrays=[
              np.array([1, 3, 7, make_examples_core.END_OF_PARTITION]),
              np.array([
                  9,
                  11,
                  make_examples_core.END_OF_PARTITION,
                  make_examples_core.END_OF_REGION,
              ]),
          ],
          expected=[1, 3, 7, 9, 11, make_examples_core.END_OF_REGION],
          check_failure=False,
      ),
      dict(
          testcase_name='one_shard',
          position_arrays=[
              np.array([
                  1,
                  2,
                  3,
                  4,
                  7,
                  make_examples_core.END_OF_PARTITION,
                  make_examples_core.END_OF_REGION,
              ])
          ],
          expected=[1, 2, 3, 4, 7, make_examples_core.END_OF_REGION],
          check_failure=False,
      ),
      dict(
          testcase_name='empty_shard',
          position_arrays=[
              np.array([
                  1,
                  2,
                  3,
                  4,
                  7,
                  make_examples_core.END_OF_PARTITION,
                  make_examples_core.END_OF_REGION,
              ]),
              np.array([]),
          ],
          expected=[1, 2, 3, 4, 7, make_examples_core.END_OF_REGION],
          check_failure=False,
      ),
      dict(
          testcase_name='unordered_input',
          position_arrays=[
              np.array([1, 7, 3, make_examples_core.END_OF_PARTITION]),
              np.array([
                  4,
                  5,
                  make_examples_core.END_OF_PARTITION,
                  make_examples_core.END_OF_REGION,
              ]),
          ],
          expected=None,
          check_failure=True,
      ),
      dict(
          testcase_name='unordered_input_2',
          position_arrays=[
              np.array([1, 3, 7, make_examples_core.END_OF_PARTITION]),
              np.array([
                  4,
                  5,
                  make_examples_core.END_OF_PARTITION,
                  make_examples_core.END_OF_REGION,
              ]),
          ],
          expected=None,
          check_failure=True,
      ),
  )
  def test_merge_ranges_from_files_sequential(
      self, position_arrays, expected, check_failure
  ):
    if check_failure:
      with self.assertRaises(AssertionError):
        make_examples_core.merge_ranges_from_files_sequential(position_arrays)
    else:
      self.assertSequenceEqual(
          list(expected),
          make_examples_core.merge_ranges_from_files_sequential(
              position_arrays
          ),
      )

  @parameterized.parameters(
      # Note that these tests aren't so comprehensive as we are trusting that
      # the intersection code logic itself is good and well-tested elsewhere.
      # Here we are focusing on some basic tests and handling of missing
      # calling_region and confident_region data.
      (['1:1-10'], ['1:1-10']),
      (['1:1-100'], ['1:1-100']),
      (['1:50-150'], ['1:50-100']),
      (None, ['1:1-100', '2:1-200']),
      (['1:20-50'], ['1:20-50']),
      # Chr3 isn't part of our contigs; make sure we tolerate it.
      (['1:20-30', '1:40-60', '3:10-50'], ['1:20-30', '1:40-60']),
      # Check that we handle overlapping calling or confident regions.
      (['1:25-30', '1:20-40'], ['1:20-40']),
  )
  def test_regions_to_process(self, calling_regions, expected):
    contigs = _make_contigs([('1', 100), ('2', 200)])
    self.assertCountEqual(
        _from_literals_list(expected),
        make_examples_core.regions_to_process(
            contigs, 1000, calling_regions=_from_literals(calling_regions)
        ),
    )

  @parameterized.parameters(
      (
          50,
          None,
          [
              '1:1-50',
              '1:51-100',
              '2:1-50',
              '2:51-76',
              '3:1-50',
              '3:51-100',
              '3:101-121',
          ],
      ),
      (120, None, ['1:1-100', '2:1-76', '3:1-120', '3:121']),
      (500, None, ['1:1-100', '2:1-76', '3:1-121']),
      (10, ['1:1-20', '1:30-35'], ['1:1-10', '1:11-20', '1:30-35']),
      (8, ['1:1-20', '1:30-35'], ['1:1-8', '1:9-16', '1:17-20', '1:30-35']),
  )
  def test_regions_to_process_partition(
      self, max_size, calling_regions, expected
  ):
    contigs = _make_contigs([('1', 100), ('2', 76), ('3', 121)])
    self.assertCountEqual(
        _from_literals_list(expected),
        make_examples_core.regions_to_process(
            contigs, max_size, calling_regions=_from_literals(calling_regions)
        ),
    )

  @parameterized.parameters(
      dict(
          ref_reader=fasta.InMemoryFastaReader([('chr1', 0, 'GATACA')]),
          expected=[],
          min_region_len=3,
      ),
      dict(
          ref_reader=fasta.InMemoryFastaReader([('chr1', 0, 'NNNGATACA')]),
          expected=[ranges.make_range('chr1', 0, 3)],
          min_region_len=3,
      ),
      dict(
          ref_reader=fasta.InMemoryFastaReader([('chr1', 0, 'GATACANNN')]),
          expected=[ranges.make_range('chr1', 6, 9)],
          min_region_len=3,
      ),
      dict(
          ref_reader=fasta.InMemoryFastaReader([('chr1', 0, 'GATACANNNTTT')]),
          expected=[ranges.make_range('chr1', 6, 9)],
          min_region_len=3,
      ),
      dict(
          ref_reader=fasta.InMemoryFastaReader(
              [('chr1', 0, 'GATACANNNAAAAANNN')]
          ),
          expected=[
              ranges.make_range('chr1', 6, 9),
              ranges.make_range('chr1', 14, 17),
          ],
          min_region_len=3,
      ),
  )
  def test_find_ref_n_regions(self, ref_reader, expected, min_region_len):
    self.assertCountEqual(
        expected,
        make_examples_core.find_ref_n_regions(ref_reader, min_region_len),
    )

  def test_regions_to_process_sorted_within_contig(self):
    # These regions are out of order but within a single contig.
    contigs = _make_contigs([('z', 100)])
    in_regions = _from_literals(['z:15', 'z:20', 'z:6', 'z:25-30', 'z:3-4'])
    sorted_regions = _from_literals_list(
        ['z:3-4', 'z:6', 'z:15', 'z:20', 'z:25-30']
    )
    actual_regions = list(
        make_examples_core.regions_to_process(
            contigs, 100, calling_regions=in_regions
        )
    )
    # The assertEqual here is checking the order is exactly what we expect.
    self.assertEqual(sorted_regions, actual_regions)

  def test_regions_to_process_sorted_contigs(self):
    # These contig names are out of order lexicographically.
    contigs = _make_contigs([('z', 100), ('a', 100), ('n', 100)])
    in_regions = _from_literals(['a:10', 'n:1', 'z:20', 'z:5'])
    sorted_regions = _from_literals_list(['z:5', 'z:20', 'a:10', 'n:1'])
    actual_regions = list(
        make_examples_core.regions_to_process(
            contigs, 100, calling_regions=in_regions
        )
    )
    # The assertEqual here is checking the order is exactly what we expect.
    self.assertEqual(sorted_regions, actual_regions)

  @parameterized.parameters([2, 3, 4, 5, 50])
  def test_regions_to_process_sharding(self, num_shards):
    """Makes sure we deterministically split up regions."""

    def get_regions(task_id, num_shards):
      return make_examples_core.regions_to_process(
          contigs=_make_contigs([('z', 100), ('a', 100), ('n', 100)]),
          partition_size=5,
          task_id=task_id,
          num_shards=num_shards,
      )

    # Check that the regions are the same unsharded vs. sharded.
    unsharded_regions = get_regions(0, 0)
    sharded_regions = []
    for task_id in range(num_shards):
      task_regions = get_regions(task_id, num_shards)
      sharded_regions.extend(task_regions)
    self.assertCountEqual(unsharded_regions, sharded_regions)

  @parameterized.parameters(
      # Providing one of task id and num_shards but not the other is bad.
      (None, 0),
      (None, 2),
      (2, None),
      (0, None),
      # Negative values are illegal.
      (-1, 2),
      (0, -2),
      # task_id >= num_shards is bad.
      (2, 2),
      (3, 2),
  )
  def test_regions_to_process_fails_with_bad_shard_args(self, task, num_shards):
    with self.assertRaises(ValueError):
      make_examples_core.regions_to_process(
          contigs=_make_contigs([('z', 100), ('a', 100), ('n', 100)]),
          partition_size=10,
          task_id=task,
          num_shards=num_shards,
      )

  @parameterized.parameters(
      # Fetch all positions
      (['chr20:1-20000000'], 221),
      # Fetch subset of positions
      (['chr20:1-10003021'], 20),
  )
  def test_fetch_vcf_positions(self, calling_regions, expected_count):
    contigs = _make_contigs([('chr20', 20000000)])
    calling_regions = _from_literals(calling_regions)
    variant_positions = make_examples_core.fetch_vcf_positions(
        [testdata.TRUTH_VARIANTS_VCF], contigs, calling_regions
    )
    self.assertLen(variant_positions, expected_count)

  @parameterized.parameters(
      # One variant in region.
      (['x:100-200'], ['x:150-151'], [0]),
      # Different chromosomes.
      (['x:100-200'], ['y:150-151'], []),
      # A variant at the beginning of a region.
      (['x:100-200', 'x:201-300'], ['x:100-101'], [0]),
      (['x:1-10', 'x:11-20', 'x:21-30'], ['x:11-12'], [1]),
      # A variant before all the regions.
      (['x:11-20', 'x:20-30'], ['x:1-2'], []),
      # A variant after all the regions.
      (['x:1-10', 'x:11-20', 'x:21-30'], ['x:40-50'], []),
      # Multiple variants in the same region.
      (
          ['x:11-20', 'x:21-30'],
          ['x:1-2', 'x:25-26', 'x:25-26', 'x:26-27', 'x:40-50'],
          [1],
      ),
      # A variant spanning multiple regions belongs where it starts.
      (
          ['x:1-10', 'x:11-20', 'x:21-30', 'x:31-40', 'x:41-50', 'x:51-60'],
          ['x:15-66'],
          [1],
      ),
  )
  def test_filter_regions_by_vcf(
      self, region_literals, variant_literals, regions_to_keep
  ):
    regions = [ranges.parse_literal(l) for l in region_literals]
    variant_positions = [ranges.parse_literal(l) for l in variant_literals]
    output = make_examples_core.filter_regions_by_vcf(
        regions, variant_positions
    )
    list_output = list(output)
    list_expected = [regions[i] for i in regions_to_keep]
    self.assertEqual(list_output, list_expected)

  @parameterized.parameters(
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2', '3'],
          vcf_names=None,
          names_to_exclude=[],
          min_coverage_fraction=1.0,
          expected_names=['1', '2', '3'],
      ),
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=None,
          names_to_exclude=[],
          min_coverage_fraction=0.66,
          expected_names=['1', '2'],
      ),
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=['1', '3'],
          names_to_exclude=[],
          min_coverage_fraction=0.33,
          expected_names=['1'],
      ),
      dict(
          ref_names=['1', '2', '3', '4', '5'],
          sam_names=['1', '2', '3'],
          vcf_names=None,
          names_to_exclude=['4', '5'],
          min_coverage_fraction=1.0,
          expected_names=['1', '2', '3'],
      ),
  )
  def test_ensure_consistent_contigs(
      self,
      ref_names,
      sam_names,
      vcf_names,
      names_to_exclude,
      min_coverage_fraction,
      expected_names,
  ):
    ref_contigs = _make_contigs([(name, 100) for name in ref_names])
    sam_contigs = _make_contigs([(name, 100) for name in sam_names])
    if vcf_names is not None:
      vcf_contigs = _make_contigs([(name, 100) for name in vcf_names])
    else:
      vcf_contigs = None
    actual = make_examples_core._ensure_consistent_contigs(
        ref_contigs,
        sam_contigs,
        vcf_contigs,
        names_to_exclude,
        min_coverage_fraction,
    )
    self.assertEqual([a.name for a in actual], expected_names)

  @parameterized.parameters(
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=None,
          names_to_exclude=[],
          min_coverage_fraction=0.67,
      ),
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=['1', '3'],
          names_to_exclude=[],
          min_coverage_fraction=0.34,
      ),
  )
  def test_ensure_inconsistent_contigs(
      self,
      ref_names,
      sam_names,
      vcf_names,
      names_to_exclude,
      min_coverage_fraction,
  ):
    ref_contigs = _make_contigs([(name, 100) for name in ref_names])
    sam_contigs = _make_contigs([(name, 100) for name in sam_names])
    if vcf_names is not None:
      vcf_contigs = _make_contigs([(name, 100) for name in vcf_names])
    else:
      vcf_contigs = None
    with self.assertRaisesRegex(ValueError, 'Reference contigs span'):
      make_examples_core._ensure_consistent_contigs(
          ref_contigs,
          sam_contigs,
          vcf_contigs,
          names_to_exclude,
          min_coverage_fraction,
      )

  @flagsaver.flagsaver
  def test_regions_and_exclude_regions_flags(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.regions = 'chr20:10,000,000-11,000,000'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.exclude_regions = 'chr20:10,010,000-10,100,000'

    options = make_examples.default_options(add_flags=True)
    _, regions_from_options = (
        make_examples_core.processing_regions_from_options(options)
    )
    self.assertCountEqual(
        list(ranges.RangeSet(regions_from_options)),
        _from_literals_list(
            ['chr20:10,000,000-10,009,999', 'chr20:10,100,001-11,000,000']
        ),
    )

  @flagsaver.flagsaver
  def test_mixed_exclude_regions_flags(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    # NOTE: This is equivalent as when we write "chr20:10,010,001-10,100,000"
    bed_path = test_utils.test_tmpfile(
        'test_mixed_exclude_regions_flags.bed',
        contents='\t'.join(['chr20', '20010000', '20100000']) + '\n',
    )
    FLAGS.exclude_regions = (
        # The following is equivalent as:
        # 'chr20:10,010,001-10,100,000 chr20:20,010,001-20,100,000'
        'chr20:10,010,001-10,100,000 '
        + bed_path
    )
    options = make_examples.default_options(add_flags=True)
    _, regions_from_options = (
        make_examples_core.processing_regions_from_options(options)
    )
    self.assertCountEqual(
        list(ranges.RangeSet(regions_from_options)),
        _from_literals_list([
            'chr20:1-10010000',
            'chr20:10100001-20010000',
            'chr20:20100001-63025520',
        ]),
    )

  @flagsaver.flagsaver
  def test_regions_exclude_n_reference(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.discard_non_dna_regions = True

    options = make_examples.default_options(add_flags=True)
    _, regions_from_options = (
        make_examples_core.processing_regions_from_options(options)
    )
    self.assertCountEqual(
        list(ranges.RangeSet(regions_from_options)),
        _from_literals_list(
            ['chr20:9,995,001-11,095,000', 'chr20:59,776,001-60,001,000']
        ),
    )

  @flagsaver.flagsaver
  def test_incorrect_empty_regions(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    # Deliberately incorrect contig name.
    FLAGS.regions = '20:10,000,000-11,000,000'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)

    options = make_examples.default_options(add_flags=True)
    with self.assertRaisesRegex(ValueError, 'The regions to call is empty.'):
      make_examples_core.processing_regions_from_options(options)


class RegionProcessorTest(parameterized.TestCase):

  def setUp(self):
    super(RegionProcessorTest, self).setUp()
    self._saved_flags = flagsaver.save_flag_values()
    self.region = ranges.parse_literal('chr20:10,000,000-10,000,100')

    FLAGS.reads = ''
    self.options = make_examples.default_options(add_flags=False)
    self.options.reference_filename = testdata.CHR20_FASTA
    main_sample = self.options.sample_options[0]
    if not main_sample.reads_filenames:
      main_sample.reads_filenames.append(testdata.CHR20_BAM)
    main_sample.variant_caller_options.sample_name = 'sample_id'
    main_sample.name = 'sample_id'
    self.options.truth_variants_filename = testdata.TRUTH_VARIANTS_VCF
    self.options.mode = deepvariant_pb2.MakeExamplesOptions.TRAINING
    self.processor = make_examples_core.RegionProcessor(self.options)
    self.ref_reader = fasta.IndexedFastaReader(self.options.reference_filename)
    self.mock_init = self.add_mock('initialize')
    for sample in self.processor.samples:
      sample.in_memory_sam_reader = mock.Mock()
      sample.small_model_variant_caller = mock.Mock()
    self.default_shape = [5, 5, 7]

  def tearDown(self):
    super(RegionProcessorTest, self).tearDown()
    flagsaver.restore_flag_values(self._saved_flags)

  def add_mock(self, name, retval='dontadd', side_effect='dontadd'):
    patcher = mock.patch.object(self.processor, name, autospec=True)
    self.addCleanup(patcher.stop)
    mocked = patcher.start()
    if retval != 'dontadd':
      mocked.return_value = retval
    if side_effect != 'dontadd':
      mocked.side_effect = side_effect
    return mocked

  def test_on_demand_initialization_called_if_not_initialized(self):
    candidates = ['Candidates']
    self.assertFalse(self.processor.initialized)
    main_sample = self.processor.samples[0]

    mock_rrna = self.add_mock('region_reads_norealign', retval=[])
    mock_rr = self.add_mock('realign_reads', retval=[])
    mock_cir = self.add_mock(
        'candidates_in_region',
        retval=(
            {'main_sample': candidates},
            {'main_sample': []},
            {'main_sample': {}},
            0,
        ),
    )
    self.processor.process(self.region)
    test_utils.assert_called_once_workaround(self.mock_init)
    mock_rrna.assert_called_once_with(
        region=self.region,
        sam_readers=None,
        reads_filenames=main_sample.options.reads_filenames,
    )
    mock_rr.assert_called_once_with(reads=[], region=self.region)
    main_sample.in_memory_sam_reader.replace_reads.assert_called_once_with([])
    mock_cir.assert_called_once_with(self.region, None)

  def test_on_demand_initialization_not_called_if_initialized(self):
    self.processor.initialized = True
    self.assertTrue(self.processor.initialized)
    main_sample = self.processor.samples[0]
    mock_rrna = self.add_mock('region_reads_norealign', retval=[])
    mock_rr = self.add_mock('realign_reads', retval=[])
    mock_cir = self.add_mock(
        'candidates_in_region',
        retval=(
            {'main_sample': []},
            {'main_sample': []},
            {'main_sample': {}},
            0,
        ),
    )
    self.processor.process(self.region)
    test_utils.assert_not_called_workaround(self.mock_init)
    mock_rrna.assert_called_once_with(
        region=self.region,
        sam_readers=None,
        reads_filenames=main_sample.options.reads_filenames,
    )
    mock_rr.assert_called_once_with(reads=[], region=self.region)
    mock_cir.assert_called_once_with(self.region, None)

  def test_process_calls_no_candidates(self):
    main_sample = self.processor.samples[0]
    mock_rrna = self.add_mock('region_reads_norealign', retval=[])
    mock_rr = self.add_mock('realign_reads', retval=[])
    mock_cir = self.add_mock(
        'candidates_in_region',
        retval=(
            {'main_sample': []},
            {'main_sample': []},
            {'main_sample': {}},
            0,
        ),
    )
    candidates, gvcfs, runtimes, read_phases, phased_reads_count = (
        self.processor.process(self.region)
    )
    self.assertEmpty(candidates['main_sample'])
    self.assertEmpty(gvcfs['main_sample'])
    self.assertEmpty(read_phases['main_sample'])
    self.assertIsInstance(runtimes, dict)
    self.assertEqual(phased_reads_count, 0)
    mock_rrna.assert_called_once_with(
        region=self.region,
        sam_readers=None,
        reads_filenames=main_sample.options.reads_filenames,
    )
    mock_rr.assert_called_once_with(reads=[], region=self.region)
    mock_cir.assert_called_once_with(self.region, None)

  @parameterized.parameters([
      deepvariant_pb2.MakeExamplesOptions.TRAINING,
      deepvariant_pb2.MakeExamplesOptions.CALLING,
  ])
  def test_process_calls_with_candidates(self, mode):
    self.processor.options.mode = mode

    main_sample = self.processor.samples[0]
    mock_read = mock.MagicMock()
    mock_candidate = mock.MagicMock()
    mock_rrna = self.add_mock('region_reads_norealign', retval=[mock_read])
    mock_rr = self.add_mock('realign_reads', retval=[mock_read])
    mock_cir = self.add_mock(
        'candidates_in_region',
        retval=(
            {'main_sample': [mock_candidate]},
            {'main_sample': []},
            {'main_sample': {}},
            0,
        ),
    )
    candidates, gvcfs, runtimes, read_phases, phased_reads_count = (
        self.processor.process(self.region)
    )
    self.assertEqual(candidates['main_sample'], [mock_candidate])
    self.assertEmpty(gvcfs['main_sample'])
    self.assertEmpty(read_phases['main_sample'])
    self.assertIsInstance(runtimes, dict)
    self.assertEqual(phased_reads_count, 0)
    mock_rrna.assert_called_once_with(
        region=self.region,
        sam_readers=None,
        reads_filenames=main_sample.options.reads_filenames,
    )
    mock_rr.assert_called_once_with(reads=[mock_read], region=self.region)
    mock_cir.assert_called_once_with(self.region, None)

  @parameterized.parameters([
      deepvariant_pb2.MakeExamplesOptions.TRAINING,
      deepvariant_pb2.MakeExamplesOptions.CALLING,
  ])
  def test_process_keeps_ordering_of_candidates_and_examples(self, mode):
    self.processor.options.mode = mode

    r1, r2 = mock.Mock(), mock.Mock()
    c1, c2 = mock.Mock(), mock.Mock()
    main_sample = self.processor.samples[0]
    self.add_mock('region_reads_norealign', retval=[r1, r2])
    self.add_mock('realign_reads', retval=[r1, r2])
    self.add_mock(
        'candidates_in_region',
        retval=(
            {'main_sample': [c1, c2]},
            {'main_sample': []},
            {'main_sample': {}},
            0,
        ),
    )
    candidates, gvcfs, runtimes, read_phases, phased_reads_count = (
        self.processor.process(self.region)
    )
    self.assertEqual(candidates['main_sample'], [c1, c2])
    self.assertEmpty(gvcfs['main_sample'])
    self.assertEmpty(read_phases['main_sample'])
    self.assertEqual(phased_reads_count, 0)
    self.assertIsInstance(runtimes, dict)
    main_sample.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [r1, r2]
    )

  def test_process_with_realigner(self):
    self.processor.options.mode = deepvariant_pb2.MakeExamplesOptions.CALLING
    self.processor.options.realigner_enabled = True
    self.processor.options.realigner_options.CopyFrom(
        realigner_pb2.RealignerOptions()
    )
    self.processor.realigner = mock.Mock()
    self.processor.realigner.realign_reads.return_value = [], []

    main_sample = self.processor.samples[0]
    main_sample.sam_readers = [mock.Mock()]
    main_sample.sam_readers[0].query.return_value = []

    c1, c2 = mock.Mock(), mock.Mock()
    self.add_mock(
        'candidates_in_region',
        retval=(
            {'main_sample': [c1, c2]},
            {'main_sample': []},
            {'main_sample': {}},
            2,
        ),
    )

    candidates, gvcfs, runtimes, read_phases, phased_reads_count = (
        self.processor.process(self.region)
    )
    self.assertEqual(candidates['main_sample'], [c1, c2])
    self.assertEmpty(gvcfs['main_sample'])
    self.assertEmpty(read_phases['main_sample'])
    self.assertIsInstance(runtimes, dict)
    self.assertEqual(phased_reads_count, 2)
    main_sample.sam_readers[0].query.assert_called_once_with(self.region)
    self.processor.realigner.realign_reads.assert_called_once_with(
        [], self.region
    )
    main_sample.in_memory_sam_reader.replace_reads.assert_called_once_with([])

  def test_call_small_model_examples(self):
    self.processor.options.mode = deepvariant_pb2.MakeExamplesOptions.CALLING
    self.processor.options.call_small_model_examples = True
    candidate1, candidate2 = mock.Mock(), mock.Mock()
    n_stats = {'n_small_model_calls': 0}
    mock_writer = mock.Mock()
    self.processor.small_model_example_factory = mock.Mock()
    fake_inference_example_set = make_small_model_examples.InferenceExampleSet(
        skipped_candidates=[],
        candidates_with_alt_allele_indices=[
            (candidate1, (0,)),
            (candidate2, (0,)),
        ],
        inference_examples=[[], []],
    )
    self.processor.small_model_example_factory.encode_inference_examples.return_value = (
        fake_inference_example_set
    )
    cvo_candidate1, cvo_candidate2 = mock.Mock(), mock.Mock()
    main_sample = self.processor.samples[0]
    main_sample.small_model_variant_caller.call_variants.return_value = (
        [cvo_candidate1, cvo_candidate2],
        [],
    )

    candidates = [candidate1, candidate2]
    self.processor.call_small_model_examples_in_region(
        candidates=candidates,
        read_phases={},
        sample=main_sample,
        writer=mock_writer,
        n_stats=n_stats,
        runtimes={},
    )

    self.processor.small_model_example_factory.encode_inference_examples.assert_called_once_with(
        candidates, {}
    )
    main_sample.small_model_variant_caller.call_variants.assert_called_once_with(
        fake_inference_example_set.candidates_with_alt_allele_indices,
        fake_inference_example_set.inference_examples,
    )
    mock_writer.write_call_variant_outputs.assert_called_once_with(
        cvo_candidate1, cvo_candidate2
    )
    self.assertEqual(n_stats['n_small_model_calls'], 2)

  def test_candidates_in_region_no_reads(self):
    main_sample = self.processor.samples[0]
    main_sample.in_memory_sam_reader.query.return_value = []
    mock_ac = self.add_mock('_make_allele_counter_for_region')

    self.assertEqual(
        ({}, {}, {}, 0), self.processor.candidates_in_region(self.region)
    )

    main_sample.in_memory_sam_reader.query.assert_called_once_with(self.region)
    # A region with no reads should return out without making an AlleleCounter.
    test_utils.assert_not_called_workaround(mock_ac)

  @parameterized.parameters(True, False)
  def test_candidates_in_region(self, include_gvcfs):
    self.options.gvcf_filename = 'foo.vcf' if include_gvcfs else ''
    main_sample = self.processor.samples[0]
    reads = ['read1', 'read2']
    main_sample.in_memory_sam_reader.query.return_value = reads

    # Setup our make_allele_counter and other mocks.
    mock_ac = mock.Mock()
    mock_make_ac = self.add_mock(
        '_make_allele_counter_for_region', retval=mock_ac
    )
    # Setup our make_variant_caller and downstream mocks.
    mock_vc = mock.Mock()
    mock_vc.calls_and_gvcfs.return_value = (
        ['variant'],
        ['gvcf'] if include_gvcfs else [],
    )
    main_sample.variant_caller = mock_vc

    actual = self.processor.candidates_in_region(self.region)

    # Make sure we're getting our reads for the region.
    main_sample.in_memory_sam_reader.query.assert_called_once_with(self.region)

    # Make sure we're creating an AlleleCounter once and adding each of our
    # reads to it.
    mock_make_ac.assert_called_once_with(self.region, [])
    self.assertEqual(
        [mock.call(r, 'sample_id') for r in reads], mock_ac.add.call_args_list
    )

    # Make sure we call CallVariant for each of the counts returned by the
    # allele counter.
    include_med_dp = False
    mock_vc.calls_and_gvcfs.assert_called_once_with(
        allele_counters={'sample_id': mock_ac},
        target_sample='sample_id',
        include_gvcfs=include_gvcfs,
        include_med_dp=include_med_dp,
        left_padding=0,
        right_padding=0,
    )

    # Finally, our actual result should be the single 'variant' and potentially
    # the gvcf records, each organized by sample.
    expected_output = (
        {'main_sample': ['variant']},
        {'main_sample': ['gvcf'] if include_gvcfs else []},
        {'main_sample': {}},
        0,
    )
    self.assertEqual(expected_output, actual)

  # TODO Uncomment sitelist tests once the bug is fixed.
  # def test_output_sitelist_calling(self):
  #   FLAGS.mode = 'calling'
  #   FLAGS.ref = testdata.CHR20_FASTA
  #   FLAGS.reads = testdata.CHR20_BAM
  #   FLAGS.regions = 'chr20:10006000-10007612'
  #   FLAGS.examples = os.path.join(self.create_tempdir(),
  # 'examples.tfrecord.gz')
  #   FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
  #   FLAGS.output_sitelist = True
  #   options = make_examples.default_options(add_flags=True)
  #   make_examples_core.make_examples_runner(options)

  #   # Check that sitelist exists
  #   with open(FLAGS.examples + '.sitelist.tsv', 'r') as f:
  #     sitelist = f.readlines()

  #   self.assertNotEmpty(sitelist)
  #   # Check that last column is -1 in calling mode, indicating no label.
  #   self.assertSameElements(
  #       [x.strip().split('\t')[-1] for x in sitelist], ['-1']
  #   )

  # def test_output_sitelist_training(self):
  #   FLAGS.mode = 'training'
  #   FLAGS.reads = testdata.CHR20_BAM
  #   FLAGS.regions = 'chr20:10006000-10007612'
  #   FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
  #   FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
  #   FLAGS.ref = testdata.CHR20_FASTA
  #   FLAGS.examples = self.create_tempfile('examples.tfrecord.gz').full_path
  #   FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
  #   FLAGS.output_sitelist = True
  #   options = make_examples.default_options(add_flags=True)
  #   make_examples_core.make_examples_runner(options)

  #   # Check that sitelist exists
  #   with open(FLAGS.examples + '.sitelist.tsv', 'r') as f:
  #     sitelist = f.readlines()

  #   self.assertNotEmpty(sitelist)
  #   # Check that last column is -1 in calling mode, indicating no label.
  #   self.assertNoCommonElements(
  #       [x.strip().split('\t')[-1] for x in sitelist], ['-1']
  #   )

  @parameterized.parameters(
      (
          {
              'sort_by_haplotypes': True,
              'phase_reads': False,
          },
          r'Parsing HP AUX tag because [^\n]+ --sort_by_haplotypes',
      ),
      (
          {
              'reverse_haplotypes': True,
              'phase_reads': False,
          },
          r'Parsing HP AUX tag because [^\n]+ --reverse_haplotypes',
      ),
      (
          {
              'hp_tag_for_assembly_polishing': True,
              'sort_by_haplotypes': True,
              'phase_reads': False,
          },
          r'Parsing HP AUX tag because [^\n]+ --hp_tag_for_assembly_polishing',
      ),
      (
          {
              'channel_list': 'BASE_CHANNELS,haplotype',
              'phase_reads': False,
          },
          r'Parsing HP AUX tag because [^\n]+ haplotype channel is present',
      ),
      (
          # When phase_reads=True, we don't need to parse HP.
          {
              'channel_list': 'BASE_CHANNELS,haplotype',
              'phase_reads': True,
              'track_ref_reads': True,
          },
          r'Parsing AUX Fields: \[\]',
      ),
      (
          {
              'use_original_quality_scores': True,
          },
          'Parsing OQ AUX tag because --use_original_quality_scores is set.',
      ),
      (
          {
              'channel_list': 'BASE_CHANNELS,base_methylation',
          },
          (
              'Parsing MM, ML, and MN AUX tags because of base modification'
              ' channel.'
          ),
      ),
  )
  def test_aux_field_handling(
      self,
      test_flags,
      expected_log_message,
  ):
    flags_to_set = {
        'mode': 'calling',
        'ref': testdata.CHR20_FASTA,
        'reads': testdata.CHR20_BAM,
        'examples': 'examples.tfrecord',
        'channel_list': ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS),
    }
    flags_to_set.update(test_flags)
    with flagsaver.flagsaver(**flags_to_set):
      with self.assertLogs() as logs:
        make_examples.default_options(add_flags=True)

    self.assertRegex(
        '\n'.join(logs.output),
        expected_log_message,
    )


if __name__ == '__main__':
  absltest.main()
