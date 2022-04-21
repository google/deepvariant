# Copyright 2019 Google LLC.
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

"""Tests for third_party.nucleus.util.vis."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import glob
import os
from absl.testing import absltest
from absl.testing import parameterized
import numpy as np

from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import vis
# pylint: disable=g-direct-tensorflow-import
from tensorflow.core.example import example_pb2
from tensorflow.core.example import feature_pb2


def _bytes_feature(list_of_strings):
  """Returns a bytes_list from a list of string / byte."""
  return feature_pb2.Feature(
      bytes_list=feature_pb2.BytesList(value=list_of_strings))


def _int_feature(list_of_ints):
  """Returns a int64_list from a list of int / bool."""
  return feature_pb2.Feature(
      int64_list=feature_pb2.Int64List(value=list_of_ints))


def _image_array(shape):
  return np.random.randint(255, size=shape, dtype=np.uint8)


def _mock_example_with_image(shape):
  arr = _image_array(shape)
  feature = {
      'image/encoded': _bytes_feature([arr.tobytes()]),
      'image/shape': _int_feature(shape)
  }
  return arr, example_pb2.Example(
      features=feature_pb2.Features(feature=feature))


def _mock_example_with_variant_and_alt_allele_indices(
    encoded_indices=b'\n\x01\x00', alleles=('A', 'C')):
  variant = test_utils.make_variant(chrom='X', alleles=alleles, start=10)
  feature = {
      'variant/encoded': _bytes_feature([variant.SerializeToString()]),
      'alt_allele_indices/encoded': _bytes_feature([encoded_indices])
  }
  return example_pb2.Example(features=feature_pb2.Features(feature=feature))


def _mock_pileup_array_with_reads():
  shape = (10, 15)  # (height, width)
  pileup = np.zeros(shape)
  # Top 5 rows and all other non-read space is left as 0.
  # Like an actual pileup with 4 reads, each 8 bases long.
  pileup[5, 0:8] = 254
  pileup[6, 1:9] = 254
  pileup[7, 4:12] = 100  # One read with low value.
  pileup[8, 6:14] = 254
  pileup[8, 8:10] = 100  # Two bases of another read with low value.
  return pileup


class VisTest(parameterized.TestCase):

  def test_get_image_array_from_example(self):
    shape = (3, 2, 4)
    arr, example = _mock_example_with_image(shape)
    decoded_image_array = vis.get_image_array_from_example(example)
    self.assertTrue((arr == decoded_image_array).all())

  @parameterized.parameters(((5, 4, 3),), ((10, 7, 5),))
  def test_split_3d_array_into_channels(self, input_shape):
    arr = np.random.random(input_shape)
    output = vis.split_3d_array_into_channels(arr)
    self.assertLen(output, input_shape[2])
    for i in range(input_shape[2]):
      self.assertEqual(output[i].shape, arr.shape[0:2])
      self.assertTrue((output[i] == arr[:, :, i]).all())

  def test_channels_from_example(self):
    shape = (3, 2, 4)
    arr, example = _mock_example_with_image(shape)
    channels = vis.channels_from_example(example)
    self.assertLen(channels, shape[2])
    self.assertTrue((channels[0] == arr[:, :, 0]).all())

  @parameterized.parameters(((4, 8), (4, 8, 3)), ((100, 20), (100, 20, 3)))
  def test_convert_6_channels_to_rgb(self, input_shape, expected_output_shape):
    channels = [np.random.random(input_shape) for _ in range(6)]
    rgb = vis.convert_6_channels_to_rgb(channels)
    self.assertEqual(rgb.shape, expected_output_shape)

  @parameterized.parameters((None,), ('RGB',))
  def test_draw_deepvariant_pileup_with_example_input(self, composite_type):
    _, example = _mock_example_with_image((100, 10, 7))
    # Test that it runs without error.
    vis.draw_deepvariant_pileup(example=example, composite_type=composite_type)

  @parameterized.parameters((None,), ('RGB',))
  def test_draw_deepvariant_pileup_with_channels_input(self, composite_type):
    channels = [_image_array((100, 221)) for _ in range(6)]
    # Test that it runs without error.
    vis.draw_deepvariant_pileup(
        channels=channels, composite_type=composite_type)

  @parameterized.parameters(
      ([[0.0, 1], [5, 10]], 0, 10, [[0, 25], [127, 255]]),
      ([[0.0, 0.1], [0.5, 1]], 0, 1, [[0, 25], [127, 255]]),
      ([[0.0, 0.1], [0.5, 1]], 0, 0.5, [[0, 51], [255, 255]]),
      ([[0.0, 0.1], [0.5, 1]], 0.5, 1, [[0, 0], [0, 255]]),
      ([[0.0, 0.1], [0.5, 1]], -1, 1, [[127, 140], [191, 255]]),
      ([[0.0, 0.1], [0.5, 1]], -1, 2, [[85, 93], [127, 170]]))
  def test_scale_colors_for_png(self, arr, vmin, vmax, expected):
    arr = np.array(arr)
    scaled = vis.scale_colors_for_png(arr, vmin=vmin, vmax=vmax)
    self.assertTrue((scaled == expected).all())

  @parameterized.parameters(
      ((100, 200), 'L'),
      ((100, 200, 3), 'RGB'),
  )
  def test_autoscale_colors_for_png(self, shape, expected_image_mode):
    arr = np.random.random(shape)
    scaled, image_mode = vis.autoscale_colors_for_png(arr)
    # Original array should be unchanged.
    self.assertLess(np.max(arr), 1)
    self.assertNotEqual(arr.dtype, np.uint8)
    # Output values have been scaled up and the array's data type changed.
    self.assertGreater(np.max(scaled), 1)
    self.assertEqual(scaled.dtype, np.uint8)
    self.assertEqual(image_mode, expected_image_mode)

  @parameterized.parameters(
      ((100, 200), 'L'),
      ((10, 1), 'L'),
      ((100, 200, 3), 'RGB'),
      ((10, 1, 3), 'RGB'),
      ((100, 200, 6), None),
      ((100, 200, 3, 1), None),
      ((100), None),
  )
  def test_get_image_type_from_array(self, shape, expected):
    arr = _image_array(shape)
    if expected is not None:
      self.assertEqual(vis._get_image_type_from_array(arr), expected)
    else:
      self.assertRaisesWithPredicateMatch(
          ValueError, lambda x: str(x).index('dimensions') != -1,
          vis.save_to_png, arr)

  @parameterized.parameters(
      ((100, 200, 3), True),
      ((100, 200), True),
      ((100, 200, 6), False),
      ((100, 200, 3, 1), False),
      ((100), False),
  )
  def test_save_to_png(self, shape, should_succeed):
    arr = _image_array(shape)

    if should_succeed:
      temp_dir = self.create_tempdir().full_path
      output_path = os.path.join(temp_dir, 'test.png')
      # check the file doesn't already exist before function runs
      self.assertEmpty(glob.glob(output_path))
      vis.save_to_png(arr, path=output_path)
      self.assertLen(glob.glob(output_path), 1)
    else:
      self.assertRaisesWithPredicateMatch(
          ValueError, lambda x: str(x).index('dimensions') != -1,
          vis.save_to_png, arr)

  @parameterized.parameters(
      ((100, 200, 3), True),
      ((100, 200), True),
      ((100, 200, 6), False),
      ((100, 200, 3, 1), False),
      ((100), False),
  )
  def test_array_to_png_works_with_floats(self, shape, should_succeed):
    arr = np.random.random(shape)

    if should_succeed:
      temp_dir = self.create_tempdir().full_path
      output_path = os.path.join(temp_dir, 'test.png')
      # Check the file doesn't already exist before function runs.
      self.assertEmpty(glob.glob(output_path))
      vis.array_to_png(arr, path=output_path)
      self.assertLen(glob.glob(output_path), 1)
    else:
      self.assertRaisesWithPredicateMatch(
          ValueError, lambda x: str(x).index('dimensions') != -1,
          vis.array_to_png, arr)

  def test_variant_from_example(self):
    example = _mock_example_with_variant_and_alt_allele_indices()
    variant = vis.variant_from_example(example)
    self.assertIsInstance(variant, variants_pb2.Variant)

  @parameterized.parameters(
      (b'\n\x01\x00', [0]),
      (b'\n\x02\x00\x01', [0, 1]),
  )
  def test_alt_allele_indices_from_example(self, encoded_indices, expected):
    example = _mock_example_with_variant_and_alt_allele_indices(encoded_indices)
    indices = vis.alt_allele_indices_from_example(example)
    self.assertEqual(indices, expected)

  @parameterized.parameters(
      ('chr1', 100, 'G', 'chr1:100_G'),
      ('X', 0, 'GACGT', 'X:0_GACGT'),
  )
  def test_locus_id_from_variant(self, chrom, pos, ref, expected):
    variant = test_utils.make_variant(
        chrom=chrom, alleles=[ref, 'A'], start=pos)
    locus_id = vis.locus_id_from_variant(variant)
    self.assertEqual(locus_id, expected)

  @parameterized.parameters(
      (b'\n\x01\x00', ['A', 'G', 'GA', 'AG'], 'G'),
      (b'\n\x02\x00\x01', ['C', 'CA', 'T', 'TA'], 'CA-T'),
      (b'\n\x02\x01\x02', ['C', 'CA', 'T', 'TA'], 'T-TA'),
  )
  def test_alt_from_example(self, encoded_indices, alleles, expected):
    example = _mock_example_with_variant_and_alt_allele_indices(
        encoded_indices=encoded_indices, alleles=alleles)
    alt = vis.alt_from_example(example)
    self.assertEqual(alt, expected)

  @parameterized.parameters(
      (b'\n\x01\x00', ['A', 'G', 'GA', 'AG'], 'X:10_A_G'),
      (b'\n\x02\x00\x01', ['C', 'CA', 'T', 'TA'], 'X:10_C_CA-T'),
      (b'\n\x02\x01\x02', ['C', 'CA', 'T', 'TA'], 'X:10_C_T-TA'),
  )
  def test_locus_id_with_alt(self, encoded_indices, alleles, expected):
    example = _mock_example_with_variant_and_alt_allele_indices(
        encoded_indices=encoded_indices, alleles=alleles)
    locus_id_with_alt = vis.locus_id_with_alt(example)
    self.assertEqual(locus_id_with_alt, expected)

  @parameterized.parameters(
      ([0], ['C'], 'C'),
      ([0, 1], ['C', 'TT'], 'C-TT'),
      ([3, 4], ['C', 'TT', 'T', 'G', 'A'], 'G-A'),
  )
  def test_alt_bases_from_indices(self, indices, alternate_bases, expected):
    alt = vis.alt_bases_from_indices(indices, alternate_bases)
    self.assertEqual(alt, expected)

  @parameterized.parameters([(0), (1), (2)])
  def test_label_from_example(self, truth_label):
    feature = {'label': _int_feature([truth_label])}
    example = example_pb2.Example(
        features=feature_pb2.Features(feature=feature))
    output = vis.label_from_example(example)
    self.assertEqual(truth_label, output)

  @parameterized.parameters([(0), (1), (2), (8), (9), (20)])
  def test_deepvariant_channel_names(self, num_channels):
    output = vis._deepvariant_channel_names(num_channels)
    self.assertLen(output, num_channels)

  def test_remove_ref_band(self):
    pileup = _mock_pileup_array_with_reads()
    bottom_part = vis.remove_ref_band(pileup)
    self.assertEqual((pileup.shape[0] - 5, pileup.shape[1]),
                     bottom_part.shape,
                     msg='Checking output shape is correct.')
    # Since the ref band is all zero, the sum should stay the same.
    self.assertEqual(
        np.sum(pileup),
        np.sum(bottom_part),
        msg='Checking bottom part of pileup is intact.')

    test_pileup = np.zeros((100, 200))
    self.assertEqual((95, 200), vis.remove_ref_band(test_pileup).shape)

    too_small = np.zeros((4, 10)) + 254
    with self.assertRaises(AssertionError):
      vis.remove_ref_band(too_small)

  def test_fraction_low_base_quality(self):
    shape = (10, 15)
    high_quality = [[], np.zeros(shape) + 254]
    low_quality = [[], np.zeros(shape) + 100]
    empty = [[], np.zeros(shape)]
    golden_pileup = [[], _mock_pileup_array_with_reads()]
    self.assertEqual(
        0, vis.fraction_low_base_quality(high_quality), msg='All high quality')
    self.assertEqual(
        1, vis.fraction_low_base_quality(low_quality), msg='All low quality')
    self.assertEqual(
        0, vis.fraction_low_base_quality(empty), msg='Empty pileup, no reads')
    self.assertEqual(
        0.3125,
        vis.fraction_low_base_quality(golden_pileup),
        msg='Mixed high and low quality')

  def test_fraction_reads_with_low_mapq(self):
    shape = (10, 15)
    filler_channels = [0] * 2
    high_quality = filler_channels + [np.zeros(shape) + 254]
    low_quality = filler_channels + [np.zeros(shape) + 100]
    empty = filler_channels + [np.zeros(shape)]
    golden_pileup = filler_channels + [_mock_pileup_array_with_reads()]
    self.assertEqual(
        0,
        vis.fraction_reads_with_low_mapq(high_quality),
        msg='All high quality')
    self.assertEqual(
        1, vis.fraction_reads_with_low_mapq(low_quality), msg='All low quality')
    self.assertEqual(
        0,
        vis.fraction_reads_with_low_mapq(empty),
        msg='Empty pileup, no reads')
    self.assertEqual(
        0.25,
        vis.fraction_reads_with_low_mapq(golden_pileup),
        msg='Mixed high and low quality')

  def test_fraction_read_support_and_describer(self):
    shape = (10, 15)
    filler_channels = [0] * 4
    all_support = filler_channels + [np.zeros(shape) + 254]
    no_support = filler_channels + [np.zeros(shape) + 100]
    empty = filler_channels + [np.zeros(shape)]
    golden_pileup = filler_channels + [_mock_pileup_array_with_reads()]

    self.assertEqual(1, vis.fraction_read_support(all_support))
    self.assertEqual(vis.ReadSupport.ALL,
                     vis.describe_read_support(all_support))

    self.assertEqual(0, vis.fraction_read_support(no_support))
    self.assertEqual(vis.ReadSupport.LOW, vis.describe_read_support(no_support))

    self.assertEqual(0, vis.fraction_read_support(empty))
    self.assertEqual(vis.ReadSupport.LOW, vis.describe_read_support(empty))

    self.assertEqual(0.75, vis.fraction_read_support(golden_pileup))
    self.assertEqual(vis.ReadSupport.HALF,
                     vis.describe_read_support(golden_pileup))

  @parameterized.parameters([
      dict(k=12, n=24, expected_p=1.0),
      dict(k=1, n=4, expected_p=0.625),
      dict(k=3, n=4, expected_p=0.625),
      dict(k=0, n=4, expected_p=0.125),
      dict(k=4, n=4, expected_p=0.125),
      dict(k=0, n=8, expected_p=0.0078125),
      dict(k=8, n=8, expected_p=0.0078125)
  ])
  def test_binomial_test(self, k, n, expected_p):
    observed_p = vis.binomial_test(k=k, n=n)
    self.assertEqual(expected_p, observed_p)

  @parameterized.parameters([
      dict(test_case='support = forward', expected=0.0625),
      dict(test_case='support = reverse', expected=0.0625),
      dict(test_case='support = 5+/5-', expected=1.0),
      dict(test_case='support = 2+/2-', expected=1.0),
      # From scipy.stats.binom_test(x=1, n=6):
      dict(test_case='support = 1+/5-', expected=0.21875),
      # For two-tailed, this must match the previous:
      dict(test_case='support = 5+/1-', expected=0.21875)
  ])
  def test_pvalue_for_strand_bias(self, test_case, expected):
    shape = (15, 4)
    strand = np.zeros(shape)
    strand[5:10, :] = 240  # Forward.
    strand[10:15, :] = 70  # Reverse.

    read_support = np.zeros(shape)
    if test_case == 'support = forward':
      read_support[5:10, :] = 254  # Supporting.
      read_support[10:15, :] = 100  # Anything not 254 means not supporting.
    elif test_case == 'support = reverse':
      read_support[5:10, :] = 100  # Not supporting.
      read_support[10:15, :] = 254  # Supporting.
    elif test_case == 'support = 5+/5-':
      read_support[5:15, :] = 254  # All support: five forward, five reverse.
    elif test_case == 'support = 2+/2-':
      read_support[5:15, :] = 100  # Most not supporting.
      read_support[8:12, :] = 254  # Two supporting from each strand.
    elif test_case == 'support = 1+/5-':
      read_support[5:15, :] = 100  # Most not supporting.
      read_support[5:6, :] = 254  # One forward support.
      read_support[10:15, :] = 254  # Five reverse support.
    elif test_case == 'support = 5+/1-':
      read_support[5:15, :] = 100  # Most not supporting.
      read_support[5:10, :] = 254  # Five forward support.
      read_support[10:11, :] = 254  # One reverse support.
    else:
      raise ValueError('test_case not recognized')

    filler_channels = [0] * 3
    channels = filler_channels + [strand, read_support]
    self.assertEqual(expected, vis.pvalue_for_strand_bias(channels))

  @parameterized.parameters([
      dict(
          test_case='nearby_variants',
          expected_description=vis.Diff.NEARBY_VARIANTS,
          expected_diff_fraction=0.0,
          expected_nearby_variants=5),
      dict(
          test_case='few_diffs',
          expected_description=vis.Diff.FEW_DIFFS,
          expected_diff_fraction=0.0,
          expected_nearby_variants=2),
      dict(
          test_case='many_diffs',
          expected_description=vis.Diff.MANY_DIFFS,
          expected_diff_fraction=0.1,
          expected_nearby_variants=0),
      dict(
          test_case='empty',
          expected_description=vis.Diff.FEW_DIFFS,
          expected_diff_fraction=0.0,
          expected_nearby_variants=0)
  ])
  def test_analyze_diff_and_nearby_variants_and_describe_diff(
      self, test_case, expected_description, expected_diff_fraction,
      expected_nearby_variants):
    shape = (15, 8)
    diff_channel = np.zeros(shape) + 100
    if test_case == 'nearby_variants':
      # Five columns with homozygous variants:
      diff_channel[5:, [0, 1, 2, 4, 6]] = 254
    elif test_case == 'few_diffs':
      # Less than five columns with homozygous variants:
      diff_channel[5:, [2, 5]] = 254
    elif test_case == 'many_diffs':
      # One read full of differences:
      diff_channel[5, 0:8] = 254
    elif test_case == 'empty':
      # No reads:
      diff_channel = np.zeros(shape)
    else:
      raise ValueError('test_case not recognized')
    filler_channels = [0] * 5
    channels = filler_channels + [diff_channel]
    diff_fraction, nearby_variants = vis.analyze_diff_and_nearby_variants(
        channels)
    self.assertEqual(diff_fraction, expected_diff_fraction)
    self.assertEqual(nearby_variants, expected_nearby_variants)
    self.assertEqual(expected_description, vis.describe_diff(channels))

  def test_curate_pileup(self):
    # Use the same pileup array for all of the channels.
    # It has 4 reads, 2 of which are high values throughout the read, one read
    # is all low values, and one read is high values except 2 lower-value bases.
    channels = [_mock_pileup_array_with_reads() for _ in range(6)]
    tags = vis.curate_pileup(channels)

    # One read plus a few bases are low quality:
    self.assertEqual(tags.base_quality, vis.BaseQuality.BAD)
    # One fully low quality read out of four is enough to be "bad" mapq:
    self.assertEqual(tags.mapping_quality, vis.MappingQuality.BAD)
    # Not enough reads to get a p-value below 0.05 for strand bias:
    self.assertEqual(tags.strand_bias, vis.StrandBias.GOOD)
    # Many differences (large fraction of high values):
    self.assertEqual(tags.diff_category, vis.Diff.MANY_DIFFS)
    # One of four reads supporting is interpreted as roughly heterozygous:
    self.assertEqual(tags.read_support, vis.ReadSupport.HALF)


if __name__ == '__main__':
  absltest.main()
