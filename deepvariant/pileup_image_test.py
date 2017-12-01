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
"""Tests for deepvariant.pileup_image."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized
import mock
import numpy as np
import numpy.testing as npt

from deepvariant import pileup_image
from deepvariant import test_utils
from deepvariant.core import ranges
from deepvariant.core.genomics import variants_pb2
from deepvariant.core.python import reference_fai
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import pileup_image_native


def _supporting_reads(*names):
  return deepvariant_pb2.DeepVariantCall.SupportingReads(read_names=names)


def _make_dv_call(ref_bases='A', alt_bases='C'):
  return deepvariant_pb2.DeepVariantCall(
      variant=variants_pb2.Variant(
          reference_name='chr1',
          start=10,
          end=11,
          reference_bases=ref_bases,
          alternate_bases=[alt_bases]),
      allele_support={
          'C': _supporting_reads('read1/1', 'read2/1')
      })


def _make_encoder(**kwargs):
  """Make a PileupImageEncoderNative with overrideable default options."""
  options = pileup_image.default_options()
  options.MergeFrom(deepvariant_pb2.PileupImageOptions(**kwargs))
  return pileup_image_native.PileupImageEncoderNative(options)


def _make_image_creator(ref_reader, sam_reader_obj, **kwargs):
  options = pileup_image.default_options()
  options.MergeFrom(deepvariant_pb2.PileupImageOptions(**kwargs))
  return pileup_image.PileupImageCreator(options, ref_reader, sam_reader_obj)


class PileupImageEncoderTest(parameterized.TestCase):

  @parameterized.parameters(('A', 250), ('G', 180), ('T', 100), ('C', 30),
                            ('N', 0), ('X', 0))
  def test_base_color(self, base, expected_color):
    pie = _make_encoder(
        base_color_offset_a_and_g=40,
        base_color_offset_t_and_c=30,
        base_color_stride=70)
    self.assertAlmostEqual(pie.base_color(base), expected_color)

  @parameterized.parameters(
      (0, 0),
      (5, int(254 * 0.25)),
      (10, int(254 * 0.5)),
      (15, int(254 * 0.75)),
      (19, int(254 * 0.95)),
      (20, 254),
      (21, 254),
      (25, 254),
      (40, 254),
  )
  def test_base_quality_color(self, base_qual, expected_color):
    pie = _make_encoder(base_quality_cap=20)
    self.assertAlmostEqual(pie.base_quality_color(base_qual), expected_color)

  @parameterized.parameters(
      (0, 0),
      (5, int(254 * 0.25)),
      (10, int(254 * 0.5)),
      (15, int(254 * 0.75)),
      (19, int(254 * 0.95)),
      (20, 254),
      (21, 254),
      (25, 254),
      (40, 254),
  )
  def test_mapping_quality_color(self, mapping_qual, expected_color):
    pie = _make_encoder(mapping_quality_cap=20)
    self.assertAlmostEqual(
        pie.mapping_quality_color(mapping_qual), expected_color)

  @parameterized.parameters((True, 1), (False, 2))
  def test_strand_color(self, on_positive_strand, expected_color):
    pie = _make_encoder(positive_strand_color=1, negative_strand_color=2)
    self.assertAlmostEqual(pie.strand_color(on_positive_strand), expected_color)

  @parameterized.parameters(
      (False, int(254.0 * 0.2)),
      (True, int(254.0 * 0.1)),
  )
  def test_supports_alt_color(self, supports_alt, expected_color):
    pie = _make_encoder(
        allele_supporting_read_alpha=0.1, allele_unsupporting_read_alpha=0.2)
    self.assertAlmostEqual(pie.supports_alt_color(supports_alt), expected_color)

  @parameterized.parameters(
      (False, int(254.0 * 0.4)),
      (True, int(254.0 * 0.3)),
  )
  def test_matches_ref_color(self, matches_ref, expected_color):
    pie = _make_encoder(
        reference_matching_read_alpha=0.3, reference_mismatching_read_alpha=0.4)
    self.assertAlmostEqual(pie.matches_ref_color(matches_ref), expected_color)

  def test_reference_encoding(self):
    self.assertImageRowEquals(
        _make_encoder().encode_reference('ACGTN'),
        np.dstack([
            # Base.
            (250, 30, 180, 100, 0),
            # Base quality.
            (254, 254, 254, 254, 254),
            # Mapping quality.
            (254, 254, 254, 254, 254),
            # Strand channel (forward or reverse)
            (70, 70, 70, 70, 70),
            # Supports alt or not.
            (152, 152, 152, 152, 152),
            # Matches ref or not.
            (50, 50, 50, 50, 50),
            # Cigar operation length. 0 for reference.
            (0, 0, 0, 0, 0)
        ]).astype(np.uint8))

  def assertImageRowEquals(self, image_row, expected):
    npt.assert_equal(image_row, expected.astype(np.uint8))

  def test_encode_read_matches(self):
    start = 10
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases[0]
    read = test_utils.make_read(
        'ACCGT', start=start, cigar='5M', quals=range(10, 15), name='read1')
    full_expected = np.dstack([
        # Base.
        (250, 30, 30, 180, 100),
        # Base quality.
        (63, 69, 76, 82, 88),
        # Mapping quality.
        (211, 211, 211, 211, 211),
        # Strand channel (forward or reverse)
        (70, 70, 70, 70, 70),
        # Supports alt or not.
        (254, 254, 254, 254, 254),
        # Matches ref or not.
        (50, 50, 254, 50, 50),
        # Cigar operation length.
        (5, 5, 5, 5, 5)
    ]).astype(np.uint8)

    self.assertImageRowEquals(_make_encoder().encode_read(
        dv_call, 'ACAGT', read, start, alt_allele), full_expected)

  @parameterized.parameters((bases_start, bases_end)
                            for bases_start in range(0, 5)
                            for bases_end in range(6, 12))
  def test_encode_read_spans2(self, bases_start, bases_end):
    bases = 'AAAACCGTCCC'
    quals = [9, 9, 9, 10, 11, 12, 13, 14, 8, 8, 8]
    bases_start_offset = 7
    ref_start = 10
    ref_size = 5
    read_bases = bases[bases_start:bases_end]
    read_quals = quals[bases_start:bases_end]
    read_start = bases_start_offset + bases_start

    # Create our expected image row encoding.
    op_len = bases_end - bases_start
    full_expected = np.dstack([
        # Base.
        (250, 30, 30, 180, 100),
        # Base quality.
        (63, 69, 76, 82, 88),
        # Mapping quality.
        (211, 211, 211, 211, 211),
        # Strand channel (forward or reverse)
        (70, 70, 70, 70, 70),
        # Supports alt or not.
        (254, 254, 254, 254, 254),
        # Matches ref or not.
        (50, 50, 254, 50, 50),
        # Cigar operation length.
        [op_len] * 5
    ]).astype(np.uint8)
    expected = np.zeros(
        (1, ref_size, pileup_image.DEFAULT_NUM_CHANNEL), dtype=np.uint8)
    for i in range(read_start, read_start + len(read_bases)):
      if ref_start <= i < ref_start + ref_size:
        expected[0, i - ref_start] = full_expected[0, i - ref_start]

    read = test_utils.make_read(
        read_bases,
        start=read_start,
        cigar=str(len(read_bases)) + 'M',
        quals=read_quals,
        name='read1')
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases[0]
    self.assertImageRowEquals(_make_encoder().encode_read(
        dv_call, 'ACAGT', read, ref_start, alt_allele), expected)

  def test_encode_read_deletion(self):
    # ref:  AACAG
    # read: AA--G
    start = 2
    read = test_utils.make_read(
        'AAG', start=start, cigar='2M2D1M', quals=range(10, 13), name='read1')
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases[0]
    full_expected = np.dstack([
        # Base. The second A is 0 because it's the anchor of the deletion.
        (250, 0, 0, 0, 180),
        # Base quality.
        (63, 69, 0, 0, 76),
        # Mapping quality.
        (211, 211, 0, 0, 211),
        # Strand channel (forward or reverse)
        (70, 70, 0, 0, 70),
        # Supports alt or not.
        (254, 254, 0, 0, 254),
        # Matches ref or not.
        (50, 254, 0, 0, 50),
        # Cigar operation length.
        (2, 2, 0, 0, 1)
    ]).astype(np.uint8)
    self.assertImageRowEquals(_make_encoder().encode_read(
        dv_call, 'AACAG', read, start, alt_allele), full_expected)

  def test_encode_read_insertion(self):
    # ref:  AA-CAG
    # read: AAACAG
    start = 2
    read = test_utils.make_read(
        'AAACAG',
        start=start,
        cigar='2M1I3M',
        quals=range(10, 16),
        name='read1')
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases[0]
    full_expected = np.dstack([
        # Base.
        (250, 0, 30, 250, 180),
        # Base quality.
        (63, 76, 82, 88, 95),
        # Mapping quality.
        (211, 211, 211, 211, 211),
        # Strand channel (forward or reverse)
        (70, 70, 70, 70, 70),
        # Supports alt or not.
        (254, 254, 254, 254, 254),
        # Matches ref or not.
        (50, 254, 50, 50, 50),
        # Cigar operation length.
        (2, 1, 3, 3, 3)
    ]).astype(np.uint8)
    self.assertImageRowEquals(_make_encoder().encode_read(
        dv_call, 'AACAG', read, start, alt_allele), full_expected)

  def test_ignores_reads_with_low_quality_bases(self):
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=2,
            end=3,
            reference_bases='A',
            alternate_bases=['C']))
    pie = _make_encoder()

    # Get the threshold the encoder uses.
    min_qual = pileup_image.DEFAULT_MIN_BASE_QUALITY

    for qual in range(0, min_qual + 5):
      quals = [min_qual - 1, qual, min_qual + 1]
      read = test_utils.make_read('AAA', start=1, cigar='3M', quals=quals)
      actual = pie.encode_read(dv_call, 'AACAG', read, 1, 'C')
      if qual < min_qual:
        self.assertIsNone(actual)
      else:
        self.assertIsNotNone(actual)

  @parameterized.parameters(
      # alt_allele is C, supported by only frag 1 of read1 and frag 2 of read 3.
      ('read1', 1, 'C', 'C', True),
      ('read1', 2, 'C', 'C', False),
      ('read2', 1, 'C', 'G', False),
      ('read2', 2, 'C', 'G', False),
      ('read3', 1, 'C', 'C', False),
      ('read3', 2, 'C', 'C', True),
      # alt_allele is now G, so only read2 (both frags) support alt.
      ('read1', 1, 'G', 'C', False),
      ('read1', 2, 'G', 'C', False),
      ('read2', 1, 'G', 'G', True),
      ('read2', 2, 'G', 'G', True),
      ('read3', 1, 'G', 'C', False),
      ('read3', 2, 'G', 'C', False),
  )
  def test_read_support_is_respected(self, read_name, read_number, alt_allele,
                                     read_base, supports_alt):
    """supports_alt is encoded as the 5th channel out of the 7 channels."""
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=10,
            end=11,
            reference_bases='A',
            alternate_bases=[alt_allele]),
        allele_support={
            'C': _supporting_reads('read1/1', 'read3/2'),
            'G': _supporting_reads('read2/1', 'read2/2'),
        })
    read = test_utils.make_read(
        read_base,
        start=dv_call.variant.start,
        cigar='1M',
        quals=[50],
        name=read_name)
    read.read_number = read_number
    actual = _make_encoder().encode_read(dv_call, 'TAT', read,
                                         dv_call.variant.start - 1, alt_allele)
    expected_base_values = {'C': 30, 'G': 180}
    expected_supports_alt_channel = [152, 254]
    expected = [
        expected_base_values[read_base], 254, 211, 70,
        expected_supports_alt_channel[supports_alt], 254, 1
    ]

    self.assertEqual(list(actual[0, 1]), expected)


class PileupImageCreatorEncodePileupTest(parameterized.TestCase):
  """Tests of PileupImageCreator build_pileup routine."""

  def setUp(self):
    self.alt_allele = 'C'
    self.dv_call = _make_dv_call(ref_bases='G', alt_bases=self.alt_allele)
    self.pic = _make_image_creator(
        None, None, width=3, height=4, reference_band_height=2)
    self.ref = 'AGC'
    self.read1 = test_utils.make_read('AGC', start=0, cigar='3M', name='read1')
    self.read2 = test_utils.make_read('AGC', start=1, cigar='3M', name='read2')
    self.read3 = test_utils.make_read('AGC', start=2, cigar='3M', name='read3')
    self.read4 = test_utils.make_read('AGC', start=3, cigar='3M', name='read4')

    self.expected_rows = {
        'ref':
            np.asarray(range(0, 21), np.uint8)
            .reshape(1, 3, pileup_image.DEFAULT_NUM_CHANNEL),
        'empty':
            np.zeros((1, 3, pileup_image.DEFAULT_NUM_CHANNEL), dtype=np.uint8),
        'read1':
            np.full(
                (1, 3, pileup_image.DEFAULT_NUM_CHANNEL), 1, dtype=np.uint8),
        'read2':
            np.full(
                (1, 3, pileup_image.DEFAULT_NUM_CHANNEL), 2, dtype=np.uint8),
        'read3':
            None,
        'read4':
            np.full(
                (1, 3, pileup_image.DEFAULT_NUM_CHANNEL), 3, dtype=np.uint8),
    }

    # Setup our shared mocks.
    mock_encoder = mock.Mock(spec=['encode_read', 'encode_reference'])
    mock_encoder.encode_reference.return_value = self.expected_rows['ref']

    # pylint: disable=unused-argument
    def get_read_row(dv_call, refbases, read, pos, alt_allele):
      return self.expected_rows[read.fragment_name]

    mock_encoder.encode_read.side_effect = get_read_row

    self.mock_enc_ref = mock_encoder.encode_reference
    self.mock_enc_read = mock_encoder.encode_read

    self.pic._encoder = mock_encoder

  def assertImageMatches(self, actual_image, *row_names):
    """Checks that actual_image matches an image from constructed row_names."""
    self.assertEqual(
        actual_image.shape,
        (self.pic.height, self.pic.width, pileup_image.DEFAULT_NUM_CHANNEL))
    expected_image = np.vstack([self.expected_rows[name] for name in row_names])
    npt.assert_equal(actual_image, expected_image)

  def test_image_no_reads(self):
    # This image is created just from reference and no reads. Checks that the
    # function is listening to all of our image creation parameters (e.g.,
    # reference_band_height, width, height, etc) and is filling the image with
    # empty rows when it runs out of reads.
    image = self.pic.build_pileup(self.dv_call, self.ref, [], {self.alt_allele})
    self.mock_enc_ref.assert_called_once_with(self.ref)
    test_utils.assert_not_called_workaround(self.mock_enc_read)
    self.assertImageMatches(image, 'ref', 'ref', 'empty', 'empty')

  def test_image_one_read(self):
    # We add a single read to our image.
    image = self.pic.build_pileup(self.dv_call, self.ref, [self.read1],
                                  {self.alt_allele})
    self.mock_enc_ref.assert_called_once_with(self.ref)
    self.mock_enc_read.assert_called_once_with(self.dv_call, self.ref,
                                               self.read1, 9, {self.alt_allele})
    self.assertImageMatches(image, 'ref', 'ref', 'read1', 'empty')

  def test_image_creation_with_more_reads_than_rows(self):
    # Read1 should be dropped because there's only space for Read2 and Read4.
    # If there are more reads than rows, a deterministic random subset is used.
    image = self.pic.build_pileup(self.dv_call, self.ref,
                                  [self.read1, self.read2, self.read4],
                                  {self.alt_allele})
    self.mock_enc_ref.assert_called_once_with(self.ref)
    self.assertEqual(self.mock_enc_read.call_args_list, [
        mock.call(self.dv_call, self.ref, self.read1, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read2, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read4, 9, {self.alt_allele}),
    ])
    self.assertImageMatches(image, 'ref', 'ref', 'read2', 'read4')

  def test_image_creation_with_bad_read(self):
    # Read 3 is bad (return value is None) so it should be skipped.
    image = self.pic.build_pileup(self.dv_call, self.ref,
                                  [self.read1, self.read3, self.read2],
                                  {self.alt_allele})
    self.mock_enc_ref.assert_called_once_with(self.ref)
    self.assertEqual(self.mock_enc_read.call_args_list, [
        mock.call(self.dv_call, self.ref, self.read1, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read3, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read2, 9, {self.alt_allele}),
    ])
    self.assertImageMatches(image, 'ref', 'ref', 'read1', 'read2')

  def test_image_creation_with_all_reads_in_new_order(self):
    # Read 3 is bad (return value is None) so it should be skipped. Read2 should
    # also be dropped because there's only space for Read1 and Read4. If there
    # are more reads than rows, a deterministic random subset is used.
    image = self.pic.build_pileup(
        self.dv_call, self.ref,
        [self.read2, self.read3, self.read4, self.read1], {self.alt_allele})
    self.mock_enc_ref.assert_called_once_with(self.ref)
    self.assertEqual(self.mock_enc_read.call_args_list, [
        mock.call(self.dv_call, self.ref, self.read2, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read3, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read4, 9, {self.alt_allele}),
        mock.call(self.dv_call, self.ref, self.read1, 9, {self.alt_allele}),
    ])
    self.assertImageMatches(image, 'ref', 'ref', 'read1', 'read4')


class PileupImageCreatorTest(parameterized.TestCase):

  def setUp(self):
    self.options = pileup_image.default_options()
    self.options.width = 5
    self.mock_ref_reader = mock.MagicMock(spec=reference_fai.GenomeReferenceFai)
    self.mock_ref_reader.bases.return_value = 'ACGT'
    self.mock_ref_reader.is_valid_interval.return_value = True
    self.mock_sam_reader = mock.MagicMock()
    self.mock_sam_reader.query.return_value = ['read1', 'read2']
    self.dv_call = _make_dv_call()
    self.variant = self.dv_call.variant
    self.pic = self._make_pic()

  def _make_pic(self, **kwargs):
    return pileup_image.PileupImageCreator(self.options, self.mock_ref_reader,
                                           self.mock_sam_reader, **kwargs)

  @parameterized.parameters(
      ('A', ['C'], [{'C'}]),
      ('A', ['C', 'G'], [{'C'}, {'G'}, {'C', 'G'}]),
  )
  def test_alt_combinations(self, ref, alts, expected):
    variant = variants_pb2.Variant(reference_bases=ref, alternate_bases=alts)
    self.assertEqual(expected, list(self.pic._alt_allele_combinations(variant)))

  @parameterized.parameters(
      ('A', ['C'], [{'C'}]),
      ('A', ['C', 'G'], [{'C'}, {'G'}]),
  )
  def test_alt_combinations_no_het_alt(self, ref, alts, expected):
    options = pileup_image.default_options()
    options.multi_allelic_mode = (
        deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES)
    pic = pileup_image.PileupImageCreator(options, self.mock_ref_reader,
                                          self.mock_sam_reader)
    variant = variants_pb2.Variant(reference_bases=ref, alternate_bases=alts)
    self.assertEqual(expected, list(pic._alt_allele_combinations(variant)))

  def test_get_reference_bases_good_region(self):
    self.dv_call.variant.start = 10
    region = ranges.make_range(self.variant.reference_name, 8, 13)

    actual = self.pic.get_reference_bases(self.variant)
    self.assertEqual('ACGT', actual)
    self.mock_ref_reader.is_valid_interval.assert_called_once_with(region)
    self.mock_ref_reader.bases.assert_called_once_with(region)

  def test_get_reference_bases_bad_region_returns_none(self):
    self.mock_ref_reader.is_valid_interval.return_value = False
    self.dv_call.variant.start = 3

    self.assertIsNone(self.pic.get_reference_bases(self.variant))
    test_utils.assert_called_once_workaround(
        self.mock_ref_reader.is_valid_interval)
    self.mock_ref_reader.bases.assert_not_called()

  def test_create_pileup_image_returns_none_for_bad_region(self):
    self.mock_ref_reader.is_valid_interval.return_value = False
    self.dv_call.variant.start = 3
    self.assertIsNone(self.pic.create_pileup_images(self.dv_call))
    test_utils.assert_called_once_workaround(
        self.mock_ref_reader.is_valid_interval)
    self.mock_ref_reader.bases.assert_not_called()

  def test_create_pileup_image(self):
    self.dv_call.variant.alternate_bases[:] = ['C', 'T']

    with mock.patch.object(
        self.pic, 'build_pileup', autospec=True) as mock_encoder:
      mock_encoder.side_effect = ['mi1', 'mi2', 'mi3']

      self.assertEqual([
          ({'C'}, 'mi1'),
          ({'T'}, 'mi2'),
          ({'C', 'T'}, 'mi3'),
      ], self.pic.create_pileup_images(self.dv_call))

      def _expected_call(alts):
        return mock.call(self.dv_call, self.mock_ref_reader.bases.return_value,
                         self.mock_sam_reader.query.return_value, alts)

      self.assertEqual(mock_encoder.call_count, 3)
      mock_encoder.assert_has_calls([
          _expected_call({'C'}),
          _expected_call({'T'}),
          _expected_call({'C', 'T'}),
      ])


if __name__ == '__main__':
  absltest.main()
