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
"""Tests for deepvariant.pileup_image."""

import itertools



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import numpy.testing as npt

from deepvariant import dv_constants
from deepvariant import pileup_image
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import pileup_image_native
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils


MAX_PIXEL_FLOAT = 254.0
MAX_HOMOPOLYMER_WEIGHTED = 30


def _supporting_reads(*names):
  return deepvariant_pb2.DeepVariantCall.SupportingReads(read_names=names)


def _make_dv_call(ref_bases='A', alt_bases='C'):
  return deepvariant_pb2.DeepVariantCall(
      variant=variants_pb2.Variant(
          reference_name='chr1',
          start=10,
          end=11,
          reference_bases=ref_bases,
          alternate_bases=[alt_bases],
      ),
      allele_support={'C': _supporting_reads('read1/1', 'read2/1')},
  )


def _make_dv_call_with_allele_frequency(
    ref_bases='A', alt_bases='C', alt_frequency=0.1
):
  return deepvariant_pb2.DeepVariantCall(
      variant=variants_pb2.Variant(
          reference_name='chr1',
          start=10,
          end=11,
          reference_bases=ref_bases,
          alternate_bases=[alt_bases],
      ),
      allele_support={'C': _supporting_reads('read1/1', 'read2/1')},
      allele_frequency=dict(A=1 - alt_frequency, C=alt_frequency),
  )


def _make_encoder(read_requirements=None, **kwargs):
  """Make a PileupImageEncoderNative with overrideable default options."""
  options = pileup_image.default_options(read_requirements)
  for channel in dv_constants.PILEUP_DEFAULT_CHANNELS:
    options.channels.append(channel)
    options.num_channels += 1
  options.MergeFrom(deepvariant_pb2.PileupImageOptions(**kwargs))
  return pileup_image_native.PileupImageEncoderNative(options)


def _make_encoder_with_channels(channel_set):
  options = pileup_image.default_options()
  for channel in channel_set:
    assert channel in dv_constants.CHANNELS
    options.channels.append(channel)
    options.num_channels += 1
  return pileup_image_native.PileupImageEncoderNative(options)


def _make_encoder_with_allele_frequency(read_requirements=None, **kwargs):
  """Make a PileupImageEncoderNative with overrideable default options."""
  options = pileup_image.default_options(read_requirements)
  for channel in dv_constants.PILEUP_DEFAULT_CHANNELS:
    options.channels.append(channel)
    options.num_channels += 1
  options.channels.append('allele_frequency')
  options.num_channels += 1
  options.MergeFrom(deepvariant_pb2.PileupImageOptions(**kwargs))
  return pileup_image_native.PileupImageEncoderNative(options)


def _make_encoder_with_hp_channel(
    read_requirements=None, hp_tag_for_assembly_polishing=None, **kwargs
):
  """Make a PileupImageEncoderNative with overrideable default options."""
  options = pileup_image.default_options(read_requirements)
  for channel in dv_constants.PILEUP_DEFAULT_CHANNELS:
    options.channels.append(channel)
    options.num_channels += 1
  options.channels.append('haplotype')
  options.num_channels += 1
  if hp_tag_for_assembly_polishing is not None:
    options.hp_tag_for_assembly_polishing = hp_tag_for_assembly_polishing
  options.MergeFrom(deepvariant_pb2.PileupImageOptions(**kwargs))
  return pileup_image_native.PileupImageEncoderNative(options)


class PileupImageEncoderTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.options = pileup_image.default_options()
    for channel in dv_constants.PILEUP_DEFAULT_CHANNELS:
      self.options.channels.append(channel)
      self.options.num_channels += 1

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
        ]).astype(np.uint8),
    )

  def assertImageRowEquals(self, image_row, expected):
    npt.assert_equal(image_row, expected.astype(np.uint8))

  def test_encode_read_matches(self):
    start = 10
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases
    read = test_utils.make_read(
        'ACCGT', start=start, cigar='5M', quals=range(10, 15), name='read1'
    )
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
    ]).astype(np.uint8)

    self.assertImageRowEquals(
        _make_encoder().encode_read(dv_call, 'ACAGT', read, start, alt_allele),
        full_expected,
    )

  @parameterized.parameters(
      # The corresponding colors here are defined in HPValueColor in
      # pileup_image_native.cc.
      # Reads with no HP values will resulted in a value of 0.
      (None, 0, None),
      (0, 0, None),
      (1, 254 / 2, None),
      (2, 254, None),
      # And, the color of 1 and 2 should swap when
      # hp_tag_for_assembly_polishing is set to 2.
      (None, 0, 2),
      (0, 0, 2),
      (1, 254, 2),
      (2, 254 / 2, 2),
  )
  def test_encode_read_matches_with_hp_channel(
      self, hp_value, hp_color, hp_tag_for_assembly_polishing
  ):
    # Same test case as test_encode_read_matches(), with
    # --channel_list=haplotype.
    start = 10
    dv_call = _make_dv_call_with_allele_frequency()
    alt_allele = dv_call.variant.alternate_bases
    read = test_utils.make_read(
        'ACCGT', start=start, cigar='5M', quals=range(10, 15), name='read1'
    )

    if hp_value is not None:
      read.info['HP'].values.add().int_value = hp_value

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
        # HP value
        (hp_color, hp_color, hp_color, hp_color, hp_color),
    ]).astype(np.uint8)

    self.assertImageRowEquals(
        _make_encoder_with_hp_channel(
            hp_tag_for_assembly_polishing=hp_tag_for_assembly_polishing
        ).encode_read(dv_call, 'ACAGT', read, start, alt_allele),
        full_expected,
    )

  def test_encode_read_matches_with_allele_frequency(self):
    # Same test case as test_encode_read_matches(), with allele frequency.
    start = 10
    dv_call = _make_dv_call_with_allele_frequency()
    alt_allele = dv_call.variant.alternate_bases
    read = test_utils.make_read(
        'ACCGT', start=start, cigar='5M', quals=range(10, 15), name='read1'
    )
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
        # Allele frequencies
        (203, 203, 203, 203, 203),
    ]).astype(np.uint8)

    self.assertImageRowEquals(
        _make_encoder_with_allele_frequency().encode_read(
            dv_call, 'ACAGT', read, start, alt_allele
        ),
        full_expected,
    )

  # pylint:disable=g-complex-comprehension
  @parameterized.parameters(
      (bases_start, bases_end)
      for bases_start in range(0, 5)
      for bases_end in range(6, 12)
  )
  # pylint:enable=g-complex-comprehension
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
    ]).astype(np.uint8)
    expected = np.zeros(
        (1, ref_size, self.options.num_channels), dtype=np.uint8
    )
    for i in range(read_start, read_start + len(read_bases)):
      if ref_start <= i < ref_start + ref_size:
        expected[0, i - ref_start] = full_expected[0, i - ref_start]

    read = test_utils.make_read(
        read_bases,
        start=read_start,
        cigar=str(len(read_bases)) + 'M',
        quals=read_quals,
        name='read1',
    )
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases
    self.assertImageRowEquals(
        _make_encoder().encode_read(
            dv_call, 'ACAGT', read, ref_start, alt_allele
        ),
        expected,
    )

  def test_encode_read_deletion(self):
    # ref:  AACAG
    # read: AA--G
    start = 2
    read = test_utils.make_read(
        'AAG', start=start, cigar='2M2D1M', quals=range(10, 13), name='read1'
    )
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases
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
    ]).astype(np.uint8)
    self.assertImageRowEquals(
        _make_encoder().encode_read(dv_call, 'AACAG', read, start, alt_allele),
        full_expected,
    )

  def test_encode_read_insertion(self):
    # ref:  AA-CAG
    # read: AAACAG
    start = 2
    read = test_utils.make_read(
        'AAACAG', start=start, cigar='2M1I3M', quals=range(10, 16), name='read1'
    )
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases
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
    ]).astype(np.uint8)
    self.assertImageRowEquals(
        _make_encoder().encode_read(dv_call, 'AACAG', read, start, alt_allele),
        full_expected,
    )

  @parameterized.parameters(
      (min_base_qual, min_mapping_qual)
      for min_base_qual, min_mapping_qual in itertools.product(
          range(0, 5), range(0, 5)
      )
  )
  def test_ignores_reads_with_low_quality_bases(
      self, min_base_qual, min_mapping_qual
  ):
    """Check that we discard reads with low quality bases at variant start site.

    We have the following scenario:

    position    0    1    2    3    4    5
    reference        A    A    C    A    G
    read             A    A    A
    variant               C

    We set the base quality of the middle base in the read to different values
    of `base_qual`. Since the middle position of the read is where the variant
    starts, the read should only be kept if `base_qual` >= `min_base_qual`.

    Args:
      min_base_qual: Reads are discarded if the base at a variant start position
        does not meet this base quality requirement.
      min_mapping_qual: Reads are discarded if they do not meet this mapping
        quality requirement.
    """
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=2,
            end=3,
            reference_bases='A',
            alternate_bases=['C'],
        )
    )

    read_requirements = reads_pb2.ReadRequirements(
        min_base_quality=min_base_qual,
        min_mapping_quality=min_mapping_qual,
        min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT,
    )
    pie = _make_encoder(read_requirements=read_requirements)

    for base_qual in range(min_base_qual + 5):
      quals = [min_base_qual, base_qual, min_base_qual]
      read = test_utils.make_read(
          'AAA', start=1, cigar='3M', quals=quals, mapq=min_mapping_qual
      )
      actual = pie.encode_read(dv_call, 'AACAG', read, 1, ['C'])
      if base_qual < min_base_qual:
        self.assertIsNone(actual)
      else:
        self.assertIsNotNone(actual)

  @parameterized.parameters(
      (min_base_qual, min_mapping_qual)
      for min_base_qual, min_mapping_qual in itertools.product(
          range(0, 5), range(0, 5)
      )
  )
  def test_keeps_reads_with_low_quality_bases(
      self, min_base_qual, min_mapping_qual
  ):
    """Check that we keep reads with adequate quality at variant start position.

    We have the following scenario:

    position    0    1    2    3    4    5
    reference        A    A    C    A    G
    read             A    A    A
    variant               C

    We set the base quality of the first and third bases in the read to
    different functions of `base_qual`. The middle position of the read is
    where the variant starts, and this position always has base quality greater
    than `min_base_qual`. Thus, the read should always be kept.

    Args:
      min_base_qual: Reads are discarded if the base at a variant start position
        does not meet this base quality requirement.
      min_mapping_qual: Reads are discarded if they do not meet this mapping
        quality requirement.
    """
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=2,
            end=3,
            reference_bases='A',
            alternate_bases=['C'],
        )
    )

    read_requirements = reads_pb2.ReadRequirements(
        min_base_quality=min_base_qual,
        min_mapping_quality=min_mapping_qual,
        min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT,
    )
    pie = _make_encoder(read_requirements=read_requirements)

    for base_qual in range(min_base_qual + 5):
      quals = [base_qual - 1, min_base_qual, base_qual + 1]
      read = test_utils.make_read(
          'AAA', start=1, cigar='3M', quals=quals, mapq=min_mapping_qual
      )
      actual = pie.encode_read(dv_call, 'AACAG', read, 1, ['C'])
      self.assertIsNotNone(actual)

  @parameterized.parameters(
      (min_base_qual, min_mapping_qual)
      for min_base_qual, min_mapping_qual in itertools.product(
          range(0, 5), range(0, 5)
      )
  )
  def test_ignores_reads_with_low_mapping_quality(
      self, min_base_qual, min_mapping_qual
  ):
    """Check that we discard reads with low mapping quality.

    We have the following scenario:

    position    0    1    2    3    4    5
    reference        A    A    C    A    G
    read             A    A    A
    variant               C

    We set the mapping quality of the read to different values of
    `mapping_qual`. All bases in the read have base quality greater than
    `min_base_qual`. The read should only be kept if
    `mapping_qual` > `min_mapping_qual`.

    Args:
      min_base_qual: Reads are discarded if the base at a variant start position
        does not meet this base quality requirement.
      min_mapping_qual: Reads are discarded if they do not meet this mapping
        quality requirement.
    """
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=2,
            end=3,
            reference_bases='A',
            alternate_bases=['C'],
        )
    )

    read_requirements = reads_pb2.ReadRequirements(
        min_base_quality=min_base_qual,
        min_mapping_quality=min_mapping_qual,
        min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT,
    )
    pie = _make_encoder(read_requirements=read_requirements)

    for mapping_qual in range(min_mapping_qual + 5):
      quals = [min_base_qual, min_base_qual, min_base_qual]
      read = test_utils.make_read(
          'AAA', start=1, cigar='3M', quals=quals, mapq=mapping_qual
      )
      actual = pie.encode_read(dv_call, 'AACAG', read, 1, ['C'])
      if mapping_qual < min_mapping_qual:
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
  def test_read_support_is_respected(
      self, read_name, read_number, alt_allele, read_base, supports_alt
  ):
    """supports_alt is encoded as the 5th channel out of the 7 channels."""
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=10,
            end=11,
            reference_bases='A',
            alternate_bases=['C', 'G'],
        ),
        allele_support={
            'C': _supporting_reads('read1/1', 'read3/2'),
            'G': _supporting_reads('read2/1', 'read2/2'),
        },
    )
    read = test_utils.make_read(
        read_base,
        start=dv_call.variant.start,
        cigar='1M',
        quals=[50],
        name=read_name,
    )
    read.read_number = read_number
    actual = _make_encoder().encode_read(
        dv_call, 'TAT', read, dv_call.variant.start - 1, [alt_allele]
    )
    expected_base_values = {'C': 30, 'G': 180}
    expected_supports_alt_channel = [152, 254]
    expected = [
        expected_base_values[read_base],
        254,
        211,
        70,
        expected_supports_alt_channel[supports_alt],
        254,
    ]

    self.assertEqual(list(actual[0, 1]), expected)

  @parameterized.parameters(
      ('read1', 1, 'C', 'C', True, int(254.0 * 1.0)),
      # This read isn't present in allele support for 'C'.
      ('read1', 2, 'C', 'C', True, int(254.0 * 0.6)),
      ('read2', 1, 'C', 'G', True, int(254.0 * 0.3)),
      ('read1', 1, 'C', 'C', False, int(254.0 * 1.0)),
      # This read isn't present in allele support for 'C'.
      ('read1', 2, 'C', 'C', False, int(254.0 * 0.6)),
      ('read2', 1, 'C', 'G', False, int(254.0 * 0.6)),
  )
  def test_read_support_multiallelic(
      self,
      read_name,
      read_number,
      alt_allele,
      read_base,
      add_supporting_other_alt_color,
      expected_color,
  ):
    """supports_alt is encoded as the 5th channel out of the 7 channels."""
    dv_call = deepvariant_pb2.DeepVariantCall(
        variant=variants_pb2.Variant(
            reference_name='chr1',
            start=10,
            end=11,
            reference_bases='A',
            alternate_bases=['C', 'G'],
        ),
        allele_support={
            'C': _supporting_reads('read1/1'),
            'G': _supporting_reads('read2/1', 'read2/2'),
        },
    )
    read = test_utils.make_read(
        read_base,
        start=dv_call.variant.start,
        cigar='1M',
        quals=[50],
        name=read_name,
    )
    read.read_number = read_number

    if add_supporting_other_alt_color:
      other_allele_supporting_read_alpha = 0.3
    else:
      other_allele_supporting_read_alpha = 0.6

    pie = _make_encoder(
        other_allele_supporting_read_alpha=other_allele_supporting_read_alpha
    )
    actual = pie.encode_read(
        dv_call, 'TAT', read, dv_call.variant.start - 1, [alt_allele]
    )
    self.assertEqual(actual[0, 1, 4], expected_color)


class PileupCustomChannels(absltest.TestCase):

  def get_encoded_read(
      self, channel_set, cigar='20M5D20M5S', fragment_length=10
  ):
    start = 500
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases[0]

    read = test_utils.make_read(
        'TTTTATGACAAAAAAGATGCGACGGTTCCGTAACCCATAAGAAAGAACGT',
        start=start,
        cigar=cigar,
        quals=range(1, 51),
        name='read1',
        fragment_length=fragment_length,
    )
    return _make_encoder_with_channels(channel_set).encode_read(
        dv_call,
        'TTTTATGACAAAAAAGATGCGACGGTTCCGTAACCCATAAGAAAGAACGT',
        read,
        start,
        [alt_allele],
    )

  def test_read_mapping_percent(self):
    result = self.get_encoded_read(['read_mapping_percent'])
    ch1 = result[:, :, 0]
    self.assertSetEqual(set(np.unique(ch1)), set([0, 203]))

  def test_avg_mapping_quality(self):
    result = self.get_encoded_read(['avg_base_quality'])
    ch1 = result[:, :, 0]
    self.assertSetEqual(set(np.unique(ch1)), set([0, 68]))

  def test_identity(self):
    result = self.get_encoded_read(['identity'], cigar='5M20D20M5S')
    ch1 = result[:, :, 0]
    self.assertSetEqual(set(np.unique(ch1)), set([0, 127]))

  def test_gap_compressed_identity(self):
    result = self.get_encoded_read(
        ['gap_compressed_identity'], cigar='5M20D20M5S'
    )
    ch1 = result[:, :, 0]
    self.assertSetEqual(set(np.unique(ch1)), set([0, 243]))

  def test_blank(self):
    result = self.get_encoded_read(['blank'])
    ch1 = result[:, :, 0]
    self.assertSetEqual(set(np.unique(ch1)), set([0]))

  def test_insert_size(self):
    result = self.get_encoded_read(['insert_size'], fragment_length=22)
    ch1 = result[:, :, 0]
    self.assertSetEqual(set(np.unique(ch1)), set([0, 5]))

  def test_multi(self):
    result = self.get_encoded_read(
        ['read_mapping_percent', 'gap_compressed_identity', 'blank']
    )
    self.assertEqual(result.shape, (1, 50, 3))


class PileupCustomChannelsParam(parameterized.TestCase):

  def make_pileup(self, seq, channels):
    start = 500
    dv_call = _make_dv_call()
    alt_allele = dv_call.variant.alternate_bases[0]
    read = test_utils.make_read(
        seq,
        start=start,
        cigar=f'{len(seq)}M',
        quals=range(1, len(seq) + 1),
        name='read1',
    )
    result = _make_encoder_with_channels(channels).encode_read(
        dv_call, seq, read, start, [alt_allele]
    )
    return result

  @parameterized.parameters(
      ('GC', 1.0),
      ('GAC', 0.66),
      ('GGAA', 0.50),
      ('ATTCTGTTAA', 0.20),
      ('TTTTTTTTTT', 0.00),
  )
  def test_gc_content(self, seq, exp_gc_content):
    result = self.make_pileup(seq, ['gc_content'])
    ch7 = result[:, :, 0][0]
    gc_content = ch7.max() / MAX_PIXEL_FLOAT
    self.assertAlmostEqual(gc_content, exp_gc_content, 2)

  @parameterized.parameters(
      ('AAATTCCC', [1, 1, 1, 0, 0, 1, 1, 1]),
      ('ATCGTTCCC', [0, 0, 0, 0, 0, 0, 1, 1, 1]),
      ('ATTCCCTTA', [0, 0, 0, 1, 1, 1, 0, 0, 0]),
      ('ATCG', [0, 0, 0, 0]),
      ('AATTCCGG', [0, 0, 0, 0, 0, 0, 0, 0]),
      ('AAAAAAAA', [1, 1, 1, 1, 1, 1, 1, 1]),
  )
  def test_is_homopolymer(self, seq, expected_result):
    result = self.make_pileup(seq, ['is_homopolymer'])
    ch7 = (result[:, :, 0][0] / MAX_PIXEL_FLOAT).astype(int)
    self.assertTrue((ch7 == expected_result).all())

  @parameterized.parameters(
      ('AAATTCCC', [3, 3, 3, 2, 2, 3, 3, 3]),
      ('ATCGTTCCC', [1, 1, 1, 1, 2, 2, 3, 3, 3]),
      ('ATTCCCTTA', [1, 2, 2, 3, 3, 3, 2, 2, 1]),
  )
  def test_weighted_homopolymer(self, seq, expected_result):
    result = self.make_pileup(seq, ['homopolymer_weighted'])
    ch7 = np.round(
        (result[:, :, 0][0] / MAX_PIXEL_FLOAT) * MAX_HOMOPOLYMER_WEIGHTED
    ).astype(int)
    self.assertTrue((ch7 == expected_result).all())


if __name__ == '__main__':
  absltest.main()
