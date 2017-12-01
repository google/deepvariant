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
"""Encodes reference and read data into a PileupImage for DeepVariant."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools



import numpy as np

from deepvariant.core import ranges
from deepvariant.core import utils
from deepvariant.core.protos import core_pb2
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import pileup_image_native

DEFAULT_MIN_BASE_QUALITY = 10
DEFAULT_MIN_MAPPING_QUALITY = 10
# redacted
DEFAULT_NUM_CHANNEL = 7


def default_options(read_requirements=None):
  """Creates a PileupImageOptions populated with good default values."""
  if not read_requirements:
    read_requirements = core_pb2.ReadRequirements(
        min_base_quality=DEFAULT_MIN_BASE_QUALITY,
        min_mapping_quality=DEFAULT_MIN_MAPPING_QUALITY,
        min_base_quality_mode=core_pb2.ReadRequirements.ENFORCED_BY_CLIENT)

  return deepvariant_pb2.PileupImageOptions(
      reference_band_height=5,
      base_color_offset_a_and_g=40,
      base_color_offset_t_and_c=30,
      base_color_stride=70,
      allele_supporting_read_alpha=1.0,
      allele_unsupporting_read_alpha=0.6,
      reference_matching_read_alpha=0.2,
      reference_mismatching_read_alpha=1.0,
      indel_anchoring_base_char='*',
      reference_alpha=0.4,
      reference_base_quality=60,
      positive_strand_color=70,
      negative_strand_color=240,
      base_quality_cap=40,
      mapping_quality_cap=60,
      height=100,
      width=221,
      read_overlap_buffer_bp=5,
      read_requirements=read_requirements,
      multi_allelic_mode=deepvariant_pb2.PileupImageOptions.ADD_HET_ALT_IMAGES,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=2101079370)


def _empty_image_row(width):
  """Creates an empty image row as an uint8 np.array."""
  return np.zeros((1, width, DEFAULT_NUM_CHANNEL), dtype=np.uint8)


def _compute_half_width(width):
  return int((width - 1) / 2)


class PileupImageCreator(object):
  """High-level API for creating images of pileups of reads and reference bases.

  This class provides a higher-level and more natural API for constructing
  images at a candidate variant call site. Given a DeepVariantCall, which
  contains the candidate variant call along with key supplementary information,
  this class provides create_pileup_images() that will do all of the necessary
  fetching of reads and reference bases from readers and pass those off to the
  lower-level PileupImageEncoder to construct the image Tensor.

  for dv_call in candidates:
    allele_and_images = pic.create_pileup_images(dv_call)
    ...

  A quick note on how we deal with multiple alt alleles:

  Suppose variant has ref and two alt alleles. Assuming the sample is diploid,
  we have the following six possible genotypes:

    ref/ref   => 0/0
    ref/alt1  => 0/1
    alt1/alt1 => 1/1
    ref/alt2  => 0/2
    alt1/alt2 => 1/2
    alt2/alt2 => 2/2

  In DeepVariant we predict the genotype count (0, 1, 2) for a specific set of
  alternate alleles. If we only had a single alt, we'd construct an image for
  ref vs. alt1:

    image1 => ref vs. alt1 => determine if we are 0/0, 0/1, 1/1

  If we add a second image for alt2, we get:

    image2 => ref vs. alt2 => determine if we are 0/0, 0/2, 2/2

  but the problem here is that we don't have a good estimate for the het-alt
  state 1/2. So we construct a third image contrasting ref vs. either alt1 or
  alt2:

    image3 => ref vs. alt1 or alt2 => determines 0/0, 0/{1,2}, {1,2}/{1,2}

  Given the predictions for each image:

    image1 => p00, p01, p11
    image2 => p00, p02, p22
    image3 => p00, p0x, pxx where x is {1,2}

  we calculate our six genotype likelihoods as:

    0/0 => p00 [from any image]
    0/1 => p01 [image1]
    1/1 => p11 [image1]
    0/2 => p02 [image2]
    2/2 => p22 [image2]
    1/2 => pxx [image3]

  The function create_pileup_images() returns all of the necessary images, along
  with the alt alleles used for each image.
  """

  def __init__(self, options, ref_reader, sam_reader):
    self._options = options
    self._encoder = pileup_image_native.PileupImageEncoderNative(self._options)
    self._ref_reader = ref_reader
    self._sam_reader = sam_reader
    self._random = np.random.RandomState(self._options.random_seed)

  def __getattr__(self, attr):
    """Gets attributes from self._options as though they are our attributes."""
    return self._options.__getattribute__(attr)

  @property
  def half_width(self):
    return _compute_half_width(self._options.width)

  @property
  def max_reads(self):
    return self.height - self.reference_band_height

  def get_reads(self, variant):
    """Gets the reads used to construct the pileup image around variant.

    Args:
      variant: A learning.genomics.deepvariant.core.genomics.Variant proto
        describing the variant we are creating the pileup image of.

    Returns:
      A list of learning.genomics.deepvariant.core.genomics.Read protos.
    """
    query_start = variant.start - self._options.read_overlap_buffer_bp
    query_end = variant.end + self._options.read_overlap_buffer_bp
    region = ranges.make_range(variant.reference_name, query_start, query_end)
    return list(self._sam_reader.query(region))

  def get_reference_bases(self, variant):
    """Gets the reference bases used to make the pileup image around variant.

    Args:
      variant: A learning.genomics.deepvariant.core.genomics.Variant proto
        describing the variant we are creating the pileup image of.

    Returns:
      A string of reference bases or None. Returns None if the reference
      interval for variant isn't valid for some reason.
    """
    start = variant.start - self.half_width
    end = start + self._options.width
    region = ranges.make_range(variant.reference_name, start, end)
    if self._ref_reader.is_valid_interval(region):
      return self._ref_reader.bases(region)
    else:
      return None

  def _alt_allele_combinations(self, variant):
    """Yields the set of all alt_alleles for variant.

    This function computes the sets of alt_alleles we want to use to cover all
    genotype likelihood calculations we need for n alt alleles (see class docs
    for background). The easiest way to do this is to calculate all combinations
    of 2 alleles from ref + alts and then strip away the reference alleles,
    leaving us with the set of alts for the pileup image encoder.

    Args:
      variant: learning.genomics.deepvariant.core.genomics.Variant to
        generate the alt allele combinations for.

    Yields:
      A series of sets containing the alt alleles we want to use for a single
      pileup image. The entire series covers all combinations of alt alleles
      needed for variant.

    Raises:
      ValueError: if options.multi_allelic_mode is UNSPECIFIED.
    """
    ref = variant.reference_bases
    alts = list(variant.alternate_bases)

    if (self.multi_allelic_mode ==
        deepvariant_pb2.PileupImageOptions.UNSPECIFIED):
      raise ValueError('multi_allelic_mode cannot be UNSPECIFIED')
    elif (self.multi_allelic_mode ==
          deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES):
      for alt in alts:
        yield set([alt])
    else:
      for combination in itertools.combinations([ref] + alts, 2):
        yield set(combination) - {ref}

  def build_pileup(self, dv_call, refbases, reads, alt_alleles):
    """Creates a pileup tensor for dv_call.

    Args:
      dv_call: learning.genomics.deepvariant.DeepVariantCall object with
        information on our candidate call and allele support information.
      refbases: A string options.width in length containing the reference base
        sequence to encode. The middle base of this string should be at the
        start of the variant in dv_call.
      reads: Iterable of learning.genomics.deepvariant.core.genomics.Read
        objects that we'll use to
        encode the read information supporting our call. Assumes each read is
        aligned and is well-formed (e.g., has bases and quality scores, cigar).
        Rows of the image are encoded in the same order as reads.
      alt_alleles: A collection of alternative_bases from dv_call.variant that
        we are treating as "alt" when constructing this pileup image. A read
        will be considered supporting the "alt" allele if it occurs in the
        support list for any alt_allele in this collection.

    Returns:
      A [self.width, self.height, DEFAULT_NUM_CHANNEL] uint8 Tensor image.

    Raises:
      ValueError: if any arguments are invalid.
    """
    if len(refbases) != self.width:
      raise ValueError('refbases is {} long but width is {}'.format(
          len(refbases), self.width))

    if not alt_alleles:
      raise ValueError('alt_alleles cannot be empty')
    if any(alt not in dv_call.variant.alternate_bases for alt in alt_alleles):
      raise ValueError('all elements of alt_alleles must be the alternate bases'
                       ' of dv_call.variant', alt_alleles, dv_call.variant)

    image_start_pos = dv_call.variant.start - self.half_width
    if (len(dv_call.variant.reference_bases) == 1 and
        refbases[self.half_width] != dv_call.variant.reference_bases):
      raise ValueError('center of refbases doesnt match variant.refbases',
                       self.half_width, refbases[self.half_width],
                       dv_call.variant)

    # We start with n copies of our encoded reference bases.
    rows = (
        [self._encoder.encode_reference(refbases)] * self.reference_band_height)

    # A generator that yields tuples of the form (position, row), iff the read
    # can be encoded as a valid row to be used in the pileup image.
    def _row_generator():
      for read in reads:
        read_row = self._encoder.encode_read(dv_call, refbases, read,
                                             image_start_pos, alt_alleles)
        if read_row is not None:
          yield read.alignment.position.position, read_row

    # We add a row for each read in order, down-sampling if the number of reads
    # is greater than self.max_reads. Sort the reads by their alignment
    # position.
    sample = sorted(
        utils.reservoir_sample(
            _row_generator(), self.max_reads, random=self._random),
        key=lambda x: x[0])

    rows += [read_row for _, read_row in sample]

    # Finally, fill in any missing rows to bring our image to self.height rows
    # with empty (all black) pixels.
    n_missing_rows = self.height - len(rows)
    if n_missing_rows > 0:
      # Add values to rows to fill it out with zeros.
      rows += [_empty_image_row(len(refbases))] * n_missing_rows

    # Vertically stack the image rows to create a single
    # h x w x DEFAULT_NUM_CHANNEL image.
    return np.vstack(rows)

  def create_pileup_images(self, dv_call):
    """Creates a DeepVariant TF.Example for the DeepVariant call dv_call.

    See class documents for more details.

    Args:
      dv_call: A learning.genomics.deepvariant.DeepVariantCall proto
        that we want to create a TF.Example pileup image of.

    Returns:
      A list of tuples. The first element of the tuple is a set of alternate
      alleles used as 'alt' when encoding this image. The second element is a
      [w, h, DEFAULT_NUM_CHANNEL] uint8 Tensor of the pileup image for those
      alt alleles.
    """
    variant = dv_call.variant
    ref = self.get_reference_bases(variant)
    if not ref:
      # This interval isn't valid => we off the edge of the chromosome so return
      # None to indicate we couldn't process this variant.
      return None
    reads = self.get_reads(variant)

    def _make_one(alt_alleles):
      image = self.build_pileup(dv_call, ref, reads, alt_alleles)
      return alt_alleles, image

    return [_make_one(alts) for alts in self._alt_allele_combinations(variant)]
