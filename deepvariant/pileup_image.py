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
"""Encodes reference and read data into a PileupImage for DeepVariant."""

import itertools
from typing import Iterable, List, Optional



from absl import logging
import numpy as np

from deepvariant import dv_constants
from deepvariant import sample as sample_lib
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import pileup_image_native
from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges


def default_options(read_requirements=None):
  """Creates a PileupImageOptions populated with good default values."""
  if not read_requirements:
    read_requirements = reads_pb2.ReadRequirements(
        min_base_quality=10,
        min_mapping_quality=10,
        min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT,
    )

  return deepvariant_pb2.PileupImageOptions(
      reference_band_height=5,
      base_color_offset_a_and_g=40,
      base_color_offset_t_and_c=30,
      base_color_stride=70,
      allele_supporting_read_alpha=1.0,
      allele_unsupporting_read_alpha=0.6,
      other_allele_supporting_read_alpha=0.6,
      reference_matching_read_alpha=0.2,
      reference_mismatching_read_alpha=1.0,
      indel_anchoring_base_char='*',
      reference_alpha=0.4,
      reference_base_quality=60,
      positive_strand_color=70,
      negative_strand_color=240,
      base_quality_cap=40,
      mapping_quality_cap=60,
      height=dv_constants.PILEUP_DEFAULT_HEIGHT,
      width=dv_constants.PILEUP_DEFAULT_WIDTH,
      num_channels=dv_constants.PILEUP_NUM_CHANNELS,
      read_overlap_buffer_bp=5,
      read_requirements=read_requirements,
      multi_allelic_mode=deepvariant_pb2.PileupImageOptions.ADD_HET_ALT_IMAGES,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=2101079370,
      sequencing_type=deepvariant_pb2.PileupImageOptions.UNSPECIFIED_SEQ_TYPE,
      alt_aligned_pileup='none',
      types_to_alt_align='indels',
      min_non_zero_allele_frequency=0.00001,
      use_allele_frequency=False,
  )


def _compute_half_width(width):
  return int((width - 1) / 2)


def _represent_alt_aligned_pileups(representation, ref_image, alt_images):
  """Combines ref and alt-aligned pileup images according to the representation.

  Args:
    representation: string, one of "rows", "base_channels", "diff_channels".
    ref_image: 3D numpy array. The original pileup image.
    alt_images: list of either one or two 3D numpy arrays, both of the same
      dimensions as ref_image. Pileup image(s) of the same reads aligned to the
      alternate haplotype(s).

  Returns:
    One 3D numpy array containing a selection of data from the input arrays.
  """

  # If there is only one alt, duplicate it to make all pileups the same size.
  if len(alt_images) == 1:
    alt_images = alt_images + alt_images
  if len(alt_images) != 2:
    raise ValueError('alt_images must contain exactly one or two arrays.')

  # Ensure that all three pileups have the same shape.
  if not ref_image.shape == alt_images[0].shape == alt_images[1].shape:
    raise ValueError(
        'Pileup images must be the same shape to be combined. '
        'ref_image.shape is {}. alt_images[0].shape is {}. '
        'alt_images[1].shape is {}.'.format(
            ref_image.shape, alt_images[0].shape, alt_images[1].shape
        )
    )

  if representation == 'rows':
    # Combine all images: [ref, alt1, alt2].
    return np.concatenate([ref_image] + alt_images, axis=0)
  elif representation == 'base_channels':
    channels = [ref_image[:, :, c] for c in range(ref_image.shape[2])]
    # Add channel 0 (bases ATCG) of both alts as channels.
    channels.append(alt_images[0][:, :, 0])
    channels.append(alt_images[1][:, :, 0])
    return np.stack(channels, axis=2)
  elif representation == 'diff_channels':
    channels = [ref_image[:, :, c] for c in range(ref_image.shape[2])]
    # Add channel 5 (base differs from ref) of both alts as channels.
    channels.append(alt_images[0][:, :, 5])
    channels.append(alt_images[1][:, :, 5])
    return np.stack(channels, axis=2)
  else:
    raise ValueError(
        '_represent_alt_aligned_pileups received invalid value for '
        'representation: "{}". Must be one of '
        'rows, base_channels, or diff_channels.'.format(representation)
    )


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

  def __init__(
      self,
      options: deepvariant_pb2.PileupImageOptions,
      ref_reader: fasta.IndexedFastaReader,
      samples: List[sample_lib.Sample],
  ):
    self._options = options
    self._encoder = pileup_image_native.PileupImageEncoderNative(self._options)
    self._channels_enum = self._encoder.all_channels_enum(
        options.alt_aligned_pileup
    )
    self._ref_reader = ref_reader
    self._samples = samples

  def __getattr__(self, attr):
    """Gets attributes from self._options as though they are our attributes."""
    return self._options.__getattribute__(attr)

  @property
  def half_width(self):
    return _compute_half_width(self._options.width)

  def get_channels(self):
    return self._channels_enum

  def get_reads(
      self, variant: variants_pb2.Variant, sam_reader: sam.SamReader
  ) -> List[reads_pb2.Read]:
    """Gets the reads used to construct the pileup image around variant.

    Args:
      variant: A third_party.nucleus.protos.Variant proto describing the variant
        we are creating the pileup image of.
      sam_reader: Nucleus sam_reader from which to query.

    Returns:
      A list of third_party.nucleus.protos.Read protos.
    """
    query_start = variant.start - self._options.read_overlap_buffer_bp
    query_end = variant.end + self._options.read_overlap_buffer_bp
    region = ranges.make_range(variant.reference_name, query_start, query_end)
    return list(sam_reader.query(region))

  def get_reference_bases(self, variant: variants_pb2.Variant) -> Optional[str]:
    """Gets the reference bases used to make the pileup image around variant.

    Args:
      variant: A third_party.nucleus.protos.Variant proto describing the variant
        we are creating the pileup image of.

    Returns:
      A string of reference bases or None. Returns None if the reference
      interval for variant isn't valid for some reason.
    """
    start = variant.start - self.half_width
    end = start + self._options.width
    region = ranges.make_range(variant.reference_name, start, end)
    if self._ref_reader.is_valid(region):
      return self._ref_reader.query(region)
    else:
      return None

  def _alt_allele_combinations(
      self, variant: variants_pb2.Variant
  ) -> Iterable[List[str]]:
    """Yields the set of all alt_alleles for variant.

    This function computes the sets of alt_alleles we want to use to cover all
    genotype likelihood calculations we need for n alt alleles (see class docs
    for background). The easiest way to do this is to calculate all combinations
    of 2 alleles from ref + alts and then strip away the reference alleles,
    leaving us with the set of alts for the pileup image encoder. The sets are
    converted to sorted lists at the end for downstream consistency.

    Args:
      variant: third_party.nucleus.protos.Variant to generate the alt allele
        combinations for.

    Yields:
      A series of lists containing the alt alleles we want to use for a single
      pileup image. The entire series covers all combinations of alt alleles
      needed for variant.

    Raises:
      ValueError: if options.multi_allelic_mode is UNSPECIFIED.
    """
    ref = variant.reference_bases
    alts = list(variant.alternate_bases)

    if (
        self.multi_allelic_mode
        == deepvariant_pb2.PileupImageOptions.UNSPECIFIED
    ):
      raise ValueError('multi_allelic_mode cannot be UNSPECIFIED')
    elif (
        self.multi_allelic_mode
        == deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES
    ):
      for alt in alts:
        yield sorted([alt])
    else:
      for combination in itertools.combinations([ref] + alts, 2):
        yield sorted(list(set(combination) - {ref}))

  def build_pileup(
      self,
      dv_call: deepvariant_pb2.DeepVariantCall,
      refbases: str,
      reads_for_samples: List[List[reads_pb2.Read]],
      alt_alleles: List[str],
      sample_order: Optional[List[int]] = None,
      custom_ref: bool = False,
  ):
    """Creates a pileup tensor for dv_call.

    Args:
      dv_call: learning.genomics.deepvariant.DeepVariantCall object with
        information on our candidate call and allele support information.
      refbases: A string options.width in length containing the reference base
        sequence to encode. The middle base of this string should be at the
        start of the variant in dv_call.
      reads_for_samples: list by sample of Iterable of
        third_party.nucleus.protos.Read objects that we'll use to encode the
        read information supporting our call. Assumes each read is aligned and
        is well-formed (e.g., has bases and quality scores, cigar). Rows of the
        image are encoded in the same order as reads.
      alt_alleles: A collection of alternative_bases from dv_call.variant that
        we are treating as "alt" when constructing this pileup image. A read
        will be considered supporting the "alt" allele if it occurs in the
        support list for any alt_allele in this collection.
      sample_order: A list of indices representing the order in which samples
        should be represented in the pileup image. Example: [1,0,2] to swap the
        first two samples out of three. This is None by default which puts the
        samples in order.
      custom_ref: True if refbases should not be checked for matching against
        variant's reference_bases.

    Returns:
      A uint8 Tensor image of shape
        [self.width, <sum of sample pileup heights>, DEFAULT_NUM_CHANNEL]

    Raises:
      ValueError: if any arguments are invalid.
    """
    if len(refbases) != self.width:
      raise ValueError(
          'refbases is {} long but width is {}'.format(
              len(refbases), self.width
          )
      )

    if not alt_alleles:
      raise ValueError('alt_alleles cannot be empty')
    if any(alt not in dv_call.variant.alternate_bases for alt in alt_alleles):
      raise ValueError(
          (
              'all elements of alt_alleles must be the alternate bases'
              ' of dv_call.variant'
          ),
          alt_alleles,
          dv_call.variant,
      )
    if len(self._samples) != len(reads_for_samples):
      raise ValueError(
          'The number of self._samples ({}) must be the same as the number of '
          'reads_for_samples ({}).'.format(
              len(self._samples), len(reads_for_samples)
          )
      )

    image_start_pos = dv_call.variant.start - self.half_width
    if not custom_ref and (
        refbases[self.half_width] != dv_call.variant.reference_bases[0]
    ):
      raise ValueError(
          'The middle base of reference sequence in the window '
          "({} at base {}) doesn't match first "
          'character of variant.reference_bases ({}).'.format(
              refbases[self.half_width],
              self.half_width,
              dv_call.variant.reference_bases,
          )
      )

    sample_sections = []
    if sample_order is None:
      sample_order = range(len(self._samples))
    for i in sample_order:
      sample = self._samples[i]
      sample_sections.extend(
          self._encoder.build_pileup_for_one_sample(
              dv_call,
              refbases,
              reads_for_samples[i],
              image_start_pos,
              alt_alleles,
              sample.options,
          )
      )

    # Vertically stack the image rows to create a single
    # h x w x DEFAULT_NUM_CHANNEL image.
    return np.vstack(sample_sections)

  def _empty_image_row(self) -> np.ndarray:
    """Creates an empty image row as an uint8 np.array."""
    return np.zeros((1, self.width, self.num_channels), dtype=np.uint8)

  def create_pileup_images(
      self,
      dv_call,
      reads_for_samples,
      sample_order=None,
      haplotype_alignments_for_samples=None,
      haplotype_sequences=None,
  ):
    """Creates a DeepVariant TF.Example for the DeepVariant call dv_call.

    See class documents for more details.

    Args:
      dv_call: A learning.genomics.deepvariant.DeepVariantCall proto that we
        want to create a TF.Example pileup image of.
      reads_for_samples: list of read generators, one for each sample.
      sample_order: A list of indices representing the order in which samples
        should be represented in the pileup image. Example: [1,0,2] to swap the
        first two samples out of three. This is None by default which puts the
        samples in order.
      haplotype_alignments_for_samples: list with a dict for each sample of read
        alignments keyed by haplotype.
      haplotype_sequences: dict of sequences keyed by haplotype.

    Returns:
      A list of tuples. The first element of the tuple is a set of alternate
      alleles used as 'alt' when encoding this image. The second element is a
      [w, h, DEFAULT_NUM_CHANNEL] uint8 Tensor of the pileup image for those
      alt alleles.
    """
    variant = dv_call.variant
    # Ref bases to show at the top of the pileup:
    ref_bases = self.get_reference_bases(variant)
    if not ref_bases:
      # This interval isn't valid => we are off the edge of the chromosome, so
      # return None to indicate we couldn't process this variant.
      return None

    alt_aligned_representation = self._options.alt_aligned_pileup

    def _pileup_for_pair_of_alts(alt_alleles):
      """Create pileup image for one combination of alt alleles."""
      # Always create the ref-aligned pileup image.
      ref_image = self.build_pileup(
          dv_call=dv_call,
          refbases=ref_bases,
          reads_for_samples=reads_for_samples,
          alt_alleles=alt_alleles,
          sample_order=sample_order,
      )
      # Optionally also create pileup images with reads aligned to alts.
      if alt_aligned_representation != 'none':
        if (
            haplotype_alignments_for_samples is None
            or haplotype_sequences is None
        ):
          # Use sample height or default to pic height.
          sample_heights = [
              sample.options.pileup_height for sample in self._samples
          ]
          if None not in sample_heights:
            pileup_height = sum(sample_heights)
          else:
            pileup_height = self.height
          pileup_shape = (pileup_height, self.width, self.num_channels)
          alt_images = [
              np.zeros(pileup_shape, dtype=np.uint8) for alt in alt_alleles
          ]
        else:
          alt_images = []
          for alt in alt_alleles:
            if len(haplotype_sequences[alt]) != self.width:
              logging.warning(
                  (
                      'haplotype_sequences[alt] is %d long but pileup '
                      'image width is %d. Giving up on this image'
                  ),
                  len(haplotype_sequences[alt]),
                  self.width,
              )
              return None
            alt_image = self.build_pileup(
                dv_call=dv_call,
                refbases=haplotype_sequences[alt],
                reads_for_samples=[
                    sample[alt] for sample in haplotype_alignments_for_samples
                ],
                alt_alleles=alt_alleles,
                sample_order=sample_order,
                custom_ref=True,
            )
            alt_images.append(alt_image)
        composite_image = _represent_alt_aligned_pileups(
            alt_aligned_representation, ref_image, alt_images
        )
        return composite_image
      else:
        return ref_image

    retval = []
    for alts in self._alt_allele_combinations(variant):
      pileup = _pileup_for_pair_of_alts(alts)
      # If the pileup is None, this can mean that we're near the edge of the
      # contig, so one pileup width is invalid.
      # Return None to indicate we couldn't process this variant.
      if pileup is None:
        return None
      retval.append((alts, pileup))
    return retval
