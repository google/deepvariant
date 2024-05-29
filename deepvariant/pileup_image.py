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

from typing import List



from deepvariant import dv_constants
from deepvariant import sample as sample_lib
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import pileup_image_native
from third_party.nucleus.io import fasta
from third_party.nucleus.protos import reads_pb2


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

  def get_channels(self):
    return self._channels_enum
