# Copyright 2023 Google LLC.
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
"""Utility functions that uses dependencies with CLIF under the hood."""

import enum



from deepvariant import dv_utils
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from tensorflow.core.example import example_pb2


class EncodedVariantType(enum.Enum):
  """Enum capturing the int64 values we encode for different variant types.

  TPUs really like fixed length features, which makes it very difficult to use
  extract the type of a variant for an example using an encoded Variant
  protobufs or even a string value like "snp". The current best option appears
  to be to encode the type of a variant directly in an example as an int64. This
  enum provides a mapping between those raw int64 values in the example and a
  human-meaningful name for that type.
  """

  UNKNOWN = 0  # A variant of unknown type.
  SNP = 1  # The variant is a SNP.
  INDEL = 2  # The variant is an indel.


def encoded_variant_type(variant):
  """Gets the EncodedVariantType for variant.

  This function examines variant and returns the EncodedVariantType that best
  describes the variation type of variant. For example, if variant has
  `reference_bases = "A"` and `alternative_bases = ["C"]` this function would
  return EncodedVariantType.SNP.

  Args:
    variant: nucleus.Variant proto. The variant whose EncodedVariantType we want
      to get.

  Returns:
    EncodedVariantType enum value.
  """
  if variant_utils.is_snp(variant):
    return EncodedVariantType.SNP
  elif variant_utils.is_indel(variant):
    return EncodedVariantType.INDEL
  else:
    return EncodedVariantType.UNKNOWN


def make_example(
    variant,
    alt_alleles,
    encoded_image,
    shape,
    second_image=None,
    sequencing_type=0,
):
  """Creates a new tf.Example suitable for use with DeepVariant.

  Args:
    variant: third_party.nucleus.protos.Variant protobuf containing information
      about a candidate variant call.
    alt_alleles: A set of strings. Indicates the alternate alleles used as "alt"
      when constructing the image.
    encoded_image: a Tensor of type tf.string. Should contain an image encoding
      the reference and read data supporting variant. The encoding should be
      consistent with the image_format argument.
    shape: a list of (width, height, channel).
    second_image: a Tensor of type tf.string or None. Contains second image that
      encodes read data from another DNA sample. Must satisfy the same
      requirements as encoded_image.
    sequencing_type: int. The sequencing type of the input image.

  Returns:
    A tf.Example proto containing the standard DeepVariant features.
  """
  example = example_pb2.Example()
  features = example.features
  features.feature['locus'].bytes_list.value.append(
      ranges.to_literal(
          ranges.make_range(variant.reference_name, variant.start, variant.end)
      ).encode('latin-1')
  )
  dv_utils.example_set_variant(example, variant)
  variant_type = encoded_variant_type(variant).value
  features.feature['variant_type'].int64_list.value.append(variant_type)
  all_alts = list(variant.alternate_bases)
  alt_indices = sorted(all_alts.index(alt) for alt in alt_alleles)

  features.feature['alt_allele_indices/encoded'].bytes_list.value.append(
      deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
          indices=alt_indices
      ).SerializeToString()
  )

  features.feature['image/encoded'].bytes_list.value.append(encoded_image)
  features.feature['image/shape'].int64_list.value.extend(shape)
  if second_image is not None:
    features.feature['second_image/encoded'].bytes_list.value.append(
        second_image.encode('latin-1')
    )
    features.feature['second_image/shape'].int64_list.value.extend(shape)
  features.feature['sequencing_type'].int64_list.value.append(sequencing_type)
  return example
