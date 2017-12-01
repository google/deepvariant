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
"""Utility functions for working with TensorFlow in DeepVariant.

A collection of utilities for working with the TF models and evaluations we use
in DeepVariant.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



import tensorflow as tf

from tensorflow.core.example import example_pb2

from deepvariant.core import io_utils
from deepvariant.core import ranges
from deepvariant.core.genomics import variants_pb2
from deepvariant.protos import deepvariant_pb2


def example_locus(example):
  """Gets the locus field from example as a string."""
  return example.features.feature['locus'].bytes_list.value[0]


def example_alt_alleles_indices(example):
  """Gets an iterable of the alt allele indices in example."""
  return deepvariant_pb2.CallVariantsOutput.AltAlleleIndices.FromString(
      example.features.feature['alt_allele_indices/encoded']
      .bytes_list.value[0]).indices


def example_alt_alleles(example, variant=None):
  """Gets a list of the alt alleles in example."""
  variant = variant if variant else example_variant(example)
  return [
      variant.alternate_bases[i] for i in example_alt_alleles_indices(example)
  ]


def example_encoded_image(example):
  """Gets image field from example as a string."""
  return example.features.feature['image/encoded'].bytes_list.value[0]


def example_variant(example):
  """Gets and decodes the variant field from example as a Variant."""
  encoded = example.features.feature['variant/encoded'].bytes_list.value[0]
  return variants_pb2.Variant.FromString(encoded)


def example_label(example):
  """Gets the label field from example as a string."""
  return example.features.feature['label'].int64_list.value[0]


def example_image_format(example):
  """Gets the image format field from example as a string."""
  return example.features.feature['image/format'].bytes_list.value[0]


def example_image_shape(example):
  """Gets the image shape field from example as a list of int64."""
  if len(example.features.feature['image/shape'].int64_list.value) != 3:
    raise ValueError('Invalid image/shape: we expect to find an image/shape '
                     'field with length 3.')
  return example.features.feature['image/shape'].int64_list.value[0:3]


def example_truth_variant(example):
  """Gets and decodes the truth_variant field from example as a Variant."""
  return variants_pb2.Variant.FromString(
      example.features.feature['truth_variant/encoded'].bytes_list.value[0])


def example_key(example):
  """Constructs a key for example based on its position and alleles."""
  variant = example_variant(example)
  alts = example_alt_alleles(example)
  return '{}:{}:{}->{}'.format(variant.reference_name, variant.start + 1,
                               variant.reference_bases, '/'.join(alts))


def example_set_label(example, numeric_label):
  """Sets the label features of example.

  Sets the label feature of example to numeric_label.

  Args:
    example: a tf.Example proto.
    numeric_label: A numeric (int64 compatible) label for example.
  """
  example.features.feature['label'].int64_list.value[:] = [numeric_label]


def example_set_truth_variant(example, truth_variant, simplify=True):
  """Sets the truth_variant feature of example to truth_variant.

  The feature 'truth_variant' of example.features gets set to the serialized
  version of truth_variant. If simplify is True, we strip out the fields of
  truth_variant that aren't useful for training purposes (such as the INFO
  field map). truth_variant isn't modified directly, regardless.

  Args:
    example: a tf.Example proto.
    truth_variant: a
      learning.genomics.deepvariant.core.genomics.Variant proto.
    simplify: boolean. If True, we will simplify truth_variant before encoding
      it. If False, truth_variant will be written out as is.
  """
  if simplify:
    truth_variant = _simplify_variant(truth_variant)
  example.features.feature['truth_variant/encoded'].bytes_list.value[:] = [
      truth_variant.SerializeToString()
  ]


def _get_one_example_from_examples_path(source):
  """Reads one record from source."""
  # redacted
  # io_utils.read_tfrecord can read wildcard file patterns.
  # The source can be a comma-separated list.
  source_paths = source.split(',')
  for source_path in source_paths:
    files = tf.gfile.Glob(io_utils.NormalizeToShardedFilePattern(source_path))
    if not files:
      raise ValueError('Unable to read shape from source {}'.format(source))
    for f in files:
      try:
        return io_utils.read_tfrecords(f).next()
      except StopIteration:
        # Getting a StopIteration from one next() means source_path is empty.
        # Move on to the next one to try to get one example.
        pass
  return None


def get_shape_from_examples_path(source):
  """Reads one record from source to determine the tensor shape for all."""
  one_example = _get_one_example_from_examples_path(source)
  if one_example:
    return example_image_shape(one_example)
  return None


def get_format_from_examples_path(source):
  """Reads one record from source to determine the format for all."""
  one_example = _get_one_example_from_examples_path(source)
  if one_example:
    return example_image_format(one_example)
  return None


def _simplify_variant(variant):
  """Returns a new Variant with only the basic fields of variant."""

  def _simplify_variant_call(call):
    """Returns a new VariantCall with the basic fields of call."""
    return variants_pb2.VariantCall(
        call_set_name=call.call_set_name,
        genotype=call.genotype,
        info=dict(call.info))  # dict() is necessary to actually set info.

  return variants_pb2.Variant(
      reference_name=variant.reference_name,
      start=variant.start,
      end=variant.end,
      reference_bases=variant.reference_bases,
      alternate_bases=variant.alternate_bases,
      filter=variant.filter,
      quality=variant.quality,
      calls=[_simplify_variant_call(call) for call in variant.calls])


def make_example(variant, alt_alleles, encoded_image, shape, image_format):
  """Creates a new tf.Example suitable for use with DeepVariant.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant protobuf
      containing information about a candidate variant call.
    alt_alleles: A set of strings. Indicates the alternate alleles used as "alt"
      when constructing the image.
    encoded_image: a Tensor of type tf.string. Should contain an image encoding
      the reference and read data supporting variant. The encoding should be
      consistent with the image_format argument.
    shape: a list of (width, height, channel).
    image_format: string. The scheme used to encode our image.

  Returns:
    A tf.Example proto containing the standard DeepVariant features.
  """
  example = example_pb2.Example()
  features = example.features
  features.feature['locus'].bytes_list.value.append(
      ranges.to_literal(
          ranges.make_range(variant.reference_name, variant.start,
                            variant.end)))
  features.feature['variant/encoded'].bytes_list.value.append(
      variant.SerializeToString())
  all_alts = list(variant.alternate_bases)
  alt_indices = sorted(all_alts.index(alt) for alt in alt_alleles)

  features.feature['alt_allele_indices/encoded'].bytes_list.value.append(
      deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
          indices=alt_indices).SerializeToString())

  features.feature['image/encoded'].bytes_list.value.append(encoded_image)
  features.feature['image/format'].bytes_list.value.append(image_format)
  features.feature['image/shape'].int64_list.value.extend(shape)
  return example


def model_shapes(checkpoint_path, variables_to_get=None):
  """Returns the shape of each tensor in the model at checkpoint_path.

  Args:
    checkpoint_path: string. The path to a tensorflow checkpoint containing a
      model whose tensor shapes we want to get.
    variables_to_get: options, list of strings. If provided, only returns the
      shapes of tensors in variables whose name is present in this list. If
      None, the default, gets all of the tensors. A KeyError will be raised if
      any variable name in variables_to_get isn't present in the checkpointed
      model.

  Returns:
    A dictionary mapping variable names [string] to tensor shapes [tuple].
  """
  reader = tf.train.NewCheckpointReader(checkpoint_path)
  var_to_shape_map = reader.get_variable_to_shape_map()
  keys = variables_to_get if variables_to_get else var_to_shape_map.keys()
  return {key: tuple(var_to_shape_map[key]) for key in keys}
