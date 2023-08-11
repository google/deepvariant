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
"""Utility functions for DeepVariant.

Started with a collection of utilities for working with the TF models. Now this
file includes broader utilities we use in DeepVariant.
"""


import enum
import numpy as np
import tensorflow as tf

from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from tensorflow.core.example import example_pb2
from typing import Optional

# Convert strings up to this length, then clip.  We picked a number that
# was less than 1K, with a bit of extra space for the length element,
# to give enough space without overflowing to a larger multiple of 128.
STRING_TO_INT_MAX_CONTENTS_LEN = 1020
# This is the length of the resulting buffer (including the length entry).
STRING_TO_INT_BUFFER_LENGTH = STRING_TO_INT_MAX_CONTENTS_LEN + 1


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


def example_variant_type(example):
  """Gets the locus field from example as a string."""
  return example.features.feature['variant_type'].int64_list.value[0]


def example_locus(example):
  """Gets the locus field from example as a string."""
  return example.features.feature['locus'].bytes_list.value[0]


def example_alt_alleles_indices(example):
  """Gets an iterable of the alt allele indices in example."""
  return deepvariant_pb2.CallVariantsOutput.AltAlleleIndices.FromString(
      example.features.feature['alt_allele_indices/encoded'].bytes_list.value[0]
  ).indices


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


def example_label(example: example_pb2.Example) -> Optional[int]:
  """Gets the label field from example as a string."""
  if 'label' not in example.features.feature:
    return None
  return int(example.features.feature['label'].int64_list.value[0])


def example_denovo_label(example: example_pb2.Example) -> Optional[int]:
  """Gets the label field from example as a string.

  Args:
    example: A tf.Example containing DeepVariant example.

  Returns:
    De novo label for the example.
  """
  if 'denovo_label' not in example.features.feature:
    return None
  return int(example.features.feature['denovo_label'].int64_list.value[0])


def example_image_shape(example):
  """Gets the image shape field from example as a list of int64."""
  if len(example.features.feature['image/shape'].int64_list.value) != 3:
    raise ValueError(
        'Invalid image/shape: we expect to find an image/shape '
        'field with length 3.'
    )
  return example.features.feature['image/shape'].int64_list.value[0:3]


def example_key(example):
  """Constructs a key for example based on its position and alleles."""
  variant = example_variant(example)
  alts = example_alt_alleles(example)
  return '{}:{}:{}->{}'.format(
      variant.reference_name,
      variant.start + 1,
      variant.reference_bases,
      '/'.join(alts),
  )


def example_set_label(example, numeric_label):
  """Sets the label features of example.

  Sets the label feature of example to numeric_label.

  Args:
    example: A tf.Example proto.
    numeric_label: A numeric (int64 compatible) label for example.
  """
  example.features.feature['label'].int64_list.value[:] = [numeric_label]


def example_set_denovo_label(
    example: example_pb2.Example, numeric_label: int
) -> None:
  """Sets the denovo label features of example.

  Sets the label feature of example to numeric_label.

  Args:
    example: a tf.Example proto.
    numeric_label: A numeric (int64 compatible) label for example.
  """
  example.features.feature['denovo_label'].int64_list.value[:] = [numeric_label]


def example_set_variant(example, variant):
  """Sets the variant/encoded feature of example to variant.SerializeToString().

  Args:
    example: a tf.Example proto.
    variant: third_party.nucleus.protos.Variant protobuf containing information
      about a candidate variant call.
  """
  example.features.feature['variant/encoded'].bytes_list.value[:] = [
      variant.SerializeToString()
  ]


def example_sequencing_type(example):
  return example.features.feature['sequencing_type'].int64_list.value[0]


def get_one_example_from_examples_path(source, proto=None):
  """Get the first record from `source`.

  Args:
    source: str. A pattern or a comma-separated list of patterns that represent
      file names.
    proto: A proto class. proto.FromString() will be called on each serialized
      record in path to parse it.

  Returns:
    The first record, or None.
  """
  files = sharded_file_utils.glob_list_sharded_file_patterns(source)
  if not files:
    raise ValueError(
        'Cannot find matching files with the pattern "{}"'.format(source)
    )
  for f in files:
    try:
      return next(tfrecord.read_tfrecords(f, proto=proto))
    except StopIteration:
      # Getting a StopIteration from one next() means source_path is empty.
      # Move on to the next one to try to get one example.
      pass
  return None


def get_shape_from_examples_path(source):
  """Reads one record from source to determine the tensor shape for all."""
  one_example = get_one_example_from_examples_path(source)
  if one_example:
    return example_image_shape(one_example)
  return None


def _simplify_variant(variant):
  """Returns a new Variant with only the basic fields of variant."""

  def _simplify_variant_call(call):
    """Returns a new VariantCall with the basic fields of call."""
    return variants_pb2.VariantCall(
        call_set_name=call.call_set_name,
        genotype=call.genotype,
        info=dict(call.info),
    )  # dict() is necessary to actually set info.

  return variants_pb2.Variant(
      reference_name=variant.reference_name,
      start=variant.start,
      end=variant.end,
      reference_bases=variant.reference_bases,
      alternate_bases=variant.alternate_bases,
      filter=variant.filter,
      quality=variant.quality,
      calls=[_simplify_variant_call(call) for call in variant.calls],
  )


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
  example_set_variant(example, variant)
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
  reader = tf.compat.v1.train.NewCheckpointReader(checkpoint_path)
  var_to_shape_map = reader.get_variable_to_shape_map()
  keys = variables_to_get if variables_to_get else var_to_shape_map.keys()
  return {key: tuple(var_to_shape_map[key]) for key in keys}


def model_num_classes(checkpoint_path, n_classes_model_variable):
  """Returns the number of classes in the checkpoint."""
  if not checkpoint_path:
    return None

  # Figure out how many classes this inception model was trained to predict.
  try:
    shapes = model_shapes(checkpoint_path, [n_classes_model_variable])
  except KeyError:
    return None
  if n_classes_model_variable not in shapes:
    return None
  return shapes[n_classes_model_variable][-1]


def string_to_int_tensor(x):
  """Graph operations decode a string into a fixed-size tensor of ints."""
  decoded = tf.compat.v1.decode_raw(x, tf.uint8)
  clipped = decoded[:STRING_TO_INT_MAX_CONTENTS_LEN]  # clip to allowed max_len
  shape = tf.shape(input=clipped)
  slen = shape[0]
  # pad to desired max_len
  padded = tf.pad(
      tensor=clipped, paddings=[[0, STRING_TO_INT_MAX_CONTENTS_LEN - slen]]
  )
  casted = tf.cast(padded, tf.int32)
  casted.set_shape([STRING_TO_INT_MAX_CONTENTS_LEN])
  return tf.concat([[slen], casted], 0)


def int_tensor_to_string(x):
  """Python operations to encode a tensor of ints into string of bytes."""
  slen = x[0]
  v = x[1 : slen + 1]
  return np.array(v, dtype=np.uint8).tostring()


def compression_type_of_files(files):
  """Return GZIP or None for the compression type of the files."""
  return 'GZIP' if all(f.endswith('.gz') for f in files) else None


def tpu_available(sess=None):
  """Return true if a TPU device is available to the default session."""
  if sess is None:
    init_op = tf.group(
        tf.compat.v1.global_variables_initializer(),
        tf.compat.v1.local_variables_initializer(),
    )
    with tf.compat.v1.Session() as sess:
      sess.run(init_op)
      devices = sess.list_devices()
  else:
    devices = sess.list_devices()
  return any(dev.device_type == 'TPU' for dev in devices)


def resolve_master(master, tpu_name, tpu_zone, gcp_project):
  """Resolve the master's URL given standard flags."""
  if master is not None:
    return master
  elif tpu_name is not None:
    return tf.distribute.cluster_resolver.TPUClusterResolver(
        tpu=[tpu_name], zone=tpu_zone, project=gcp_project
    ).get_master()
  else:
    # For k8s TPU we do not have/need tpu_name. See
    # https://cloud.google.com/tpu/docs/kubernetes-engine-setup#tensorflow-code
    return tf.distribute.cluster_resolver.TPUClusterResolver().get_master()


def get_example_info_json_filename(
    examples_filename: str, task_id: Optional[int]
) -> str:
  """Returns corresponding example_info.json filename for examples_filename."""
  if sharded_file_utils.is_sharded_file_spec(examples_filename):
    assert task_id is not None
    # If examples_filename has the @shards representation, resolve it into
    # the first shard. We only write .example_info.json to the first shard.
    example_info_prefix = sharded_file_utils.sharded_filename(
        examples_filename, task_id
    )
  else:
    # In all other cases, including non-sharded files,
    # or sharded filenames with -ddddd-of-ddddd, just append.
    example_info_prefix = examples_filename
  return example_info_prefix + '.example_info.json'


def get_tf_record_writer(output_filename: str) -> tf.io.TFRecordWriter:
  tf_options = None
  if output_filename.endswith('.gz'):
    tf_options = tf.io.TFRecordOptions(compression_type='GZIP')
  return tf.io.TFRecordWriter(output_filename, options=tf_options)


def preprocess_images(images):
  """Applies preprocessing operations for Inception images.

  Because this will run in model_fn, on the accelerator, we use operations
  that efficiently execute there.

  Args:
    images: A Tensor of shape [batch_size height, width, channel] with uint8
      values.

  Returns:
    A tensor of images of shape [batch_size height, width, channel]
    containing floating point values, with all points rescaled between
    -1 and 1 and possibly resized.
  """
  images = tf.cast(images, dtype=tf.float32)
  images = tf.subtract(images, 128.0)
  images = tf.math.divide(images, 128.0)
  return images
