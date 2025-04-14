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
"""Data providers for deepvariant images.

tf.data.Dataset and data providers for standard DeepVariant datasets for
training and evaluating germline calling accuracy.
"""

import itertools
import os
from typing import Callable, Dict, Tuple, Union

from absl import logging
import ml_collections
import tensorflow as tf

from deepvariant import dv_utils
from deepvariant.protos import deepvariant_pb2
from google.protobuf import text_format
from third_party.nucleus.io import sharded_file_utils


_PROTO_FEATURES = {
    'locus': tf.io.FixedLenFeature((), tf.string),
    'variant/encoded': tf.io.FixedLenFeature((), tf.string),
    'variant_type': tf.io.FixedLenFeature((), tf.int64),
    'alt_allele_indices/encoded': tf.io.FixedLenFeature((), tf.string),
    'image/encoded': tf.io.FixedLenFeature((), tf.string),
    'sequencing_type': tf.io.FixedLenFeature((), tf.int64, default_value=0),
    'label': tf.io.FixedLenFeature((1), tf.int64),
}


def _select_features(proto_features, features):
  return {x: y for x, y in proto_features.items() if x in features}


def create_parse_example_fn(
    config: ml_collections.ConfigDict,
    mode: str,
    input_shape: Tuple[int, int, int],
) -> Callable[
    [tf.train.Example],
    Dict[str, tf.Tensor],
]:
  """Generates a function for parsing tf.train.Examples.

  If mode=debug, all fields are returned including the label.

  Args:
    config: A DeepVariant configuration.
    mode: One of 'train', 'tune', 'predict' or 'debug'.
    input_shape: The shape of the input image.

  Returns:
    Callable that returns dicts of DeepVariant tf.Examples.
  """
  if mode == 'debug':
    # Load all features in debug mode.
    proto_features = _PROTO_FEATURES
  elif mode == 'predict':
    proto_features = _select_features(
        _PROTO_FEATURES,
        ['image/encoded', 'variant/encoded', 'alt_allele_indices/encoded'],
    )
  elif mode in ['train', 'tune']:
    proto_features = _select_features(
        _PROTO_FEATURES,
        ['image/encoded', 'label', 'variant_type'],
    )
  else:
    raise ValueError(
        'Mode must be set to one of train, tune, predict, or debug'
    )

  if config.denovo_enabled:
    proto_features['denovo_label'] = tf.io.FixedLenFeature((1), tf.int64)

  # class_weights is a comma-delimited string of weights for each class.
  # e.g. 1,1,10 or 1,1,1
  if config.class_weights:
    class_weights = tf.constant(
        list(map(float, config.class_weights.split(','))), dtype=tf.float16
    )

  def parse_example(
      example: tf.train.Example,
  ) -> Union[
      Dict[str, tf.Tensor], Tuple[tf.Tensor, tf.Tensor, tf.Tensor, tf.Tensor]
  ]:
    """Parses a serialized tf.Example, preprocesses the image, and one-hot encodes the label."""

    result = tf.io.parse_single_example(
        serialized=example, features=proto_features
    )

    image = tf.io.decode_raw(result['image/encoded'], tf.uint8)
    result['image'] = tf.reshape(image, input_shape)

    # Image preprocessing and one-hot labeling takes place on accelerators to
    # reduce memory usage and speed up training. Do not move those ops here.

    if config.denovo_enabled and result['denovo_label'][0] == 1:
      # If the example is denovo then set the denovo example weights for this.
      result['sample_weight'] = tf.constant(
          config.denovo_weight, dtype=tf.float16
      )
    elif 'label' in result and config.class_weights:
      result['sample_weight'] = tf.squeeze(
          tf.gather(class_weights, result['label'])
      )
    else:
      result['sample_weight'] = tf.constant(1, dtype=tf.int8)

    if 'label' in result:
      result['label'] = tf.cast(result['label'], dtype=tf.int8)
    if mode in ['train', 'tune']:
      result['variant_type'] = tf.cast(result['variant_type'], dtype=tf.int8)
      return (
          result['image'],
          result['label'],
          result['sample_weight'],
          result['variant_type'],
      )
    else:
      del result['image/encoded']
      return result

  return parse_example


def input_fn(
    path: str,
    config: ml_collections.ConfigDict,
    mode: str = 'train',
    strategy: tf.distribute.Strategy = tf.distribute.get_strategy(),
) -> tf.data.Dataset:
  """tf.data.Dataset loading function.

  If mode == 'train' or 'tune', the output will be a tuple of
  (image, labels, weights).

  If mode == 'predict' or mode == 'debug', output will be a dictionary
  containing additional variant information.

  Args:
    path: the input filename for a tfrecord[.gz] file containing examples. Can
      contain sharding designators.
    config: A configuration file.
    mode: One of ['train', 'tune', 'predict', 'debug']. If in train/tune, a
      tuple is returned. If predict is used, fields required for inference are
      returned. If debug is used, all proto fields including label are returned.
    strategy: A tf.distribute.Strategy.

  Returns:
    tf.data.Dataset
  """

  if mode not in ['train', 'tune', 'predict', 'debug']:
    raise ValueError(
        'Mode must be set to one of train, tune, predict, or debug'
    )
  is_training = mode in ['train', 'tune']

  # Get input shape from input path.
  # First, try finding example_info.json in the same directory as the input
  # path.
  config_dir = os.path.dirname(path)
  example_info_json = os.path.join(config_dir, 'example_info.json')
  if tf.io.gfile.exists(example_info_json):
    logging.info('Reading example_info.json: %s', example_info_json)
    input_shape, _, _ = dv_utils.get_shape_and_channels_from_json(
        example_info_json
    )
  else:
    logging.info(
        'example_info.json does not exist in the directory of %s. Reading the'
        ' shape from the examples',
        path,
    )
    input_shape = dv_utils.get_shape_from_examples_path(path)

  file_list = [
      tf.io.gfile.glob(sharded_file_utils.normalize_to_sharded_file_pattern(x))
      for x in path.split(',')
  ]
  file_list = list(itertools.chain(*file_list))

  if is_training:
    file_list = tf.random.shuffle(file_list)

  ds = tf.data.Dataset.from_tensor_slices(file_list)

  def load_dataset(filename: str) -> tf.data.Dataset:
    return tf.data.TFRecordDataset(
        filename,
        buffer_size=config.prefetch_buffer_bytes,
        compression_type='GZIP',
    )

  ds = ds.interleave(
      load_dataset,
      cycle_length=config.input_read_threads,
      num_parallel_calls=tf.data.AUTOTUNE,
      deterministic=False,
  )

  ds = ds.repeat()

  # Retrieve preprocess function
  parse_example = create_parse_example_fn(
      config,
      mode,
      input_shape,
  )

  ds = ds.map(
      map_func=parse_example,
      num_parallel_calls=tf.data.AUTOTUNE,
      deterministic=False,
  )

  if is_training and config.shuffle_buffer_elements:
    ds = ds.shuffle(
        config.shuffle_buffer_elements, reshuffle_each_iteration=True
    )

  ds = ds.batch(
      batch_size=config.batch_size,
      deterministic=False,
      drop_remainder=True,
  )

  # Prefetch overlaps in-feed with training
  ds = ds.prefetch(tf.data.AUTOTUNE)

  # Distribute the dataset
  ds = strategy.experimental_distribute_dataset(ds)

  return ds


def read_dataset_config(dataset_config_filename):
  """Returns a DeepVariantDatasetConfig proto read from the dataset config file.

  Args:
    dataset_config_filename: String. Path to the dataset config pbtxt file.

  Returns:
    A DeepVariantDatasetConfig proto from the dataset_config file.

  Raises:
    ValueError: if the dataset config doesn't have the necessary information.
  """
  with tf.io.gfile.GFile(dataset_config_filename) as f:
    dataset_config = text_format.Parse(
        f.read(), deepvariant_pb2.DeepVariantDatasetConfig()
    )

  if not dataset_config.name:
    raise ValueError('dataset_config needs to have a name')

  if not dataset_config.tfrecord_path:
    raise ValueError(
        'The dataset in the config {} does not have a tfrecord_path.'.format(
            dataset_config_filename
        )
    )

  # TODO: remove this check once we're able to deal with absence
  # of num_examples.
  if not dataset_config.num_examples:
    raise ValueError(
        'The dataset in the config {} does not have a num_examples.'.format(
            dataset_config_filename
        )
    )

  return dataset_config


def write_dataset_config_to_pbtxt(dataset_config, dataset_config_filename):
  """Writes the dataset_config to a human-readable text format.

  Args:
    dataset_config: DeepVariantDatasetConfig. The config to be written out.
    dataset_config_filename: String. Path to the output pbtxt file.
  """
  with tf.io.gfile.GFile(dataset_config_filename, mode='w') as writer:
    writer.write(text_format.MessageToString(dataset_config))
