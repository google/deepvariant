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
from typing import Callable, Dict, Optional, Tuple



import ml_collections
import tensorflow as tf

from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant.protos import deepvariant_pb2
from google.protobuf import text_format
from third_party.nucleus.io import sharded_file_utils


# These are empirically determined to work well on TPU with our data sets,
# where lots of buffering and concurrency is necessary to keep the device
# busy.
# These are settable in the constructor.
_DEFAULT_INPUT_READ_THREADS = 32
_DEFAULT_SHUFFLE_BUFFER_ELEMENTS = 50_000
_DEFAULT_BUFFER_SIZE = 16 * 1000 * 1000

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
    mode: str,
    config: Optional[ml_collections.ConfigDict] = None,
    debugging_true_label_mode: bool = False,
) -> Callable[
    [tf.train.Example, Tuple[int, int, int]],
    Dict[str, tf.Tensor],
]:
  """Generates a function for parsing tf.train.Examples.

  If mode=debug, all fields are returned including the label.

  Args:
    mode: One of 'train', 'tune', 'predict' or 'debug'.
    config: A DeepVariant configuration.
    debugging_true_label_mode: If true, adds label info when mode=debug or
      mode=predict. An error is thrown if the label field is not available.

  Returns:
    Callable that returns dicts of DeepVariant tf.Examples.
  """
  proto_features = _PROTO_FEATURES.copy()
  if mode in ['debug', 'predict'] and not debugging_true_label_mode:
    del proto_features['label']
  if mode == 'debug':
    # Use all features in debug mode.
    pass
  elif mode == 'predict':
    features = [
        'image/encoded',
        'variant/encoded',
        'alt_allele_indices/encoded',
    ]
    if debugging_true_label_mode:
      features += ['label']
    proto_features = _select_features(
        proto_features,
        features,
    )
  elif mode in ['train', 'tune']:
    proto_features = _select_features(
        proto_features,
        ['image/encoded', 'label'],
    )
  else:
    raise ValueError(
        'Mode must be set to one of train, tune, predict, or debug'
    )

  use_denovo_label = config and config.denovo_enabled
  if use_denovo_label:
    proto_features['denovo_label'] = tf.io.FixedLenFeature((1), tf.int64)

  # class_weights is a comma-delimited string of weights for each class.
  # e.g. 1,1,10 or 1,1,1
  use_class_weights = (
      config and config.class_weights and proto_features.get('label')
  )
  if use_class_weights:
    class_weights = tf.constant(
        list(map(float, config.class_weights.split(','))), dtype=tf.float32
    )

  def parse_example(
      example: tf.train.Example,
      input_shape: Tuple[int, int, int],
  ) -> Dict[str, tf.Tensor]:
    """Parses a serialized tf.Example, preprocesses the image, and one-hot encodes the label."""

    result = tf.io.parse_single_example(
        serialized=example, features=proto_features
    )

    image = tf.io.decode_raw(result['image/encoded'], tf.uint8)
    image = tf.reshape(image, input_shape)
    image = tf.cast(image, tf.float32)

    # Preprocess image
    result['image'] = dv_utils.preprocess_images(image)

    if 'label' in result:
      result['label'] = tf.keras.layers.CategoryEncoding(
          num_tokens=dv_constants.NUM_CLASSES, output_mode='one_hot'
      )(result['label'])

    if use_denovo_label and result['denovo_label'][0] == 1:
      # If the example is denovo then set the denovo example weights for this.
      result['sample_weight'] = tf.constant(
          config.denovo_weight, dtype=tf.float32
      )
    elif use_class_weights:
      result['sample_weight'] = tf.gather(
          class_weights,
          tf.cast(tf.math.argmax(result['label'], axis=-1), dtype=tf.int32),
      )
    else:
      result['sample_weight'] = tf.constant(1.0, dtype=tf.float32)

    return result

  return parse_example


def input_fn(
    path: str,
    batch_size: int,
    mode: str = 'train',
    debugging_true_label_mode: bool = False,
    strategy: tf.distribute.Strategy = tf.distribute.get_strategy(),
    config: Optional[ml_collections.ConfigDict] = None,
    buffer_size: int = _DEFAULT_BUFFER_SIZE,
    input_reader_threads: int = _DEFAULT_INPUT_READ_THREADS,
    shuffle_buffer_elements: int = _DEFAULT_SHUFFLE_BUFFER_ELEMENTS,
    limit: int = -1,
) -> tf.data.Dataset:
  """tf.data.Dataset loading function.

  If mode == 'train' or 'tune', the output will be a tuple of
  (image, labels, weights).

  If mode == 'predict' or mode == 'debug', output will be a dictionary
  containing additional variant information.

  Args:
    path: the input filename for a tfrecord[.gz] file containing examples. Can
      contain sharding designators.
    batch_size: Examples per batch.
    mode: One of ['train', 'tune', 'predict', 'debug']. If in train/tune, a
      tuple is returned. If predict is used, fields required for inference are
      returned. If debug is used, all proto fields including label are returned.
    debugging_true_label_mode: If mode=predict or mode=debug, this can be used
      to include label information. An error will be thrown if no label field is
      present.
    strategy: A tf.distribute.Strategy.
    config: A configuration file.
    buffer_size: Number of bytes to prefetch when loading data.
    input_reader_threads: Number of threads to use when reading in examples.
    shuffle_buffer_elements: Number of shuffle buffer elements during training.
    limit: Limit the number of examples. If -1, then output is repeated. If 0,
      no limit is applied. If > 0, limit to that many variants.

  Returns:
    tf.data.Dataset
  """

  if mode not in ['train', 'tune', 'predict', 'debug']:
    raise ValueError(
        'Mode must be set to one of train, tune, predict, or debug'
    )
  if mode in ['train', 'tune'] and config is None:
    raise ValueError('A config must be provided when mode=train')

  is_training = mode == 'train'

  # Get input shape from input path.
  input_shape = dv_utils.get_shape_from_examples_path(path)

  def load_dataset(filename: str) -> tf.data.Dataset:
    return tf.data.TFRecordDataset(
        filename,
        buffer_size=buffer_size,
        compression_type='GZIP',
    )

  file_list = [
      tf.io.gfile.glob(sharded_file_utils.normalize_to_sharded_file_pattern(x))
      for x in path.split(',')
  ]
  file_list = list(itertools.chain(*file_list))
  ds = tf.data.Dataset.from_tensor_slices(file_list)

  if is_training:
    ds = ds.shuffle(ds.cardinality(), reshuffle_each_iteration=True)

  ds = ds.interleave(
      load_dataset,
      cycle_length=input_reader_threads,
      num_parallel_calls=tf.data.AUTOTUNE,
      deterministic=False,
  )

  if limit == -1:
    ds = ds.repeat()
  elif limit > 0:
    ds = ds.take(limit)

  if is_training:
    ds = ds.shuffle(shuffle_buffer_elements, reshuffle_each_iteration=True)

  # Retrieve preprocess function
  parse_example = create_parse_example_fn(
      mode,
      config,
      debugging_true_label_mode,
  )

  ds = ds.map(
      map_func=lambda example: parse_example(example, input_shape),
      num_parallel_calls=tf.data.AUTOTUNE,
      deterministic=False,
  )
  # Convert dict output to training tuple.
  if mode in ['train', 'tune']:
    ds = ds.map(
        map_func=lambda x: (
            x['image'],
            x['label'],
            x['sample_weight'],
        ),
        num_parallel_calls=tf.data.AUTOTUNE,
        deterministic=False,
    )
  ds = ds.batch(batch_size=batch_size, drop_remainder=is_training)

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
