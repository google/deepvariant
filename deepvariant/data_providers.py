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
"""Data providers for deepvariant images.

tf.data.Dataset and data providers for standard DeepVariant datasets for
training and evaluating germline calling accuracy.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



import tensorflow as tf

from google.protobuf import text_format
from third_party.nucleus.util import io_utils
from deepvariant import dv_constants
from deepvariant import tf_utils
from deepvariant.protos import deepvariant_pb2

slim = tf.contrib.slim

# These are emperically determined to work well on TPU with our data sets,
# where lots of buffering and concurrency is necessary to keep the device
# busy.
# These are settable in the constructor.
_DEFAULT_INPUT_READ_THREADS = 32
_DEFAULT_INPUT_MAP_THREADS = 48
_DEFAULT_SHUFFLE_BUFFER_ELEMENTS = 100
_DEFAULT_INITIAL_SHUFFLE_BUFFER_ELEMENTS = 1024
_DEFAULT_PREFETCH_BUFFER_BYTES = 16 * 1000 * 1000
# This one doesn't seem useful to adjust at runtime.
_PREFETCH_BATCHES = 4


class DeepVariantInput(object):
  """This class serves as an `input_fn` for the `tf.estimator` framework."""

  # Calling this object like a function returns a stream of variadic tuples.
  # Essentially it is a buffered io library, that handles concurrently
  # reading and possibly shuffling input records from a set of files. It
  # knows how to parse features we care about from tf.examples. It records
  # some extra information about the source of the input, such as the name
  # and number of classes.

  def __init__(
      self,
      mode,
      input_file_spec,
      num_examples=None,
      num_classes=dv_constants.NUM_CLASSES,
      tensor_shape=None,
      name=None,
      use_tpu=False,
      input_read_threads=_DEFAULT_INPUT_READ_THREADS,
      input_map_threads=_DEFAULT_INPUT_MAP_THREADS,
      shuffle_buffer_size=_DEFAULT_SHUFFLE_BUFFER_ELEMENTS,
      initial_shuffle_buffer_size=_DEFAULT_INITIAL_SHUFFLE_BUFFER_ELEMENTS,
      prefetch_dataset_buffer_size=_DEFAULT_PREFETCH_BUFFER_BYTES,
      sloppy=True,
      list_files_shuffle=True,
      debugging_true_label_mode=False):
    """Create an DeepVariantInput object, usable as an `input_fn`.

    Args:
      mode: the mode string (from `tf.estimator.ModeKeys`).
      input_file_spec: the input filename for a tfrecord[.gz] file containing
        examples.  Can contain sharding designators.
      num_examples: the number of examples contained in the input file.
        Required for setting learning rate schedule in train/eval only.
      num_classes: The number of classes in the labels of
        this dataset. Currently defaults to DEFAULT_NUM_CLASSES.
      tensor_shape: None (which means we get the shape from the first example in
        source), or list of int [height, width, channel] for testing.
      name: string, name of the dataset.
      use_tpu: use code paths tuned for TPU, in particular protobuf encoding.
        Default False.
      input_read_threads: number of threads for reading data.  Default 32.
      input_map_threads: number of threads for mapping data.  Default 48.
      shuffle_buffer_size: size of the final shuffle buffer, in elements.
        Default 100.
      initial_shuffle_buffer_size: int; the size of the dataset.shuffle buffer
        in elements.  Default is 1024.
      prefetch_dataset_buffer_size: int; the size of the TFRecordDataset buffer
        in bytes.  Default is 16 * 1000 * 1000.
      sloppy: boolean, allow parallel_interleave to be sloppy.  Default True.
      list_files_shuffle: boolean, allow list_files to shuffle.  Default True.
      debugging_true_label_mode: boolean. If true, the input examples are
                                 created with "training" mode. We'll parse the
                                 'label' field even if the `mode` is PREDICT.
    Raises:
      ValueError: if `num_examples` not provided, in a context requiring it.
    """
    self.mode = mode
    self.input_file_spec = input_file_spec
    self.name = name
    self.num_examples = num_examples
    self.num_classes = num_classes

    self.use_tpu = use_tpu
    self.sloppy = sloppy
    self.list_files_shuffle = list_files_shuffle
    self.input_read_threads = input_read_threads
    self.input_map_threads = input_map_threads
    self.shuffle_buffer_size = shuffle_buffer_size
    self.initial_shuffle_buffer_size = initial_shuffle_buffer_size
    self.prefetch_dataset_buffer_size = prefetch_dataset_buffer_size
    self.debugging_true_label_mode = debugging_true_label_mode
    self.feature_extraction_spec = self.features_extraction_spec_for_mode(
        mode in (tf.estimator.ModeKeys.TRAIN, tf.estimator.ModeKeys.EVAL) or
        debugging_true_label_mode)

    if num_examples is None and mode != tf.estimator.ModeKeys.PREDICT:
      raise ValueError('num_examples argument required for DeepVariantInput'
                       'in TRAIN/EVAL modes.')
    if tensor_shape:
      self.tensor_shape = tensor_shape
    else:
      self.tensor_shape = tf_utils.get_shape_from_examples_path(input_file_spec)

    self.input_files = tf.gfile.Glob(
        io_utils.NormalizeToShardedFilePattern(self.input_file_spec))

  def features_extraction_spec_for_mode(self, include_label_and_locus):
    """Returns a dict describing features from a TF.example."""
    spec = {
        'image/encoded': tf.FixedLenFeature((), tf.string),
        'variant/encoded': tf.FixedLenFeature((), tf.string),
        'alt_allele_indices/encoded': tf.FixedLenFeature((), tf.string),
    }
    if include_label_and_locus:
      # N.B. int32 fails here on TPU.
      spec['label'] = tf.FixedLenFeature((), tf.int64)
      spec['locus'] = tf.FixedLenFeature((), tf.string)
    return spec

  def parse_tfexample(self, tf_example):
    """Parse a DeepVariant pileup tf.Example to features and labels.

    This potentially stores parsed strings as fixed length tensors of integers,
    as required by TPU.  They have to be handled properly by consumers.

    Args:
      tf_example: a serialized tf.Example for a DeepVariant "pileup".
    Returns:
      If (mode is EVAL or TRAIN) or debugging_true_label_mode:
        (features, label) ...
      If mode is PREDICT,
        features ...
    """
    with tf.name_scope('input'):
      parsed = tf.parse_single_example(tf_example, self.feature_extraction_spec)
      image = parsed['image/encoded']
      if self.tensor_shape:
        # If the input is empty there won't be a tensor_shape.
        image = tf.reshape(tf.decode_raw(image, tf.uint8), self.tensor_shape)
        if self.use_tpu:
          # Cast to int32 for loading onto the TPU
          image = tf.cast(image, tf.int32)

      variant = parsed['variant/encoded']
      alt_allele_indices = parsed['alt_allele_indices/encoded']
      if self.use_tpu:
        # Passing a string to a TPU draws this error: TypeError: <dtype:
        # 'string'> is not a supported TPU infeed type. Supported types are:
        # [tf.float32, tf.int32, tf.complex64, tf.int64, tf.bool, tf.bfloat16]
        # Thus, we must encode the string as a tensor of int.
        variant = tf_utils.string_to_int_tensor(variant)
        alt_allele_indices = tf_utils.string_to_int_tensor(alt_allele_indices)

      features = {
          'image': image,
          'variant': variant,
          'alt_allele_indices': alt_allele_indices,
      }

      if (self.mode in (tf.estimator.ModeKeys.TRAIN, tf.estimator.ModeKeys.EVAL)
          or self.debugging_true_label_mode):
        if self.use_tpu:
          features['locus'] = tf_utils.string_to_int_tensor(parsed['locus'])
        else:
          features['locus'] = parsed['locus']

        if self.mode in (tf.estimator.ModeKeys.TRAIN,
                         tf.estimator.ModeKeys.EVAL):
          label = parsed['label']
          return features, label
        features['label'] = parsed['label']

      # For predict model, label is not present. So, returns features only.
      return features

  def __call__(self, params):
    """Interface to get a data batch, fulfilling `input_fn` contract.

    Args:
      params: a dict containing an integer value for key 'batch_size'.

    Returns:
      the tuple (features, labels), where:
        - features is a dict of Tensor-valued input features; keys populated
          are:
            'image'
            'variant'
            'alt_allele_indices'
          and, if not PREDICT mode, also:
            'locus'

          Aside from 'image', these may be encoded specially for TPU.

        - label is the Tensor-valued prediction label; in train/eval
          mode the label value is is populated from the data source; in
          inference mode, the value is a constant empty Tensor value "()".
    """
    # See https://cloud.google.com/tpu/docs/tutorials/inception-v3-advanced
    # for some background on tuning this on TPU.

    batch_size = params['batch_size']

    compression_type = tf_utils.compression_type_of_files(self.input_files)

    # NOTE: The order of the file names returned can be non-deterministic,
    # even if shuffle is false.  See b/73959787 and the note in cl/187434282.
    # We need the shuffle flag to be able to disable reordering in EVAL mode.
    dataset = tf.data.Dataset.list_files(
        io_utils.NormalizeToShardedFilePattern(self.input_file_spec),
        shuffle=self.mode == tf.estimator.ModeKeys.TRAIN,
    )

    # This shuffle applies to the set of files.
    if (self.mode == tf.estimator.ModeKeys.TRAIN and
        self.initial_shuffle_buffer_size > 0):
      dataset = dataset.shuffle(self.initial_shuffle_buffer_size)

    def load_dataset(filename):
      dataset = tf.data.TFRecordDataset(
          filename,
          buffer_size=self.prefetch_dataset_buffer_size,
          compression_type=compression_type)
      return dataset

    if self.mode == tf.estimator.ModeKeys.EVAL:
      # When EVAL, avoid parallel reads for the sake of reproducibility.
      dataset = dataset.interleave(
          load_dataset, cycle_length=self.input_read_threads, block_length=1)
    else:
      dataset = dataset.apply(
          # parallel_interleave requires tf 1.5 or later; this is
          # necessary for good performance.
          tf.contrib.data.parallel_interleave(
              load_dataset,
              cycle_length=self.input_read_threads,
              sloppy=self.sloppy))

    if self.mode == tf.estimator.ModeKeys.TRAIN:
      dataset = dataset.repeat()

    dataset = dataset.map(
        self.parse_tfexample, num_parallel_calls=self.input_map_threads)

    dataset = dataset.prefetch(_PREFETCH_BATCHES * batch_size)

    # This shuffle applies to the set of records.
    if self.mode == tf.estimator.ModeKeys.TRAIN:
      if self.shuffle_buffer_size > 0:
        dataset = dataset.shuffle(self.shuffle_buffer_size)

    if self.mode == tf.estimator.ModeKeys.PREDICT:
      dataset = dataset.batch(batch_size)
    else:
      # N.B.: we drop the final partial batch in eval mode, to work around TPU
      # static shape limitations in the current TPUEstimator.
      dataset = dataset.apply(
          tf.contrib.data.batch_and_drop_remainder(batch_size))

    return dataset

  def __str__(self):
    return (
        'DeepVariantInput(name={}, input_file_spec={}, num_examples={}, mode={}'
    ).format(self.name, self.input_file_spec, self.num_examples, self.mode)


# This is the entry point to get a DeepVariantInput when you start with
# a dataset configuration file name.
def get_input_fn_from_dataset(dataset_config_filename,
                              mode,
                              tensor_shape=None,
                              use_tpu=False):
  """Creates an input_fn from the dataset config file.

  Args:
    dataset_config_filename: str. Path to the dataset config pbtxt file.
    mode: one of tf.estimator.ModeKeys.{TRAIN,EVAL,PREDICT}
    tensor_shape: None, or list of int [height, width, channel] for testing.
    use_tpu: use the tpu code path in the input_fn.

  Returns:
    An input_fn from the specified split in the dataset_config file.

  Raises:
    ValueError: if the dataset config doesn't have the necessary information.
  """
  # Get the metadata.
  dataset_config = read_dataset_config(dataset_config_filename)
  # Return a reader for the data.
  return get_input_fn_from_filespec(
      input_file_spec=dataset_config.tfrecord_path,
      num_examples=dataset_config.num_examples,
      name=dataset_config.name,
      mode=mode,
      tensor_shape=tensor_shape,
      use_tpu=use_tpu)


# This is the entry point to get a DeepVariantInput when you start with
# a tf.example file specification, and associated metadata.
def get_input_fn_from_filespec(input_file_spec,
                               mode,
                               num_examples=None,
                               name=None,
                               tensor_shape=None,
                               use_tpu=False,
                               input_read_threads=_DEFAULT_INPUT_READ_THREADS,
                               debugging_true_label_mode=False):
  """Create a DeepVariantInput function object from a file spec.

  Args:
    input_file_spec: the tf.example input file specification, possibly sharded.
    mode: tf.estimator.ModeKeys.
    num_examples: the number of examples in the files (or None).
    name: the name for this data set (or None).
    tensor_shape: None, or list of int [height, width, channel] for testing.
    use_tpu: use the tpu code path in the input_fn.
    input_read_threads: number of threads reading the input files.
    debugging_true_label_mode: boolean. If true, the input examples are created
                               with "training" mode. We'll parse the 'label'
                               field even if the `mode` is PREDICT.

  Returns:
    A DeepVariantInput object usable as an input_fn.
  """

  return DeepVariantInput(
      mode=mode,
      input_file_spec=input_file_spec,
      num_examples=num_examples,
      tensor_shape=tensor_shape,
      name=name,
      use_tpu=use_tpu,
      input_read_threads=input_read_threads,
      debugging_true_label_mode=debugging_true_label_mode)


# Return the stream of batched images from a dataset.
def get_batches(tf_dataset, model, batch_size):
  """Provides batches of pileup images from this dataset.

  Creates a DeepVariantInput for tf_dataset. It instantiates an iterator
  on the dataset, and returns the images, labels, encoded_variant
  features in batches. This calls model.preprocess_images on the images
  (but note that we will be moving that step into model_fn for the
  Estimator api).

  Args:
    tf_dataset: a DeepVariantInput object
    model: a model object
    batch_size: int batch size

  Returns:
    (images, labels, encoded_variant)

  Raises:
    ValueError: if the dataset has the wrong mode.
  """
  if tf_dataset.mode not in (tf.estimator.ModeKeys.TRAIN,
                             tf.estimator.ModeKeys.EVAL):
    raise ValueError(
        'tf_dataset.mode is {} but must be one of TRAIN or EVAL.'.format(
            tf_dataset.mode))

  params = dict(batch_size=batch_size)
  features, labels = tf_dataset(params).make_one_shot_iterator().get_next()

  images = features['image']
  encoded_variant = features['variant']

  images = model.preprocess_images(images)
  return images, labels, encoded_variant


# Return the stream of batched images from a dataset.
def get_infer_batches(tf_dataset, model, batch_size):
  """Provides batches of pileup images from this dataset.

  This instantiates an iterator on the dataset, and returns the
  image, variant, alt_allele_indices, features in batches. It calls
  model.preprocess_images on the images (but note that we will be moving
  that step into model_fn for the Estimator api).

  Args:
    tf_dataset: DeepVariantInput.
    model: DeepVariantModel.
    batch_size: int.  The batch size.

  Returns:
    (image, variant, alt_allele_indices)

  Raises:
    ValueError: if the dataset has the wrong mode.
  """
  if tf_dataset.mode != tf.estimator.ModeKeys.PREDICT:
    raise ValueError('tf_dataset.mode is {} but must be PREDICT.'.format(
        tf_dataset.mode))

  params = dict(batch_size=batch_size)
  features = tf_dataset(params).make_one_shot_iterator().get_next()

  images = features['image']
  if tf_dataset.tensor_shape:
    # tensor_shape will be None if the input was an empty file.
    images = model.preprocess_images(images)
  variant = features['variant']
  alt_allele_indices = features['alt_allele_indices']

  return images, variant, alt_allele_indices


# This reads a pbtxt file and returns the config proto.
def read_dataset_config(dataset_config_filename):
  """Returns a DeepVariantDatasetConfig proto read from the dataset config file.

  Args:
    dataset_config_filename: String. Path to the dataset config pbtxt file.

  Returns:
    A DeepVariantDatasetConfig proto from the dataset_config file.

  Raises:
    ValueError: if the dataset config doesn't have the necessary information.
  """
  with tf.gfile.GFile(dataset_config_filename) as f:
    dataset_config = text_format.Parse(
        f.read(), deepvariant_pb2.DeepVariantDatasetConfig())

  if not dataset_config.name:
    raise ValueError('dataset_config needs to have a name')

  if not dataset_config.tfrecord_path:
    raise ValueError('The dataset in the config {} does not have a '
                     'tfrecord_path.'.format(dataset_config_filename))

  # redacted
  # of num_examples.
  if not dataset_config.num_examples:
    raise ValueError('The dataset in the config {} does not have a '
                     'num_examples.'.format(dataset_config_filename))

  return dataset_config


def write_dataset_config_to_pbtxt(dataset_config, dataset_config_filename):
  """Writes the dataset_config to a human-readable text format.

  Args:
    dataset_config: DeepVariantDatasetConfig. The config to be written out.
    dataset_config_filename: String. Path to the output pbtxt file.
  """
  with tf.gfile.GFile(dataset_config_filename, mode='w') as writer:
    writer.write(text_format.MessageToString(dataset_config))
