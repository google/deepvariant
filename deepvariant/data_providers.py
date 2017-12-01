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

tf.slim datasets and data providers for standard DeepVariant datasets for
training and evaluating germline calling accuracy.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



import tensorflow as tf

from google.protobuf import text_format
from deepvariant import tf_utils
from deepvariant.core import io_utils
from deepvariant.protos import deepvariant_pb2

slim = tf.contrib.slim

# Number of classes represented in the data set. The three classes are
# homozygous reference (0), heterozygous (1) and homozygous alternative (2).
DEFAULT_NUM_CLASSES = 3


def make_training_batches(dataset, model, batch_size):
  """Provides batches of pileup images from this dataset.

  Creates a DataSetProvider for dataset, extracts image, label, and
  truth_variant from it, preprocesses each image with model.preprocess_image()
  and finally batches these up.

  Args:
    dataset: a slim DataSet we want to turn into batches. Must provide data
      items "image", "label", and "truth_variant".
    model: a DeepVariantModel to use for preprocessing each image before
      batching.
    batch_size: the number of images in each batch.

  Returns:
    images: 4-D float Tensor of a batch of images with shape
      (batch_size, height, width, 3).
    labels: a 1-D integer Tensor shape (batch_size,) containing the labels for
      each image, in the same order.
    encoded_truth_variants: Tensor of strings with shape (batch_size,).
      Each element of this tensor is a byte-encoded learning.genomics.v1.Variant
      protobuf in the same order as images and one_hot_labels.
  """
  data_provider = slim.dataset_data_provider.DatasetDataProvider(
      dataset,
      common_queue_capacity=2 * batch_size,
      common_queue_min=batch_size,
      reader_kwargs={
          'options': io_utils.make_tfrecord_options(dataset.data_sources)
      })
  # Load the data.
  image, label, truth_variant = data_provider.get(
      ['image', 'label', 'truth_variant'])
  image = model.preprocess_image(image)
  return tf.train.shuffle_batch(
      [image, label, truth_variant],
      batch_size=batch_size,
      num_threads=4,
      capacity=5000,
      # redacted
      min_after_dequeue=min(1000, dataset.num_samples))


# redacted
class DeepVariantDataSet(object):
  """All of the information needed to create and use a DeepVariant dataset."""

  def __init__(self,
               name,
               source,
               num_examples,
               num_classes=DEFAULT_NUM_CLASSES,
               tensor_shape=None):
    """Creates a dataset.

    Args:
      name: str. The name of this dataset. Used to refer to this dataset on
        the command line.
      source: str or list[str]. A file path pattern or a comma-separated list of
        file path patterns pointing to TF.Example PIC images containing the data
        for this dataset.
      num_examples: A positive integer. The number of examples in this dataset.
      num_classes: A positive integer. The number of classes in the labels of
        this dataset. Currently defaults to DEFAULT_NUM_CLASSES.
      tensor_shape: None (whihc means we get the shape from the first example in
        source), or list of int [height, width, channel] for testing.
    """
    self.name = name
    self.source = source
    self.num_examples = num_examples
    self.num_classes = num_classes
    if tensor_shape:
      self.tensor_shape = tensor_shape
    else:
      self.tensor_shape = tf_utils.get_shape_from_examples_path(source)

  def __str__(self):
    return ('DeepVariantDataSet(name={}, source={}, num_examples={}, '
            'num_classes={}').format(self.name, self.source, self.num_examples,
                                     self.num_classes)

  def get_slim_dataset(self):
    """Returns a Slim dataset for this dataset.

    Returns:
      A tf.slim.dataset.Dataset with the data in this DataSet.
    """
    keys_to_features = {
        'image/encoded': tf.FixedLenFeature((), tf.string),
        'variant/encoded': tf.FixedLenFeature((), tf.string),
        'truth_variant/encoded': tf.FixedLenFeature((), tf.string),
        'image/format': tf.FixedLenFeature((), tf.string),
        'label': tf.FixedLenFeature((1,), tf.int64),
        'locus': tf.FixedLenFeature((), tf.string),
    }
    items_to_handlers = {
        'image':
            slim.tfexample_decoder.Image(
                'image/encoded', 'image/format', shape=self.tensor_shape),
        'label':
            slim.tfexample_decoder.Tensor('label', shape=[]),
        'locus':
            slim.tfexample_decoder.Tensor('locus', shape=[]),
        'variant':
            slim.tfexample_decoder.Tensor('variant/encoded', shape=[]),
        'truth_variant':
            slim.tfexample_decoder.Tensor('truth_variant/encoded', shape=[]),
    }

    # redacted
    # shuffled correctly in training.
    return slim.dataset.Dataset(
        data_sources=self.source.split(','),
        reader=tf.TFRecordReader,
        decoder=slim.tfexample_decoder.TFExampleDecoder(keys_to_features,
                                                        items_to_handlers),
        num_samples=self.num_examples,
        items_to_descriptions=None)


def _get_dataset(name, path, num_examples, tensor_shape=None):
  """Creates a dataset with a specified name a path to the source file.

  Args:
    name: String. The name of the dataset.
    path: String. The path to the source file of the dataset.
    num_examples: int. The number of examples in the dataset.
    tensor_shape: None, or list of int [height, width, channel] for testing.

  Returns:
    A DeepVariantDataSet where name is the specified name, and the source is the
    path.

  Raises:
    ValueError: if name and path are not specified.
  """
  if not name:
    raise ValueError('Name must not be None', name)
  if not path:
    raise ValueError('Path must not be None', path)

  return DeepVariantDataSet(
      name=name,
      source=path,
      num_examples=num_examples,
      tensor_shape=tensor_shape)


def get_dataset(dataset_config_filename, tensor_shape=None):
  """Creates a DeepVariantDataSet from the dataset config file.

  Args:
    dataset_config_filename: String. Path to the dataset config pbtxt file.
    tensor_shape: None, or list of int [height, width, channel] for testing.

  Returns:
    A DeepVariantDataSet from the specified split in the dataset_config file.

  Raises:
    ValueError: if the dataset config doesn't have the necessary information.
  """

  def read_dataset_config(dataset_config_pbtxt_filename):
    with tf.gfile.GFile(dataset_config_pbtxt_filename) as f:
      return text_format.Parse(f.read(),
                               deepvariant_pb2.DeepVariantDatasetConfig())

  dataset_config = read_dataset_config(dataset_config_filename)

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

  return _get_dataset(
      dataset_config.name,
      dataset_config.tfrecord_path,
      dataset_config.num_examples,
      tensor_shape=tensor_shape)


def write_dataset_config_to_pbtxt(dataset_config, dataset_config_filename):
  """Writes the dataset_config to a human-readable text format.

  Args:
    dataset_config: DeepVariantDatasetConfig. The config to be written out.
    dataset_config_filename: String. Path to the output pbtxt file.
  """
  with tf.gfile.GFile(dataset_config_filename, mode='w') as writer:
    writer.write(text_format.MessageToString(dataset_config))
