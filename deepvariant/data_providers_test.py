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
"""Tests for learning.genomics.deepvariant.data_provider."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import functools



from absl.testing import absltest
from absl.testing import parameterized
import mock
import tensorflow as tf

from tensorflow.core.example import example_pb2
from deepvariant import data_providers
from deepvariant import modeling
from deepvariant import test_utils
from deepvariant.core import io_utils
from deepvariant.core import variantutils
from deepvariant.protos import deepvariant_pb2

slim = tf.contrib.slim


def setUpModule():
  test_utils.init()


def make_golden_dataset(compressed_inputs=False):
  if compressed_inputs:
    source_path = test_utils.test_tmpfile('make_golden_dataset.tfrecord.gz')
    io_utils.write_tfrecords(
        io_utils.read_tfrecords(test_utils.GOLDEN_TRAINING_EXAMPLES),
        source_path)
  else:
    source_path = test_utils.GOLDEN_TRAINING_EXAMPLES
  return data_providers.DeepVariantDataSet(
      name='labeled_golden', source=source_path, num_examples=49)


def _test_dataset_config(filename, **kwargs):
  """Creates a DeepVariantDatasetConfig(**kwargs) and writes it to filename."""
  dataset_config_pbtext_filename = test_utils.test_tmpfile(filename)
  dataset_config = deepvariant_pb2.DeepVariantDatasetConfig(**kwargs)
  data_providers.write_dataset_config_to_pbtxt(dataset_config,
                                               dataset_config_pbtext_filename)
  return dataset_config_pbtext_filename


class DataProviderTest(parameterized.TestCase):

  def test_get_dataset(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'golden.dataset_config.pbtxt',
        name='some_dataset_name',
        tfrecord_path='/path/to/dataset',
        num_examples=1000)
    ds = data_providers.get_dataset(
        dataset_config_pbtext_filename, tensor_shape=[3, 4, 7])

    self.assertEqual('some_dataset_name', ds.name)
    self.assertEqual('/path/to/dataset', ds.source)
    self.assertEqual(1000, ds.num_examples)
    self.assertEqual([3, 4, 7], ds.tensor_shape)

  def test_get_dataset_raises_error_for_empty_name(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_get_dataset_raises_error_for_empty_name.pbtxt')
    with self.assertRaisesRegexp(ValueError,
                                 'dataset_config needs to have a name'):
      data_providers.get_dataset(dataset_config_pbtext_filename)

  def test_get_dataset_raises_error_for_empty_data_split(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_get_dataset_raises_error_for_empty_data_split.pbtxt',
        name='some_dataset_name')
    expected_exception_message = ('The dataset in the config {} does not '
                                  'have a tfrecord_path.'
                                  .format(dataset_config_pbtext_filename))
    with self.assertRaisesRegexp(ValueError, expected_exception_message):
      data_providers.get_dataset(dataset_config_pbtext_filename)

  def test_get_dataset_raises_error_for_empty_num_examples(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_get_dataset_raises_error_for_empty_num_examples.pbtxt',
        name='some_dataset_name',
        tfrecord_path='/path/to/dataset')
    expected_exception_message = ('The dataset in the config {} does not have '
                                  'a num_examples.'
                                  .format(dataset_config_pbtext_filename))
    with self.assertRaisesRegexp(ValueError, expected_exception_message):
      data_providers.get_dataset(dataset_config_pbtext_filename)

  def test_dataset_definition(self):
    ds = data_providers.DeepVariantDataSet(
        name='name',
        source='test.tfrecord',
        num_examples=10,
        num_classes=2,
        tensor_shape=[11, 13, 7])
    self.assertEqual('name', ds.name)
    self.assertEqual('test.tfrecord', ds.source)
    self.assertEqual(10, ds.num_examples)
    self.assertEqual(2, ds.num_classes)
    self.assertEqual([11, 13, 7], ds.tensor_shape)

  def test_good_dataset(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_good_dataset.pbtxt',
        name='some_dataset_name',
        tfrecord_path='/path/to/dataset',
        num_examples=1000)
    ds = data_providers.get_dataset(
        dataset_config_pbtext_filename, tensor_shape=[100, 221, 7])
    # Test that the slim.DataSet we create from the dataset has the values
    # and fields we expect.
    with tf.Session():
      slim_ds = ds.get_slim_dataset()
      self.assertEqual(ds.num_examples, slim_ds.num_samples)
      self.assertItemsEqual(
          ['image', 'label', 'locus', 'variant', 'truth_variant'],
          slim_ds.decoder.list_items())
      self.assertEqual([100, 221, 7], ds.tensor_shape)

  def assertDataSetExamplesMatchExpected(self, dataset, expected_dataset):
    with tf.Session() as sess:
      provider = slim.dataset_data_provider.DatasetDataProvider(
          expected_dataset.get_slim_dataset(),
          shuffle=False,
          reader_kwargs={
              'options': io_utils.make_tfrecord_options(expected_dataset.source)
          })
      sess.run(tf.global_variables_initializer())
      coord = tf.train.Coordinator()
      threads = tf.train.start_queue_runners(coord=coord, sess=sess)
      image, label, locus = provider.get(['image', 'label', 'locus'])
      seen = [
          sess.run([image, label, locus])[2]
          for _ in range(expected_dataset.num_examples)
      ]
      coord.request_stop()
      coord.join(threads)

    expected_loci = [
        example.features.feature['locus'].bytes_list.value[0]
        for example in io_utils.read_tfrecords(expected_dataset.source)
    ]
    self.assertEqual(len(expected_loci), expected_dataset.num_examples)
    self.assertEqual(expected_loci, seen)
    # Note that this expected shape comes from the golden dataset. If the data
    # is remade in the future, the values might need to be modified accordingly.
    self.assertEqual([100, 221, 7], expected_dataset.tensor_shape)

  @parameterized.parameters(True, False)
  def test_reading_dataset(self, compressed_inputs):
    golden_dataset = make_golden_dataset(compressed_inputs)
    self.assertDataSetExamplesMatchExpected(golden_dataset.get_slim_dataset(),
                                            golden_dataset)

  @parameterized.parameters(True, False)
  def test_reading_sharded_dataset(self, compressed_inputs):
    golden_dataset = make_golden_dataset(compressed_inputs)
    n_shards = 3
    sharded_path = test_utils.test_tmpfile('sharded@{}'.format(n_shards))
    io_utils.write_tfrecords(
        io_utils.read_tfrecords(golden_dataset.source), sharded_path)

    config_file = _test_dataset_config(
        'test_sharded.pbtxt',
        name='sharded_test',
        tfrecord_path=sharded_path,
        num_examples=golden_dataset.num_examples)

    self.assertDataSetExamplesMatchExpected(
        data_providers.get_dataset(config_file).get_slim_dataset(),
        golden_dataset)

  @parameterized.parameters(True, False)
  def test_get_training_batches(self, compressed_inputs):
    golden_dataset = make_golden_dataset(compressed_inputs)
    batch_size = 16
    with tf.Session() as sess:
      mock_model = mock.MagicMock(autospec=modeling.DeepVariantModel)
      mock_model.preprocess_image.side_effect = functools.partial(
          tf.image.resize_image_with_crop_or_pad,
          target_height=107,
          target_width=221)
      batch = data_providers.make_training_batches(
          golden_dataset.get_slim_dataset(), mock_model, batch_size)

      # We should have called our preprocess_image exactly once. We don't have
      # the actual objects to test for the call, though.
      test_utils.assert_called_once_workaround(mock_model.preprocess_image)

      # Get our images, labels, and variants for further testing.
      sess.run(tf.global_variables_initializer())
      coord = tf.train.Coordinator()
      threads = tf.train.start_queue_runners(coord=coord, sess=sess)
      images, labels, variants = sess.run(batch)

      # Checks that our labels are the right shape and are one-hot encoded.
      self.assertEqual((batch_size, 107, 221, 7), images.shape)
      self.assertEqual((batch_size,), labels.shape)
      for label in labels:
        self.assertTrue(0 <= label <= 2)

      # Check that our variants has the shape we expect and actually contain
      # variants by decoding them and checking the reference_name.
      self.assertEqual((batch_size,), variants.shape)
      for variant in variantutils.decode_variants(variants):
        self.assertEqual(variant.reference_name, 'chr20')

      # Shutdown tensorflow
      coord.request_stop()
      coord.join(threads)

  @parameterized.parameters(
      ('test_shape.gz', 'test_shape.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape@1.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-?????-of-00001.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-*.gz'), ('output', 'output'),
      ('test_shape-00000-of-00001', 'test_shape@1'),
      ('test_shape-00000-of-00001', 'test_shape-?????-of-00001'),
      ('test_shape-00000-of-00001', 'test_shape-*'))
  def test_get_shape_from_examples_path(self, file_name_to_write,
                                        tfrecord_path_to_match):
    example = example_pb2.Example()
    valid_shape = [1, 2, 3]
    example.features.feature['image/shape'].int64_list.value.extend(valid_shape)
    output_file = test_utils.test_tmpfile(file_name_to_write)
    io_utils.write_tfrecords([example], output_file)
    ds = data_providers.DeepVariantDataSet(
        name='test_shape',
        source=test_utils.test_tmpfile(tfrecord_path_to_match),
        num_examples=1)
    self.assertEqual(valid_shape, ds.tensor_shape)

  def test_get_shape_from_examples_path_invalid_path(self):
    with self.assertRaisesRegexp(Exception, '/this/path/does/not'):
      data_providers.DeepVariantDataSet(
          name='test_invalid_path',
          source='/this/path/does/not/exist',
          num_examples=1)


if __name__ == '__main__':
  absltest.main()
