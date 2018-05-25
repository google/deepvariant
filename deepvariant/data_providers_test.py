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

import math



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import six
import tensorflow as tf

from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import io_utils
from third_party.nucleus.util import variant_utils
from tensorflow.core.example import example_pb2
from deepvariant import data_providers
from deepvariant import dv_constants
from deepvariant import testdata
from deepvariant import tf_utils
from deepvariant.protos import deepvariant_pb2


def setUpModule():
  testdata.init()


# Return a DeepVariantInput attached to the golden training data.
# Run with shuffling off, and in eval mode.
def make_golden_dataset(compressed_inputs=False,
                        mode=tf.estimator.ModeKeys.EVAL,
                        use_tpu=False):
  if compressed_inputs:
    source_path = test_utils.test_tmpfile('make_golden_dataset.tfrecord.gz')
    io_utils.write_tfrecords(
        io_utils.read_tfrecords(testdata.GOLDEN_TRAINING_EXAMPLES), source_path)
  else:
    source_path = testdata.GOLDEN_TRAINING_EXAMPLES
  return data_providers.get_input_fn_from_filespec(
      input_file_spec=source_path,
      num_examples=testdata.N_GOLDEN_TRAINING_EXAMPLES,
      name='labeled_golden',
      mode=mode,
      tensor_shape=None,
      use_tpu=use_tpu)


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
        tfrecord_path='/dev/null',
        num_examples=1000)
    ds = data_providers.get_input_fn_from_dataset(
        dataset_config_pbtext_filename,
        mode=tf.estimator.ModeKeys.EVAL,
        tensor_shape=[3, 4, dv_constants.PILEUP_NUM_CHANNELS])

    self.assertEqual('some_dataset_name', ds.name)
    self.assertEqual('/dev/null', ds.input_file_spec)
    self.assertEqual(1000, ds.num_examples)
    self.assertEqual([3, 4, dv_constants.PILEUP_NUM_CHANNELS], ds.tensor_shape)

  def test_get_dataset_raises_error_for_empty_name(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_get_dataset_raises_error_for_empty_name.pbtxt')
    with self.assertRaisesRegexp(ValueError,
                                 'dataset_config needs to have a name'):
      data_providers.get_input_fn_from_dataset(
          dataset_config_pbtext_filename, mode=tf.estimator.ModeKeys.EVAL)

  def test_get_dataset_raises_error_for_empty_data_split(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_get_dataset_raises_error_for_empty_data_split.pbtxt',
        name='some_dataset_name')
    expected_exception_message = ('The dataset in the config {} does not '
                                  'have a tfrecord_path.'
                                  .format(dataset_config_pbtext_filename))
    with self.assertRaisesRegexp(ValueError, expected_exception_message):
      data_providers.get_input_fn_from_dataset(
          dataset_config_pbtext_filename, mode=tf.estimator.ModeKeys.EVAL)

  def test_get_dataset_raises_error_for_empty_num_examples(self):
    dataset_config_pbtext_filename = _test_dataset_config(
        'test_get_dataset_raises_error_for_empty_num_examples.pbtxt',
        name='some_dataset_name',
        tfrecord_path='/path/to/dataset')
    expected_exception_message = ('The dataset in the config {} does not have '
                                  'a num_examples.'
                                  .format(dataset_config_pbtext_filename))
    with self.assertRaisesRegexp(ValueError, expected_exception_message):
      data_providers.get_input_fn_from_dataset(
          dataset_config_pbtext_filename, mode=tf.estimator.ModeKeys.EVAL)

  def test_dataset_definition(self):
    ds = data_providers.DeepVariantInput(
        mode=tf.estimator.ModeKeys.PREDICT,
        name='name',
        input_file_spec='test.tfrecord',
        num_examples=10,
        num_classes=dv_constants.NUM_CLASSES,
        tensor_shape=[11, 13, dv_constants.PILEUP_NUM_CHANNELS])
    self.assertEqual('name', ds.name)
    self.assertEqual('test.tfrecord', ds.input_file_spec)
    self.assertEqual(10, ds.num_examples)
    self.assertEqual(dv_constants.NUM_CLASSES, ds.num_classes)
    self.assertEqual([11, 13, dv_constants.PILEUP_NUM_CHANNELS],
                     ds.tensor_shape)

  def assertTfDataSetExamplesMatchExpected(self,
                                           input_fn,
                                           expected_dataset,
                                           use_tpu=False,
                                           workaround_list_files=False):
    # Note that we use input_fn to get an iterator, while we use
    # expected_dataset to get a filename, even though they are the same
    # type (DeepVariantInput), and may even be the same object.
    with tf.Session() as sess:
      params = {'batch_size': 1}
      batch_feed = input_fn(params).make_one_shot_iterator().get_next()

      sess.run(tf.global_variables_initializer())
      sess.run(tf.local_variables_initializer())
      seen = []
      while True:
        try:
          features, _ = sess.run(batch_feed)
        except tf.errors.OutOfRangeError:
          break
        locus = features['locus'][0]
        if use_tpu:
          locus = tf_utils.int_tensor_to_string(locus)
        # NB, this looks like: array(['chr20:10001019-10001019'], dtype=object)
        seen.append(locus)

    if workaround_list_files:
      # This really only works for loci, because those are string valued and
      # are expected to show up in sorted order.  For arbitrary data that's
      # not true.  In prod we have the version of tf that lets us turn off
      # shuffling so this path is skipped, but kokoro hits this.
      seen = sorted(seen)

    expected_loci = [
        example.features.feature['locus'].bytes_list.value[0]
        for example in io_utils.read_tfrecords(expected_dataset.input_file_spec)
    ]
    self.assertEqual(len(expected_loci), expected_dataset.num_examples)
    if seen != expected_loci:
      print('\n\nlen expected seen', len(expected_loci), len(seen))
      print('\n\nexpected=', expected_loci)
      print('\n\nseen=', seen)
    self.assertEqual(expected_loci, seen)
    # Note that this expected shape comes from the golden dataset. If the data
    # is remade in the future, the values might need to be modified accordingly.
    self.assertEqual(dv_constants.PILEUP_DEFAULT_DIMS,
                     expected_dataset.tensor_shape)

  @parameterized.parameters(
      dict(compressed_inputs=compressed_inputs, use_tpu=use_tpu)
      for compressed_inputs in [True, False]
      for use_tpu in [True, False])
  def test_reading_dataset(self, compressed_inputs, use_tpu):
    golden_dataset = make_golden_dataset(compressed_inputs, use_tpu=use_tpu)
    self.assertTfDataSetExamplesMatchExpected(
        input_fn=golden_dataset,
        expected_dataset=golden_dataset,
        use_tpu=use_tpu)

  # It looks like tf.data.Dataset.list_files is potentially nondeterministic.
  # There's no guaranteed way to get around that (yet, b/73959787).
  # A list_files() flag I want is only available in tf 1.7,
  # so for the short term, work around the problem by asking
  # self.assertTfDataSetExamplesMatchExpected to sort the
  # loci it sees.  That doesn't generalize well, but we should
  # be able to fix this soon.
  @parameterized.parameters(
      dict(compressed_inputs=compressed_inputs, use_tpu=use_tpu)
      for compressed_inputs in [True, False]
      for use_tpu in [True, False])
  def test_reading_sharded_dataset(self, compressed_inputs, use_tpu):
    golden_dataset = make_golden_dataset(compressed_inputs, use_tpu=use_tpu)
    n_shards = 3
    sharded_path = test_utils.test_tmpfile('sharded@{}'.format(n_shards))
    io_utils.write_tfrecords(
        io_utils.read_tfrecords(golden_dataset.input_file_spec), sharded_path)

    config_file = _test_dataset_config(
        'test_sharded.pbtxt',
        name='sharded_test',
        tfrecord_path=sharded_path,
        num_examples=golden_dataset.num_examples)

    self.assertTfDataSetExamplesMatchExpected(
        data_providers.get_input_fn_from_dataset(
            config_file, mode=tf.estimator.ModeKeys.EVAL),
        golden_dataset,
        # workaround_list_files is needed because wildcards, and so sharded
        # files, are nondeterministicly ordered (for now).
        workaround_list_files=True,
    )

  @parameterized.parameters(
      dict(compressed_inputs=compressed_inputs, mode=mode, use_tpu=use_tpu)
      for compressed_inputs in [True, False] for use_tpu in [True, False]
      for mode in ['TRAIN', 'EVAL'])
  def test_get_batches(self, compressed_inputs, mode, use_tpu):
    mode = (
        tf.estimator.ModeKeys.EVAL
        if mode == 'EVAL' else tf.estimator.ModeKeys.TRAIN)
    input_fn = make_golden_dataset(
        compressed_inputs, mode=mode, use_tpu=use_tpu)
    batch_size = 16
    with tf.Session() as sess:
      batch = input_fn(
          dict(batch_size=batch_size)).make_one_shot_iterator().get_next()

      # Get our images, labels, and variants for further testing.
      sess.run(tf.global_variables_initializer())
      features, labels = sess.run(batch)
      variants = features['variant']
      images = features['image']

      # Checks that our labels are the right shape and are one-hot encoded.
      # Note that the shape is 100, not 107, because we only adjust the image
      # in the model_fn now, where previously it was done in the input_fn.
      self.assertEqual([batch_size] + dv_constants.PILEUP_DEFAULT_DIMS,
                       list(images.shape))
      self.assertEqual((batch_size,), labels.shape)
      for label in labels:
        # pylint: disable=g-generic-assert
        self.assertTrue(0 <= label < dv_constants.NUM_CLASSES)

      # Check that our variants has the shape we expect and actually contain
      # variants by decoding them and checking the reference_name.
      self.assertEqual(batch_size, variants.shape[0])
      for variant in variants:
        if use_tpu:
          variant = tf_utils.int_tensor_to_string(variant)
        for v in variant_utils.decode_variants([variant]):
          self.assertEqual(v.reference_name, 'chr20')

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
    ds = data_providers.DeepVariantInput(
        mode=tf.estimator.ModeKeys.PREDICT,
        name='test_shape',
        input_file_spec=test_utils.test_tmpfile(tfrecord_path_to_match),
        num_examples=1)
    self.assertEqual(valid_shape, ds.tensor_shape)

  def test_get_shape_from_examples_path_invalid_path(self):
    with self.assertRaisesRegexp(Exception, '/this/path/does/not'):
      data_providers.DeepVariantInput(
          mode=tf.estimator.ModeKeys.PREDICT,
          name='test_invalid_path',
          input_file_spec='/this/path/does/not/exist',
          num_examples=1)


class InputTest(
    six.with_metaclass(parameterized.TestGeneratorMetaclass, tf.test.TestCase)):
  """Tests of input_fn, doing end-to-end I/O.

  These tests instantiate an input stream and then check it in various ways,
  in increasing complexity.
  """

  def get_batch_feed(self, batch_size=1, use_tpu=False):
    # This is an input_fn reading test_utils.N_GOLDEN_CALLING_EXAMPLES records.
    # Use PREDICT mode so we get finite input.
    dvi = data_providers.DeepVariantInput(
        mode=tf.estimator.ModeKeys.PREDICT,
        input_file_spec=testdata.GOLDEN_CALLING_EXAMPLES,
        num_examples=testdata.N_GOLDEN_CALLING_EXAMPLES,
        tensor_shape=None,
        use_tpu=use_tpu)
    params = {'batch_size': batch_size}
    batch_feed = dvi(params).make_one_shot_iterator().get_next()
    return batch_feed

  def check_batch_feed(self, batch_feed, use_tpu, expected_batch_size,
                       expected_n_batches):
    # Consume batch_feed, check that the right number of things is seen.
    with self.test_session() as sess:
      sess.run(tf.local_variables_initializer())
      sess.run(tf.global_variables_initializer())

      n = 0
      n_valid_entries = 0
      while True:
        try:
          features = sess.run(batch_feed)
        except tf.errors.OutOfRangeError:
          break
        n += 1
        a = features['image']  # np.ndarray
        self.assertTrue(a is not None)
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
        else:
          self.assertEqual(a.dtype, np.dtype('uint8'))
        current_batch_size = a.shape[0]
        self.assertLessEqual(current_batch_size, expected_batch_size)
        self.assertEqual(
            list(a.shape),
            [current_batch_size] + dv_constants.PILEUP_DEFAULT_DIMS)
        n_valid_entries += current_batch_size

      self.assertEqual(expected_n_batches, n)
      self.assertEqual(testdata.N_GOLDEN_CALLING_EXAMPLES, n_valid_entries)

  @parameterized.parameters(False, True)
  def testInputStream(self, use_tpu):
    # Read batch_feed one at a time, check the shape of each, and the
    # total count.
    batch_size = 1
    batch_feed = self.get_batch_feed(batch_size=batch_size, use_tpu=use_tpu)
    expected_n_batches = math.ceil(
        float(testdata.N_GOLDEN_CALLING_EXAMPLES) / batch_size)
    self.check_batch_feed(batch_feed, use_tpu, batch_size, expected_n_batches)

  @parameterized.parameters(False, True)
  def testBatching(self, use_tpu):
    # Test reading with a larger batch size.  Similar to testInputStream,
    # but note that the last batch may be truncated when not in predict mode,
    # so current_batch_size has to be recovered from the actual output.
    batch_size = 1024
    batch_feed = self.get_batch_feed(batch_size=batch_size, use_tpu=use_tpu)
    expected_n_batches = math.ceil(
        float(testdata.N_GOLDEN_CALLING_EXAMPLES) / batch_size)
    self.check_batch_feed(batch_feed, use_tpu, batch_size, expected_n_batches)

  @parameterized.parameters(False, True)
  def testGoldenCallingExamples(self, use_tpu):
    # Read the golden calling examples, and read the batch_feed instantiated
    # from the golden calling examples, and ensure that we get the same
    # parsed records in both cases.

    # Read and parse the canonical data.
    expected_decoded_records = list(
        io_utils.read_tfrecords(
            testdata.GOLDEN_CALLING_EXAMPLES, proto=example_pb2.Example))

    # Read and parse the data using tf.  This is the function under test,
    # although we indirectly check parse_tfexample as well.
    batch_feed = self.get_batch_feed(batch_size=1, use_tpu=use_tpu)

    with self.test_session() as sess:
      sess.run(tf.local_variables_initializer())
      sess.run(tf.global_variables_initializer())

      n = 0
      while True:
        # Read from batch.
        try:
          features = sess.run(batch_feed)
        except tf.errors.OutOfRangeError:
          break

        # Get the corresponding parsed golden example.
        example = expected_decoded_records[n]
        expected_alt_allele_indices_encoded = example.features.feature[
            'alt_allele_indices/encoded'].bytes_list.value[0]
        expected_variant_encoded = example.features.feature[
            'variant/encoded'].bytes_list.value[0]

        # Compare against the parsed batch feed.

        a = features['image'][0]  # np.ndarray
        self.assertEqual(list(a.shape), dv_constants.PILEUP_DEFAULT_DIMS)
        self.assertIsNotNone(a)
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
        else:
          self.assertEqual(a.dtype, np.dtype('uint8'))

        a = features['alt_allele_indices'][0]
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
          self.assertEqual(a.shape, (tf_utils.STRING_TO_INT_BUFFER_LENGTH,))
          actual_alt_allele_indices_encoded = tf_utils.int_tensor_to_string(a)
        else:
          self.assertIsInstance(a, six.string_types)
          actual_alt_allele_indices_encoded = a
        self.assertEqual(expected_alt_allele_indices_encoded,
                         actual_alt_allele_indices_encoded)

        a = features['variant'][0]
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
          self.assertEqual(a.shape, (tf_utils.STRING_TO_INT_BUFFER_LENGTH,))
          actual_variant_encoded = tf_utils.int_tensor_to_string(a)
        else:
          self.assertIsInstance(a, six.string_types)
          actual_variant_encoded = a
        self.assertEqual(expected_variant_encoded, actual_variant_encoded)

        n += 1

      self.assertEqual(n, testdata.N_GOLDEN_CALLING_EXAMPLES)


if __name__ == '__main__':
  absltest.main()
