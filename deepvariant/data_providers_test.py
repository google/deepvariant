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
"""Tests for learning.genomics.deepvariant.data_provider."""

import math



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tensorflow as tf
from tensorflow import estimator as tf_estimator

from deepvariant import data_providers
from deepvariant import dv_config
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io import tfrecord
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import variant_utils
from tensorflow.core.example import example_pb2

WGS_PILEUP_DEFAULT_DIMS = [
    dv_constants.PILEUP_DEFAULT_HEIGHT,
    dv_constants.PILEUP_DEFAULT_WIDTH,
    dv_constants.PILEUP_NUM_CHANNELS + 1,  # insert_size
]


def setUpModule():
  testdata.init()


# Return a DeepVariantInput attached to the golden training data.
# Run with shuffling off, and in eval mode.
def make_golden_dataset(
    compressed_inputs=False, mode=tf_estimator.ModeKeys.EVAL, use_tpu=False
):
  if compressed_inputs:
    source_path = test_utils.test_tmpfile('make_golden_dataset.tfrecord.gz')
    tfrecord.write_tfrecords(
        tfrecord.read_tfrecords(testdata.GOLDEN_TRAINING_EXAMPLES), source_path
    )
  else:
    source_path = testdata.GOLDEN_TRAINING_EXAMPLES
  return data_providers.get_input_fn_from_filespec(
      input_file_spec=source_path,
      num_examples=testdata.N_GOLDEN_TRAINING_EXAMPLES,
      name='labeled_golden',
      mode=mode,
      tensor_shape=None,
      use_tpu=use_tpu,
  )


def _test_dataset_config(filename, **kwargs):
  """Creates a DeepVariantDatasetConfig(**kwargs) and writes it to filename."""
  dataset_config_pbtext_filename = test_utils.test_tmpfile(filename)
  dataset_config = deepvariant_pb2.DeepVariantDatasetConfig(**kwargs)
  data_providers.write_dataset_config_to_pbtxt(
      dataset_config, dataset_config_pbtext_filename
  )
  return dataset_config_pbtext_filename


class ParseExampleTest(absltest.TestCase):

  def test_parse_example(self):
    path = testdata.GOLDEN_TRAINING_EXAMPLES
    ds = tf.data.TFRecordDataset(path, compression_type='GZIP')
    item = ds.take(1).get_single_element()
    input_shape = dv_utils.get_shape_from_examples_path(path)
    config = dv_config.get_config('exome')
    parse_example = data_providers.create_parse_example_fn(config)
    output = parse_example(item, input_shape)
    self.assertIsInstance(output, tuple)
    self.assertIsInstance(output[0], tf.Tensor)


class CreateExamplesTest(absltest.TestCase):

  def setUp(self):
    super().setUp()
    self.config = dv_config.get_config('exome')

  def test_invalid_mode(self):
    with self.assertRaisesRegex(ValueError, 'Mode must be set to'):
      _ = data_providers.input_fn(
          path=testdata.GOLDEN_TRAINING_EXAMPLES,
          config=self.config,
          mode='invalid_mode',
      )

  def test_create_input_dataset(self):
    ds = data_providers.input_fn(
        path=testdata.GOLDEN_TRAINING_EXAMPLES,
        config=self.config,
        mode='train',
    )
    item = ds.take(1).get_single_element()
    self.assertIsInstance(item, tuple)


class DataProviderTest(parameterized.TestCase):

  def test_dataset_definition(self):
    ds = data_providers.DeepVariantInput(
        mode=tf_estimator.ModeKeys.PREDICT,
        name='name',
        input_file_spec='test.tfrecord',
        num_examples=10,
        num_classes=dv_constants.NUM_CLASSES,
        tensor_shape=[11, 13, dv_constants.PILEUP_NUM_CHANNELS],
    )
    self.assertEqual('name', ds.name)
    self.assertEqual('test.tfrecord', ds.input_file_spec)
    self.assertEqual(10, ds.num_examples)
    self.assertEqual(dv_constants.NUM_CLASSES, ds.num_classes)
    self.assertEqual(
        [11, 13, dv_constants.PILEUP_NUM_CHANNELS], ds.tensor_shape
    )

  def assertTfDataSetExamplesMatchExpected(
      self,
      input_fn,
      expected_dataset,
      use_tpu=False,
      workaround_list_files=False,
  ):
    # Note that we use input_fn to get an iterator, while we use
    # expected_dataset to get a filename, even though they are the same
    # type (DeepVariantInput), and may even be the same object.
    with tf.compat.v1.Session() as sess:
      params = {'batch_size': 1}
      batch_feed = tf.compat.v1.data.make_one_shot_iterator(
          input_fn(params)
      ).get_next()

      sess.run(tf.compat.v1.global_variables_initializer())
      sess.run(tf.compat.v1.local_variables_initializer())
      seen = []
      while True:
        try:
          features, _ = sess.run(batch_feed)
        except tf.errors.OutOfRangeError:
          break
        locus = features['locus'][0]
        if use_tpu:
          locus = dv_utils.int_tensor_to_string(locus)
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
        for example in tfrecord.read_tfrecords(expected_dataset.input_file_spec)
    ]
    self.assertLen(expected_loci, expected_dataset.num_examples)
    if seen != expected_loci:
      print('\n\nlen expected seen', len(expected_loci), len(seen))
      print('\n\nexpected=', expected_loci)
      print('\n\nseen=', seen)
    self.assertEqual(expected_loci, seen)
    # Note that this expected shape comes from the golden dataset. If the data
    # is remade in the future, the values might need to be modified accordingly.
    self.assertEqual(WGS_PILEUP_DEFAULT_DIMS, expected_dataset.tensor_shape)

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(
      dict(compressed_inputs=compressed_inputs, use_tpu=use_tpu)
      for compressed_inputs in [True, False]
      for use_tpu in [True, False]
  )
  # pylint: enable=g-complex-comprehension
  def test_reading_dataset(self, compressed_inputs, use_tpu):
    golden_dataset = make_golden_dataset(compressed_inputs, use_tpu=use_tpu)
    self.assertTfDataSetExamplesMatchExpected(
        input_fn=golden_dataset,
        expected_dataset=golden_dataset,
        use_tpu=use_tpu,
    )

  @parameterized.parameters(
      dict(compressed_inputs=compressed_inputs, mode=mode, use_tpu=use_tpu)
      for compressed_inputs in [True, False]
      for use_tpu in [True, False]
      for mode in ['TRAIN', 'EVAL']
  )
  def test_get_batches(self, compressed_inputs, mode, use_tpu):
    mode = (
        tf_estimator.ModeKeys.EVAL
        if mode == 'EVAL'
        else tf_estimator.ModeKeys.TRAIN
    )
    input_fn = make_golden_dataset(
        compressed_inputs, mode=mode, use_tpu=use_tpu
    )
    batch_size = 16
    with tf.compat.v1.Session() as sess:
      batch = tf.compat.v1.data.make_one_shot_iterator(
          input_fn(dict(batch_size=batch_size))
      ).get_next()

      # Get our images, labels, and variants for further testing.
      sess.run(tf.compat.v1.global_variables_initializer())
      features, labels = sess.run(batch)
      variants = features['variant']
      images = features['image']

      # Checks that our labels are the right shape and are one-hot encoded.
      # Note that the shape is 100, not 107, because we only adjust the image
      # in the model_fn now, where previously it was done in the input_fn.
      self.assertEqual(
          [batch_size] + WGS_PILEUP_DEFAULT_DIMS, list(images.shape)
      )
      self.assertEqual((batch_size,), labels.shape)
      for label in labels:
        # pylint: disable=g-generic-assert
        self.assertTrue(0 <= label < dv_constants.NUM_CLASSES)

      # Check that our variants has the shape we expect and actually contain
      # variants by decoding them and checking the reference_name.
      self.assertEqual(batch_size, variants.shape[0])
      for variant in variants:
        if use_tpu:
          variant = dv_utils.int_tensor_to_string(variant)
        for v in variant_utils.decode_variants([variant]):
          self.assertEqual(v.reference_name, 'chr20')

  @parameterized.parameters(
      ('test_shape.gz', 'test_shape.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape@1.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-?????-of-00001.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-*.gz'),
      ('output', 'output'),
      ('test_shape-00000-of-00001', 'test_shape@1'),
      ('test_shape-00000-of-00001', 'test_shape-?????-of-00001'),
      ('test_shape-00000-of-00001', 'test_shape-*'),
  )
  def test_get_shape_from_examples_path(
      self, file_name_to_write, tfrecord_path_to_match
  ):
    example = example_pb2.Example()
    valid_shape = [1, 2, 3]
    example.features.feature['image/shape'].int64_list.value.extend(valid_shape)
    output_file = test_utils.test_tmpfile(file_name_to_write)
    tfrecord.write_tfrecords([example], output_file)
    ds = data_providers.DeepVariantInput(
        mode=tf_estimator.ModeKeys.PREDICT,
        name='test_shape',
        input_file_spec=test_utils.test_tmpfile(tfrecord_path_to_match),
        num_examples=1,
    )
    self.assertEqual(valid_shape, ds.tensor_shape)

  def test_get_shape_from_examples_path_invalid_path(self):
    with self.assertRaisesRegex(Exception, '/this/path/does/not'):
      data_providers.DeepVariantInput(
          mode=tf_estimator.ModeKeys.PREDICT,
          name='test_invalid_path',
          input_file_spec='/this/path/does/not/exist',
          num_examples=1,
      )

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(
      dict(max_examples=max_examples, batch_size=batch_size)
      for max_examples in [2, 4, 8]
      for batch_size in [4, 8, 16]
  )
  # pylint: enable=g-complex-comprehension
  def test_max_examples(self, max_examples, batch_size):
    input_fn = data_providers.get_input_fn_from_filespec(
        input_file_spec=testdata.GOLDEN_TRAINING_EXAMPLES,
        num_examples=testdata.N_GOLDEN_TRAINING_EXAMPLES,
        name='labeled_golden',
        max_examples=max_examples,
        mode=tf_estimator.ModeKeys.TRAIN,
    )

    n_batches_to_read = 100
    with tf.compat.v1.Session() as sess:
      sess.run(tf.compat.v1.global_variables_initializer())

      iterator = tf.compat.v1.data.make_one_shot_iterator(
          input_fn(dict(batch_size=batch_size))
      )
      next_element = iterator.get_next()

      def read_loci_in_batches():
        features, _ = sess.run(next_element)
        return features['locus']

      batches = [read_loci_in_batches() for _ in range(n_batches_to_read)]
      # pylint: disable=g-complex-comprehension
      unique_loci = {locus for batch in batches for locus in batch}
      # pylint: enable=g-complex-comprehension
      # assertLen not available OSS.
      # pylint: disable=g-generic-assert
      self.assertEqual(len(unique_loci), max_examples)
      # pylint: enable=g-generic-assert

  @parameterized.parameters(
      # When max_examples is None, dataset.num_examples will equal num_examples
      # arg.
      dict(num_examples=10, max_examples=None, expected=10),
      # When max_examples is larger than num_examples, dataset.num_examples will
      # equal the smaller value.
      dict(num_examples=10, max_examples=100, expected=10),
      # When max_examples is smaller than num_examples, dataset.num_examples
      # will equal the smaller max_examples value.
      dict(num_examples=10, max_examples=5, expected=5),
      # When num_examples isn't provided (None), but max_examples is, we don't
      # update num_examples so it remains None.
      dict(num_examples=None, max_examples=5, expected=None),
  )
  def test_max_examples_overrides_num_examples(
      self, num_examples, max_examples, expected
  ):
    dataset = data_providers.DeepVariantInput(
        # Use predict mode so we can have num_examples == None.
        mode=tf_estimator.ModeKeys.PREDICT,
        input_file_spec=testdata.GOLDEN_TRAINING_EXAMPLES,
        num_examples=num_examples,
        max_examples=max_examples,
    )
    self.assertEqual(expected, dataset.num_examples)

  def test_features_extraction_spec_for_mode(self):
    dataset = make_golden_dataset()

    shared_feature_names = {
        'image/encoded',
        'variant/encoded',
        'alt_allele_indices/encoded',
        'variant_type',
        'sequencing_type',
    }
    self.assertEqual(
        shared_feature_names,
        set(
            dataset.features_extraction_spec_for_mode(
                include_label_and_locus=False
            ).keys()
        ),
    )
    self.assertEqual(
        shared_feature_names.union({'label', 'locus'}),
        set(
            dataset.features_extraction_spec_for_mode(
                include_label_and_locus=True
            ).keys()
        ),
    )


class InputTest(
    tf.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
):
  """Tests of input_fn, doing end-to-end I/O.

  These tests instantiate an input stream and then check it in various ways,
  in increasing complexity.
  """

  def get_batch_feed(self, batch_size=1, use_tpu=False):
    # This is an input_fn reading test_utils.N_GOLDEN_CALLING_EXAMPLES records.
    # Use PREDICT mode so we get finite input.
    dvi = data_providers.DeepVariantInput(
        mode=tf_estimator.ModeKeys.PREDICT,
        input_file_spec=testdata.GOLDEN_CALLING_EXAMPLES,
        num_examples=testdata.N_GOLDEN_CALLING_EXAMPLES,
        tensor_shape=None,
        use_tpu=use_tpu,
    )
    params = {'batch_size': batch_size}
    batch_feed = tf.compat.v1.data.make_one_shot_iterator(
        dvi(params)
    ).get_next()
    return batch_feed

  def check_batch_feed(
      self, batch_feed, use_tpu, expected_batch_size, expected_n_batches
  ):
    # Consume batch_feed, check that the right number of things is seen.
    with tf.compat.v1.Session() as sess:
      sess.run(tf.compat.v1.local_variables_initializer())
      sess.run(tf.compat.v1.global_variables_initializer())

      n = 0
      n_valid_entries = 0
      while True:
        try:
          features = sess.run(batch_feed)
        except tf.errors.OutOfRangeError:
          break
        n += 1
        a = features['image']  # np.ndarray
        self.assertIsNot(a, None)
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
        else:
          self.assertEqual(a.dtype, np.dtype('uint8'))
        current_batch_size = a.shape[0]
        self.assertLessEqual(current_batch_size, expected_batch_size)
        self.assertEqual(
            list(a.shape),
            [current_batch_size] + WGS_PILEUP_DEFAULT_DIMS,
        )
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
        float(testdata.N_GOLDEN_CALLING_EXAMPLES) / batch_size
    )
    self.check_batch_feed(batch_feed, use_tpu, batch_size, expected_n_batches)

  @parameterized.parameters(False, True)
  def testBatching(self, use_tpu):
    # Test reading with a larger batch size.  Similar to testInputStream,
    # but note that the last batch may be truncated when not in predict mode,
    # so current_batch_size has to be recovered from the actual output.
    batch_size = 1024
    batch_feed = self.get_batch_feed(batch_size=batch_size, use_tpu=use_tpu)
    expected_n_batches = math.ceil(
        float(testdata.N_GOLDEN_CALLING_EXAMPLES) / batch_size
    )
    self.check_batch_feed(batch_feed, use_tpu, batch_size, expected_n_batches)

  @parameterized.parameters(False, True)
  def testGoldenCallingExamples(self, use_tpu):
    # Read the golden calling examples, and read the batch_feed instantiated
    # from the golden calling examples, and ensure that we get the same
    # parsed records in both cases.

    # Read and parse the canonical data.
    expected_decoded_records = list(
        tfrecord.read_tfrecords(
            testdata.GOLDEN_CALLING_EXAMPLES, proto=example_pb2.Example
        )
    )

    # Read and parse the data using tf.  This is the function under test,
    # although we indirectly check parse_tfexample as well.
    batch_feed = self.get_batch_feed(batch_size=1, use_tpu=use_tpu)

    with tf.compat.v1.Session() as sess:
      sess.run(tf.compat.v1.local_variables_initializer())
      sess.run(tf.compat.v1.global_variables_initializer())

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
            'alt_allele_indices/encoded'
        ].bytes_list.value[0]
        expected_variant_encoded = example.features.feature[
            'variant/encoded'
        ].bytes_list.value[0]
        expected_sequencing_type = example.features.feature[
            'sequencing_type'
        ].int64_list.value[0]

        # Compare against the parsed batch feed.

        a = features['image'][0]  # np.ndarray
        self.assertEqual(list(a.shape), WGS_PILEUP_DEFAULT_DIMS)
        self.assertIsNotNone(a)
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
        else:
          self.assertEqual(a.dtype, np.dtype('uint8'))

        a = features['alt_allele_indices'][0]
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
          self.assertEqual(a.shape, (dv_utils.STRING_TO_INT_BUFFER_LENGTH,))
          actual_alt_allele_indices_encoded = dv_utils.int_tensor_to_string(a)
        else:
          self.assertIsInstance(a, bytes)
          actual_alt_allele_indices_encoded = a
        self.assertEqual(
            expected_alt_allele_indices_encoded,
            actual_alt_allele_indices_encoded,
        )

        a = features['variant'][0]
        if use_tpu:
          self.assertEqual(a.dtype, np.dtype('int32'))
          self.assertEqual(a.shape, (dv_utils.STRING_TO_INT_BUFFER_LENGTH,))
          actual_variant_encoded = dv_utils.int_tensor_to_string(a)
        else:
          self.assertIsInstance(a, bytes)
          actual_variant_encoded = a
        self.assertEqual(expected_variant_encoded, actual_variant_encoded)

        self.assertEqual(features['sequencing_type'], expected_sequencing_type)

        n += 1

      self.assertEqual(n, testdata.N_GOLDEN_CALLING_EXAMPLES)


if __name__ == '__main__':
  tf.compat.v1.disable_eager_execution()
  absltest.main()
