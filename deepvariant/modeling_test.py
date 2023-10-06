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
"""Tests for learning.genomics.deepvariant.modeling."""



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tensorflow as tf
import tf_slim

from deepvariant import dv_constants
from deepvariant import dv_utils_using_clif
from deepvariant import modeling

tf.compat.v1.disable_eager_execution()
slim = tf_slim


class ModelingTest(
    tf.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
):

  @parameterized.parameters(
      (model_class().name, type(model_class()))
      for model_class in modeling.all_models()
  )
  def test_get_model_existing_models(self, model_name, expected_class):
    self.assertIsInstance(modeling.get_model(model_name), expected_class)

  def test_get_model_unknown_model_signals_error(self):
    with self.assertRaisesRegex(ValueError, 'Unknown model'):
      modeling.get_model('unknown_model_1234')

  def test_make_deepvariant_slim_model(self):
    model = modeling.DeepVariantSlimModel(
        name='foo',
        n_classes_model_variable=['n_classes'],
        excluded_scopes_for_incompatible_classes=['logits'],
        excluded_scopes_for_incompatible_channels=['logits'],
        pretrained_model_path='path',
    )

    self.assertEqual('foo', model.name)
    self.assertEqual(['n_classes'], model.n_classes_model_variable)
    self.assertEqual(['logits'], model.excluded_scopes_for_incompatible_classes)
    self.assertEqual(
        ['logits'], model.excluded_scopes_for_incompatible_channels
    )
    self.assertEqual('path', model.pretrained_model_path)

  def test_is_encoded_variant_type(self):
    types = [
        dv_utils_using_clif.EncodedVariantType.SNP.value,
        dv_utils_using_clif.EncodedVariantType.INDEL.value,
    ]
    tensor = tf.constant(types * 4, dtype=tf.int64)

    def _run(tensor_to_run):
      with self.test_session() as sess:
        return list(sess.run(tensor_to_run))

    self.assertEqual(
        _run(
            modeling.is_encoded_variant_type(
                tensor, dv_utils_using_clif.EncodedVariantType.SNP
            )
        ),
        [True, False] * 4,
    )
    self.assertEqual(
        _run(
            modeling.is_encoded_variant_type(
                tensor, dv_utils_using_clif.EncodedVariantType.INDEL
            )
        ),
        [False, True] * 4,
    )

  @parameterized.parameters(
      dict(labels=[0, 2, 1, 0], target_class=0, expected=[0, 1, 1, 0]),
      dict(labels=[0, 2, 1, 0], target_class=1, expected=[1, 1, 0, 1]),
      dict(labels=[0, 2, 1, 0], target_class=2, expected=[1, 0, 1, 1]),
  )
  def test_binarize(self, labels, target_class, expected):
    with self.test_session() as sess:
      result = sess.run(
          modeling.binarize(np.array(labels), np.array(target_class))
      )
      self.assertListEqual(result.tolist(), expected)

  @parameterized.parameters([True, False])
  def test_eval_metric_fn(self, include_variant_types):
    labels = tf.constant([1, 0], dtype=tf.int64)
    predictions = tf.constant([[1, 0], [0, 1]], dtype=tf.int64)
    if include_variant_types:
      variant_types = tf.constant([0, 1], dtype=tf.int64)
    else:
      variant_types = None

    expected = modeling.eval_function_metrics(
        has_variant_types=include_variant_types
    )
    actual = modeling.eval_metric_fn(labels, predictions, variant_types)
    self.assertEqual(set(expected.keys()), set(actual.keys()))

  def test_variables_to_restore_from_model(self):
    model = modeling.DeepVariantModel('test', 'path')
    # We haven't created a slim model, so the variables_to_restore_from_model
    # should be returning an empty list.
    self.assertEqual([], model.variables_to_restore_from_model())

    # Create two model variable and one regular variables.
    with tf.compat.v1.variable_scope('model'):
      with tf.compat.v1.variable_scope('l1'):
        w1 = slim.model_variable('w1', shape=[10, 3, 3])
      with tf.compat.v1.variable_scope('l2'):
        w2 = slim.model_variable('w2', shape=[10, 3, 3])
        w3 = slim.model_variable('w3', shape=[10, 3, 3])
    v1 = slim.variable('my_var', shape=[20, 1])

    # The only variables in the system are the three we've created.
    self.assertCountEqual([w1, w2, w3, v1], slim.get_variables())

    # We get just the three model variables without any excludes.
    self.assertCountEqual([w1, w2, w3], model.variables_to_restore_from_model())
    # As well as when exclude_scopes is an empty list.
    self.assertCountEqual(
        [w1, w2, w3], model.variables_to_restore_from_model(exclude_scopes=[])
    )

    # Excluding model/l1 variables gives us w2 and w3.
    self.assertCountEqual(
        [w2, w3],
        model.variables_to_restore_from_model(exclude_scopes=['model/l1']),
    )
    # Excluding model/l2 gives us just w1 back.
    self.assertCountEqual(
        [w1], model.variables_to_restore_from_model(exclude_scopes=['model/l2'])
    )
    # Excluding multiple scopes works as expected.
    self.assertCountEqual(
        [],
        model.variables_to_restore_from_model(
            exclude_scopes=['model/l1', 'model/l2']
        ),
    )
    # Excluding the root model scope also produces no variables..
    self.assertCountEqual(
        [], model.variables_to_restore_from_model(exclude_scopes=['model'])
    )


# Hide the baseclass inside an enclosing scope so that unittest doesn't try to
# run our baseclass tests directly. http://stackoverflow.com/a/1323554.
class HiddenFromUnitTest(object):

  class SlimModelBaseTest(
      tf.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
  ):

    @parameterized.parameters(
        dict(is_training=True),
        dict(is_training=False),
    )
    def test_create(self, is_training):
      # Creates a training=False model.
      self.assertEqual(
          len(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)), 0
      )
      images = tf.compat.v1.placeholder(
          tf.float32,
          (
              4,
              dv_constants.PILEUP_DEFAULT_HEIGHT,
              dv_constants.PILEUP_DEFAULT_WIDTH,
              dv_constants.PILEUP_NUM_CHANNELS,
          ),
      )
      endpoints = self.model.create(
          images, dv_constants.NUM_CLASSES, is_training=is_training
      )
      if is_training:
        self.assertNotEqual(
            len(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)),
            0,
        )
      else:
        self.assertEqual(
            len(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)),
            0,
        )
      self.assertIn('Predictions', endpoints)
      self.assertIn('Logits', endpoints)
      self.assertEqual(
          endpoints['Predictions'].shape, (4, dv_constants.NUM_CLASSES)
      )

    def test_preprocess_images(self):
      with self.test_session() as sess:
        batch_size = 3
        values = range(91, 91 + 2 * 1 * dv_constants.PILEUP_NUM_CHANNELS)
        all_values = list(np.tile(values, batch_size))
        raw = np.array(all_values, dtype='uint8').reshape(
            (batch_size, 2, 1, dv_constants.PILEUP_NUM_CHANNELS)
        )
        images = sess.run(self.model.preprocess_images(raw))
        for i in range(batch_size):
          image = images[i]

          # Check that our image has the right shape and that all values are
          # floats between between -1 and 1.
          self.assertEqual(tf.float32, image.dtype)
          self.assertTrue((image >= -1).all() and (image <= 1).all())
          self.assertEqual(
              (2, 1, dv_constants.PILEUP_NUM_CHANNELS), image.shape
          )

          # The preprocess step resizes the image to h x w as needed by
          # inception. We don't really care where it goes in the image (and the
          # calculation is complex. So we are simply checking here that all
          # values are zero except for the transformed values we see in values.
          # We are relying here on the tf operations to be correct and to not
          # change their behavior over time. Because we are doing assertEqual
          # we are also testing the order of the values, which means that we
          # are sure that the pixels have been translated in the right order in
          # the image, wherever the actual translation might be.
          self.assertEqual(
              [(x - 128.0) / 128.0 for x in values],
              [x for x in np.nditer(image) if x != 0.0],
          )


class InceptionV3ModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    super(InceptionV3ModelTest, cls).setUpClass()
    cls.model = modeling.get_model('inception_v3')

  # Note this test is only applied to inception_v3.
  @parameterized.parameters(
      dict(width=221, height=100),
      dict(width=221, height=200),
      dict(width=75, height=362),
  )
  def test_image_dimensions(self, width, height):
    with self.test_session():
      images = tf.compat.v1.placeholder(tf.float32, (4, height, width, 3))
      # We shouldn't get an exception creating images with these sizes.
      _ = self.model.create(images, 3, is_training=True)

  @parameterized.parameters(
      dict(width=73, height=100),
      dict(width=221, height=2000),
      dict(width=73, height=2000),
  )
  def test_bad_inception_v3_image_dimensions_get_custom_exception(
      self, width, height
  ):
    with self.test_session():
      images = tf.compat.v1.placeholder(tf.float32, (4, height, width, 3))
      expected_message = (
          'Unsupported image dimensions.* model inception_v3.*w={} x h={}.*'
      ).format(width, height)
      with self.assertRaisesRegex(
          modeling.UnsupportedImageDimensionsError, expected_message
      ):
        self.model.create(images, 3, is_training=True)


class InceptionV3EmbeddingModelTest(
    tf.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
):

  @classmethod
  def setUpClass(cls):
    super(InceptionV3EmbeddingModelTest, cls).setUpClass()
    cls.model = modeling.get_model('inception_v3_embedding')

  @parameterized.parameters(
      dict(is_training=True),
      dict(is_training=False),
  )
  def test_create(self, is_training):
    self.assertEqual(
        len(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)), 0
    )
    images = tf.compat.v1.placeholder(
        tf.float32,
        (
            4,
            dv_constants.PILEUP_DEFAULT_HEIGHT,
            dv_constants.PILEUP_DEFAULT_WIDTH,
            dv_constants.PILEUP_NUM_CHANNELS,
        ),
    )
    seq_type = tf.compat.v1.placeholder(tf.int64, (4,))
    endpoints = self.model._create(
        (images, seq_type), dv_constants.NUM_CLASSES, is_training=is_training
    )
    if is_training:
      self.assertNotEqual(
          len(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)), 0
      )
    else:
      self.assertEqual(
          len(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)), 0
      )
    self.assertIn('Predictions', endpoints)
    self.assertIn('Logits', endpoints)
    self.assertEqual(
        endpoints['Predictions'].shape, (4, dv_constants.NUM_CLASSES)
    )
    self.assertIn('Embeddings', endpoints)
    self.assertEqual(
        endpoints['Embeddings'].shape, (4, 2048 + self.model.embedding_size)
    )

  def test_create_embeddings(self):
    indices = tf.compat.v1.placeholder(tf.int64, (4,))
    embeddings = self.model._create_embeddings(indices)
    self.assertEqual(embeddings.shape, (4, self.model.embedding_size))

  def test_embedding_lookup(self):
    indices = tf.compat.v1.placeholder(tf.int64, (4,))
    embeddings = self.model._embedding_lookup(indices)
    self.assertEqual(embeddings.shape, (4, self.model.embedding_size))


if __name__ == '__main__':
  absltest.main()
