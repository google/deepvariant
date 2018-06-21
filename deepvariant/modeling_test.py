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
"""Tests for learning.genomics.deepvariant.modeling."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import six
import tensorflow as tf

from deepvariant import dv_constants
from deepvariant import modeling

slim = tf.contrib.slim


class ModelingTest(
    six.with_metaclass(parameterized.TestGeneratorMetaclass, tf.test.TestCase)):

  @parameterized.parameters(
      (model.name, type(model)) for model in modeling.all_models())
  def test_get_model_existing_models(self, model_name, expected_class):
    self.assertIsInstance(modeling.get_model(model_name), expected_class)

  def test_get_model_unknown_model_signals_error(self):
    with self.assertRaisesRegexp(ValueError, 'Unknown model'):
      modeling.get_model('unknown_model_1234')

  def test_make_deepvariant_slim_model(self):
    model = modeling.DeepVariantSlimModel(
        name='foo',
        n_classes_model_variable=['n_classes'],
        excluded_scopes_for_incompatible_shapes=['logits'],
        pretrained_model_path='path')

    self.assertEqual('foo', model.name)
    self.assertEqual(['n_classes'], model.n_classes_model_variable)
    self.assertEqual(['logits'], model.excluded_scopes_for_incompatible_shapes)
    self.assertEqual('path', model.pretrained_model_path)

  def test_variables_to_restore_from_model(self):
    model = modeling.DeepVariantModel('test', 'path')
    # We haven't created a slim model, so the variables_to_restore_from_model
    # should be returning an empty list.
    self.assertEqual([], model.variables_to_restore_from_model())

    # Create two model variable and one regular variables.
    with tf.variable_scope('model'):
      with tf.variable_scope('l1'):
        w1 = slim.model_variable('w1', shape=[10, 3, 3])
      with tf.variable_scope('l2'):
        w2 = slim.model_variable('w2', shape=[10, 3, 3])
        w3 = slim.model_variable('w3', shape=[10, 3, 3])
    v1 = slim.variable('my_var', shape=[20, 1])

    # The only variables in the system are the three we've created.
    self.assertItemsEqual([w1, w2, w3, v1], slim.get_variables())

    # We get just the three model variables without any excludes.
    self.assertItemsEqual([w1, w2, w3], model.variables_to_restore_from_model())
    # As well as when exclude_scopes is an empty list.
    self.assertItemsEqual(
        [w1, w2, w3], model.variables_to_restore_from_model(exclude_scopes=[]))

    # Excluding model/l1 variables gives us w2 and w3.
    self.assertItemsEqual(
        [w2, w3],
        model.variables_to_restore_from_model(exclude_scopes=['model/l1']))
    # Excluding model/l2 gives us just w1 back.
    self.assertItemsEqual(
        [w1],
        model.variables_to_restore_from_model(exclude_scopes=['model/l2']))
    # Excluding multiple scopes works as expected.
    self.assertItemsEqual(
        [],
        model.variables_to_restore_from_model(
            exclude_scopes=['model/l1', 'model/l2']))
    # Excluding the root model scope also produces no variables..
    self.assertItemsEqual(
        [], model.variables_to_restore_from_model(exclude_scopes=['model']))


# Hide the baseclass inside an enclosing scope so that unittest doesn't try to
# run our baseclass tests directly. http://stackoverflow.com/a/1323554.
class HiddenFromUnitTest(object):

  class SlimModelBaseTest(
      six.with_metaclass(parameterized.TestGeneratorMetaclass,
                         tf.test.TestCase)):

    @parameterized.parameters(
        dict(is_training=True),
        dict(is_training=False),
    )
    def test_create(self, is_training):
      # Creates a training=False model.
      self.assertEqual(len(tf.get_collection(tf.GraphKeys.UPDATE_OPS)), 0)
      images = tf.placeholder(
          tf.float32,
          (4, dv_constants.PILEUP_DEFAULT_HEIGHT,
           dv_constants.PILEUP_DEFAULT_WIDTH, dv_constants.PILEUP_NUM_CHANNELS))
      endpoints = self.model.create(
          images, dv_constants.NUM_CLASSES, is_training=is_training)
      if is_training:
        self.assertNotEqual(len(tf.get_collection(tf.GraphKeys.UPDATE_OPS)), 0)
      else:
        self.assertEqual(len(tf.get_collection(tf.GraphKeys.UPDATE_OPS)), 0)
      self.assertIn('Predictions', endpoints)
      self.assertIn('Logits', endpoints)
      self.assertEqual(endpoints['Predictions'].shape,
                       (4, dv_constants.NUM_CLASSES))

    def test_preprocess_images(self):
      with self.test_session() as sess:
        batch_size = 3
        values = range(91, 91 + 2 * 1 * dv_constants.PILEUP_NUM_CHANNELS)
        all_values = values * batch_size
        raw = np.array(
            all_values, dtype='uint8').reshape(
                (batch_size, 2, 1, dv_constants.PILEUP_NUM_CHANNELS))
        images = sess.run(self.model.preprocess_images(raw))
        for i in range(batch_size):
          image = images[i]

          # Check that our image has the right shape and that all values are
          # floats between between -1 and 1.
          self.assertEqual(tf.float32, image.dtype)
          self.assertTrue((image >= -1).all() and (image <= 1).all())
          self.assertEqual((2, 1, dv_constants.PILEUP_NUM_CHANNELS),
                           image.shape)

          # The preprocess step resizes the image to h x w as needed by
          # inception. We don't really care where it goes in the image (and the
          # calculation is complex. So we are simply checking here that all
          # values are zero except for the transformed values we see in values.
          # We are relying here on the tf operations to be correct and to not
          # change their behavior over time. Because we are doing assertEqual
          # we are also testing the order of the values, which means that we
          # are sure that the pixels have been translated in the right order in
          # the image, wherever the actual translation might be.
          self.assertEqual([(x - 128.0) / 128.0 for x in values],
                           [x for x in np.nditer(image) if x != 0.0])


class InceptionV3ModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    cls.model = modeling.get_model('inception_v3')

  # Note this test is only applied to inception_v3 since v2 and mobilenet don't
  # support some of these dimensions.
  @parameterized.parameters(
      dict(width=221, height=100),
      dict(width=221, height=200),
      dict(width=75, height=362),
  )
  def test_image_dimensions(self, width, height):
    with self.test_session():
      images = tf.placeholder(tf.float32, (4, height, width, 3))
      # We shouldn't get an exception creating images with these sizes.
      _ = self.model.create(images, 3, is_training=True)

  @parameterized.parameters(
      dict(width=73, height=100),
      dict(width=221, height=2000),
      dict(width=73, height=2000),
  )
  def test_bad_inception_v3_image_dimensions_get_custom_exception(
      self, width, height):
    with self.test_session():
      images = tf.placeholder(tf.float32, (4, height, width, 3))
      expected_message = ('Unsupported image dimensions.* model '
                          'inception_v3.*w={} x h={}.*').format(width, height)
      with self.assertRaisesRegexp(modeling.UnsupportedImageDimensions,
                                   expected_message):
        self.model.create(images, 3, is_training=True)


class InceptionV2ModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    cls.model = modeling.get_model('inception_v2')


class MobileNetModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    cls.model = modeling.get_model('mobilenet_v1')


class RandomGuessModelTest(tf.test.TestCase):

  def test_deterministic_predictions_for_fixed_seed(self):

    def predictions(seed):
      with self.test_session() as sess:
        model = modeling.DeepVariantRandomGuessModel(seed=seed)
        images = tf.placeholder(tf.float32, (4, 10, 10, 3))
        predictions = sess.run(model.create(images, 3, False)['Predictions'])
        return predictions

    # Note we do not use assertAllClose here as there's no assertNotAllClose().
    self.assertTrue((predictions(seed=123) == predictions(seed=123)).all())
    self.assertFalse((predictions(seed=123) == predictions(seed=456)).all())


if __name__ == '__main__':
  absltest.main()
