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

import os


from absl.testing import parameterized
import mock
import numpy as np
import six
import tensorflow as tf

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
        excluded_scopes=['logits'],
        pretrained_model_path='path')

    self.assertEqual('foo', model.name)
    self.assertEqual(['n_classes'], model.n_classes_model_variable)
    self.assertEqual(['logits'], model.excluded_scopes)
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
        {
            'is_training': True
        },
        {'is_training': False},
    )
    def test_create(self, is_training):
      # Creates a training=False model.
      self.assertEqual(len(tf.get_collection(tf.GraphKeys.UPDATE_OPS)), 0)
      h, w = (100, 221)
      images = tf.placeholder(tf.float32, (4, h, w, 3))
      endpoints = self.model.create(images, 3, is_training=is_training)
      if is_training:
        self.assertNotEqual(len(tf.get_collection(tf.GraphKeys.UPDATE_OPS)), 0)
      else:
        self.assertEqual(len(tf.get_collection(tf.GraphKeys.UPDATE_OPS)), 0)
      self.assertIn('Predictions', endpoints)
      self.assertEqual(endpoints['Predictions'].shape, (4, 3))

    def test_preprocess_image(self):
      with self.test_session() as sess:
        values = [91, 92, 93, 94, 95, 96]
        raw = np.array(values, dtype='uint8').reshape((2, 1, 3))
        image = sess.run(self.model.preprocess_image(raw))

        # Check that our image has the right shape and that all values are
        # floats between between -1 and 1.
        self.assertEqual(tf.float32, image.dtype)
        self.assertTrue((image >= -1).all() and (image <= 1).all())
        # Check that the shape is the same, except for inception_v3 model we had
        # to resize the height to at least 107.
        if self.model.name == 'inception_v3':
          self.assertEqual((107, 1, 3), image.shape)
        else:
          self.assertEqual((2, 1, 3), image.shape)

        # The preprocess step resizes the image to h x w as needed by
        # inception. We don't really care where it goes in the image (and the
        # calculation is complex. So we are simply checking here that all values
        # are zero except for the transformed values we see in values. We are
        # relying here on the tf operations to be correct and to not change
        # their behavior over time. Because we are doing assertEqual we are also
        # testing the order of the values, which means that we are sure that the
        # pixels have been translated in the right order in the image, wherever
        # the actual translation might be.
        self.assertEqual([(x - 128.0) / 128.0 for x in values],
                         [x for x in np.nditer(image) if x != 0.0])

    @parameterized.parameters([(3, True), (1000, True), (3, False)])
    @mock.patch('deepvariant'
                '.tf_utils.model_shapes')
    def test_initialize_from_checkpoint(self, n_checkpoint_classes, is_training,
                                        mock_model_shapes):

      def _create_checkpoint(checkpoint_dir, decay_factor):
        with self.test_session() as sess:
          checkpoint_prefix = os.path.join(checkpoint_dir, 'model')
          checkpoint_state_name = 'checkpoint'
          v1_data = [10.0] * 12
          v1 = slim.variables.Variable(v1_data, name='v1')
          slim.add_model_variable(v1)
          assign_to_v1 = v1.assign([20.0] * 12)
          variable_averages = tf.train.ExponentialMovingAverage(
              decay_factor, num_updates=999, zero_debias=False)
          average_op = variable_averages.apply([v1])
          v1_average = variable_averages.average(v1)
          self.assertItemsEqual([v1], slim.variables.moving_average_variables())
          self.assertEqual('v1/ExponentialMovingAverage:0', v1_average.name)
          sess.run(slim.variables.global_variables_initializer())
          sess.run(assign_to_v1)
          sess.run(average_op)
          saver = tf.train.Saver()
          saver.save(
              sess,
              checkpoint_prefix,
              global_step=3,
              latest_filename=checkpoint_state_name)
          return tf.train.latest_checkpoint(checkpoint_dir)

      checkpoint_dir = self.get_temp_dir()
      decay_factor = 0.9753
      checkpoint_file = _create_checkpoint(checkpoint_dir, decay_factor)
      n_prediction_classes = 3

      mock_model_shapes.return_value = {
          self.model.n_classes_model_variable: (n_checkpoint_classes,),
      }
      with self.test_session(graph=tf.Graph()) as sess:
        self.assertItemsEqual([], slim.variables.moving_average_variables())
        v1 = slim.variable_scope.get_variable('v1', [12])
        slim.add_model_variable(v1)
        sess.run(slim.variables.global_variables_initializer())
        self.model.initialize_from_checkpoint(
            checkpoint_file, n_prediction_classes, is_training)(
                sess)
        mock_model_shapes.assert_called_once_with(
            checkpoint_file, [self.model.n_classes_model_variable])
        v1_val = v1.eval(sess)
        if is_training:
          v1_expected = [20.0] * 12
        else:
          v1_expected = [10.0 * decay_factor + 20.0 * (1 - decay_factor)] * 12
        self.assertAllClose(v1_expected, v1_val)

    def test_initialize_from_checkpoint_fails_with_bad_path(self):
      with self.assertRaisesRegexp(ValueError, 'Checkpoint cannot be None'):
        self.model.initialize_from_checkpoint(None, 3, True)

    @mock.patch('deepvariant.tf_utils.model_shapes')
    def test_initialize_raises_inference_shape_mismatch(self,
                                                        mock_model_shapes):
      mock_model_shapes.return_value = {
          self.model.n_classes_model_variable: (100,)
      }
      with self.assertRaisesRegexp(
          ValueError, ('Checkpoint has 100 classes but we want to use 3 and '
                       'is_training=False')):
        self.model.initialize_from_checkpoint('path', 3, False)

    @mock.patch('deepvariant'
                '.modeling.slim.losses.softmax_cross_entropy')
    @mock.patch('deepvariant'
                '.modeling.slim.losses.get_total_loss')
    def test_loss(self, mock_total_loss, mock_cross):
      endpoints = {'Logits': 'Logits'}
      labels = [[0, 1, 0], [1, 0, 0]]
      actual = self.model.loss(endpoints, labels)
      mock_total_loss.assert_called_once_with()
      self.assertEqual(actual, mock_total_loss.return_value)
      # We really only want to test that the endpoints and labels are being
      # passed in correctly, without specifying all of the additional keywords
      # again here like weight that can change without breaking the correctness
      # of loss.
      self.assertEqual([('Logits', labels)],
                       [args for args, _ in mock_cross.call_args_list])


class InceptionV3ModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    cls.model = modeling.get_model('inception_v3')


class InceptionV2ModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    cls.model = modeling.get_model('inception_v2')


class Resnet50ModelTest(HiddenFromUnitTest.SlimModelBaseTest):

  @classmethod
  def setUpClass(cls):
    cls.model = modeling.get_model('resnet_v2_50')


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
  tf.test.main()
