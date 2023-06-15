# Copyright 2023 Google LLC.
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
"""Tests for keras_modeling."""

from typing import Tuple

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tensorflow as tf

from deepvariant import dv_constants
from deepvariant import keras_modeling
from third_party.nucleus.testing import test_utils


def _filepath_for_weights_ckpt_from_shape(shape: Tuple[int, int, int]) -> str:
  model = keras_modeling.inceptionv3(input_shape=shape)
  tmp_model_name = test_utils.test_tmpfile(
      'filepath_for_weights_ckpt_from_shape_{}.weights'.format(
          '_'.join(str(i) for i in shape)
      )
  )
  model.save_weights(tmp_model_name)
  return tmp_model_name


class KerasModelingTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name='Model with 3 channels',
          model_num_channel=3,
      ),
  )
  def test_prediction(self, model_num_channel):
    model = keras_modeling.inceptionv3(
        input_shape=(100, 221, model_num_channel)
    )
    input_shape = (1, 100, 221, model_num_channel)
    x = np.random.uniform(size=input_shape)
    y = model.predict(tf.keras.applications.inception_v3.preprocess_input(x))
    expected_shape = (1, dv_constants.NUM_CLASSES)
    self.assertEqual(y.shape, expected_shape)
    self.assertTrue(np.all(y >= 0))
    self.assertTrue(np.all(y <= 1))
    self.assertAlmostEqual(np.sum(y), 1.0, delta=1e-4)

  @parameterized.named_parameters(
      dict(
          testcase_name='Model and input have different number of channels',
          model_num_channel=7,
          input_num_channel=6,
      ),
  )
  def test_prediction_failure(self, model_num_channel, input_num_channel):
    # Confirm that imcompatible shape causes an issue.
    model = keras_modeling.inceptionv3(
        input_shape=(100, 221, model_num_channel)
    )
    input_shape = (1, 100, 221, input_num_channel)
    x = np.random.uniform(size=input_shape)
    with self.assertRaisesRegex(
        ValueError,
        'Input 0 of layer "inceptionv3" is incompatible with the layer: ',
    ):
      _ = model.predict(tf.keras.applications.inception_v3.preprocess_input(x))

  @parameterized.named_parameters(
      dict(
          testcase_name='Model with 3 channels',
          model_num_channel=3,
      ),
  )
  def test_model_training(self, model_num_channel):
    # Define a model.
    input_shape = (100, 221, model_num_channel)
    model = keras_modeling.inceptionv3(input_shape=input_shape)

    # Generate random input and target data.
    x_train = np.random.rand(32, *input_shape)
    y_train = np.random.randint(dv_constants.NUM_CLASSES, size=(32,))

    # Compile the model.
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy')

    # Train the model for a few steps.
    model.train_on_batch(x_train, y_train)
    model.train_on_batch(x_train, y_train)

    # Check that the loss has decreased after training.
    loss_before = model.evaluate(x_train, y_train, verbose=0)
    model.train_on_batch(x_train, y_train)
    loss_after = model.evaluate(x_train, y_train, verbose=0)
    self.assertNotEqual(
        loss_after, loss_before, 'Training is expected to change loss'
    )

  @parameterized.parameters(
      dict(input_shape=(75, 75, 8)),
      dict(input_shape=(100, 221, 4)),
      dict(input_shape=(300, 100, 5)),
  )
  def test_num_channels_from_checkpoint(self, input_shape):
    # Create a model and save it to a checkpoint. Then test whether we can
    # detect its number of channels correctly.
    weights_ckpt_path = _filepath_for_weights_ckpt_from_shape(input_shape)

    # Load it back and determine the num_channels.
    detected_num_channels = keras_modeling.num_channels_from_checkpoint(
        weights_ckpt_path
    )
    self.assertEqual(detected_num_channels, input_shape[-1])

  @parameterized.named_parameters(
      dict(
          testcase_name='Weights ckpt has the same shape as the model.',
          checkpoint_weights_shape=(75, 75, 3),
          input_shape=(75, 75, 3),
      ),
      dict(
          testcase_name='Weights ckpt has fewer #channels than the model.',
          checkpoint_weights_shape=(75, 75, 3),
          input_shape=(75, 75, 4),
      ),
      dict(
          testcase_name='Weights ckpt has more #channels than the model.',
          checkpoint_weights_shape=(75, 75, 4),
          input_shape=(75, 75, 3),
      ),
      dict(
          testcase_name='Larger-than-input height/width do not break.',
          checkpoint_weights_shape=(100, 100, 3),
          input_shape=(75, 75, 4),
      ),
      dict(
          testcase_name='Smaller-than-input height/width do not break.',
          checkpoint_weights_shape=(75, 75, 3),
          input_shape=(80, 80, 4),
      ),
  )
  def test_inceptionv3_with_init_weights(
      self, checkpoint_weights_shape, input_shape
  ):
    """keras_modeling.inceptionv3 can load weights (even different #channels).

    Args:
      checkpoint_weights_shape: The shape of the weights (checkpoint) file.
      input_shape: The shape of the model we're training now.
    """
    # Create a model and save it to a checkpoint. Then test whether we can
    # detect its number of channels correctly.
    weights_ckpt_path = _filepath_for_weights_ckpt_from_shape(
        checkpoint_weights_shape
    )

    model = keras_modeling.inceptionv3(
        input_shape=input_shape, weights=weights_ckpt_path
    )
    x = np.random.uniform(size=(1,) + input_shape)
    y = model.predict(tf.keras.applications.inception_v3.preprocess_input(x))
    expected_shape = (1, dv_constants.NUM_CLASSES)
    self.assertEqual(y.shape, expected_shape)
    self.assertTrue(np.all(y >= 0))
    self.assertTrue(np.all(y <= 1))
    self.assertAlmostEqual(np.sum(y), 1.0, delta=1e-4)


if __name__ == '__main__':
  absltest.main()
