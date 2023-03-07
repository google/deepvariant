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

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tensorflow as tf

from deepvariant import dv_constants
from deepvariant import keras_modeling


class KerasModelingTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name='Model with 7 channels',
          model_num_channel=7,
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


if __name__ == '__main__':
  absltest.main()
