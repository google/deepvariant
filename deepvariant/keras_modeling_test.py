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

from deepvariant import dv_config
from deepvariant import dv_constants
from deepvariant import keras_modeling
from third_party.nucleus.testing import test_utils


def _filepath_for_weights_ckpt_from_shape(shape: Tuple[int, int, int]) -> str:
  model = keras_modeling.inceptionv3(
      input_shape=shape, init_backbone_with_imagenet=False
  )
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
        input_shape=(100, 221, model_num_channel),
        init_backbone_with_imagenet=False,
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
        input_shape=(100, 221, model_num_channel),
        init_backbone_with_imagenet=False,
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
    model = keras_modeling.inceptionv3(
        input_shape=input_shape, init_backbone_with_imagenet=False
    )

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
        input_shape=input_shape,
        weights=weights_ckpt_path,
        init_backbone_with_imagenet=False,
    )
    x = np.random.uniform(size=(1,) + input_shape)
    y = model.predict(tf.keras.applications.inception_v3.preprocess_input(x))
    expected_shape = (1, dv_constants.NUM_CLASSES)
    self.assertEqual(y.shape, expected_shape)
    self.assertTrue(np.all(y >= 0))
    self.assertTrue(np.all(y <= 1))
    self.assertAlmostEqual(np.sum(y), 1.0, delta=1e-4)


class GetModelTest(parameterized.TestCase):

  @parameterized.parameters(
      dict(model_type='inception_v3'),
  )
  def test_retrieve_model_and_fn(self, model_type):
    config = dv_config.get_config('exome')

    with config.unlocked():
      config.model_type = model_type

    # This test should not throw any errors when retrieving the model
    # and it's corresponding preprocess function.
    keras_modeling.get_model(config)

  def test_get_model_error(self):
    config = dv_config.get_config('exome')

    with config.unlocked():
      config.model_type = 'not_a_model'

    with self.assertRaisesRegex(ValueError, 'Unsupported model type'):
      keras_modeling.get_model(config)

  def test_get_activations_model(self):
    # If you need to run with DeepVariant WGS v1.6.0 checkpoint, it can be found
    # in gs://deepvariant/models/DeepVariant/1.6.0/checkpoints/wgs/.
    # Here, for a unit test, I'll just hardcode the shape, initialize the model,
    # but not actually loading the weights.
    example_shape = [100, 221, 7]
    model = keras_modeling.inceptionv3(
        example_shape, init_backbone_with_imagenet=False
    )
    # Usually, here we would call something like:
    #     model.load_weights(checkpoint_path).expect_partial()
    # to load the weights. But here we're just testing the architecture, so
    # we skip loading specific weights.
    layers = ['max_pooling2d', 'mixed0', 'mixed5']
    activation_model = keras_modeling.get_activations_model(model, layers)
    input_image = np.zeros(example_shape, dtype=np.float32)
    input_image = tf.expand_dims(input_image, axis=0)
    output_layers = activation_model(input_image)
    self.assertCountEqual(output_layers.keys(), layers)
    self.assertEqual(
        list(output_layers['max_pooling2d'].shape), [1, 23, 53, 64]
    )
    self.assertEqual(list(output_layers['mixed0'].shape), [1, 10, 25, 256])
    self.assertEqual(list(output_layers['mixed5'].shape), [1, 4, 12, 768])


if __name__ == '__main__':
  absltest.main()
