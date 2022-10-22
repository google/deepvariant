# Copyright 2022 Google LLC.
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
"""Provides an abstraction around deep learning Keras models in DeepVariant."""

import os
import tempfile
from typing import Optional, Tuple, Type

from absl import logging
import tensorflow as tf

from deepvariant import dv_constants

_DEFAULT_WEIGHT_DECAY = 0.00004
_DEFAULT_BACKBONE_DROP_DRATE = 0.2


def build_classification_head(inputs: tf.Tensor, l2: float = 0.0) -> tf.Tensor:
  """Returns an output head tensor configured for classification.

  In the future, this can be extended for regression, or with different params
  for different heads.

  Args:
    inputs: The backbone output tensor; used as the input to the head.
    l2: The l2 regularization factor used in `tf.keras.layers.Dense` layers.

  Returns:
    A tensor representing the output of the given head.
  """
  l2_regularizer = tf.keras.regularizers.L2(l2) if l2 else None
  head = tf.keras.layers.Dense(
      dv_constants.NUM_CLASSES,
      activation='softmax',
      dtype=tf.float32,
      name='classification',
      kernel_regularizer=l2_regularizer)
  return head(inputs)


def add_l2_regularizers(
    model: tf.keras.Model,
    layer_class: Type[tf.keras.layers.Layer],
    l2: float = _DEFAULT_WEIGHT_DECAY,
    regularizer_attr: str = 'kernel_regularizer') -> tf.keras.Model:
  """Adds L2 regularizers to all `layer_class` layers in `model`.

  Models from `tf.keras.applications` do not support specifying kernel or bias
  regularizers. However, adding regularization is important when fine tuning
  'imagenet' pretrained weights. In order to do this, this function updates the
  current model's configuration to include regularizers and reloads the model so
  that the newly created losses are registered.
  Note: this will not overwrite existing `kernel_regularizer` regularizers on
  the given layer.
  Args:
    model: The base model.
    layer_class: We add regularizers to all layers of type `layer_class`.
    l2: The l2 regularization factor.
    regularizer_attr: The layer's regularizer attribute.

  Returns:
    A model with l2 regularization added to each `layer_class` layer.
  """
  # Save the original model weights.
  tmp_weights_dir = tempfile.gettempdir()
  tmp_weights_path = os.path.join(tmp_weights_dir, 'tmp_weights.h5')
  model.save_weights(tmp_weights_path)

  # Clone the original model.
  reg_model = tf.keras.models.clone_model(model)

  # Set the L2 `regularizer_attr` on all layers of type `layer_class`. This
  # change is only reflected in the model's config file.
  num_regularizers_added = 0
  for layer in reg_model.layers:
    if not isinstance(layer, layer_class):
      continue
    if not hasattr(layer, regularizer_attr):
      continue
    if getattr(layer, regularizer_attr) is not None:
      continue
    setattr(layer, regularizer_attr, tf.keras.regularizers.l2(l2=l2))
    num_regularizers_added += 1

  # Save the updated model configuration.
  reg_model_json = reg_model.to_json()

  # Create a "new" model from the updated configuration and load the original
  # model's weights.
  reg_model = tf.keras.models.model_from_json(reg_model_json)
  reg_model.load_weights(tmp_weights_path, by_name=True)

  # Ensure model weights have not changed after adding regularization layers.
  for layer, reg_layer in zip(model.layers, reg_model.layers):
    weights = layer.weights
    reg_weights = reg_layer.weights
    if not weights:
      assert not reg_weights
    else:
      for i, weight in enumerate(weights):
        tf.debugging.assert_near(weight, reg_weights[i])

  # Ensure the newly added regularizers are registered as losses.
  assert len(reg_model.losses) == (len(model.losses) + num_regularizers_added)

  return reg_model


def inceptionv3(input_shape: Tuple[int, int, int],
                weights: Optional[str] = None) -> tf.keras.Model:
  """Returns an InceptionV3 architecture.

  See https://tensorflow.org/api_docs/python/tf/keras/applications/InceptionV3.

  Args:
    input_shape: a 3-tuple describing the input shape.
    weights: str. To initial weights from.

  Returns:
    An InceptionV3-based model.
  """
  backbone = tf.keras.applications.InceptionV3(
      include_top=False,
      # The old pretrained weights can't be loaded in
      # TF2/Keras, but could investigate converting the checkpoint file
      # (reference: internal)
      weights=None,
      input_shape=input_shape,
      classes=dv_constants.NUM_CLASSES,
      pooling='avg')

  weight_decay = _DEFAULT_WEIGHT_DECAY
  backbone = add_l2_regularizers(
      backbone, tf.keras.layers.Conv2D, l2=weight_decay)
  backbone_drop_rate = _DEFAULT_BACKBONE_DROP_DRATE

  inputs_image = tf.keras.Input(shape=input_shape, name='image')
  hid = backbone(inputs_image)
  hid = tf.keras.layers.Dropout(backbone_drop_rate)(hid)

  outputs = []
  outputs.append(build_classification_head(hid, l2=weight_decay))

  model = tf.keras.Model(
      inputs=[inputs_image], outputs=outputs, name='inceptionv3')
  model.summary()
  logging.info('Number of l2 regularizers: %s.', len(model.losses))
  if weights:
    logging.info('inceptionv3: load_weights from init_weights_file: %s',
                 weights)
    model.load_weights(weights)
  else:
    logging.info('inceptionv3: No init_weights_file.')
  return model
