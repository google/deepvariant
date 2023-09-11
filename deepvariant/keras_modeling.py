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
from typing import Callable, Optional, Tuple, Type, Union

from absl import logging
import ml_collections
import tensorflow as tf
import tensorflow_addons as tfa

from deepvariant import dv_constants
from deepvariant import dv_utils

_DEFAULT_WEIGHT_DECAY = 0.00004
_DEFAULT_BACKBONE_DROPOUT_RATE = 0.2


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
      kernel_regularizer=l2_regularizer,
  )
  return head(inputs)


def add_l2_regularizers(
    model: tf.keras.Model,
    layer_class: Type[tf.keras.layers.Layer],
    l2: float = _DEFAULT_WEIGHT_DECAY,
    regularizer_attr: str = 'kernel_regularizer',
) -> tf.keras.Model:
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


def load_weights_to_model_with_different_channels(
    model: tf.keras.Model, input_model: tf.keras.Model
) -> tf.keras.Model:
  """Initialize `model` with weights from `input_model` (different #channels).

  Args:
    model: The model we want to output.
    input_model: The input model that contains the weights to initialize from.

  Returns:
    `model` with updated weights from `input_model`
  """
  for layer_i, (input_model_layer, new_layer) in enumerate(
      zip(input_model.layers, model.layers)
  ):
    if not new_layer.weights:
      continue
    if len(new_layer.weights) != len(input_model_layer.weights):
      raise ValueError(
          'We expect input weights and the model we train to both '
          'be InceptionV3. The top level dict should be the same.'
      )

    # Create a list of ndarray, which will be used input for `set_weights`
    # at the end.
    new_weights_to_assign = new_layer.get_weights()

    for i, (input_model_layer_weights, new_layer_weights) in enumerate(
        zip(input_model_layer.get_weights(), new_layer.get_weights())
    ):
      if input_model_layer_weights.shape == new_layer_weights.shape:
        new_weights_to_assign[i] = input_model_layer_weights
      else:
        logging.info(
            (
                'input weights file layer %s:%s has shape %s, '
                'target model layer %s:%s has shape %s'
            ),
            input_model_layer.name,
            i,
            input_model_layer_weights.shape,
            new_layer.name,
            i,
            new_layer_weights.shape,
        )
        min_num_channels = min(
            input_model_layer_weights.shape[2], new_layer_weights.shape[2]
        )
        new_weights_to_assign[i][:, :, :min_num_channels, :] = (
            input_model_layer_weights[:, :, :min_num_channels, :]
        )
    # Now that the new_layer_weights list has the value we want to load with,
    # and has the right shape.
    model.layers[layer_i].set_weights(new_weights_to_assign)

  return model


def num_channels_from_checkpoint(filepath: str) -> int:
  """Determine the number of channels from a checkpoint path."""
  reader = tf.train.load_checkpoint(filepath)

  # Loop over variables in the checkpoint.
  for name in reader.get_variable_to_shape_map().keys():
    # 'layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE' seems to the main
    # variable to look at. I'm not sure if this heuristics will always work.
    # TODO
    if name.startswith(
        'layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE'
    ):
      weight_tensor = reader.get_tensor(name)
      return weight_tensor.shape[2]
    # We used to wire the model a bit differently.
    # If we see 'layer_with_weights-0/layer_with_weights-0/kernel', stop and
    # alert the user.
    if name.startswith('layer_with_weights-0/layer_with_weights-0/kernel'):
      raise ValueError(
          'You are using an older DeepVariant Keras model '
          'architecture. Please use a new model.'
      )
  raise ValueError('Unexpected model format.')


def inceptionv3_with_imagenet(
    input_shape: Tuple[int, int, int]
) -> tf.keras.Model:
  """Returns `inceptionv3` model with 3 channels; init with `weights=imagenet`.

  Our `inceptionv3` model (defined in this file as well) allows #channels other
  than 3. When the #channels is not 3, we couldn't set `weights=imagenet` to
  our backbone tf.keras.applications.InceptionV3 because it will complain the
  number of channels is not 3. We created this "inceptionv3_with_imagenet"
  model which has the same architecture, which we can then use our own
  "load_weights_to_model_with_different_channels" function to initiate an
  inceptionv3 model with any numbers of channels.

  The reason why this function is separate (instead of parameterized in the same
  implementation of inceptionv3) is to make it easier to read.

  Args:
    input_shape: a 3-tuple describing the input shape. The 3rd dimension is not
      used in this function. We always set that to 3 in this function.

  Returns:
    An InceptionV3-based model with 3 channels and init with `weights=imagenet`.
  """
  input_shape = list(input_shape)
  input_shape = [input_shape[0], input_shape[1], 3]

  backbone = tf.keras.applications.InceptionV3(
      include_top=False,
      weights='imagenet',
      input_shape=input_shape,
      classes=dv_constants.NUM_CLASSES,
      pooling='avg',
  )

  weight_decay = _DEFAULT_WEIGHT_DECAY
  backbone = add_l2_regularizers(
      backbone, tf.keras.layers.Conv2D, l2=weight_decay
  )
  backbone_drop_rate = _DEFAULT_BACKBONE_DROPOUT_RATE

  hid = tf.keras.layers.Dropout(backbone_drop_rate)(backbone.output)

  outputs = []
  outputs.append(build_classification_head(hid, l2=weight_decay))

  model = tf.keras.Model(
      inputs=backbone.input, outputs=outputs, name='inceptionv3'
  )
  return model


def inceptionv3(
    input_shape: Tuple[int, int, int],
    weights: Optional[str] = None,
    init_backbone_with_imagenet: bool = True,
    config: Optional[ml_collections.ConfigDict] = None,
) -> tf.keras.Model:
  """Returns an InceptionV3 architecture.

  See https://tensorflow.org/api_docs/python/tf/keras/applications/InceptionV3.

  Args:
    input_shape: a 3-tuple describing the input shape.
    weights: str. To initial weights from.
    init_backbone_with_imagenet: If True, get a model with InceptionV3 that has
      `weights='imagenet'` to start with. This will download a model. It should
      be set to False in unit tests, or when specific model weights will be
      loaded afterwards.
    config: a model configuration.

  Returns:
    An InceptionV3-based model.
  """
  backbone = tf.keras.applications.InceptionV3(
      include_top=False,
      weights=None,
      input_shape=input_shape,
      classes=dv_constants.NUM_CLASSES,
      pooling='avg',
  )

  if config:
    weight_decay = config.weight_decay
    backbone_dropout_rate = config.backbone_dropout_rate
  else:
    weight_decay = _DEFAULT_WEIGHT_DECAY
    backbone_dropout_rate = _DEFAULT_BACKBONE_DROPOUT_RATE

  backbone = add_l2_regularizers(
      backbone, tf.keras.layers.Conv2D, l2=weight_decay
  )

  hid = tf.keras.layers.Dropout(backbone_dropout_rate)(backbone.output)

  outputs = []
  outputs.append(build_classification_head(hid, l2=weight_decay))

  model = tf.keras.Model(
      inputs=backbone.input, outputs=outputs, name='inceptionv3'
  )
  model.summary()
  logging.info('Number of l2 regularizers: %s.', len(model.losses))

  # If no weights file is specified, initialize with `imagenet`.
  # The `init_backbone_with_imagenet` flag should be set to False for unit
  # tests to avoid loading the `imagenet` model from online.
  if not weights and init_backbone_with_imagenet:
    logging.info('inceptionv3: Initiate the model with imagenet (3 channels).')
    model = load_weights_to_model_with_different_channels(
        model, inceptionv3_with_imagenet(input_shape)
    )
    return model
  elif not weights and not init_backbone_with_imagenet:
    logging.info('inceptionv3: No init_weights_file.')
    return model

  weights_num_channels = num_channels_from_checkpoint(weights)
  # If the input weights have different number of channels, need some special
  # care:
  model_num_channels = input_shape[2]
  if weights_num_channels != model_num_channels:
    # This step is harder to do directly from `weights`, or even the Checkpoint
    # file format. So, create a `input_model` with expected #chanenls, load
    # the weights, and then post-process.
    # Improve later if possible: find a more readable alternative for this.
    weights_input_shape = list(input_shape)
    weights_input_shape[2] = weights_num_channels
    input_model = inceptionv3(
        tuple(weights_input_shape), weights, init_backbone_with_imagenet=False
    )
    logging.info(
        'inceptionv3: Assigning weights from %s channels to %s channels',
        weights_num_channels,
        model_num_channels,
    )
    model = load_weights_to_model_with_different_channels(model, input_model)
    return model
  else:
    logging.info(
        'inceptionv3: load_weights from init_weights_file: %s', weights
    )
    model.load_weights(weights)
    return model


class F1ScorePerClass(tfa.metrics.F1Score):
  """Reports F1 Score for a target class."""

  # TODO: Create custom metrics.py module.

  def __init__(self, num_classes: int, target_class: int, name: str):
    self.target_class = target_class
    super().__init__(num_classes=num_classes, name=name)

  def result(self) -> tf.Tensor:
    return super().result()[self.target_class]


def create_metrics():
  return [
      tf.keras.metrics.CategoricalAccuracy(),
      tf.keras.metrics.CategoricalCrossentropy(),
      tf.keras.metrics.TruePositives(),
      tf.keras.metrics.TrueNegatives(),
      tf.keras.metrics.FalsePositives(),
      tf.keras.metrics.FalseNegatives(),
      tf.keras.metrics.Precision(),
      tf.keras.metrics.Precision(name='precision_homref', class_id=0),
      tf.keras.metrics.Precision(name='precision_het', class_id=1),
      tf.keras.metrics.Precision(name='precision_homalt', class_id=2),
      tf.keras.metrics.Recall(),
      tf.keras.metrics.Recall(name='recall_homref', class_id=0),
      tf.keras.metrics.Recall(name='recall_het', class_id=1),
      tf.keras.metrics.Recall(name='recall_homalt', class_id=2),
      tfa.metrics.F1Score(
          num_classes=3, average='weighted', name='f1_weighted'
      ),
      tfa.metrics.F1Score(num_classes=3, average='micro', name='f1_micro'),
      tfa.metrics.F1Score(num_classes=3, average='macro', name='f1_macro'),
      F1ScorePerClass(num_classes=3, target_class=0, name='f1_homref'),
      F1ScorePerClass(num_classes=3, target_class=1, name='f1_het'),
      F1ScorePerClass(num_classes=3, target_class=2, name='f1_homalt'),
      # Leave mean loss as the last metric as it is updated differently.
      tf.keras.metrics.Mean(name='loss'),
  ]


def create_state(
    config: ml_collections.ConfigDict,
    model_dir: str,
    model: tf.keras.Model,
    optimizer: tf.optimizers.Optimizer,
    strategy: tf.distribute.Strategy,
) -> tf.train.CheckpointManager:
  """Initializes a checkpoint manager, and restores a checkpoint if one exists.

  Args:
    config: Training configuration.
    model_dir: Where model is stored.
    model: a tf Model.
    optimizer: A tf Optimizer.
    strategy: Distribution strategy.

  Returns:
    The state as `tf.train.Checkpoint`. This includes the `model` (network),
    the `optimizer`, metrics (train and tune), and the `global_step` variable.
  """
  with strategy.scope():
    # TODO: Load model and optimizer within this function.
    global_step = tf.Variable(
        0, trainable=False, name='global_step', dtype=tf.int64
    )
    early_stopping = tf.Variable(
        0, trainable=False, name='early_stopping', dtype=tf.int64
    )

    state = tf.train.Checkpoint(
        model,
        optimizer=optimizer,
        global_step=global_step,
        early_stopping=early_stopping,
        train_metrics=create_metrics(),
        tune_metrics=create_metrics(),
    )
    ckpt_manager = tf.train.CheckpointManager(
        checkpoint=state,
        directory=model_dir,
        max_to_keep=5,
    )
    ckpt_manager.restore_or_initialize()

    if ckpt_manager.latest_checkpoint:
      # Report current checkpoint state.
      best_checkpoint_metric_idx = [
          f'tune/{x.name}' for x in state.tune_metrics
      ].index(config.best_checkpoint_metric)
      logging.info(
          'Restored checkpoint %s at step=%s. %s=%s',
          os.path.basename(ckpt_manager.latest_checkpoint),
          global_step.numpy(),
          config.best_checkpoint_metric,
          state.tune_metrics[best_checkpoint_metric_idx].result(),
      )
    else:
      logging.info('Initialized Checkpoint')

    return ckpt_manager


def get_model(
    config: ml_collections.ConfigDict,
) -> Union[tf.keras.Model, Callable[..., tf.keras.Model]]:
  if config.model_type == 'inception_v3':
    return inceptionv3
  else:
    raise ValueError('Unsupported model type.')


def get_model_preprocess_fn(
    config: ml_collections.ConfigDict,
) -> Callable[[tf.train.Example], tf.train.Example]:
  """Returns the preprocess function for the model type."""

  if config.model_type == 'inception_v3':
    return dv_utils.preprocess_images
  else:
    raise ValueError('Unsupported model type.')
