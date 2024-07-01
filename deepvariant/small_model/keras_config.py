# Copyright 2024 Google LLC.
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
"""Module for configuring Keras models and training.

This module is used by the training and inference libraries.
"""
import os
from typing import Sequence
import ml_collections
import tensorflow as tf
import tensorflow_addons.metrics as tfa_metrics
from deepvariant.small_model import make_small_model_examples


@tf.keras.utils.register_keras_serializable(package="small_model")
class F1ScorePerClass(tfa_metrics.F1Score):
  """Reports F1 Score for a target class.

  # TODO: move custom metrics.py module. This module can be shared
  by both deepvariant/keras_modeling.py and small_model/train_utils.py
  """

  def __init__(
      self,
      num_classes: int = 0,
      target_class: int = 0,
      name: str = "",
      *args,
      **kwargs,
  ):
    self.target_class = target_class
    super().__init__(num_classes=num_classes, name=name)

  def result(self) -> tf.Tensor:
    return super().result()[self.target_class]


def keras_model_metrics() -> list[tf.keras.metrics.Metric]:
  """Returns a list of Keras model metrics."""
  return [
      tf.keras.metrics.CategoricalAccuracy(),
      tf.keras.metrics.CategoricalCrossentropy(),
      tf.keras.metrics.TruePositives(),
      tf.keras.metrics.TrueNegatives(),
      tf.keras.metrics.FalsePositives(),
      tf.keras.metrics.FalseNegatives(),
      tf.keras.metrics.Precision(),
      tf.keras.metrics.Precision(name="precision_homref", class_id=0),
      tf.keras.metrics.Precision(name="precision_het", class_id=1),
      tf.keras.metrics.Precision(name="precision_homalt", class_id=2),
      tf.keras.metrics.Recall(),
      tf.keras.metrics.Recall(name="recall_homref", class_id=0),
      tf.keras.metrics.Recall(name="recall_het", class_id=1),
      tf.keras.metrics.Recall(name="recall_homalt", class_id=2),
      tfa_metrics.F1Score(
          num_classes=3, average="weighted", name="f1_weighted"
      ),
      tfa_metrics.F1Score(num_classes=3, average="micro", name="f1_micro"),
      tfa_metrics.F1Score(num_classes=3, average="macro", name="f1_macro"),
      F1ScorePerClass(num_classes=3, target_class=0, name="f1_homref"),
      F1ScorePerClass(num_classes=3, target_class=1, name="f1_het"),
      F1ScorePerClass(num_classes=3, target_class=2, name="f1_homalt"),
  ]


def keras_mlp_model(model_params: ml_collections.ConfigDict) -> tf.keras.Model:
  """Creates a Keras MLP model."""
  model = tf.keras.Sequential()
  input_shape = len(make_small_model_examples.MODEL_FEATURES)
  hidden_layers = model_params.hidden_layer_sizes
  model.add(
      tf.keras.layers.Dense(
          hidden_layers[0],
          activation=model_params.activation,
          input_shape=(input_shape,),
      )
  )
  if len(hidden_layers) > 1:
    for layer_size in hidden_layers[1:]:
      model.add(
          tf.keras.layers.Dense(layer_size, activation=model_params.activation)
      )
  output_shape = len(make_small_model_examples.GenotypeEncoding)
  model.add(tf.keras.layers.Dense(output_shape, activation="softmax"))

  model.summary()
  model.compile(
      optimizer=model_params.optimizer,
      loss="categorical_crossentropy",
      metrics=keras_model_metrics(),
  )
  return model


def load_keras_model(checkpoint_path: str) -> tf.keras.Model:
  """Loads a Keras model from the given checkpoint path."""
  # return tf.saved_model.load(checkpoint_path)
  # tf.train.Checkpoint(checkpoint_path).restore(checkpoint_path)
  # return tf.saved_model.load(checkpoint_path)
  return tf.keras.models.load_model(checkpoint_path)


class LegacyFormatModelCheckpoint(tf.keras.callbacks.ModelCheckpoint):
  """Checkpoint callback that saves the model in the legacy format.

  This is required while we are on Keras 2.11.
  """

  def _save_handler(self, filepath):
    if filepath.endswith(".keras"):
      raise ValueError(
          "The filepath cannot end in .keras in Keras 2.11. Please remove the "
          "suffix and leave the file extension empty."
      )
    tf.keras.models.save_model(self.model, filepath, save_format="tf")


def get_keras_training_callbacks(
    checkpoint_filepath: str,
    is_xmanager_run: bool,
    tensorboard_log_dir: str,
    logging_frequency: int,
    batch_size: int,
    num_train_samples: int,
) -> Sequence[tf.keras.callbacks.Callback]:
  """If provisioned, provides a callback to log metrics to Xmanager/TB."""
  metric_loggers = []
  checkpoint_path_template = (
      "%s_epoch:{epoch:02d}_val_loss:{val_loss:.5f}" % checkpoint_filepath
  )
  model_checkpoint_callback = LegacyFormatModelCheckpoint(
      filepath=checkpoint_path_template,
      monitor="val_accuracy",
      mode="max",
      save_freq="epoch",
      save_best_only=False,
  )
  callbacks = [model_checkpoint_callback]

  return callbacks
