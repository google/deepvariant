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
import ml_collections
import tensorflow as tf
from deepvariant import metrics
from deepvariant.small_model import make_small_model_examples


class LearningRateMetric(tf.keras.metrics.Metric):
  """Reports the learning rate of the optimizer."""

  def __init__(
      self,
      optimizer: tf.keras.optimizers.Optimizer,
      name="learning_rate",
      **kwargs,
  ):
    super().__init__(name=name, **kwargs)
    self.learning_rate = self.add_weight(name="lr", initializer="zeros")
    self.optimizer = optimizer

  def update_state(self, y_true, y_pred, sample_weight=None):
    self.learning_rate.assign(self.optimizer.learning_rate)

  def result(self):
    return self.learning_rate

  def reset_state(self):
    self.learning_rate.assign(0.0)


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
      tf.keras.metrics.F1Score(average="weighted", name="f1_weighted"),
      tf.keras.metrics.F1Score(num_classes=3, average="micro", name="f1_micro"),
      tf.keras.metrics.F1Score(num_classes=3, average="macro", name="f1_macro"),
      metrics.F1ScorePerClass(target_class=0, name="f1_homref"),
      metrics.F1ScorePerClass(target_class=1, name="f1_het"),
      metrics.F1ScorePerClass(target_class=2, name="f1_homalt"),
  ]


def get_learning_rate(
    config: ml_collections.ConfigDict,
) -> tf.keras.optimizers.schedules.ExponentialDecay:
  """Returns an exponential decay learning rate for the model."""
  model_params = config.model_params
  steps_per_epoch = config.num_train_samples // config.batch_size
  decay_steps = int(
      steps_per_epoch * model_params.learning_rate_num_epochs_per_decay
  )
  return tf.keras.optimizers.schedules.ExponentialDecay(
      initial_learning_rate=model_params.learning_rate,
      decay_steps=decay_steps,
      decay_rate=model_params.learning_rate_decay_rate,
      staircase=True,
  )


def keras_mlp_model(config: ml_collections.ConfigDict) -> tf.keras.Model:
  """Creates a Keras MLP model."""
  model_params = config.model_params
  model = tf.keras.Sequential()
  input_shape = len(
      make_small_model_examples.SmallModelExampleFactory(
          vaf_context_window_size=model_params.vaf_context_window_size,
          expand_by_haplotype=model_params.expand_by_haplotype,
      ).model_features
  )
  if model_params.features:
    input_shape = len(model_params.features)
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

  optimizer = tf.keras.optimizers.Adam(
      learning_rate=get_learning_rate(config),
      weight_decay=model_params.weight_decay,
  )
  model.summary()
  model.compile(
      optimizer=optimizer,
      steps_per_execution=model_params.steps_per_execution,
      loss="categorical_crossentropy",
      metrics=keras_model_metrics() + [LearningRateMetric(optimizer=optimizer)],
  )
  return model


def load_keras_model(checkpoint_path: str) -> tf.keras.Model:
  """Loads a Keras model from the given checkpoint path."""
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
) -> list[tf.keras.callbacks.Callback]:
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
