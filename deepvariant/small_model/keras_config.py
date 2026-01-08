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
import keras
import ml_collections
import tensorflow as tf
from deepvariant.small_model import make_small_model_examples


@keras.saving.register_keras_serializable(package="CustomMetrics")
class LearningRateMetric(keras.metrics.Metric):
  """Reports the learning rate of the optimizer."""

  def __init__(
      self,
      optimizer: keras.optimizers.Optimizer,
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


@keras.saving.register_keras_serializable(package="CustomMetrics")
class F1ScorePerClass(keras.metrics.F1Score):
  """Reports F1 Score for a target class."""


  def __init__(self, target_class: int, name: str):
    self.target_class = target_class
    super().__init__(name=name)

  def result(self) -> tf.Tensor:
    return super().result()[self.target_class]


def keras_model_metrics() -> list[keras.metrics.Metric]:
  """Returns a list of Keras model metrics."""
  return [
      keras.metrics.CategoricalAccuracy(),
      keras.metrics.CategoricalCrossentropy(),
      keras.metrics.TruePositives(),
      keras.metrics.TrueNegatives(),
      keras.metrics.FalsePositives(),
      keras.metrics.FalseNegatives(),
      keras.metrics.Precision(),
      keras.metrics.Precision(name="precision_homref", class_id=0),
      keras.metrics.Precision(name="precision_het", class_id=1),
      keras.metrics.Precision(name="precision_homalt", class_id=2),
      keras.metrics.Recall(),
      keras.metrics.Recall(name="recall_homref", class_id=0),
      keras.metrics.Recall(name="recall_het", class_id=1),
      keras.metrics.Recall(name="recall_homalt", class_id=2),
      keras.metrics.F1Score(average="weighted", name="f1_weighted"),
      keras.metrics.F1Score(average="micro", name="f1_micro"),
      keras.metrics.F1Score(average="macro", name="f1_macro"),
      F1ScorePerClass(target_class=0, name="f1_homref"),
      F1ScorePerClass(target_class=1, name="f1_het"),
      F1ScorePerClass(target_class=2, name="f1_homalt"),
  ]


def get_learning_rate(
    config: ml_collections.ConfigDict,
) -> keras.optimizers.schedules.ExponentialDecay:
  """Returns an exponential decay learning rate for the model."""
  model_params = config.model_params
  steps_per_epoch = config.num_train_samples // config.batch_size
  decay_steps = int(
      steps_per_epoch * model_params.learning_rate_num_epochs_per_decay
  )
  return keras.optimizers.schedules.ExponentialDecay(
      initial_learning_rate=model_params.learning_rate,
      decay_steps=decay_steps,
      decay_rate=model_params.learning_rate_decay_rate,
      staircase=True,
  )


def keras_mlp_model(config: ml_collections.ConfigDict) -> keras.Model:
  """Creates a Keras MLP model."""
  model_params = config.model_params
  model = keras.Sequential()
  input_shape = len(
      make_small_model_examples.SmallModelExampleFactory(
          vaf_context_window_size=model_params.vaf_context_window_size,
          sample_names=[str(i) for i in range(model_params.num_samples)],
          expand_by_haplotype=model_params.expand_by_haplotype,
      ).model_features
  )
  if model_params.features:
    input_shape = len(model_params.features)
  hidden_layers = model_params.hidden_layer_sizes
  model.add(
      keras.layers.Dense(
          hidden_layers[0],
          activation=model_params.activation,
          input_shape=(input_shape,),
      )
  )
  if len(hidden_layers) > 1:
    for layer_size in hidden_layers[1:]:
      model.add(
          keras.layers.Dense(layer_size, activation=model_params.activation)
      )
  output_shape = len(make_small_model_examples.GenotypeEncoding)
  model.add(keras.layers.Dense(output_shape, activation="softmax"))

  optimizer = keras.optimizers.Adam(
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


def load_keras_model(
    checkpoint_path: str, compile_model: bool = True
) -> keras.Model:
  """Loads a Keras model from the given checkpoint path."""
  return keras.models.load_model(checkpoint_path, compile=compile_model)


def get_keras_training_callbacks(
    checkpoint_filepath: str,
    is_xmanager_run: bool,
    tensorboard_log_dir: str,
    logging_frequency: int,
    batch_size: int,
    num_train_samples: int,
) -> list[keras.callbacks.Callback]:
  """If provisioned, provides a callback to log metrics to Xmanager/TB."""
  metric_loggers = []
  checkpoint_path_template = (
      "%s_epoch:{epoch:02d}_val_loss:{val_loss:.5f}.keras" % checkpoint_filepath
  )
  model_checkpoint_callback = keras.callbacks.ModelCheckpoint(
      filepath=checkpoint_path_template,
      monitor="val_accuracy",
      mode="max",
      save_freq="epoch",
      save_best_only=False,
  )
  callbacks = [model_checkpoint_callback]

  return callbacks


