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
r"""Train a DeepVariant Keras Model.
"""

import os
import sys
import warnings

from absl import app
from absl import flags
from absl import logging
from clu import metric_writers
from clu import periodic_actions
import ml_collections
from ml_collections.config_flags import config_flags
import tensorflow as tf
import tensorflow_addons as tfa

from deepvariant import data_providers
from deepvariant import dv_utils
from deepvariant import keras_modeling
from tensorflow_models.official.modeling import optimization


_BEST_CHECKPOINT_METRIC = 'tune/f1_weighted'

_LEADER = flags.DEFINE_string(
    'leader',
    'local',
    (
        'The leader flag specifies the host-controller. Possible values: '
        '(1) local=runs locally. If GPUs are available they will be detected'
        ' and used.'
    ),
)

_EXPERIMENT_DIR = flags.DEFINE_string(
    'experiment_dir',
    None,
    (
        'The directory where the model weights, training/tuning summaries, '
        'and backup information are stored.'
    ),
)

_LIMIT = flags.DEFINE_integer(
    'limit', None, 'Limit the number of steps for testing.'
)

_DEBUG = flags.DEFINE_bool(
    'debug', False, 'Run tensorflow eagerly in debug mode.'
)

config_flags.DEFINE_config_file('config', None)

FLAGS = flags.FLAGS


def construct_hyperparam_str() -> str:
  return (
      f'bs{FLAGS.config.batch_size}'
      f'lr{FLAGS.config.learning_rate}'
      f'lrnepochsperdecay{FLAGS.config.learning_rate_num_epochs_per_decay}'
      f'lrdecayrate{FLAGS.config.learning_rate_decay_rate}'
      f'ad{FLAGS.config.average_decay}'
      f'ls{FLAGS.config.label_smoothing}'
      f'rho{FLAGS.config.rho}'
      f'mom{FLAGS.config.momentum}'
      f'eps{FLAGS.config.epsilon}'
  )


class F1Class(tfa.metrics.F1Score):
  """Reports F1 Score for a target class."""

  # TODO: Create custom metrics.py module.

  def __init__(self, num_classes: int, target_class: int, name: str):
    self.target_class = target_class
    super().__init__(num_classes=num_classes, name=name)

  def result(self) -> tf.Tensor:
    return super().result()[self.target_class]


def train(config: ml_collections.ConfigDict):
  """Train a model."""
  logging.info('Running with debug=%s', _DEBUG.value)
  tf.config.run_functions_eagerly(_DEBUG.value)
  if _DEBUG.value:
    tf.data.experimental.enable_debug_mode()

  experiment_dir = _EXPERIMENT_DIR.value
  experiment_dir = os.path.join(experiment_dir, construct_hyperparam_str())
  model_dir = f'{experiment_dir}/checkpoints'
  experiment_tensorboard_dir = f'{experiment_dir}/tensorboard_log'
  logging.info(
      'Use TPU at %s', _LEADER.value if _LEADER.value is not None else 'local'
  )
  logging.info('experiment_dir: %s', experiment_dir)
  tpu = _LEADER.value
  if tpu != 'local':
    resolver = tf.distribute.cluster_resolver.TPUClusterResolver(tpu=tpu)
    tf.config.experimental_connect_to_cluster(resolver)
    tf.tpu.experimental.initialize_tpu_system(resolver)
    strategy = tf.distribute.TPUStrategy(resolver)
  else:
    strategy = tf.distribute.MirroredStrategy()

  # Load config
  train_dataset_config = data_providers.read_dataset_config(
      config.train_dataset_pbtxt
  )

  tune_dataset_config = data_providers.read_dataset_config(
      config.tune_dataset_pbtxt
  )

  input_shape = dv_utils.get_shape_from_examples_path(
      train_dataset_config.tfrecord_path
  )

  steps_per_epoch = train_dataset_config.num_examples // config.batch_size
  steps_per_tune = (
      config.num_validation_examples
      or tune_dataset_config.num_examples // config.batch_size
  )

  if _LIMIT.value:
    steps_per_epoch = _LIMIT.value
    steps_per_tune = _LIMIT.value

  # =========== #
  # Setup Model #
  # =========== #

  with strategy.scope():
    model = keras_modeling.inceptionv3(
        input_shape=input_shape,
        weights=config.init_weights_file,
        init_backbone_with_imagenet=config.init_backbone_with_imagenet,
    )

    # Define Loss Function.
    # TODO: Add function for retrieving custom loss fn.
    loss_function = tf.keras.losses.CategoricalCrossentropy(
        label_smoothing=config.label_smoothing,
        reduction=tf.keras.losses.Reduction.NONE,
    )

    def compute_loss(probabilities, labels):
      per_example_loss = loss_function(y_pred=probabilities, y_true=labels)
      # We divide per-replica losses by global batch size and sum this value
      # across all replicas to compute average loss scaled by global batch size.
      return tf.nn.compute_average_loss(
          per_example_loss, global_batch_size=config.batch_size
      )

    decay_steps = int(
        steps_per_epoch * config.learning_rate_num_epochs_per_decay
    )

    # TODO: Define learning rate via config.
    learning_rate = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=config.learning_rate,
        decay_steps=decay_steps,
        decay_rate=config.learning_rate_decay_rate,
        staircase=True,
    )

    if config.warmup_steps > 0:
      warmup_steps = config.warmup_steps
      warmup_learning_rate = config.learning_rate / 10
      logging.info(
          (
              'Use LinearWarmup. warmup_steps=%s, warmup_steps=%s, '
              'warmup_learning_rate=%s'
          ),
          config.warmup_steps,
          warmup_steps,
          warmup_learning_rate,
      )
      learning_rate = optimization.LinearWarmup(
          # This is initial learning rate.
          warmup_learning_rate=warmup_learning_rate,
          after_warmup_lr_sched=learning_rate,
          warmup_steps=warmup_steps,
      )

    # Define Optimizer.
    # TODO: Add function for retrieving custom optimizer.
    optimizer = tfa.optimizers.MovingAverage(
        optimizer=tf.keras.optimizers.RMSprop(
            learning_rate=learning_rate,
            rho=config.rho,
            momentum=config.momentum,
            epsilon=config.epsilon,
        ),
        average_decay=config.average_decay,
    )

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
          F1Class(num_classes=3, target_class=0, name='f1_homref'),
          F1Class(num_classes=3, target_class=1, name='f1_het'),
          F1Class(num_classes=3, target_class=2, name='f1_homalt'),
          # Leave mean loss as the last metric as it is updated differently.
          tf.keras.metrics.Mean(name='loss'),
      ]

    train_metrics = create_metrics()
    tune_metrics = create_metrics()

  @tf.function
  def run_train_step(inputs):
    model_inputs, labels = inputs

    with tf.GradientTape() as tape:
      logits = model(model_inputs)
      probabilities = tf.nn.softmax(logits)
      train_loss = compute_loss(probabilities=probabilities, labels=labels)

    gradients = tape.gradient(train_loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))

    for metric in train_metrics[:-1]:
      metric.update_state(
          y_pred=probabilities,
          y_true=labels,
      )
    train_metrics[-1].update_state(train_loss)
    return train_loss

  @tf.function
  def run_tune_step(tune_inputs):
    """Single non-distributed tune step."""
    model_inputs, labels = tune_inputs
    logits = model(model_inputs)
    probabilities = tf.nn.softmax(logits)
    tune_loss = compute_loss(probabilities=probabilities, labels=labels)

    for metric in tune_metrics[:-1]:
      metric.update_state(
          y_pred=probabilities,
          y_true=labels,
      )
      tune_metrics[-1].update_state(tune_loss)
    return tune_loss

  @tf.function
  def distributed_train_step(iterator):
    per_replica_losses = strategy.run(run_train_step, args=(next(iterator),))
    # state.global_step.assign_add(1)
    return strategy.reduce(
        tf.distribute.ReduceOp.SUM, per_replica_losses, axis=None
    )

  @tf.function
  def distributed_tune_step(iterator):
    per_replica_losses = strategy.run(run_tune_step, args=(next(iterator),))
    return strategy.reduce(
        tf.distribute.ReduceOp.SUM, per_replica_losses, axis=None
    )

  # ============== #
  # Setup Datasets #
  # ============== #
  train_ds = data_providers.input_fn(
      train_dataset_config.tfrecord_path,
      mode='train',
      strategy=strategy,
      n_epochs=config.num_epochs,
      config=config,
      limit=_LIMIT.value,
  )
  tune_ds = data_providers.input_fn(
      tune_dataset_config.tfrecord_path,
      mode='tune',
      strategy=strategy,
      config=config,
      limit=steps_per_tune,
  )
  train_iter, tune_iter = iter(train_ds), iter(tune_ds)

  logging.info(
      (
          '1) Training Examples: %s.\n'
          '2) Batch Size: %s.\n'
          '3) Epochs: %s.\n'
          '4) Steps per epoch: %s.\n'
      ),
      train_dataset_config.num_examples,
      config.batch_size,
      config.num_epochs,
      steps_per_epoch,
  )

  # ================= #
  # Setup Checkpoint  #
  # ================= #

  with strategy.scope():
    checkpoint = tf.train.Checkpoint(
        optimizer=optimizer,
        model=model,
    )
    ckpt_manager = tf.train.CheckpointManager(
        checkpoint=checkpoint,
        directory=model_dir,
        max_to_keep=5,
    )
    starting_epoch = 0
    # TODO: Restore checkpoint using a state object.
    # if ckpt_manager.latest_checkpoint:
    #   ckpt_manager.restore(ckpt_manager.latest_checkpoint)
    #   restored_step = checkpoint.global_step.numpy().item()
    #   starting_epoch = restored_step // steps_per_epoch
    #   # Calculate the number of remaining training steps.
    #   best_tune = checkpoint.tune_metrics.result()[_BEST_CHECKPOINT_METRIC]
    #   logging.info(
    #       'Restored from checkpoint %s at epoch=%s',
    #       ckpt_manager.latest_checkpoint,
    #       starting_epoch,
    #   )
    # else:
    #   # Set n_steps train to full number of training steps.
    #   logging.info('Initializing from scratch.')
    #   starting_epoch = 0
    #   best_tune = 0

  # ============= #
  # Training Loop #
  # ============= #
  metric_writer = metric_writers.create_default_writer(
      logdir=experiment_tensorboard_dir
  )
  num_train_steps = steps_per_epoch * config.num_epochs
  report_progress = periodic_actions.ReportProgress(
      num_train_steps=num_train_steps, writer=metric_writer, every_secs=600
  )

  metric_writer = metric_writers.create_default_writer(
      logdir=experiment_tensorboard_dir
  )

  with strategy.scope():

    def get_checkpoint_metric():
      """Returns the metric we are optimizing for."""
      best_checkpoint_metric_idx = [
          f'tune/{x.name}' for x in tune_metrics
      ].index(_BEST_CHECKPOINT_METRIC)
      return tune_metrics[best_checkpoint_metric_idx].result()

    best_checkpoint_metric = get_checkpoint_metric()

    with metric_writers.ensure_flushes(metric_writer):
      for epoch in range(starting_epoch, config.num_epochs):
        logging.info('Starting epoch %s', epoch)
        for loop_train_step in range(steps_per_epoch):
          # ===== #
          # train #
          # ===== #
          # Calculate full train step.
          train_step = loop_train_step + (epoch * steps_per_epoch)

          is_last_step = (train_step == (steps_per_epoch - 1)) and (
              epoch == (config.num_epochs - 1)
          )

          with tf.profiler.experimental.Trace(
              'train', step_num=train_step, _r=1
          ):
            distributed_train_step(train_iter)

            # Quick indication that training is happening.
            logging.log_first_n(
                logging.INFO, 'Finished training step %d.', 5, train_step
            )

          # Log metrics
          report_progress(train_step)

          if (train_step % config.log_every_steps == 0) or is_last_step:
            metrics_to_write = {
                f'train/{x.name}': x.result() for x in train_metrics
            }
            metrics_to_write['train/learning_rate'] = optimizer.learning_rate(
                train_step
            )
            metrics_to_write['epoch'] = epoch
            metric_writer.write_scalars(
                train_step,
                metrics_to_write,
            )
            # Reset train metrics.
            for metric in train_metrics:
              metric.reset_states()

          # ==== #
          # tune #
          # ==== #
          # Run tune at every epoch, periodically, and at final step.
          if (
              (train_step > 0 and train_step % steps_per_epoch == 0)
              or (train_step > 0 and train_step % config.tune_every_steps == 0)
              or is_last_step
          ):
            logging.info('Running tune at step=%d epoch=%d', train_step, epoch)
            for loop_step_tune in range(steps_per_tune):
              step_tune = loop_step_tune + (epoch * steps_per_tune)
              with tf.profiler.experimental.Trace(
                  'tune', step_num=step_tune, _r=1
              ):
                distributed_tune_step(tune_iter)

            metric_writer.write_scalars(
                train_step,
                {f'tune/{x.name}': x.result() for x in tune_metrics},
            )

            if get_checkpoint_metric() > best_checkpoint_metric:
              best_checkpoint_metric = get_checkpoint_metric()
              checkpoint_path = ckpt_manager.save(train_step)
              logging.info(
                  'Saved checkpoint %s=%s (path=%s)',
                  _BEST_CHECKPOINT_METRIC,
                  get_checkpoint_metric(),
                  checkpoint_path,
              )
            else:
              logging.info(
                  'Skipping checkpoint with %s=%s < previous best %s=%s',
                  _BEST_CHECKPOINT_METRIC,
                  get_checkpoint_metric(),
                  _BEST_CHECKPOINT_METRIC,
                  best_checkpoint_metric,
              )

            # Reset tune metrics
            for metric in tune_metrics:
              metric.reset_states()


def main(unused_argv):
  keep_running = True
  while keep_running:
    try:
      train(FLAGS.config)
      keep_running = False  # Training completed successfully.
    except tf.errors.UnavailableError as error:
      logging.warning(
          'UnavailableError encountered during training: %s.', error
      )
      sys.exit(42)


if __name__ == '__main__':
  logging.set_verbosity(logging.INFO)
  warnings.filterwarnings(
      'ignore', module='tensorflow_addons.optimizers.average_wrapper'
  )
  app.run(main)
