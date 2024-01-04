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

import math
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

from deepvariant import data_providers
from deepvariant import dv_utils
from deepvariant import keras_modeling
from official.modeling import optimization


_LEADER = flags.DEFINE_string(
    'leader',
    'local',
    (
        'The leader flag specifies the host-controller. Possible values: '
        '(1) local=runs locally. If GPUs are available they will be detected'
        ' and used.'
    ),
)

_STRATEGY = flags.DEFINE_enum(
    'strategy',
    'mirrored',
    ['tpu', 'mirrored'],
    'The strategy to use.',
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
    'limit', None, 'Limit the number of steps used for train/eval.'
)

_DEBUG = flags.DEFINE_bool(
    'debug', False, 'Run tensorflow eagerly in debug mode.'
)

config_flags.DEFINE_config_file('config', None)

FLAGS = flags.FLAGS


def roundup(n: int, nearest: int) -> int:
  """Rounds up an integer to the nearest value."""
  return int(math.ceil(n / float(nearest)) * nearest)


def train(config: ml_collections.ConfigDict):
  """Train a model."""
  logging.info('Running with debug=%s', _DEBUG.value)
  tf.config.run_functions_eagerly(_DEBUG.value)
  if _DEBUG.value:
    tf.data.experimental.enable_debug_mode()

  experiment_dir = _EXPERIMENT_DIR.value

  model_dir = f'{experiment_dir}/checkpoints'
  logging.info(
      'Use TPU at %s', _LEADER.value if _LEADER.value is not None else 'local'
  )
  logging.info('experiment_dir: %s', experiment_dir)
  if _STRATEGY.value == 'tpu':
    resolver = tf.distribute.cluster_resolver.TPUClusterResolver(
        tpu=_LEADER.value
    )
    tf.config.experimental_connect_to_cluster(resolver, protocol='grpc+loas')
    tf.tpu.experimental.initialize_tpu_system(resolver)
    strategy = tf.distribute.TPUStrategy(resolver)
  elif _STRATEGY.value in ['mirrored']:
    strategy = tf.distribute.MirroredStrategy()
  else:
    raise ValueError(f'Unknown strategy: {_STRATEGY.value}')

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

  # Copy example_info.json to checkpoint path.
  example_info_json_path = os.path.join(
      os.path.dirname(train_dataset_config.tfrecord_path), 'example_info.json'
  )
  if not tf.io.gfile.exists(example_info_json_path):
    example_info_json_path = (
        train_dataset_config.tfrecord_path + '.example_info.json'
    )
  if not tf.io.gfile.exists(example_info_json_path):
    raise FileNotFoundError(
        'example_info.json not found in directory'
        f' {os.path.dirname(train_dataset_config.tfrecord_path)}'
    )
  tf.io.gfile.makedirs(experiment_dir)
  tf.io.gfile.copy(
      example_info_json_path,
      os.path.join(experiment_dir, 'example_info.json'),
      overwrite=True,
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
        weights=config.init_checkpoint,
        init_backbone_with_imagenet=config.init_backbone_with_imagenet,
        config=config,
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

    logging.info(
        'Exponential Decay:'
        ' initial_learning_rate=%s\n'
        ' decay_steps=%s\n'
        ' learning_rate_decay_rate=%s',
        config.learning_rate,
        decay_steps,
        config.learning_rate_decay_rate,
    )

    if config.warmup_steps > 0:
      warmup_learning_rate = config.learning_rate / 10
      logging.info(
          'Use LinearWarmup: \n warmup_steps=%s\n warmup_learning_rate=%s',
          config.warmup_steps,
          warmup_learning_rate,
      )
      learning_rate = optimization.LinearWarmup(
          # This is initial learning rate.
          warmup_learning_rate=warmup_learning_rate,
          after_warmup_lr_sched=learning_rate,
          warmup_steps=config.warmup_steps,
      )

    # Define Optimizer.
    # TODO: Add function for retrieving custom optimizer.
    if config.optimizer == 'nadam':
      optimizer = tf.keras.optimizers.Nadam(learning_rate=learning_rate)
    elif config.optimizer == 'rmsprop':
      optimizer = tf.keras.optimizers.RMSprop(
          learning_rate=learning_rate,
          rho=config.rho,
          momentum=config.momentum,
          epsilon=config.epsilon,
          use_ema=config.use_ema,
          weight_decay=config.optimizer_weight_decay,
      )
    else:
      raise ValueError(f'Unknown optimizer: {config.optimizer}')

  # ================= #
  # Setup Checkpoint  #
  # ================= #

  # The state object allows checkpointing of the model and associated variables
  # for the optimizer, step, and train/tune metrics.
  ckpt_manager = keras_modeling.create_state(
      config,
      model_dir,
      model,
      optimizer,
      strategy,
  )
  state = ckpt_manager.checkpoint

  @tf.function
  def run_train_step(inputs):
    model_inputs, labels = inputs

    with tf.GradientTape() as tape:
      logits = model(model_inputs, training=True)
      probabilities = tf.nn.softmax(logits)
      train_loss = compute_loss(probabilities=probabilities, labels=labels)

    gradients = tape.gradient(train_loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))

    for metric in state.train_metrics[:-1]:
      metric.update_state(
          y_pred=probabilities,
          y_true=labels,
      )
    state.train_metrics[-1].update_state(train_loss)
    return train_loss

  @tf.function
  def run_tune_step(tune_inputs):
    """Single non-distributed tune step."""
    model_inputs, labels = tune_inputs
    logits = model(model_inputs, training=False)
    probabilities = tf.nn.softmax(logits)
    tune_loss = compute_loss(probabilities=probabilities, labels=labels)

    for metric in state.tune_metrics[:-1]:
      metric.update_state(
          y_pred=probabilities,
          y_true=labels,
      )
      state.tune_metrics[-1].update_state(tune_loss)
    return tune_loss

  @tf.function
  def distributed_train_step(iterator, num_steps):
    for _ in tf.range(num_steps):
      strategy.run(run_train_step, args=(next(iterator),))
    state.global_step.assign_add(num_steps)

  @tf.function
  def distributed_tune_step(iterator, num_steps):
    for _ in tf.range(num_steps):
      strategy.run(run_tune_step, args=(next(iterator),))

  # ============== #
  # Setup Datasets #
  # ============== #
  train_ds = data_providers.input_fn(
      train_dataset_config.tfrecord_path,
      mode='train',
      strategy=strategy,
      config=config,
  )
  tune_ds = data_providers.input_fn(
      tune_dataset_config.tfrecord_path,
      mode='tune',
      strategy=strategy,
      config=config,
  )
  train_ds, tune_ds = iter(train_ds), iter(tune_ds)
  num_train_steps = steps_per_epoch * config.num_epochs

  logging.info(
      (
          '\n\n'
          'Training Examples: %s\n'
          'Tune Examples: %s\n'
          'Batch Size: %s\n'
          'Epochs: %s\n'
          'Steps per epoch: %s\n'
          'Steps per tune: %s\n'
          'Num train steps: %s\n'
          'Steps per iter: %s\n'
          '\n'
      ),
      train_dataset_config.num_examples,
      tune_dataset_config.num_examples,
      config.batch_size,
      config.num_epochs,
      steps_per_epoch,
      steps_per_tune,
      num_train_steps,
      config.steps_per_iter,
  )

  if config.log_every_steps % config.steps_per_iter != 0:
    logging.warning(
        'log_every_steps should be a multiple of steps_per_iter. '
        'log_every_steps=%s, steps_per_iter=%s',
        config.log_every_steps,
        config.steps_per_iter,
    )
  if config.tune_every_steps % config.steps_per_iter != 0:
    logging.warning(
        'tune_every_steps should be a multiple of steps_per_iter. '
        'tune_every_steps=%s, steps_per_iter=%s',
        config.tune_every_steps,
        config.steps_per_iter,
    )

  # ============= #
  # Training Loop #
  # ============= #

  metric_writer = metric_writers.create_default_writer(logdir=experiment_dir)
  report_progress = periodic_actions.ReportProgress(
      num_train_steps=num_train_steps,
      writer=metric_writer,
      every_secs=300,
      on_steps=[0, num_train_steps - 1],
  )

  with strategy.scope():

    def get_checkpoint_metric():
      """Returns the metric we are optimizing for."""
      best_checkpoint_metric_idx = [
          f'tune/{x.name}' for x in state.tune_metrics
      ].index(config.best_checkpoint_metric)
      return state.tune_metrics[best_checkpoint_metric_idx].result().numpy()

    best_checkpoint_metric_value = get_checkpoint_metric()

    with metric_writers.ensure_flushes(metric_writer):

      def run_tune(train_step, epoch, steps_per_tune):
        logging.info('Running tune at step=%d epoch=%d', train_step, epoch)
        for tune_step in range(0, steps_per_tune, config.steps_per_iter):
          with tf.profiler.experimental.Trace('tune', step_num=tune_step, _r=1):
            if roundup(tune_step, config.steps_per_iter) > steps_per_tune:
              steps_per_iter = steps_per_tune % config.steps_per_iter
            else:
              steps_per_iter = config.steps_per_iter
            if (
                tune_step
                % roundup(config.log_every_steps, config.steps_per_iter)
                == 0
            ):
              logging.info(
                  'Tune step %s / %s (%s%%)',
                  tune_step,
                  steps_per_tune,
                  round(float(tune_step) / float(steps_per_tune), 1) * 100.0,
              )
            distributed_tune_step(tune_ds, steps_per_iter)

        metric_writer.write_scalars(
            train_step,
            {f'tune/{x.name}': x.result() for x in state.tune_metrics},
        )

      for train_step in range(
          state.global_step.numpy(), num_train_steps, config.steps_per_iter
      ):
        # Calculate current epoch
        epoch = train_step // steps_per_epoch
        if train_step % roundup(steps_per_epoch, config.steps_per_iter) == 0:
          logging.info('Starting epoch %s', epoch)

        # If we are warmstarting, establish an initial best_checkpoint_metric
        # value before beginning any training.
        if train_step == 0 and (
            config.init_checkpoint or config.init_backbone_with_imagenet
        ):
          logging.info('Performing initial evaluation of warmstart model.')
          run_tune(train_step, epoch, steps_per_tune)
          best_checkpoint_metric_value = get_checkpoint_metric()
          logging.info(
              'Warmstart checkpoint best checkpoint metric: %s=%s',
              config.best_checkpoint_metric,
              best_checkpoint_metric_value,
          )
          # Reset tune metrics
          for metric in state.tune_metrics:
            metric.reset_states()

        # ===== #
        # train #
        # ===== #
        # Calculate full train step.
        is_last_step = (
            roundup(train_step, config.steps_per_iter) >= num_train_steps
        )
        if is_last_step:
          n_steps_iter = train_step % config.steps_per_iter
          logging.info(
              'Running last training step. n_steps remaining=%d', n_steps_iter
          )
        else:
          n_steps_iter = config.steps_per_iter

        with tf.profiler.experimental.Trace('train', step_num=train_step, _r=1):
          distributed_train_step(train_ds, n_steps_iter)

          # Quick indication that training is happening.
          logging.log_first_n(
              logging.INFO, 'Finished training step %d.', 5, train_step
          )

        # Log metrics
        report_progress(train_step)

        if (
            train_step % roundup(config.log_every_steps, config.steps_per_iter)
            == 0
        ) or is_last_step:
          metrics_to_write = {
              f'train/{x.name}': x.result() for x in state.train_metrics
          }
          if isinstance(
              optimizer.learning_rate, tf.distribute.DistributedValues
          ):
            current_learning_rate = optimizer.learning_rate.numpy()
          else:
            current_learning_rate = optimizer.learning_rate(train_step)
          metrics_to_write['train/learning_rate'] = current_learning_rate
          metrics_to_write['epoch'] = epoch
          metric_writer.write_scalars(
              train_step,
              metrics_to_write,
          )
          # Reset train metrics.
          for metric in state.train_metrics:
            metric.reset_states()

        # ==== #
        # tune #
        # ==== #
        # Run tune at every epoch, periodically, and at final step.
        if (
            (
                train_step > 0
                and train_step % roundup(steps_per_epoch, config.steps_per_iter)
                == 0
            )
            or (
                train_step > 0
                and train_step
                % roundup(config.tune_every_steps, config.steps_per_iter)
                == 0
            )
            or is_last_step
        ):
          run_tune(train_step, epoch, steps_per_tune)

          if get_checkpoint_metric() > best_checkpoint_metric_value:
            best_checkpoint_metric_value = get_checkpoint_metric()
            checkpoint_path = ckpt_manager.save(train_step)
            # Reset early stopping counter
            state.early_stopping.assign(0)
            logging.info(
                'Saved checkpoint %s=%s step=%s epoch=%s path=%s',
                config.best_checkpoint_metric,
                get_checkpoint_metric(),
                train_step,
                epoch,
                checkpoint_path,
            )
          else:
            if (
                config.early_stopping_patience
                and state.early_stopping.value()
                >= config.early_stopping_patience
            ):
              break
            logging.info(
                'Skipping checkpoint with %s=%s < previous best %s=%s',
                config.best_checkpoint_metric,
                get_checkpoint_metric(),
                config.best_checkpoint_metric,
                best_checkpoint_metric_value,
            )
            state.early_stopping.assign_add(1)
          if config.early_stopping_patience:
            metric_writer.write_scalars(
                train_step,
                {'tune/early_stopping': state.early_stopping.value()},
            )

          # Reset tune metrics
          for metric in state.tune_metrics:
            metric.reset_states()

    # After training completes, load the latest checkpoint and create
    # a saved model (.pb) and keras model formats.
    checkpoint_path = ckpt_manager.latest_checkpoint
    if not checkpoint_path:
      logging.info('No checkpoint found.')
      return

    # The latest checkpoint will be the best performing checkpoint.
    logging.info('Loading best checkpoint: %s', checkpoint_path)
    tf.train.Checkpoint(model).restore(checkpoint_path).expect_partial()

    logging.info('Saving model using saved_model format.')
    saved_model_dir = checkpoint_path
    model.save(saved_model_dir, save_format='tf')
    # Copy example_info.json to saved_model directory.
    tf.io.gfile.copy(
        example_info_json_path,
        os.path.join(
            saved_model_dir,
            'example_info.json',
        ),
        overwrite=True,
    )


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
