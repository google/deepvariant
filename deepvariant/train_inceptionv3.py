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
r"""Train a TF2/Keras Inception V3 model.
"""

import gc
import os
import re
from typing import Any, Dict, Optional, Tuple
import warnings

from absl import app
from absl import flags
from absl import logging
from ml_collections.config_flags import config_flags
import tensorflow as tf
import tensorflow_addons as tfa

from deepvariant import data_providers
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import keras_modeling as modeling
from deepvariant.average_model_checkpoint_patched import AverageModelCheckpointPatched
from third_party.nucleus.io import sharded_file_utils
from official.modeling import optimization


# This gets set by XManager even though we don't pass it as a coordinator_flag
_TPU = flags.DEFINE_string('tpu', None, 'Name of the TPU to use.')

_EXPERIMENT_DIR = flags.DEFINE_string(
    'experiment_dir',
    None,
    (
        'The directory where the model weights, training/evaluation summaries, '
        'and backup information are stored.'
    ),
)

_SAVE_BEST_ONLY = flags.DEFINE_bool(
    'save_best_only', True, 'If False, save all ckpts in experiment_dir.'
)

config_flags.DEFINE_config_file('train_config', None)

FLAGS = flags.FLAGS

_DEFAULT_INPUT_READ_THREADS = 32
_DEFAULT_SHUFFLE_BUFFER_ELEMENTS = 100
_DEFAULT_INITIAL_SHUFFLE_BUFFER_ELEMENTS = 1024
_DEFAULT_PREFETCH_BUFFER_BYTES = 16 * 1000 * 1000


class XManagerCallback(tf.keras.callbacks.Callback):
  """Registers metrics with XManager."""

  def __init__(self):
    super(tf.keras.callbacks.Callback, self).__init__()
    self.xm_client = xmanager_api.XManagerApi()
    self.work_unit = self.xm_client.get_current_work_unit()

  def on_epoch_end(self, epoch: int, logs: Optional[Dict[str, Any]] = None):
    if logs:
      for metric_name, metric_value in logs.items():
        measurements = self.work_unit.get_measurement_series(label=metric_name)
        measurements.create_measurement(metric_value, step=epoch)


class F1ScorePerClass(tfa.metrics.F1Score):
  """Reports F1 Score for a target class."""

  # TODO: Create custom metrics.py module.

  def __init__(self, num_classes: int, target_class: int, name: str):
    self.target_class = target_class
    super().__init__(num_classes=num_classes, name=name)

  def result(self) -> tf.Tensor:
    return super().result()[self.target_class]


def parse_example(
    example: tf.train.Example, input_shape: Tuple[int, int, int]
) -> Tuple[tf.Tensor, tf.Tensor, tf.Tensor]:
  """Parses a serialized tf.Example, preprocesses the image, and one-hot encodes the label."""
  proto_features = {
      'image/encoded': tf.io.FixedLenFeature((), tf.string),
      'variant/encoded': tf.io.FixedLenFeature((), tf.string),
      'alt_allele_indices/encoded': tf.io.FixedLenFeature((), tf.string),
      'label': tf.io.FixedLenFeature((1), tf.int64),
  }
  proto_features_denovo = {
      'denovo_label': tf.io.FixedLenFeature((1), tf.int64),
  } | proto_features

  if FLAGS.train_config.denovo_enabled:
    parsed_features = tf.io.parse_single_example(
        serialized=example, features=proto_features_denovo
    )
  else:
    parsed_features = tf.io.parse_single_example(
        serialized=example, features=proto_features
    )
  image = tf.io.decode_raw(parsed_features['image/encoded'], tf.uint8)
  image = tf.reshape(image, input_shape)
  image = tf.cast(image, tf.float32)

  # Image preprocessing and one-hot encoding were previously done inside the
  # TF Estimator API's model_fn. Though we can subclass the Keras InceptionV3
  # class and do it in the forward pass, it seems more fitting to do it during
  # dataset loading along with the above image loading steps.
  image = dv_utils.preprocess_images(image)
  label = tf.keras.layers.CategoryEncoding(
      num_tokens=dv_constants.NUM_CLASSES, output_mode='one_hot'
  )(parsed_features['label'])

  if (
      FLAGS.train_config.denovo_enabled
      and parsed_features['denovo_label'][0] == 1
  ):
    # If the example is denovo then set the denovo example weights for this.
    sample_weight = tf.constant(
        FLAGS.train_config.denovo_weight, dtype=tf.float32
    )
  else:
    sample_weight = tf.constant(1.0, dtype=tf.float32)

  return image, label, sample_weight


def load_dataset(filename: str) -> tf.data.Dataset:
  return tf.data.TFRecordDataset(
      filename,
      buffer_size=_DEFAULT_PREFETCH_BUFFER_BYTES,
      compression_type='GZIP',
  )


def input_fn(
    path: str,
    input_shape: Tuple[int, int, int],
    max_examples: Optional[int] = None,
    mode: str = 'train',
) -> Optional[tf.data.Dataset]:
  """Load the dataset and shuffle for training."""
  assert mode in ['train', 'validation']
  is_training = mode == 'train'
  dataset = None
  for pattern in path.split(','):
    one_dataset = tf.data.TFRecordDataset.list_files(
        sharded_file_utils.normalize_to_sharded_file_pattern(pattern),
        shuffle=is_training,
    )
    dataset = dataset.concatenate(one_dataset) if dataset else one_dataset

  if dataset is None:
    return None

  if is_training:
    dataset = dataset.shuffle(
        _DEFAULT_INITIAL_SHUFFLE_BUFFER_ELEMENTS, reshuffle_each_iteration=True
    )

  dataset = dataset.interleave(
      load_dataset,
      cycle_length=_DEFAULT_INPUT_READ_THREADS,
      num_parallel_calls=tf.data.AUTOTUNE,
      deterministic=False,
  )

  if max_examples is not None:
    dataset = dataset.take(max_examples)

  dataset = dataset.repeat()

  if is_training:
    dataset = dataset.shuffle(
        _DEFAULT_INITIAL_SHUFFLE_BUFFER_ELEMENTS, reshuffle_each_iteration=True
    )
  # Best practices suggest batching before mapping, but this fails so I swapped
  # the order.
  dataset = dataset.map(
      map_func=lambda example: parse_example(example, input_shape),
      num_parallel_calls=tf.data.AUTOTUNE,
      deterministic=False,
  )
  batch_size = (
      FLAGS.train_config.batch_size_train
      if mode == 'train'
      else FLAGS.train_config.batch_size_eval
  )
  dataset = dataset.batch(batch_size=batch_size, drop_remainder=True)

  # Prefetch overlaps in-feed with training
  dataset = dataset.prefetch(tf.data.AUTOTUNE)
  return dataset


def construct_hyperparam_str() -> str:
  return (
      f'mepochs{FLAGS.train_config.num_mini_epochs_per_epoch}'
      f'bstrain{FLAGS.train_config.batch_size_train}'
      f'bseval{FLAGS.train_config.batch_size_eval}'
      f'lr{FLAGS.train_config.learning_rate}'
      f'lrnepochsperdecay{FLAGS.train_config.learning_rate_num_epochs_per_decay}'
      f'lrdecayrate{FLAGS.train_config.learning_rate_decay_rate}'
      f'ad{FLAGS.train_config.average_decay}'
      f'ls{FLAGS.train_config.label_smoothing}'
      f'rho{FLAGS.train_config.rho}'
      f'mom{FLAGS.train_config.momentum}'
      f'eps{FLAGS.train_config.epsilon}'
      f'warmupmepochs{FLAGS.train_config.warmup_mini_epochs}'
  )


class ClearSessionCallback(tf.keras.callbacks.Callback):

  def on_epoch_end(self, epoch: int, logs: Optional[Dict[str, Any]] = None):
    gc.collect()


def train():
  """Train a model."""
  if FLAGS.train_config is None:
    raise ValueError('train_config must be specified')

  experiment_dir = _EXPERIMENT_DIR.value
  experiment_dir = os.path.join(experiment_dir, construct_hyperparam_str())
  experiment_model_dir = f'{experiment_dir}/checkpoints'
  experiment_tensorboard_dir = f'{experiment_dir}/tensorboard_log'
  experiment_backup_dir = f'{experiment_dir}/backup'
  logging.info(
      'Use TPU at %s', _TPU.value if _TPU.value is not None else 'local'
  )
  logging.info('experiment_dir: %s', experiment_dir)
  tpu = _TPU.value
  if tpu is not None:
    resolver = tf.distribute.cluster_resolver.TPUClusterResolver(tpu=tpu)
    tf.config.experimental_connect_to_cluster(resolver)
    tf.tpu.experimental.initialize_tpu_system(resolver)
    strategy = tf.distribute.TPUStrategy(resolver)
  else:
    strategy = tf.distribute.MirroredStrategy()

  train_dataset_config = data_providers.read_dataset_config(
      FLAGS.train_config.train_dataset_pbtxt
  )

  eval_dataset_config = data_providers.read_dataset_config(
      FLAGS.train_config.tune_dataset_pbtxt
  )
  input_shape = dv_utils.get_shape_from_examples_path(
      train_dataset_config.tfrecord_path
  )
  # Set the total number epochs. This is the number of times we will iterate
  # over the data.
  epochs = FLAGS.train_config.num_epochs
  # Calculate the total number of mini_epochs, we use it to report metrics
  # sooner. This divides a long waited epoch into smaller epochs. This is the
  # value fed into model.fit so it is assumed to be one epoch when it reaches
  # a mini epoch.
  mini_epochs = epochs * FLAGS.train_config.num_mini_epochs_per_epoch
  # Steps per epoch is the total number of steps/batches we have in each epoch.
  steps_per_epoch = (
      train_dataset_config.num_examples // FLAGS.train_config.batch_size_train
  )
  # Steps per mini-epoch is after how many batches we report a metric. This is
  # what used in model.fit so it reports metrics after steps_per_mini_epoch
  # numbers.
  steps_per_mini_epoch = (
      steps_per_epoch // FLAGS.train_config.num_mini_epochs_per_epoch
  )
  if (
      FLAGS.train_config.denovo_enabled
      and FLAGS.train_config.denovo_weight == 1.0
  ):
    logging.warning('Denovo is enabled but denovo weight is set to 1.0')

  logging.info(
      (
          '1) Total examples for training: %s.\n'
          '2) Train batch size or examples per step: %s.\n'
          '3) Total epochs: %s.\n'
          '4) Total steps per epoch: %s.\n'
          '5) Number of mini-epochs per epoch: %s.\n'
          '6) Total steps per mini-epoch: %s.\n'
          '7) Total mini-epochs to training completion (x-axis in xm): %s.\n'
          '8) Denovo weight: %s.\n'
          '9) Denovo enabled: %s. \b'
      ),
      train_dataset_config.num_examples,
      FLAGS.train_config.batch_size_train,
      epochs,
      steps_per_epoch,
      FLAGS.train_config.num_mini_epochs_per_epoch,
      steps_per_mini_epoch,
      epochs * FLAGS.train_config.num_mini_epochs_per_epoch,
      FLAGS.train_config.denovo_weight,
      FLAGS.train_config.denovo_enabled,
  )

  with strategy.scope():
    init_weights_file = None
    if FLAGS.train_config.init_weights_file:
      init_weights_file = FLAGS.train_config.init_weights_file
    model = modeling.inceptionv3(input_shape, init_weights_file)

    decay_steps = int(
        steps_per_epoch * FLAGS.train_config.learning_rate_num_epochs_per_decay
    )

    learning_rate = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=FLAGS.train_config.learning_rate,
        decay_steps=decay_steps,
        decay_rate=FLAGS.train_config.learning_rate_decay_rate,
        staircase=True,
    )

    if FLAGS.train_config.warmup_mini_epochs > 0:
      warmup_steps = int(
          steps_per_mini_epoch * FLAGS.train_config.warmup_mini_epochs
      )
      warmup_learning_rate = FLAGS.train_config.learning_rate / 10
      logging.info(
          (
              'Use LinearWarmup. warmup_mini_epochs=%s, warmup_steps=%s, '
              'warmup_learning_rate=%s'
          ),
          FLAGS.train_config.warmup_mini_epochs,
          warmup_steps,
          warmup_learning_rate,
      )
      learning_rate = optimization.LinearWarmup(
          # This is initial learning rate.
          warmup_learning_rate=warmup_learning_rate,
          after_warmup_lr_sched=learning_rate,
          warmup_steps=warmup_steps,
      )

    # Our "steps_per_epoch" in model.fit is:
    #   steps_per_epoch // FLAGS.train_config.num_mini_epochs_per_epoch.
    # I divided that by 10 just to be safe. In case that gets too small,
    # I set a max with 128 here.
    # Note that this might only make sense for TPUs. In the future we'll want to
    # check whether this works for GPU or not.
    # We should also keep an eye on this feature which might allow autotune in
    # the future: https://github.com/keras-team/keras/issues/16573.
    steps_per_execution = max(
        128,
        steps_per_epoch // FLAGS.train_config.num_mini_epochs_per_epoch // 10,
    )
    logging.info('steps_per_execution=%s', steps_per_execution)
    model.compile(
        # TODO: Try different optimizers
        # and sweep over these hyperparams
        optimizer=tfa.optimizers.MovingAverage(
            optimizer=tf.keras.optimizers.RMSprop(
                learning_rate=learning_rate,
                rho=FLAGS.train_config.rho,
                momentum=FLAGS.train_config.momentum,
                epsilon=FLAGS.train_config.epsilon,
            ),
            average_decay=FLAGS.train_config.average_decay,
        ),
        # This is from:
        # https://www.tensorflow.org/guide/tpu#train_the_model_using_keras_high-level_apis.
        # Anything between 2 and `steps_per_epoch` could help here.
        steps_per_execution=steps_per_execution,
        # TODO: Sweep over these hyperparams
        # TODO: Does computing fewer TensorBoard metrics speed up
        # training?
        loss=tf.keras.losses.CategoricalCrossentropy(
            label_smoothing=FLAGS.train_config.label_smoothing
        ),
        metrics=[
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
            tfa.metrics.F1Score(
                num_classes=3, average='micro', name='f1_micro'
            ),
            tfa.metrics.F1Score(
                num_classes=3, average='macro', name='f1_macro'
            ),
            F1ScorePerClass(num_classes=3, target_class=0, name='f1_homref'),
            F1ScorePerClass(num_classes=3, target_class=1, name='f1_het'),
            F1ScorePerClass(num_classes=3, target_class=2, name='f1_homalt'),
        ],
        weighted_metrics=[
            tf.keras.metrics.CategoricalAccuracy(
                name='weighted_categorical_accuracy'
            ),
            tf.keras.metrics.CategoricalCrossentropy(
                name='weighted_categorical_crossentropy'
            ),
            tf.keras.metrics.TruePositives(name='weighted_true_positives'),
            tf.keras.metrics.TrueNegatives(name='weighted_true_negatives'),
            tf.keras.metrics.FalsePositives(name='weighted_false_positives'),
            tf.keras.metrics.FalseNegatives(name='weighted_false_negatives'),
            tf.keras.metrics.Precision(name='weighted_precision'),
        ],
    )

  tensorboard_callback = tf.keras.callbacks.TensorBoard(
      log_dir=experiment_tensorboard_dir,
      write_images=False,
      profile_batch=0,
  )

  checkpoint_file_avg = os.path.join(
      experiment_model_dir,
      'weights-{epoch:02d}-{%s:f}.ckpt' % FLAGS.train_config.best_metrics,
  )
  logging.info(
      'Metrics used to choose best ckpt is: %s', FLAGS.train_config.best_metrics
  )
  logging.info('Best ckpt file path: %s', checkpoint_file_avg)
  logging.info('save_best_only: %s', _SAVE_BEST_ONLY.value)

  checkpoint_path = os.path.join(experiment_model_dir, 'checkpoint')
  best_metrics_from_checkpoint = None
  weights_path_from_checkpoint = None
  if tf.io.gfile.exists(checkpoint_path):
    with tf.io.gfile.GFile(checkpoint_path, 'r') as ckpt_file:
      first_line = ckpt_file.readline()
      match = re.search(r'model_checkpoint_path: "(.*)"', first_line)
      if match is None:
        raise ValueError('Unexpected checkpoint format.')
      weights_path_from_checkpoint = match.group(1)
      match = re.search(
          r'weights-\d+-(0.\d+).ckpt', weights_path_from_checkpoint
      )
      if match is None:
        raise ValueError('Unexpected checkpoint format.')
      match = re.search(
          r'model_checkpoint_path: "weights-\d+-(0.\d+).ckpt"', first_line
      )
      if match is None:
        raise ValueError('Unexpected checkpoint format.')
      best_metrics_from_checkpoint = float(match.group(1))
      logging.info(
          'Resuming: Loading best metrics from %s: %s',
          checkpoint_path,
          best_metrics_from_checkpoint,
      )
  average_model_checkpoint_callback = AverageModelCheckpointPatched(
      filepath=str(checkpoint_file_avg),
      update_weights=True,
      save_best_only=_SAVE_BEST_ONLY.value,
      save_weights_only=True,
      monitor=FLAGS.train_config.best_metrics,
      mode='max',
      save_freq='epoch',
      best_metrics_from_checkpoint=best_metrics_from_checkpoint,
  )
  if weights_path_from_checkpoint is not None:
    weights_path_from_checkpoint = os.path.join(
        experiment_model_dir, weights_path_from_checkpoint
    )
    logging.info(
        'Resuming: load_weights from from %s', weights_path_from_checkpoint
    )
    model.load_weights(weights_path_from_checkpoint)

  backup_callback = tf.keras.callbacks.BackupAndRestore(experiment_backup_dir)

  # XManager callback.
  xmanager_callback = XManagerCallback()

  training_callbacks = [
      tensorboard_callback,
      average_model_checkpoint_callback,
      backup_callback,
      xmanager_callback,
      ClearSessionCallback(),
      tf.keras.callbacks.EarlyStopping(
          monitor=FLAGS.train_config.best_metrics,
          # `patient` below might need to be tuned.
          patience=FLAGS.train_config.early_stopping.patience
          * FLAGS.train_config.num_mini_epochs_per_epoch,
          mode='max',
          verbose=1,
      ),
  ]
  # Calculate validation attributes.
  num_validation_examples = eval_dataset_config.num_examples
  if FLAGS.train_config.num_validation_examples:
    logging.info(
        (
            'train_config.num_validation_examples is set. Using '
            'min(%s, %s) as num_validation_examples.'
        ),
        FLAGS.train_config.num_validation_examples,
        eval_dataset_config.num_examples,
    )
    num_validation_examples = min(
        FLAGS.train_config.num_validation_examples,
        eval_dataset_config.num_examples,
    )
  logging.info('num_validation_examples = %s', num_validation_examples)
  validation_steps = (
      num_validation_examples // FLAGS.train_config.batch_size_eval
  )
  logging.info('validation_steps: %s', validation_steps)
  model.fit(
      input_fn(train_dataset_config.tfrecord_path, input_shape, mode='train'),
      epochs=mini_epochs,
      verbose=2,
      steps_per_epoch=steps_per_mini_epoch,
      callbacks=training_callbacks,
      validation_data=input_fn(
          eval_dataset_config.tfrecord_path,
          input_shape,
          max_examples=num_validation_examples,
          mode='validation',
      ),
      # Even with "mini epoch", we still want to evaluate the same amount of
      # validation examples per point. So, I'm not dividing this by
      # FLAGS.train_config.num_mini_epochs_per_epoch
      validation_steps=validation_steps,
  )


def main(unused_argv):
  keep_running = True
  while keep_running:
    try:
      train()
      keep_running = False  # Training completed successfully.
    except tf.errors.UnavailableError as error:
      logging.warning(
          'UnavailableError encountered during training: %s.', error
      )


if __name__ == '__main__':
  logging.set_verbosity(logging.INFO)
  warnings.filterwarnings(
      'ignore', module='tensorflow_addons.optimizers.average_wrapper'
  )
  app.run(main)
