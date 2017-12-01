# Copyright 2017 Google Inc.
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
"""Trains the DeepVariant model."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json
import os



import tensorflow as tf

from absl import logging

from deepvariant import data_providers
from deepvariant import logging_level
from deepvariant import modeling
from deepvariant.core import proto_utils

slim = tf.contrib.slim

FLAGS = tf.flags.FLAGS

# Data set selection parameters
tf.flags.DEFINE_string('dataset_config_pbtxt', None,
                       'The path to the dataset config file.')

tf.flags.DEFINE_string('model_name', 'inception_v3',
                       'The name of the model to use for predictions.')

tf.flags.DEFINE_integer('batch_size', 64,
                        'The number of samples in each batch.')

tf.flags.DEFINE_string('master', '',
                       'The TensorFlow master to use. Set to the empty string '
                       'to let TF pick a default.')

tf.flags.DEFINE_string('train_dir', '/tmp/deepvariant/',
                       'Directory where to write event logs.')

tf.flags.DEFINE_integer('worker_replicas', 1, 'Number of worker replicas.')

tf.flags.DEFINE_integer(
    'ps_tasks', 0,
    'The number of parameter servers. If the value is 0, then the parameters '
    'are handled locally by the worker.')

tf.flags.DEFINE_integer(
    'save_summaries_secs', 600,
    'The frequency with which summaries are saved, in seconds.')

tf.flags.DEFINE_integer(
    'save_interval_secs', 600,
    'The frequency with which the model is saved, in seconds.')

tf.flags.DEFINE_integer('startup_delay_steps', 15,
                        'Number of training steps between replicas startup.')

tf.flags.DEFINE_integer('task', 0,
                        'Task id of the replica running the training.')

tf.flags.DEFINE_integer('number_of_steps', 30000000,
                        'Maximum number of global steps to take when training.')

# Training parameters.
tf.flags.DEFINE_float('learning_rate', 0.001, 'Initial learning rate.')

tf.flags.DEFINE_float('rmsprop_momentum', 0.9, 'Momentum.')

tf.flags.DEFINE_float('rmsprop_decay', 0.9, 'Decay term for RMSProp.')

tf.flags.DEFINE_float('rmsprop_epsilon', 1.0, 'Epsilon term for RMSProp.')

tf.flags.DEFINE_float('learning_rate_decay_factor', 0.94,
                      'Learning rate decay factor.')

tf.flags.DEFINE_float('num_epochs_per_decay', 2.0,
                      'Number of epochs after which learning rate decays.')

tf.flags.DEFINE_integer(
    'replicas_to_aggregate', 1,
    'The Number of gradients to collect before updating params.')

tf.flags.DEFINE_float('moving_average_decay', 0.9999,
                      'The decay to use for the moving average.')

tf.flags.DEFINE_integer(
    'num_retries', 0,
    'The number of times to retry on InternalError or UnavailableError.')

# Pre-trained model parameters
tf.flags.DEFINE_string(
    'start_from_checkpoint', 'model_default',
    'A path to a checkpoint of model weights to initalize our model at the '
    'start of training. If None or "", the model will start from random weights'
    '. The special value "model_default" will use the default pretrained '
    'path for the selected model.')

tf.flags.DEFINE_integer('max_checkpoints_to_keep', 10,
                        'Number of last checkpoints to keep during traning. '
                        'Passing "0" preserves all checkpoints.')

tf.flags.DEFINE_float(
    'keep_checkpoint_every_n_hours', 1.0,
    'If specified, in addition to keeping the last "max_checkpoints_to_keep" '
    'checkpoints, an additional checkpoint will be kept for every n hours of '
    'training.')


def model_init_function(model, num_classes, checkpoint_path):
  """Creates an init_fn for slim.learning.train.

  Args:
    model: DeepVariantModel. The model we want an init_fn for.
    num_classes: int. The number of class labels we want to predict with this
      model.
    checkpoint_path: str or ''/None. A path to a model checkpoint file that we
      will load our model parameters from. If bool(checkpoint_path) == False, we
      not load a checkpoint but rather return None, indicating no initialization
      is needed.

  Returns:
    A init_fn suitable for use with slim.learning.train, or None if
    bool(checkpoint_path) == False.
  """
  # If the special value "model_default" was passed, ask the model for
  # its default.
  if checkpoint_path == 'model_default':
    checkpoint_path = model.pretrained_model_path

  # If the path is non-False, use it.
  if checkpoint_path:
    logging.info('Initializing model from checkpoint at %s', checkpoint_path)
    return model.initialize_from_checkpoint(
        checkpoint_path, num_classes, is_training=True)
  else:
    logging.info('Initializing model with random parameters')
    return None


def run(target, is_chief, device_fn):
  """Run training.

  Args:
     target: The target of the TensorFlow standard server to use. Can be the
       empty string to run locally using an inprocess server.
     is_chief: Boolean indicating whether this process is the chief.
     device_fn: Device function used to assign ops to devices.
  """
  if not FLAGS.dataset_config_pbtxt:
    logging.error('Need to specify --dataset_config_pbtxt')
    return

  g = tf.Graph()
  with g.as_default():
    model = modeling.get_model(FLAGS.model_name)
    dataset = data_providers.get_dataset(FLAGS.dataset_config_pbtxt)
    print('Running training on {} with model {}\n'.format(dataset, model))

    with tf.device(device_fn):
      # If ps_tasks is zero, the local device is used. When using multiple
      # (non-local) replicas, the ReplicaDeviceSetter distributes the variables
      # across the different devices.
      images, labels, _ = data_providers.make_training_batches(
          dataset.get_slim_dataset(), model, FLAGS.batch_size)
      endpoints = model.create(images, dataset.num_classes, is_training=True)
      labels = slim.one_hot_encoding(labels, dataset.num_classes)
      total_loss = model.loss(endpoints, labels)

      # Setup the moving averages:
      moving_average_variables = slim.get_model_variables()
      moving_average_variables.extend(slim.losses.get_losses())
      moving_average_variables.append(total_loss)

      variable_averages = tf.train.ExponentialMovingAverage(
          FLAGS.moving_average_decay, slim.get_or_create_global_step())

      tf.add_to_collection(tf.GraphKeys.UPDATE_OPS,
                           variable_averages.apply(moving_average_variables))

      # Configure the learning rate using an exponetial decay.
      decay_steps = int(((1.0 * dataset.num_examples) / FLAGS.batch_size) *
                        FLAGS.num_epochs_per_decay)

      learning_rate = tf.train.exponential_decay(
          FLAGS.learning_rate,
          slim.get_or_create_global_step(),
          decay_steps,
          FLAGS.learning_rate_decay_factor,
          staircase=True)

      opt = tf.train.RMSPropOptimizer(learning_rate, FLAGS.rmsprop_decay,
                                      FLAGS.rmsprop_momentum,
                                      FLAGS.rmsprop_epsilon)

      # Create training op
      train_tensor = slim.learning.create_train_op(
          total_loss,
          optimizer=opt,
          update_ops=tf.get_collection(tf.GraphKeys.UPDATE_OPS))

      # Summaries:
      slim.summaries.add_histogram_summaries(slim.get_model_variables())
      slim.summaries.add_scalar_summaries(slim.losses.get_losses(), 'losses')
      slim.summaries.add_scalar_summary(total_loss, 'Total_Loss', 'losses')
      slim.summaries.add_scalar_summary(learning_rate, 'Learning_Rate',
                                        'training')
      slim.summaries.add_histogram_summaries(endpoints.values())
      slim.summaries.add_zero_fraction_summaries(endpoints.values())
      # redacted

      # Set start-up delay
      startup_delay_steps = FLAGS.task * FLAGS.startup_delay_steps

      init_fn = model_init_function(model, dataset.num_classes,
                                    FLAGS.start_from_checkpoint)

      saver = tf.train.Saver(
          max_to_keep=FLAGS.max_checkpoints_to_keep,
          keep_checkpoint_every_n_hours=FLAGS.keep_checkpoint_every_n_hours)

      # Train model
      slim.learning.train(
          train_tensor,
          number_of_steps=FLAGS.number_of_steps,
          logdir=FLAGS.train_dir,
          master=target,
          init_fn=init_fn,
          is_chief=is_chief,
          saver=saver,
          startup_delay_steps=startup_delay_steps,
          save_summaries_secs=FLAGS.save_summaries_secs,
          save_interval_secs=FLAGS.save_interval_secs)


def parse_and_run():
  """Parse TF_CONFIG to cluster_spec and call run().

  TF_CONFIG environment variable is available when running using
  gcloud either locally or on cloud. It has all the information required
  to create a ClusterSpec which is important for running distributed code.

  Raises:
    ValueError: If flags are invalid.
  """
  tf_config = os.environ.get('TF_CONFIG')

  for name in ['master', 'task', 'ps_tasks']:

    if getattr(FLAGS, name) and tf_config:
      raise ValueError(
          'Either the flag --%s or the environment variable TF_CONFIG can be'
          ' set but not both.' % name)

  # If TF_CONFIG is not available we are either running locally in Cloud
  # or distributed inside Google. On Cloud the default values of
  # FLAGS.master and FLAGS.task correspond to running training locally.
  # Inside Google they will be set as needed to configure local or distributed
  # training. Inside Google we don't need to explicitly set worker_device
  # in replica_device_setter becaue this will be set automatically based
  # on various flags.
  if not tf_config:
    device_fn = tf.train.replica_device_setter(FLAGS.ps_tasks)
    return run(FLAGS.master, FLAGS.task == 0, device_fn=device_fn)

  tf_config_json = json.loads(tf_config)

  cluster = tf_config_json.get('cluster')
  job_name = tf_config_json.get('task', {}).get('type')
  task_index = tf_config_json.get('task', {}).get('index')

  # If cluster information is empty run local
  if job_name is None or task_index is None:
    device_fn = tf.train.replica_device_setter(0)
    return run('', True, device_fn=device_fn)

  ps = cluster.get('ps', [])
  num_ps = len(ps)

  cluster_spec = tf.train.ClusterSpec(cluster)
  server = tf.train.Server(
      cluster_spec, job_name=job_name, task_index=task_index)

  if job_name == 'ps':
    server.join()
    return
  elif job_name in ['master', 'worker']:
    device_fn = tf.train.replica_device_setter(
        num_ps,
        worker_device='/job:%s/task:%d' % (job_name, task_index),
        cluster=cluster_spec)
    return run(server.target, job_name == 'master', device_fn=device_fn)


def main(_):
  """Run and handle retryable errors."""
  proto_utils.uses_fast_cpp_protos_or_die()

  logging_level.set_from_flag()
  for _ in range(FLAGS.num_retries + 1):
    try:
      parse_and_run()
      return
    except tf.errors.UnavailableError as e:
      # An UnavailableError indicates a gRPC error, typically this is
      # retryable.
      logging.error('Caught UnavailableError %s; will retry.', e)
    except tf.errors.InternalError as e:
      # Retry on an InternalError.
      logging.error('Caught InternalError %s; will retry.', e)


if __name__ == '__main__':
  tf.app.run()
