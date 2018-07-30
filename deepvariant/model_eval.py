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
"""Evaluates a DeepVariant model during training.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json
import os



from absl import flags
from absl import logging
import tensorflow as tf

from third_party.nucleus.util import proto_utils
from deepvariant import data_providers
from deepvariant import logging_level
from deepvariant import modeling
from deepvariant import tf_utils


slim = tf.contrib.slim
FLAGS = flags.FLAGS

flags.DEFINE_integer('batch_size', 64, 'The number of samples in each batch.')

# Cloud TPU Cluster Resolvers
flags.DEFINE_string(
    'gcp_project', None,
    'Project name for the Cloud TPU-enabled project. If not specified, we '
    'will attempt to automatically detect the GCE project from metadata.')
flags.DEFINE_string(
    'tpu_zone', None,
    'GCE zone where the Cloud TPU is located in. If not specified, we '
    'will attempt to automatically detect the GCE project from metadata.')
flags.DEFINE_string(
    'tpu_name', None,
    'Name of the Cloud TPU for Cluster Resolvers. You must specify either '
    'this flag or --master.')

flags.DEFINE_string(
    'master', None,
    'GRPC URL of the master (e.g. grpc://ip.address.of.tpu:8470). You '
    'must specify either this flag or --tpu_name.')

flags.DEFINE_string('checkpoint_dir', '/tmp/deepvariant/',
                    'Directory where the model was written to.')

flags.DEFINE_string(
    'eval_name', None,
    'Name of the evaluation if user needs to run multiple evaluations on '
    'different data sets, such as on training data vs test data. Metrics for '
    'different evaluations are saved in separate directories, and appear '
    'separately in tensorboard.  The directory will be named "eval_"+eval_name')

flags.DEFINE_string(
    'eval_dir', None,
    'This is used only to generate eval_name, if that is not provided.')

flags.DEFINE_integer('min_eval_interval_s', 180,
                     'Minimum seconds between evaluations.')

flags.DEFINE_integer(
    'eval_timeout', None,
    'Maximum seconds between checkpoints before evaluation '
    'terminates.')

flags.DEFINE_integer('max_evaluations', None,
                     'Max number of batches to evaluate')

flags.DEFINE_integer('max_examples', 64 * 1024 * 4,
                     'Maximum number of examples to evaluate.')

flags.DEFINE_string('model_name', 'inception_v3',
                    'The name of the model to use for predictions.')

flags.DEFINE_string('dataset_config_pbtxt', None,
                    'The path to the dataset config file.')

flags.DEFINE_boolean('use_tpu', False, 'use tpu if available')

flags.DEFINE_string(
    'kmp_blocktime', '0',
    'Value to set the KMP_BLOCKTIME environment variable to for efficient MKL '
    'evaluation. See https://www.tensorflow.org/performance/performance_guide '
    'for more information. The default value is 0, which provides the best '
    'performance in our tests. Set this flag to "" to not set the variable.')


def main(_):
  proto_utils.uses_fast_cpp_protos_or_die()

  if not FLAGS.dataset_config_pbtxt:
    logging.error('Need to specify --dataset_config_pbtxt')
  logging_level.set_from_flag()

  if FLAGS.kmp_blocktime:
    os.environ['KMP_BLOCKTIME'] = FLAGS.kmp_blocktime
    logging.info('Set KMP_BLOCKTIME to %s', os.environ['KMP_BLOCKTIME'])

  master = tf_utils.resolve_master(FLAGS.master, FLAGS.tpu_name, FLAGS.tpu_zone,
                                   FLAGS.gcp_project)
  eval_loop(
      master=master,
      dataset_config_pbtxt=FLAGS.dataset_config_pbtxt,
      checkpoint_dir=FLAGS.checkpoint_dir,
      model_name=FLAGS.model_name,
      batch_size=FLAGS.batch_size,
      max_examples=FLAGS.max_examples,
      eval_name=FLAGS.eval_name,
      max_evaluations=FLAGS.max_evaluations,
      use_tpu=FLAGS.use_tpu,
  )


def checkpoints_iterator(checkpoint_dir,
                         min_interval_secs=0,
                         timeout=None,
                         timeout_fn=None):
  # This is here to make it easy to mock out the iterator for tests.
  return tf.contrib.training.checkpoints_iterator(
      checkpoint_dir, min_interval_secs, timeout, timeout_fn)


def eval_loop(master,
              dataset_config_pbtxt,
              checkpoint_dir,
              model_name,
              batch_size,
              max_examples,
              eval_name,
              max_evaluations,
              use_tpu=False):
  """Evaluate incoming checkpoints, until the specified end."""
  logging.info('Running fixed eval for: %s', dataset_config_pbtxt)

  tf_dataset = data_providers.get_input_fn_from_dataset(
      dataset_config_filename=dataset_config_pbtxt,
      mode=tf.estimator.ModeKeys.EVAL,
      use_tpu=use_tpu,
  )

  model = modeling.get_model(model_name)
  logging.info('Running evaluations on %s with model %s', tf_dataset, model)

  # Compute when to stop reading, in terms of batches.
  num_batches = min(max_examples, tf_dataset.num_examples) // batch_size
  num_samples = batch_size * num_batches
  logging.info(
      'Dataset has %d samples, doing eval over %d; '
      'max_examples is %d, num_batches is %d', tf_dataset.num_examples,
      num_samples, max_examples, num_batches)
  batches_per_epoch = tf_dataset.num_examples / batch_size

  # This loads EMA variables.
  eval_hooks = [h(checkpoint_dir) for h in model.session_eval_hooks()]

  classifier = model.make_estimator(
      batch_size=batch_size,
      model_dir=checkpoint_dir,
      params={'batches_per_epoch': batches_per_epoch},
      use_tpu=use_tpu,
      master=master,
  )

  def terminate_eval():
    logging.info('Terminating eval after %d seconds of no checkpoints',
                 FLAGS.eval_timeout)
    return True

  # Run evaluation when there's a new checkpoint
  num_evaluations = 0
  for ckpt in checkpoints_iterator(
      checkpoint_dir=checkpoint_dir,
      min_interval_secs=FLAGS.min_eval_interval_s,
      timeout=FLAGS.eval_timeout,
      timeout_fn=terminate_eval):

    logging.info('Starting to evaluate.')

    # For each step, calls input_fn, which returns one batch of data.
    # Evaluates until either steps batches are processed, or input_fn raises an
    # end-of-input exception (OutOfRangeError or StopIteration).
    eval_results = classifier.evaluate(
        input_fn=tf_dataset,
        steps=num_batches,
        hooks=eval_hooks,
        checkpoint_path=ckpt,
        name=eval_name)
    logging.info('Eval results: %s', eval_results)

    _write_checkpoint_metrics(ckpt, eval_results, eval_name)

    # An alternative strategy might check step-number-of-ckpt >= train_steps.
    num_evaluations += 1
    if max_evaluations is not None and num_evaluations >= max_evaluations:
      logging.info('Evaluation finished after %d evaluations', num_evaluations)
      break

  return


def checkpoint_metrics_path(checkpoint_path, eval_name):
  """Gets a path to the JSON of eval metrics for checkpoint in eval_name."""
  checkpoint_dir = os.path.dirname(checkpoint_path)
  checkpoint_name = os.path.basename(checkpoint_path)
  if eval_name:
    # This bit of magic is defined by the estimator framework, and isn't easy
    # to change.  We only get to specify the suffix.
    d = os.path.join(checkpoint_dir, 'eval_' + eval_name)
  else:
    d = checkpoint_dir
  return os.path.join(d, checkpoint_name + '.metrics')


def read_metrics(checkpoint_path, eval_name):
  """Reads the JSON of metrics for checkpoint_path in eval_dir."""
  metrics_path = checkpoint_metrics_path(checkpoint_path, eval_name)
  with tf.gfile.GFile(metrics_path) as fin:
    return {k: float(v) for k, v in json.load(fin).iteritems()}


def _write_checkpoint_metrics(checkpoint_path, metrics_and_values, eval_name):
  """Writes a JSON of metrics for checkpoint_path in eval_name.

  This function writes out metrics to a JSON for a checkpoint into
  .../eval_name.
  The exact path of this file will be computed with:

    `checkpoint_metrics_path(checkpoint_path, eval_name)`

  and the values for metrics_and_values (a dict of strings => objects) written
  out as key: str(object) into a JSON file.

  Args:
    checkpoint_path: str; a path to the checkpoint we computed metrics on.
    metrics_and_values: dict[string,object]; a dictionary of key/value pairs
      containing our metrics. These will be converted to a JSON of key/string
      pairs and written out to disk.
    eval_name: str; the name of the eval run, which is used to derive the
      the subdirectory of checkpoint_path where the eval metrics will be
      written.
  """
  path = checkpoint_metrics_path(checkpoint_path, eval_name)
  serializable = {k: str(v) for k, v in metrics_and_values.iteritems()}
  logging.info('Writing checkpoint metrics %s', path)
  try:
    with tf.gfile.GFile(path, 'w') as fout:
      json.dump(serializable, fout, sort_keys=True, indent=4)
  except:  # pylint: disable=bare-except
    # Note we have a bare exception here as as there's no clear TF base
    # exception to catch will cover all of the potential issues that might arise
    # trying to write our metrics to our metrics file.
    logging.warning('Failed to write checkpoint metrics to path %s', path)


if __name__ == '__main__':
  tf.app.run()
