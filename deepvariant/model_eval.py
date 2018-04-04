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
import math
import os



from tensorflow import flags
import numpy as np
import tensorflow as tf

from absl import logging

from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import variant_utils
from deepvariant import data_providers
from deepvariant import logging_level
from deepvariant import modeling

slim = tf.contrib.slim
FLAGS = flags.FLAGS

flags.DEFINE_integer('batch_size', 64, 'The number of samples in each batch.')

flags.DEFINE_string('master', '',
                    'The TensorFlow master to use. Set to the empty string '
                    'to let TF pick a default.')

flags.DEFINE_string('checkpoint_dir', '/tmp/deepvariant/',
                    'Directory where the model was written to.')

flags.DEFINE_string('eval_dir', '/tmp/deepvariant/',
                    'Directory where the results are saved to.')

flags.DEFINE_integer('max_evaluations', None,
                     'Max number of batches to evaluate')

flags.DEFINE_integer('max_examples', 64 * 1024 * 4,
                     'Maximum number of examples to evaluate.')

flags.DEFINE_string('model_name', 'inception_v3',
                    'The name of the model to use for predictions.')

flags.DEFINE_string('dataset_config_pbtxt', None,
                    'The path to the dataset config file.')

flags.DEFINE_float('moving_average_decay', 0.9999,
                   'The decay to use for the moving average.')


def select_variants_weights(variant_p_func, encoded_variants, name=None):
  """Creates a Tensor with 1.0 values anywhere variant_p_func returns True.

  Creates a TensorFlow operation with tf.py_func that calls variant_p_func on
  each Variant proto in encoded_variants (after decoding it), returning a 1.0
  for each variant where variant_p_func returns True and 0.0 where it returns
  value. For example, if is_snp returns True when a Variant is a SNP, then:

    weights = select_variants_weights(is_snp, encoded_variants)

  produces a weights tensor with 1.0 for each SNP in encoded_variants and 0.0
  for any non-SNP variants.

  Args:
    variant_p_func: a unary function accepting a nucleus.genomics.v1.Variant pb
      and returning a boolean value with True indicating the variant is part of
      the set and False indicating it is not. This function should be stateless.
    encoded_variants: An iterable of elements where each element is a string
      encoded value Variant protobuf.
    name: The optional name for this TF op.

  Returns:
    A TensorFlow Op.
  """

  def _select(encoded_variants):
    weights = [
        1.0 * variant_p_func(variant)
        for variant in variant_utils.decode_variants(encoded_variants)
    ]
    return np.array(weights, dtype=np.float32)

  return tf.py_func(
      _select, [encoded_variants], tf.float32, stateful=False, name=name)


def calling_metrics(metrics_map, selectors_map, predictions, labels):
  """Creates a dictionary of name to slim.metric functions.

  This function creates a dictionary of metric names to metric ops processing
  (using the predictions / labels for evaluation). The selectors_map
  defines a mapping from a selection name (e.g., "SNPs") to a weights tensor
  that defines elements of predictions/class_lables are part of the selection
  and which are not.

  This function creates a new metric that applies the raw metric from
  metrics_map to the weights for each selector in selectors_map.

  Args:
    metrics_map: A map from string name to a function accepting predictions,
      labels, and weights
      arguments.
    selectors_map: A map from a string name to a weight tensor. The weights
      tensor
      should have dimensions compatible with predictions and labels.
    predictions: A Tensor of predictions for our labels. Should have
      dimensions compatible with labels.
    labels: A Tensor of true labels with dimensions compatible with predictions.

  Returns:
    A dictionary of string => object where object is the type returned by a call
    to a metrics_map.values() element.
  """
  return {
      mname + '/' + sname: mfunc(
          predictions=predictions, labels=labels, weights=weights)
      for sname, weights in selectors_map.items()
      for mname, mfunc in metrics_map.items()
  }


def main(_):
  proto_utils.uses_fast_cpp_protos_or_die()

  if not FLAGS.dataset_config_pbtxt:
    logging.error('Need to specify --dataset_config_pbtxt')
  logging_level.set_from_flag()
  eval_loop(
      master=FLAGS.master,
      dataset_config_pbtxt=FLAGS.dataset_config_pbtxt,
      checkpoint_dir=FLAGS.checkpoint_dir,
      model_name=FLAGS.model_name,
      batch_size=FLAGS.batch_size,
      moving_average_decay=FLAGS.moving_average_decay,
      max_examples=FLAGS.max_examples,
      eval_dir=FLAGS.eval_dir,
      max_evaluations=FLAGS.max_evaluations,
  )


def make_metrics(predictions, labels, encoded_variants):
  """Creates our evaluation metrics."""
  # Define the metrics we'll get for each variant selection:
  raw_metrics = {
      'Accuracy': tf.metrics.accuracy,
      'Precision': tf.metrics.precision,
      'Recall': tf.metrics.recall,
      'FPs': tf.metrics.false_positives,
      'FNs': tf.metrics.false_negatives,
      'TPs': tf.metrics.true_positives,
      # redacted
      # 'TNs': tf.metrics.true_negatives,
  }

  def _make_selector(func):
    return select_variants_weights(func, encoded_variants)

  selectors = {
      'All': None,
      'SNPs': _make_selector(variant_utils.is_snp),
      'Indels': _make_selector(variant_utils.is_indel),
      # These haven't proven particularly useful, but are commented out here
      # in case someone wants to do some more explorations.
      # 'Insertions': _make_selector(variant_utils.has_insertion),
      # 'Deletions': _make_selector(variant_utils.has_deletion),
      # 'BiAllelic': _make_selector(variant_utils.is_biallelic),
      # 'MultiAllelic': _make_selector(variant_utils.is_multiallelic),
      # 'HomRef': tf.equal(labels, 0),
      # 'Het': tf.equal(labels, 1),
      # 'HomAlt': tf.equal(labels, 2),
      # 'NonRef': tf.greater(labels, 0),
  }

  return tf.contrib.metrics.aggregate_metric_map(
      calling_metrics(raw_metrics, selectors, predictions, labels))


def checkpoints_iterator(checkpoint_dir):
  # This is here to make it easy to mock out the iterator for tests.
  return tf.contrib.training.checkpoints_iterator(checkpoint_dir)


def eval_loop(master, dataset_config_pbtxt, checkpoint_dir, model_name,
              batch_size, moving_average_decay, max_examples, eval_dir,
              max_evaluations):
  logging.info('Running fixed eval for: %s', dataset_config_pbtxt)

  num_evaluations = 0
  for checkpoint_path in checkpoints_iterator(checkpoint_dir):
    logging.info('Using checkpoint %s %d', checkpoint_path, num_evaluations)

    g = tf.Graph()
    with g.as_default():
      tf_global_step = tf.train.get_or_create_global_step()

      # redacted
      model = modeling.get_model(model_name)
      dataset = data_providers.get_dataset(dataset_config_pbtxt)
      logging.info('Running evaluations on %s with model %s', dataset, model)

      images, labels, encoded_variant = data_providers.make_batches(
          dataset.get_slim_dataset(), model, batch_size, mode='EVAL')
      endpoints = model.create(images, dataset.num_classes, is_training=False)
      predictions = tf.argmax(endpoints['Predictions'], 1)

      # For eval, explicitly add moving_mean and moving_variance variables to
      # the MOVING_AVERAGE_VARIABLES collection.
      variable_averages = tf.train.ExponentialMovingAverage(
          moving_average_decay, tf_global_step)

      for var in tf.get_collection('moving_vars'):
        tf.add_to_collection(tf.GraphKeys.MOVING_AVERAGE_VARIABLES, var)
      for var in slim.get_model_variables():
        tf.add_to_collection(tf.GraphKeys.MOVING_AVERAGE_VARIABLES, var)

      variables_to_restore = variable_averages.variables_to_restore()
      variables_to_restore[tf_global_step.op.name] = tf_global_step

      names_to_values, names_to_updates = make_metrics(predictions, labels,
                                                       encoded_variant)

      for name, value in names_to_values.iteritems():
        slim.summaries.add_scalar_summary(value, name, print_summary=True)

      num_batches = int(
          math.floor(
              min(max_examples, dataset.num_examples) / float(batch_size)))
      num_samples = batch_size * num_batches
      logging.info('Dataset has %d samples, doing eval over %d',
                   dataset.num_examples, num_samples)

      names_to_values = slim.evaluation.evaluate_once(
          master=master,
          checkpoint_path=checkpoint_path,
          logdir=eval_dir,
          variables_to_restore=variables_to_restore,
          num_evals=num_batches,
          initial_op=tf.group(tf.global_variables_initializer(),
                              tf.local_variables_initializer()),
          eval_op=names_to_updates.values(),
          final_op=names_to_values,
      )

      # --- LOW LEVEL [WIP], hangs, initialization seems busted ---
      # This is (marginally) nicer as it can eliminate the slim dep.
      # saver = tf.train.Saver(variables_to_restore)
      # scaffold = tf.train.Scaffold(saver=saver)
      # names_to_values = tf.contrib.training.evaluate_once(
      #     checkpoint_path=checkpoint_path,
      #     master=FLAGS.master,
      #     scaffold=scaffold,
      #     eval_ops=names_to_updates.values(),
      #     final_ops=names_to_values,
      # )

      _write_checkpoint_metrics(
          checkpoint_path, names_to_values, eval_dir=eval_dir)

    num_evaluations += 1
    if max_evaluations is not None and num_evaluations >= max_evaluations:
      return


def checkpoint_metrics_path(checkpoint_path, eval_dir):
  """Gets a path to the JSON of eval metrics for checkpoint in eval_dir."""
  return os.path.join(eval_dir, os.path.basename(checkpoint_path) + '.metrics')


def read_metrics(checkpoint_path, eval_dir):
  """Reads the JSON of metrics for checkpoint_path in eval_dir."""
  metrics_path = checkpoint_metrics_path(checkpoint_path, eval_dir)
  with tf.gfile.GFile(metrics_path) as fin:
    return {k: float(v) for k, v in json.load(fin).iteritems()}


def _write_checkpoint_metrics(checkpoint_path, metrics_and_values, eval_dir):
  """Writes a JSON of metrics for checkpoint_path in eval_dir.

  This function writes out metrics to a JSON for a checkpoint into eval_dir. The
  exact path of this file will be computed with:

    `checkpoint_metrics_path(checkpoint_path, eval_dir)`

  and the values for metrics_and_values (a dict of strings => objects) written
  out as key: str(object) into a JSON file.

  Args:
    checkpoint_path: str; a path to the checkpoint we computed metrics on.
    metrics_and_values: dict[string,object]; a dictionary of key/value pairs
      containing our metrics. These will be converted to a JSON of key/string
      pairs and written out to disk.
    eval_dir: str; a path to a directory where we will write out our checkpoint
      metrics.
  """
  path = checkpoint_metrics_path(checkpoint_path, eval_dir)
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
