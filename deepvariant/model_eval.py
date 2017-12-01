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



import numpy as np
import tensorflow as tf

from absl import logging

from deepvariant import data_providers
from deepvariant import logging_level
from deepvariant import modeling
from deepvariant.core import proto_utils
from deepvariant.core import variantutils

slim = tf.contrib.slim
FLAGS = tf.flags.FLAGS

tf.flags.DEFINE_integer('batch_size', 64,
                        'The number of samples in each batch.')

tf.flags.DEFINE_string('master', '',
                       'The TensorFlow master to use. Set to the empty string '
                       'to let TF pick a default.')

tf.flags.DEFINE_string('checkpoint_dir', '/tmp/deepvariant/',
                       'Directory where the model was written to.')

tf.flags.DEFINE_string('eval_dir', '/tmp/deepvariant/',
                       'Directory where the results are saved to.')

tf.flags.DEFINE_integer(
    'eval_interval_secs', 600,
    'The frequency, in seconds, with which evaluation is run.')

tf.flags.DEFINE_integer('batches_per_eval_step', 1000,
                        'Number of batches to evaluate in each eval step.')

tf.flags.DEFINE_integer('max_evaluations', None,
                        'Max number of batches to evaluate')

tf.flags.DEFINE_string('model_name', 'inception_v3',
                       'The name of the model to use for predictions.')

tf.flags.DEFINE_string('dataset_config_pbtxt', None,
                       'The path to the dataset config file.')

tf.flags.DEFINE_float('moving_average_decay', 0.9999,
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
    variant_p_func: a unary function accepting a learning.genomics.v1.Variant pb
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
        for variant in variantutils.decode_variants(encoded_variants)
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
      mname + '/' + sname: mfunc(predictions, labels, weights=weights)
      for sname, weights in selectors_map.items()
      for mname, mfunc in metrics_map.items()
  }


def main(_):
  proto_utils.uses_fast_cpp_protos_or_die()

  if not FLAGS.dataset_config_pbtxt:
    logging.error('Need to specify --dataset_config_pbtxt')
  logging_level.set_from_flag()

  g = tf.Graph()
  with g.as_default():
    tf_global_step = slim.get_or_create_global_step()

    model = modeling.get_model(FLAGS.model_name)
    dataset = data_providers.get_dataset(FLAGS.dataset_config_pbtxt)
    print('Running evaluations on {} with model {}\n'.format(dataset, model))

    batch = data_providers.make_training_batches(dataset.get_slim_dataset(),
                                                 model, FLAGS.batch_size)
    images, labels, encoded_truth_variants = batch
    endpoints = model.create(images, dataset.num_classes, is_training=False)
    predictions = tf.argmax(endpoints['Predictions'], 1)

    # For eval, explicitly add moving_mean and moving_variance variables to
    # the MOVING_AVERAGE_VARIABLES collection.
    variable_averages = tf.train.ExponentialMovingAverage(
        FLAGS.moving_average_decay, tf_global_step)

    for var in tf.get_collection('moving_vars'):
      tf.add_to_collection(tf.GraphKeys.MOVING_AVERAGE_VARIABLES, var)
    for var in slim.get_model_variables():
      tf.add_to_collection(tf.GraphKeys.MOVING_AVERAGE_VARIABLES, var)

    variables_to_restore = variable_averages.variables_to_restore()
    variables_to_restore[tf_global_step.op.name] = tf_global_step

    # Define the metrics:
    metrics = {
        'Accuracy':
            tf.contrib.metrics.streaming_accuracy,
        'Mean_absolute_error':
            tf.contrib.metrics.streaming_mean_absolute_error,
        'FPs':
            tf.contrib.metrics.streaming_false_positives,
        'FNs':
            tf.contrib.metrics.streaming_false_negatives,
    }

    def _make_selector(func):
      return select_variants_weights(func, encoded_truth_variants)

    selectors = {
        'All': None,
        'SNPs': _make_selector(variantutils.is_snp),
        'Indels': _make_selector(variantutils.is_indel),
        'Insertions': _make_selector(variantutils.has_insertion),
        'Deletions': _make_selector(variantutils.has_deletion),
        'BiAllelic': _make_selector(variantutils.is_biallelic),
        'MultiAllelic': _make_selector(variantutils.is_multiallelic),
        # These haven't proven particularly useful, but are commented out here
        # in case someone wants to do some more explorations.
        # 'HomRef': tf.equal(labels, 0),
        # 'Het': tf.equal(labels, 1),
        # 'HomAlt': tf.equal(labels, 2),
        # 'NonRef': tf.greater(labels, 0),
    }
    metrics = calling_metrics(metrics, selectors, predictions, labels)
    names_to_values, names_to_updates = slim.metrics.aggregate_metric_map(
        metrics)

    for name, value in names_to_values.iteritems():
      slim.summaries.add_scalar_summary(value, name, print_summary=True)

    slim.evaluation.evaluation_loop(
        FLAGS.master,
        FLAGS.checkpoint_dir,
        logdir=FLAGS.eval_dir,
        num_evals=FLAGS.batches_per_eval_step,
        eval_op=names_to_updates.values(),
        variables_to_restore=variables_to_restore,
        max_number_of_evaluations=FLAGS.max_evaluations,
        eval_interval_secs=FLAGS.eval_interval_secs)


if __name__ == '__main__':
  tf.app.run()
