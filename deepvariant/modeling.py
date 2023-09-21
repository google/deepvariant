# Copyright 2017 Google LLC.
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
"""Provides an abstraction around deep learning models in DeepVariant.

This class allows us to encapsulate all of the model management, loading,
saving, and data processing in a single place so those details don't spill over
into the more general deepvariant codebase. The key thing we are aiming for here
is to make sure we can easily play with other model architectures without
modifying the surrounding training and evaluation code.
"""

import enum
import itertools
import math



from absl import flags
from absl import logging
from tensorflow import estimator as tf_estimator
from tensorflow.compat.v1 import estimator as tf_compat_v1_estimator
import tensorflow as tf
import tf_slim
from tf_slim.nets import inception_v3

from deepvariant import dv_constants
from deepvariant import dv_utils
# pylint: disable=g-direct-tensorflow-import
from tensorflow.python.framework import ops
from tensorflow.python.tpu import tpu_config
from tensorflow.python.tpu import tpu_estimator
from tensorflow.python.tpu import tpu_optimizer
# pylint: enable=g-direct-tensorflow-import

tf.compat.v1.disable_eager_execution()

flags.DEFINE_float(
    'label_smoothing',
    1e-6,
    (
        'Amount of label smoothing to use. By default this is 0.0001% '
        'meaning that we expect a label error at a rate of 1 / 1,000,000'
    ),
)

# Training parameters.
flags.DEFINE_float('learning_rate', 0.064, 'Initial learning rate.')

flags.DEFINE_float('rmsprop_momentum', 0.9, 'Momentum.')

flags.DEFINE_float('rmsprop_decay', 0.9, 'Decay term for RMSProp.')

flags.DEFINE_float('rmsprop_epsilon', 1.0, 'Epsilon term for RMSProp.')

flags.DEFINE_float(
    'learning_rate_decay_factor', 0.94, 'Learning rate decay factor.'
)

flags.DEFINE_float(
    'num_epochs_per_decay',
    2.0,
    'Number of epochs after which learning rate decays.',
)

flags.DEFINE_float(
    'moving_average_decay', 0.9999, 'The decay to use for the moving average.'
)

flags.DEFINE_integer(
    'save_summary_steps',
    100,
    'Number of steps which must have run before showing summaries.',
)

flags.DEFINE_integer(
    'save_interval_secs',
    60 * 10,
    (
        'Interval (in seconds) at which the model data '
        'should be checkpointed. Set to 0 to disable, -1 to ignore. '
        'Exclusive with save_interval_steps.'
    ),
)

flags.DEFINE_integer(
    'save_interval_steps',
    -1,
    (
        'Interval (in steps) at which the model data '
        'should be checkpointed. Set to 0 to disable, -1 to ignore. '
        'Exclusive with save_interval_secs.'
    ),
)

flags.DEFINE_integer(
    'seq_type_embedding_size',
    200,
    (
        'Set the embedding size for the sequencing type embeddings. Default is'
        ' 200. This flag is only useful when model_name is'
        ' `inception_v3_embedding`.'
    ),
)

flags.DEFINE_bool(
    'allow_warmstart_from_different_num_channels',
    False,
    (
        'If True, always allow warmstarting from model checkpoints '
        'that has different number of channels.'
    ),
)

FLAGS = flags.FLAGS

slim = tf_slim


class UnsupportedImageDimensionsError(Exception):
  """Exception indicating the image dimensions aren't supported by our model."""


def binarize(labels, target_class):
  """Binarize labels and predictions.

  The labels that are equal to target_class parameter are set to 0, else
  set to 1.

  Args:
    labels: the ground-truth labels for the examples.
    target_class: index of the class that is left as zero.

  Returns:
    Tensor of the same shape as labels.
  """
  labels_binary = tf.compat.v1.where(
      tf.equal(labels, tf.constant(target_class, dtype=tf.int64)),
      tf.zeros_like(labels),
      tf.ones_like(labels),
  )
  return labels_binary


def get_class_recall(labels, predicted_class, target_class):
  """Compute recall from labels and predicted_class for target_class.

  Examples with label target_class are positives. Other classes are negatives.

  Args:
    labels: the ground-truth labels for the examples.
    predicted_class: the predicted labels for the examples.
    target_class: index of the class that is left as non-zero.

  Returns:
    Tensor containing the recall value.
  """
  labels_binary = binarize(labels, target_class)
  predicted_class_binary = binarize(predicted_class, target_class)
  return tf.compat.v1.metrics.recall(labels_binary, predicted_class_binary)


def get_class_precision(labels, predicted_class, target_class):
  """Compute precision from labels and predicted_class for target_class.

  Examples with label target_class are positives. Other classes are negatives.

  Args:
    labels: the ground-truth labels for the examples.
    predicted_class: the predicted labels for the examples.
    target_class: index of the class that is left as non-zero.

  Returns:
    Tensor containing the precision value.
  """
  labels_binary = binarize(labels, target_class)
  predicted_class_binary = binarize(predicted_class, target_class)
  return tf.compat.v1.metrics.precision(labels_binary, predicted_class_binary)


# TODO: Verify this F1 score is correct.
def get_f1_score(labels, predictions, target_class=None):
  """Compute F1 score of predictions with respect to the labels.

  Args:
    labels: tensor whose dimensions must match predictions. The ground-truth
      labels for the examples.
    predictions: tensor of arbitrary dimension. The predicted labels for the
      examples.
    target_class: int. Index of the class that is left as non-zero.

  Returns:
    f1_score: scalar float tensor whose dimensions match predictions. The
    calculated f1 score.
    update_op: operation that updates the f1 score streaming metric.
  """
  if target_class:
    labels = binarize(labels, target_class)
    predictions = binarize(predictions, target_class)

  precision, precision_op = tf.compat.v1.metrics.precision(labels, predictions)
  recall, recall_op = tf.compat.v1.metrics.recall(labels, predictions)

  def compute_f1_score(name):
    pr_product = tf.multiply(precision, recall)
    return tf.math.divide(
        tf.multiply(2.0, pr_product),
        tf.add(tf.add(precision, recall), 1e-12),
        name,
    )

  f1_score = compute_f1_score('value')
  with ops.control_dependencies([precision_op, recall_op]):
    update_op = compute_f1_score('update_op')

  return f1_score, update_op


def is_encoded_variant_type(variant_types_tensor, type_to_select):
  """Returns a bool tensor indicating which variant_types match type_to_select.

  Args:
    variant_types_tensor: Tensor of shape (batch_size, 1) containing
      EncodedVariantType.value int64 values. Each element of this tensor should
      be a EncodedVariantType.value int64 value indicating the type of the
      variant.
    type_to_select: EncodedVariantType. The type of variant we want to select.

  Returns:
    Tensor of shape (batch_size, 1) of type tf.bool. A True value indicates that
    the variant_type at that position matched type_to_select. Has a False
    otherwise.
  """
  return tf.equal(
      variant_types_tensor, tf.constant(type_to_select.value, dtype=tf.int64)
  )


# This dictionary contains a mapping from the human readable name of a metric
# function (e.g., Accuracy) and its associated TensorFlow metric function. All
# of the entries here will be stratified by variant_type in eval_metric_fn.
_METRICS_FUNCS_BY_VARIANT_TYPE = {
    'Accuracy': tf.compat.v1.metrics.accuracy,
    'Precision': tf.compat.v1.metrics.precision,
    'Recall': tf.compat.v1.metrics.recall,
    'FPs': tf.compat.v1.metrics.false_positives,
    'FNs': tf.compat.v1.metrics.false_negatives,
    'TPs': tf.compat.v1.metrics.true_positives,
    'TNs': tf.compat.v1.metrics.true_negatives,
}

# A set containing the names of the variant types we split our metrics by type
# by. This data structure isn't a dictionary like it's neighbors because
# eval_metric_fn requires special logic to compute the values here associated
# with each of these names.
_METRICS_BY_VARIANT_TYPE = {'All', 'SNPs', 'Indels'}

# This dictionary contains a mapping from the human readable name of a genotype
# class (e.g., Het) and its associated class label (e.g., 1). All of the entries
# here will be stratified by genotype_class in eval_metric_fn.
_METRICS_GENOTYPE_CLASSES = {
    'HomRef': 0,
    'Het': 1,
    'HomVar': 2,
}

# This dictionary contains a mapping from the human readable name of a metric
# function (e.g., Accuracy) and its associated metric function. All
# of the entries here will be stratified by genotype class (e.g., Het) in
# eval_metric_fn.
_METRICS_FUNCS_BY_GENOTYPE_CLASS = {
    'Precision': get_class_precision,
    'Recall': get_class_recall,
    'F1': get_f1_score,
}


def _eval_name(metric_name, stratification_name):
  return metric_name + '/' + stratification_name


class EvalMetricOrdering(enum.Enum):
  """Enum capturing whether a better metric should be larger or smaller."""

  BIGGER_IS_BETTER = 1
  SMALLER_IS_BETTER = 2


def eval_function_metrics(has_variant_types=True):
  """Gets the set of eval_metrics names and their directionality.

  Args:
    has_variant_types: bool. Will we be providing variant_type information
      during eval so that we'll have metrics stratified by variant_type?

  Returns:
    dict mapping from a metric name string (e.g., "F1/All") and a
    EvalMetricOrdering enum indicating whether larger metric values are better
    or worse.
  """
  names = {_eval_name('F1', 'All'): EvalMetricOrdering.BIGGER_IS_BETTER}

  if has_variant_types:
    variant_type_names = _METRICS_BY_VARIANT_TYPE
  else:
    variant_type_names = {'All'}
  for m, s in itertools.product(
      _METRICS_FUNCS_BY_VARIANT_TYPE, variant_type_names
  ):
    names[_eval_name(m, s)] = EvalMetricOrdering.BIGGER_IS_BETTER

  for m, s in itertools.product(
      _METRICS_FUNCS_BY_GENOTYPE_CLASS, _METRICS_GENOTYPE_CLASSES
  ):
    names[_eval_name(m, s)] = EvalMetricOrdering.BIGGER_IS_BETTER

  return names


# NB. This includes only a subset of our usual metrics.
# We'll add the rest back in a subsequent change.
def eval_metric_fn(labels, predictions, variant_types):
  """Calculate eval metrics from Tensors, on CPU host.

  Args:
    labels: the ground-truth labels for the examples.
    predictions: the predicted labels for the examples.
    variant_types: variant types (int64 of EncodedVariantType.value) as a tensor
      of (batch_size,) or None. The types of these variants. If None, no type
      specific evals will be performed.

  Returns:
    A dictionary of string name to metric.
  """
  predicted_classes = tf.argmax(input=predictions, axis=1)

  metrics = {}

  # Add the metrics stratified by variant_type
  weights_by_type = {'All': None}
  if variant_types is not None:
    weights_by_type['SNPs'] = is_encoded_variant_type(
        variant_types, dv_utils.EncodedVariantType.SNP
    )
    weights_by_type['Indels'] = is_encoded_variant_type(
        variant_types, dv_utils.EncodedVariantType.INDEL
    )

  for metric_name, metric_func in _METRICS_FUNCS_BY_VARIANT_TYPE.items():
    for weight_name, weights in weights_by_type.items():
      metrics[_eval_name(metric_name, weight_name)] = metric_func(
          labels, predicted_classes, weights=weights
      )

  # Add the metrics stratified by predicted class.
  for metric_name, metric_func in _METRICS_FUNCS_BY_GENOTYPE_CLASS.items():
    for class_name, class_value in _METRICS_GENOTYPE_CLASSES.items():
      metrics[_eval_name(metric_name, class_name)] = metric_func(
          labels, predicted_classes, class_value
      )

  # Special case F1/All to avoid a clash between the two different ways that we
  # can compute Precision and Recall (e.g., get_class_precision vs.
  # tf.compat.v1.metrics.precision.
  metrics[_eval_name('F1', 'All')] = get_f1_score(labels, predicted_classes)

  logging.info('Metrics are %s', metrics.keys())

  # Make sure our metrics are consistent with the expected names from
  # eval_function_metrics.
  expected_metrics = eval_function_metrics(
      has_variant_types=variant_types is not None
  )
  if set(expected_metrics) != set(metrics):
    raise AssertionError(
        'Bug: actual metrics={} not equal to expected={}'.format(
            ','.join(metrics), ','.join(expected_metrics)
        )
    )

  return metrics


# The following two classes support loading exponential moving averages into
# their corresponding variables when a checkpoint is loaded. They're called
# as hooks by the Estimators. Note for future work: this is the documented
# way, but someone on the mailing list suggested that using the scaffold_fn
# mechanism might be better.


class LoadEMAHook(tf_estimator.SessionRunHook):
  """Hook to load EMA into their corresponding variables.

  This looks for the latest checkpoint in the model dir.
  """

  def __init__(self, model_dir, ignore_missing_vars=False):
    super(LoadEMAHook, self).__init__()
    self._model_dir = model_dir
    self._ignore_missing_vars = ignore_missing_vars

  def begin(self):
    ema = tf.train.ExponentialMovingAverage(FLAGS.moving_average_decay)
    variables_to_restore = ema.variables_to_restore()
    self._load_ema = slim.assign_from_checkpoint_fn(
        tf.train.latest_checkpoint(self._model_dir),
        variables_to_restore,
        ignore_missing_vars=self._ignore_missing_vars,
    )

  def after_create_session(self, sess, coord):
    tf.compat.v1.logging.info('Reloading EMA...')
    self._load_ema(sess)


class PredictEMAHook(tf_estimator.SessionRunHook):
  """Hook to load EMA into their corresponding variables.

  This reads the specified checkpoint.
  """

  def __init__(self, checkpoint_path, ignore_missing_vars=False):
    super(PredictEMAHook, self).__init__()
    self._checkpoint_path = checkpoint_path
    self._ignore_missing_vars = ignore_missing_vars

  def begin(self):
    ema = tf.train.ExponentialMovingAverage(FLAGS.moving_average_decay)
    variables_to_restore = ema.variables_to_restore()
    self._load_ema = slim.assign_from_checkpoint_fn(
        self._checkpoint_path,
        variables_to_restore,
        ignore_missing_vars=self._ignore_missing_vars,
    )

  def after_create_session(self, sess, coord):
    tf.compat.v1.logging.info('Reloading EMA...')
    self._load_ema(sess)


class DeepVariantModel(object):
  """Base class for models that compute genotype likelihoods from an image.

  This class is intended for use anywhere in DeepVariant where we want to train
  or evaluate a model that computes genotype likelihoods from a pileup image. A
  bit of encapsulation helps us to try new models (beyond inception_v3) and unit
  test our code.

  The base class cannot be used directly; concrete subclasses actually implement
  specific models and all of the associated machinery to create/load/save
  models.

  Attributes:
    name: str. The name of this model, such as `inception_v3`.
    pretrained_model_path: str. Path to a root checkpoint where we can start
      training the model, if we are not starting from scratch.
    supported_dimensions_message: str. A human-readable string containing info
      about what image dimensions are supported by this model. E.g., "only
      widths between 42 and 189".
    use_tpu: bool or None. If True, we are executing the model on a TPU, False
      if we are using some other hardware. If None, the execution hardware is
      not yet known.
    model_dir: str or None. The path to the location where model checkpoint are
      being stored. If None, the path hasn't been set yet or is unknown.
  """

  def __init__(self, name, pretrained_model_path):
    """Creates a new DeepVariantModel with name and pretrained_model_path.

    Args:
      name: str. The name of the model. Passed to DeepVariantModel name.
      pretrained_model_path: str. A path to a pretrained model to initialize our
        network from when starting from the 'model_default'.  If None, training
        will start from randomly-initialized parameters.

    Raises:
      ValueError: if any of the arguments is invalid.
    """
    if not name:
      raise ValueError('Got an empty value for name', name)
    self.name = name
    self.pretrained_model_path = pretrained_model_path
    self.supported_dimensions_message = 'unknown'
    self.use_tpu = None

    # Set the model_dir to None by default. We capture its actual value during
    # a call to make_estimator below.
    self.model_dir = None

  def construct_scalar_host_call(
      self, metric_dict, model_dir, prefix='', record_frequency_in_steps=100
  ):
    """Construct a host call to log scalars when training on TPU.

    Args:
      metric_dict: A dict of the tensors to be logged.
      model_dir: The location to write the summary.
      prefix: The prefix (if any) to prepend to the metric names.
      record_frequency_in_steps: int; How often should we log our metrics in
        step units.

    Returns:
      A tuple of (function, args_to_be_passed_to_said_function)
    """
    # type: (dict, str) -> (function, list)
    metric_names = list(metric_dict.keys())

    def host_call_fn(global_step, *args):
      """Training host call.

      Creates scalar summaries for training metrics.

      This function is executed on the CPU and should not directly reference
      any Tensors in the rest of the `model_fn`. To pass Tensors from the
      model to the `metric_fn`, provide as part of the `host_call`. See
      https://www.tensorflow.org/api_docs/python/tf/compat/v1/estimator/tpu/TPUEstimator
      for more information.
      Arguments should match the list of `Tensor` objects passed as the second
      element in the tuple passed to `host_call`.
      Args:
        global_step: Tensor with shape `[batch]` for the global_step
        *args: Remaining tensors to log.

      Returns:
        List of summary ops to run on the CPU host.
      """
      # TODO: When updating TF to v2.9.1, I had to remove
      # create_file_writer because it was giving me:
      # ValueError: Invalid argument to flush(): <tf.Tensor 'create_file_writer/SummaryWriter:0' shape=() dtype=resource>
      # In the future, try to add create_file_writer back for learning_rate.
      return tf.compat.v1.summary.all_v2_summary_ops()

    # To log the current learning rate, and gradient norm for Tensorboard, the
    # summary op needs to be run on the host CPU via host_call. host_call
    # expects [batch_size, ...] Tensors, thus reshape to introduce a batch
    # dimension. These Tensors are implicitly concatenated to
    # [params['batch_size']].
    global_step_tensor = tf.reshape(
        tf.compat.v1.train.get_or_create_global_step(), [1]
    )
    other_tensors = [tf.reshape(metric_dict[key], [1]) for key in metric_names]

    return host_call_fn, [global_step_tensor] + other_tensors

  def _create_warm_start_settings(self, start_from_checkpoint):
    """Create a proper WarmStartSettings based on start_from_checkpoint."""
    # If the special value "model_default" was passed, ask the model for
    # its default.
    if start_from_checkpoint == 'model_default':
      start_from_checkpoint = self.pretrained_model_path

    # If the path is non-False, use it.
    if start_from_checkpoint:
      logging.info(
          'Initializing model from checkpoint at %s', start_from_checkpoint
      )
      excluded_scopes = set()
      reader = tf.compat.v1.train.NewCheckpointReader(start_from_checkpoint)
      var_to_shape_map = reader.get_variable_to_shape_map()
      if (
          dv_utils.model_num_classes(
              start_from_checkpoint, self.n_classes_model_variable
          )
          != dv_constants.NUM_CLASSES
      ):
        excluded_scopes.update(self.excluded_scopes_for_incompatible_classes)
      if FLAGS.allow_warmstart_from_different_num_channels:
        excluded_scopes.update(self.excluded_scopes_for_incompatible_channels)
      if excluded_scopes:
        logging.info(
            (
                'The model checkpoint to warm start from has different '
                'shapes. If this is in training, we will '
                'exclude: %s'
            ),
            excluded_scopes,
        )
        vars_to_include = [
            v
            for v in var_to_shape_map.keys()
            if not v.startswith(tuple(excluded_scopes))
        ]
      else:
        logging.info(
            'The model checkpoint to warm start from should have the '
            'same number of classes and same numbers of channels.'
            'If this is in training, we will include everything for '
            'warm starting....'
        )
        vars_to_include = var_to_shape_map.keys()
      return tf_estimator.WarmStartSettings(
          ckpt_to_initialize_from=start_from_checkpoint,
          vars_to_warm_start='|'.join(vars_to_include),
      )
    else:
      # If warm_start_from is an empty string, specifically set it to None.
      logging.vlog(3, 'Initializing model with random parameters')
      return None

  def make_estimator(
      self,
      batch_size,
      model_dir=None,
      max_checkpoints_to_keep=100000,
      iterations_per_loop=100,
      params=None,
      unused_device_fn=None,
      master='',
      use_tpu=False,
      start_from_checkpoint=None,
      session_config=None,
      include_debug_info=False,
  ):
    """Returns a new tf.estimator.Estimator object for training or prediction.

    The estimator needs to know batch_size. We use the same value for all
    of eval, train, and predict. The estimator will automatically save
    checkpoints to model_dir and keep the specified number of them. The value
    of iterations_per_loop is not critical, and we default to the recommended
    value. Some optional arguments are only required for use with TPU.

    This function will use self.model_fn and self.use_tpu when constructing the
    model specific Estimator object.

    Estimators are also sometimes called classifiers.

    Args:
      batch_size: the batch size to use (for TRAIN, EVAL, and PREDICT modes).
      model_dir: an (optional) string directory to use as the model directory.
      max_checkpoints_to_keep: an (optional) integer count of saved checkpoints.
      iterations_per_loop: an (optional) integer count of log_step_count_steps.
      params: an (optional) dictionary of parameters to pass to the Estimator
        constructor.
      unused_device_fn: a device_fn to pass to RunConfig, if not use_tpu.
      master: a string necessary for TPU, pass FLAGS.master through.
      use_tpu: boolean.  set self.use_tpu if not None.
      start_from_checkpoint: string. If not None, initialize model from this
        path. According to the current implementation of Estimator, this will
        only be used in training. The inference checkpoint is loaded in a
        different place.
      session_config: a tf.ConfigProto to pass to RunConfig, if not use_tpu.
      include_debug_info: from call_variants. If True, PREDICT mode will include
        extra info such as logits and prelogits.

    Returns:
      an object implementing the tf.estimator.Estimator interface (will be a
      TPUEstimator if self.use_tpu is True).
    """
    if use_tpu is not None:
      self.use_tpu = use_tpu

    self.include_debug_info = include_debug_info

    # Set the model dir of this class to the model_dir passed in here. It's not
    # so clean but it appears to be necessary due to the way estimators are
    # constructed (i.e., model_dir is set late).
    self.model_dir = model_dir

    # These flags are exclusive if not None, and 0 means disable.
    save_checkpoints_secs = None
    save_checkpoints_steps = None
    if FLAGS.save_interval_secs >= 0:
      save_checkpoints_secs = FLAGS.save_interval_secs
    if FLAGS.save_interval_steps >= 0:
      save_checkpoints_steps = FLAGS.save_interval_steps

    params = params if params is not None else {}
    warm_start_from = self._create_warm_start_settings(start_from_checkpoint)
    if self.use_tpu:
      tpu_cfg=tpu_config.TPUConfig(
          iterations_per_loop=iterations_per_loop)
      config = tpu_config.RunConfig(
          master=master,
          evaluation_master=master,
          model_dir=model_dir,
          log_step_count_steps=iterations_per_loop,
          keep_checkpoint_max=max_checkpoints_to_keep,
          save_checkpoints_secs=save_checkpoints_secs,
          save_checkpoints_steps=save_checkpoints_steps,
          save_summary_steps=FLAGS.save_summary_steps,
          tpu_config=tpu_cfg,
      )

      classifier = tpu_estimator.TPUEstimator(
          use_tpu=self.use_tpu,
          model_fn=self.model_fn,
          config=config,
          # TODO: enable setting these independently.
          train_batch_size=batch_size,
          eval_batch_size=batch_size,
          predict_batch_size=batch_size,
          params=params,
          warm_start_from=warm_start_from,
      )
    else:
      config = tf_estimator.RunConfig(
          model_dir=model_dir,
          log_step_count_steps=iterations_per_loop,
          keep_checkpoint_max=max_checkpoints_to_keep,
          # device_fn=device_fn,  # Not in tf1.8?
          save_checkpoints_secs=save_checkpoints_secs,
          save_checkpoints_steps=save_checkpoints_steps,
          save_summary_steps=FLAGS.save_summary_steps,
          session_config=session_config,
      )
      # The TPUEstimator interface implicitly adds batch_size to the params
      # dict. Do so explicitly here, so that we can use the same model_fn.
      params_with_batch_size = {'batch_size': batch_size}
      params_with_batch_size.update(params)

      classifier = tf_estimator.Estimator(
          model_fn=self.model_fn,
          config=config,
          params=params_with_batch_size,
          warm_start_from=warm_start_from,
      )
    return classifier

  def model_fn(self, features, labels, mode, params):
    """A model_fn satisfying the Estimator API.

    Args:
      features: a dictionary supplying features.
      labels: a tensor of labels.
      mode: one of tf.estimator.ModeKeys.{EVAL,TRAIN}
      params: a dictionary of parameters.

    Returns:
      a tf.estimator.EstimatorSpec or tpu_estimator.TPUEstimatorSpec,
      depending on self.use_tpu.
    """
    raise NotImplementedError

  def session_eval_hooks(self):
    """Returns a list of tf.train.SessionRunHook classes.

    A typical use case is to provide a hook to load the EMA variables.

    These will be instantiated and invoked by
      eval_hooks = [
          h(model_dir) for h in model.session_eval_hooks()
      ]
      estimator.evaluate(hooks=...).

    Note that this is done according to the instructions in
    cloud_tpu/models/inception/inception_v3.py. A newer idea is in
    tpuestimator-scaffold, but we haven't tried that approach.
    """
    return []

  def session_predict_hooks(self):
    """Returns a list of tf.train.SessionRunHook classes.

    A typical use case is to provide a hook to load the EMA variables.

    These will be instantiated and invoked by
      predict_hooks = [
          h(checkpoint_path) for h in model.session_predict_hooks()
      ]
      estimator.predict(hooks=...).

    Note that this is done according to the instructions in
    cloud_tpu/models/inception/inception_v3.py. A newer idea is in
    tpuestimator-scaffold, but we haven't tried that approach.
    """
    return []

  def create(self, images, num_classes, is_training):
    """Creates a new model.

    Args:
      images: A 4-D tensor of (batch_size, height, width, channels) of pileup
        images.
      num_classes: integer. How many prediction classes are we expecting in
        model?
      is_training: boolean. Should we setup model for training (True) or for
        inference (False).

    Returns:
      A dictionary, containing string keys mapped to endpoint tensors of this
      model. The dictionary must contain a key 'Predictions' that contains the
      probability of having each of 'num_classes' classes.
    """
    try:
      return self._create(images, num_classes, is_training)
    except (ValueError, tf.errors.OpError) as e:
      if self._is_bad_image_dimension_exception(e):
        _, height, width, _ = images.get_shape().as_list()
        message = (
            'Unsupported image dimensions detected: model {} was given images '
            'of w={} x h={} but a TensorFlow exception occurred while building '
            'the model, which typically indicates those dimensions are not '
            'supported by the model. The supported dimensions for {} are {}'
        ).format(
            self.name,
            width,
            height,
            self.name,
            self.supported_dimensions_message,
        )
        raise UnsupportedImageDimensionsError(message)
      else:
        raise

  def _is_bad_image_dimension_exception(self, exception):
    return any(
        x in str(exception) for x in ['Negative dimension', 'SpatialSqueeze']
    )

  def _create(self, images, num_classes, is_training):
    """To be overloaded by subclasses to actually create the model."""
    raise NotImplementedError

  def preprocess_images(self, images):
    """Preprocessing steps needed for this model to process a batch of images.

    Args:
      images: A (batch_size, height, width, channels) 4-D Tensor of type uint8.

    Returns:
      A new batch of images, potentially with different dimensions, based on the
      input but transformed as necessary to use with this model.
    """
    raise NotImplementedError

  @property
  def is_trainable(self):
    """Returns True if this model can be trained."""
    return True

  # TODO: Add export to save representation suitable for inference.

  def __str__(self):
    return 'DeepVariantModel(name={})'.format(self.name)

  def variables_to_restore_from_model(self, exclude_scopes=None):
    """Gets the list of model variables that should be restored.

    The primary use of this function is to get a subset of tf.Variables from a
    slim-defined model that we'd like to restore from a checkpoint. The
    checkpoint generally contains all of the variables in the graph during
    training, including things like the backprop variables, moving averages for
    visualization, etc. Simply restoring all of those variables is brittle, as
    we often want to start a new training run, maybe using a different
    optimizer, different visualization variables, or replacing part of the model
    with a new classification layer, as unneeded variables from the checkpoint
    get loaded into the graph and/or new TF variables not present in the graph
    cannot be found, raising exceptions. This function allows a clean API to get
    just the *model* variables from a graph, excluding all of those non-model
    variables, along with optionally removing parts of the model graph via
    exclude scopes.

    This function calls slim.get_model_variables() to get the raw list of all
    variables associated with the MODEL_VARIABLES collection. It then filters
    away all variables that match any of the scopes in exclude_scopes. For
    example, suppose we have a model with three variables with names:

      w1 = model/l1/weight1
      w2 = model/l2/weight2
      w3 = model/l2/weight3

    Without any exclude scopes, we would return these three variables [w1, w2,
    and w3]. Providing exclude_scopes=['model/l2'] would return only [w1], while
    exclude_scopes=['model/l1'] would return [w2, w3].

    Args:
      exclude_scopes: None, or a list of strings. Each string is a scope
        specification, such as "model/l1" to match all variables whose name
        starts with "model/l1".

    Returns:
      A list of tf.Variable objects.
    """
    vars_to_include = slim.get_model_variables()
    # We aren't excluding any variables, so just return vars_to_include.
    if not exclude_scopes:
      return vars_to_include

    vars_to_exclude = set()
    for scope in exclude_scopes:
      vars_to_exclude |= set(slim.get_variables(scope))
    return [v for v in vars_to_include if v not in vars_to_exclude]


class DeepVariantSlimModel(DeepVariantModel):
  """Baseclass for DeepVariant models based on Slim networks."""

  def __init__(
      self,
      name,
      pretrained_model_path,
      n_classes_model_variable,
      excluded_scopes_for_incompatible_classes,
      excluded_scopes_for_incompatible_channels,
  ):
    """Creates an DeepVariant CNN network based on a tf.slim model.

    Args:
      name: see baseclass.
      pretrained_model_path: see baseclass.
      n_classes_model_variable: str. A fully-qualitified TF variable name in the
        model that we can use to determine the shape of the output
        classification layer of the model. For example, in inception-v3 from
        slim this is 'InceptionV3/Logits/Conv2d_1c_1x1/weights'.
      excluded_scopes_for_incompatible_classes: set of str. A set of scopes that
        will be excluded when restoring from a checkpoint to avoid loading
        incompatible #classes.
      excluded_scopes_for_incompatible_channels: set of str. A set of scopes
        that will be excluded when restoring from a checkpoint to avoid loading
        incompatible #channels.

    Raises:
      ValueError: If any of the arguments are invalid.
    """
    super(DeepVariantSlimModel, self).__init__(
        name=name, pretrained_model_path=pretrained_model_path
    )
    self.n_classes_model_variable = n_classes_model_variable
    self.excluded_scopes_for_incompatible_classes = (
        excluded_scopes_for_incompatible_classes
    )
    self.excluded_scopes_for_incompatible_channels = (
        excluded_scopes_for_incompatible_channels
    )

  def preprocess_images(self, images):
    """Applies preprocessing operations for Inception images.

    Because this will run in model_fn, on the accelerator, we use operations
    that efficiently execute there.

    Args:
      images: An Tensor of shape [batch_size height, width, channel] with uint8
        values.

    Returns:
      A tensor of images of shape [batch_size height, width, channel]
      containing floating point values, with all points rescaled between
      -1 and 1 and possibly resized.
    """
    images = tf.cast(images, dtype=tf.float32)
    images = tf.subtract(images, 128.0)
    images = tf.math.divide(images, 128.0)
    return images

  def model_fn(self, features, labels, mode, params):
    """A model_fn for slim (really inception_v3), satisfying the Estimator API.

    Args:
      features: a single Tensor or dict of same (from input_fn).
      labels: a single Tensor or dict of same (from input_fn).
      mode: tf.estimator.ModeKeys.
      params: dict.

    Returns:
      EstimatorSpec or TPUEstimatorSpec depending on self.use_tpu.
    """
    # NB. The basic structure of this started from
    # //third_party/cloud_tpu/models/inception/inception_v3.py

    # TODO: get this from the model.
    num_classes = dv_constants.NUM_CLASSES

    images = features['image']
    images = self.preprocess_images(images)

    endpoints = self.create(
        images=images,
        num_classes=num_classes,
        is_training=mode == tf_estimator.ModeKeys.TRAIN,
    )

    logits = endpoints['Logits']

    predictions = endpoints
    predictions.update({
        'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor'),
    })
    prelogits = endpoints['PreLogits'] if self.include_debug_info else None
    if mode == tf_estimator.ModeKeys.PREDICT:
      return self._model_fn_predict(mode, features, logits, prelogits=prelogits)

    # Compute loss.
    one_hot_labels = tf.one_hot(labels, num_classes, dtype=tf.int32)
    tf.compat.v1.losses.softmax_cross_entropy(
        onehot_labels=one_hot_labels,
        logits=logits,
        weights=1.0,
        label_smoothing=FLAGS.label_smoothing,
    )
    total_loss = tf.compat.v1.losses.get_total_loss(
        add_regularization_losses=True
    )
    return self.make_ops_and_estimator(
        features,
        endpoints,
        labels,
        logits,
        predictions,
        total_loss,
        mode,
        params,
    )

  def make_ops_and_estimator(
      self,
      features,
      endpoints,
      labels,
      logits,
      predictions,
      total_loss,
      mode,
      params,
  ):
    """Make EstimatorSpec for the current model.

    Args:
      features: a single Tensor or dict of same (from input_fn).
      endpoints:  a dictionary, containing string keys mapped to endpoint
        tensors of this model. The dictionary must contain a key 'Predictions'
        that contains the probability of having each of 'num_classes' classes.
      labels: a single Tensor or dict of same (from input_fn).
      logits: a single Tensor with logits
      predictions: A dictionaty that must contain the following keys: 'Logits'
        and 'Predictions'.
      total_loss:  a single Tensor with a loss
      mode: tf.estimator.ModeKeys.
      params: dict.

    Returns:
      EstimatorSpec or TPUEstimatorSpec depending on self.use_tpu.
    """
    # Note, below, one of train_op or eval_metrics will be None, and the other
    # will be populated, depending on mode.

    # There are a lot of arguments here; that's to avoid referencing flags in
    # leaf functions.
    train_op, host_call = self._model_fn_train(
        mode=mode,
        total_loss=total_loss,
        # get() here to be robust when we are in eval mode and batches_per_epoch
        # hasn't been provided. In eval mode, model_fn_train will return without
        # doing anything.
        batches_per_epoch=params.get('batches_per_epoch', None),
        num_epochs_per_decay=FLAGS.num_epochs_per_decay,
        initial_learning_rate=FLAGS.learning_rate,
        learning_rate_decay_factor=FLAGS.learning_rate_decay_factor,
        rmsprop_decay=FLAGS.rmsprop_decay,
        rmsprop_momentum=FLAGS.rmsprop_momentum,
        rmsprop_epsilon=FLAGS.rmsprop_epsilon,
        moving_average_decay=FLAGS.moving_average_decay,
    )

    eval_metrics = self._model_fn_eval(
        mode=mode,
        features=features,
        labels=labels,
        endpoints=endpoints,
        logits=logits,
        use_logits=False,
    )

    spec = tpu_estimator.TPUEstimatorSpec(
        mode=mode,
        loss=total_loss,
        train_op=train_op,
        host_call=host_call,
        eval_metrics=eval_metrics,
        predictions=predictions,
    )
    if self.use_tpu:
      return spec
    else:
      return spec.as_estimator_spec()

  def _model_fn_predict(self, mode, features, logits, prelogits=None):
    """This is the PREDICT part of model_fn."""
    assert mode == tf_estimator.ModeKeys.PREDICT
    predictions = {
        # We don't actually use classes downstream right now.
        # 'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor'),
        # DV2 call_variants wants these passed through.
        'variant': features['variant'],
        'alt_allele_indices': features['alt_allele_indices'],
    }
    if self.include_debug_info:
      if logits is not None:
        predictions.update({'logits': logits})
      if prelogits is not None:
        predictions.update({'prelogits': prelogits})
    if 'label' in features:
      predictions['label'] = features['label']
    if 'locus' in features:
      predictions['locus'] = features['locus']
    if self.use_tpu:
      return tpu_estimator.TPUEstimatorSpec(mode=mode, predictions=predictions)
    else:
      return tf_estimator.EstimatorSpec(mode=mode, predictions=predictions)

  def _model_fn_eval(
      self, mode, features, labels, endpoints, logits, use_logits
  ):
    """This is the EVAL part of model_fn."""
    if mode != tf_estimator.ModeKeys.EVAL:
      return None
    if use_logits:
      eval_predictions = logits
    else:
      eval_predictions = endpoints['Predictions']
    variant_type = features['variant_type']
    eval_metrics = (eval_metric_fn, [labels, eval_predictions, variant_type])
    if not self.use_tpu:
      for name, value in eval_metrics[0](*eval_metrics[1]).items():
        tf.compat.v1.summary.scalar(tensor=value, name=name)
    return eval_metrics

  def _model_fn_train(
      self,
      mode,
      total_loss,
      batches_per_epoch,
      num_epochs_per_decay,
      initial_learning_rate,
      learning_rate_decay_factor,
      rmsprop_decay,
      rmsprop_momentum,
      rmsprop_epsilon,
      moving_average_decay,
  ):
    """This is the TRAIN part of model_fn."""
    if mode != tf_estimator.ModeKeys.TRAIN:
      return None, None

    # Configure the learning rate using an exponetial decay.
    global_step = tf.compat.v1.train.get_or_create_global_step()
    current_epoch = tf.cast(global_step, tf.float32) / batches_per_epoch
    decay_steps = int(1.0 * batches_per_epoch * num_epochs_per_decay)

    learning_rate = tf.compat.v1.train.exponential_decay(
        learning_rate=initial_learning_rate,
        global_step=global_step,
        decay_steps=decay_steps,
        decay_rate=learning_rate_decay_factor,
        staircase=True,
    )

    # Set a minimum boundary for the learning rate to be a fixed value of 1e-9.
    # It's common to see these tf.max(...) operations when training inception,
    # with a max of 1e-4 * initial_learning_rate but this makes it hard to
    # explore learning rate schedules that decay quickly or by a lot of each
    # step. Here we just use a very small constant 1e-9 as the minimum value.
    learning_rate = tf.maximum(learning_rate, 1e-9, name='learning_rate')

    optimizer = tf.compat.v1.train.RMSPropOptimizer(
        learning_rate,
        rmsprop_decay,
        momentum=rmsprop_momentum,
        epsilon=rmsprop_epsilon,
    )
    if self.use_tpu:
      optimizer = tpu_optimizer.CrossShardOptimizer(optimizer)
    update_ops = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)
    with tf.control_dependencies(update_ops):
      train_op = optimizer.minimize(total_loss, global_step=global_step)

    # NB. In the inception code this was "tf.trainable_variables()
    # + tf.moving_average_variables()", but we've settled on just
    # tf.model_variables() in the existing production DV2.
    variables_to_average = tf.compat.v1.model_variables()
    variable_averages = tf.train.ExponentialMovingAverage(
        decay=moving_average_decay, num_updates=global_step
    )
    with tf.control_dependencies([train_op]), tf.compat.v1.name_scope(
        'moving_average'
    ):
      train_op = variable_averages.apply(variables_to_average)
    tf.compat.v1.add_to_collection(tf.compat.v1.GraphKeys.UPDATE_OPS, train_op)

    # Compute the current epoch and associated learning rate from global_step.
    metric_dict = {
        'current_epoch': current_epoch,
        'total_loss': total_loss,
        'learning_rate': learning_rate,
    }
    host_call = self.construct_scalar_host_call(
        metric_dict=metric_dict, model_dir=self.model_dir, prefix='training/'
    )

    return train_op, host_call

  def session_eval_hooks(self):
    return [LoadEMAHook]

  def session_predict_hooks(self):
    return [PredictEMAHook]


class DeepVariantInceptionV3(DeepVariantSlimModel):
  """DeepVariant inception_v3 network."""

  def __init__(self):
    """Creates an inception-v3 network for DeepVariant."""
    super(DeepVariantInceptionV3, self).__init__(
        name='inception_v3',
        n_classes_model_variable='InceptionV3/Logits/Conv2d_1c_1x1/weights',
        excluded_scopes_for_incompatible_classes=[
            'InceptionV3/Logits',
            'InceptionV3/Conv2d_1a_3x3',
        ],
        excluded_scopes_for_incompatible_channels=['InceptionV3/Conv2d_1a_3x3'],
        pretrained_model_path=(
            '/namespace/vale-project/models/classification/'
            'imagenet/inception_v3/model.ckpt-9591376'
        ),
    )
    self.supported_dimensions_message = (
        'odd widths between 75-361 and any heights between 75-362'
    )

  def _create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(inception_v3.inception_v3_arg_scope()):
      _, endpoints = inception_v3.inception_v3(
          images, num_classes, create_aux_logits=False, is_training=is_training
      )
      return endpoints


class DeepVariantInceptionV3Embedding(DeepVariantInceptionV3):
  """DeepVariant inception_v3_embedding network."""

  def __init__(self):
    """Creates an inception_v3_embedding network for DeepVariant."""
    super(DeepVariantInceptionV3Embedding, self).__init__()
    self.name = 'inception_v3_embedding'
    # vocab_size should be a number larger than the number of sequencing types
    self.vocab_size = 5
    self.embedding_size = 200
    self.dropout_keep_prob = 0.8

  def _create(self, inputs, num_classes, is_training):
    """Creates a new inception_v3_embedding model.

    Args:
      inputs: A tuple of two elements (images, sequencing_types). images is a
        4-D tensor of (batch_size, height, width, channels) of pileup images.
        sequencing_types is a 1-D tensor of (batch_size) of example sequencing
        types.
      num_classes: integer. How many prediction classes are we expecting in
        model?
      is_training: boolean. Should we setup model for training (True) or for
        inference (False).

    Returns:
      A dictionary, containing string keys mapped to endpoint tensors of this
      model.
    """
    images, sequencing_type = inputs
    endpoints = super(DeepVariantInceptionV3Embedding, self)._create(
        images, num_classes, is_training
    )

    with tf.compat.v1.variable_scope('Embeddings'):
      # Take the graph all the way till PreLogits
      net = endpoints['PreLogits']
      net = slim.flatten(net)
      embeddings = self._create_embeddings(sequencing_type)
      net = tf.concat([net, embeddings], 1)

      endpoints['Embeddings'] = net

    with tf.compat.v1.variable_scope('Logits'):
      if isinstance(net.shape[1], int):
        hidden_size = net.shape[1] // 2
      else:
        hidden_size = net.shape[1].value // 2

      net = slim.fully_connected(net, hidden_size, activation_fn=None)
      # TODO: Explore using ReLU before norm
      net = slim.layer_norm(net, scale=False, activation_fn=tf.nn.relu)
      net = slim.dropout(net, self.dropout_keep_prob, is_training=is_training)
      net = slim.fully_connected(net, num_classes, activation_fn=None)

      endpoints.update({'Logits': net, 'Predictions': tf.nn.softmax(net)})

    return endpoints

  def _create_embeddings(self, indices):
    """Create word embeddings."""
    embeddings = self._embedding_lookup(indices)
    embeddings = slim.fully_connected(
        embeddings, self.embedding_size, activation_fn=None
    )
    return embeddings

  def _embedding_lookup(self, input_ids, word_embedding_name='seq_type_emb'):
    """Looks up words embeddings for id tensor.

    Args:
      input_ids: int64 Tensor of shape [batch_size, ] containing word ids.
      word_embedding_name: string. Name of the embedding table.

    Returns:
      float Tensor of shape [batch_size, embedding_size].
    """
    embedding_table = tf.compat.v1.get_variable(
        name=word_embedding_name,
        shape=[self.vocab_size, self.embedding_size],
        initializer=tf.compat.v1.keras.initializers.VarianceScaling(
            scale=1.0, mode='fan_avg', distribution='uniform'
        ),
        collections=[
            tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES,
            tf.compat.v1.GraphKeys.MODEL_VARIABLES,
            tf.compat.v1.GraphKeys.GLOBAL_VARIABLES,
        ],
    )

    return tf.nn.embedding_lookup(params=embedding_table, ids=input_ids)

  def model_fn(self, features, labels, mode, params):
    """A model_fn for slim, satisfying the Estimator API.

    Args:
      features: a single Tensor or dict of same (from input_fn).
      labels: a single Tensor or dict of same (from input_fn).
      mode: tf.estimator.ModeKeys.
      params: dict.

    Returns:
      EstimatorSpec or TPUEstimatorSpec depending on self.use_tpu.

    Raises:
      ValueError: if FLAGS.seq_type_embedding_size is not positive.
    """
    # NB. The basic structure of this started from
    # //third_party/cloud_tpu/models/inception/inception_v3.py

    # TODO: get this from the model.
    num_classes = dv_constants.NUM_CLASSES

    if FLAGS.seq_type_embedding_size <= 0:
      raise ValueError(
          'Expected seq_type_embedding_size to be a positive number but saw %i '
          'instead.'
          % FLAGS.seq_type_embedding_size
      )
    self.embedding_size = FLAGS.seq_type_embedding_size

    images = features['image']
    images = self.preprocess_images(images)
    sequencing_type = features['sequencing_type']

    endpoints = self.create(
        images=(images, sequencing_type),
        num_classes=num_classes,
        is_training=mode == tf_estimator.ModeKeys.TRAIN,
    )

    logits = endpoints['Logits']

    predictions = endpoints
    predictions.update({
        'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor'),
    })

    prelogits = endpoints['PreLogits'] if self.include_debug_info else None
    if mode == tf_estimator.ModeKeys.PREDICT:
      return self._model_fn_predict(mode, features, logits, prelogits=prelogits)

    # Compute loss.
    one_hot_labels = tf.one_hot(labels, num_classes, dtype=tf.int32)
    tf.compat.v1.losses.softmax_cross_entropy(
        onehot_labels=one_hot_labels,
        logits=logits,
        weights=1.0,
        label_smoothing=FLAGS.label_smoothing,
    )
    total_loss = tf.compat.v1.losses.get_total_loss(
        add_regularization_losses=True
    )

    return self.make_ops_and_estimator(
        features,
        endpoints,
        labels,
        logits,
        predictions,
        total_loss,
        mode,
        params,
    )


class DeepVariantPlaceholderModel(DeepVariantModel):
  """BaseClass for placeholder models that are useful for testing and benchmarking."""

  def __init__(self, name):
    """Creates a Placeholder model."""
    # Note the pretrained model path isn't used but we must return a valid
    # string so here we just return "UNUSED".
    super(DeepVariantPlaceholderModel, self).__init__(
        name=name, pretrained_model_path='UNUSED'
    )

  def preprocess_images(self, images):
    """Preprocess images for placeholder model."""
    # Note these calculations aren't necessary, but they are included here to
    # mimic the data processing pipeline used by inception. We may consider
    # removing them in a future CL, or making them optional, to reduce CPU cost
    # of this model.
    images = tf.cast(images, dtype=tf.float32)
    images = tf.subtract(images, 128.0)
    images = tf.math.divide(images, 128.0)
    return images

  @property
  def is_trainable(self):
    """A placeholder model cannot be trained."""
    return False


class DeepVariantConstantModel(DeepVariantPlaceholderModel):
  """Returns a constant probability distribution for each example."""

  def __init__(self, predictions=None):
    """Creates a constant model.

    Args:
      predictions: list[float]. Values to return for Predictions, which should
        be a floatting point value between 0 and 1 for each class, normalized so
        the sum of the values is 1. Predictions should have dimension
        [num_classes].

    Raises:
      ValueError: if sum(predictions) is not close to 1.
    """
    # Note the pretrained model path isn't used but we must return a valid
    # string so here we just return "UNUSED".
    super(DeepVariantConstantModel, self).__init__(name='constant')
    if predictions is None:
      self.predictions = [0.0, 1.0, 0.0]
    elif math.abs(sum(predictions) - 1) > 1e-6:
      raise ValueError('Sum of predictions should be ~1', predictions)
    else:
      self.predictions = predictions

  @staticmethod
  def _predictions(pred_const, batch_size):
    return {
        'Predictions': tf.reshape(
            tf.tile(pred_const, [batch_size]),
            shape=(batch_size, tf.shape(input=pred_const)[0]),
        )
    }

  def _create(self, images, num_classes, is_training):
    assert num_classes == len(self.predictions)
    batch_size = tf.shape(input=images)[0]
    pred_const = tf.constant(self.predictions)
    return self._predictions(pred_const, batch_size)

  def model_fn(self, features, labels, mode, params):
    """A model_fn for the constant model."""
    if mode == tf_estimator.ModeKeys.PREDICT:
      batch_size = tf.shape(input=features['image'])[0]
      logging.info('actual_batch_size %s', batch_size)
    else:
      batch_size = params['batch_size']
      logging.info('batch_size %s', batch_size)
    pred_const = tf.constant(self.predictions)
    endpoints = self._predictions(pred_const, batch_size)
    encoded_variants = features['variant']
    # For the constant model, which is for testing only, we just set the
    # variant_types to 0s. This is needed because it doesn't work to fetch
    # 'variant_type' from either features or endpoints here. Annoying.
    # variant_types = features['variant_type']     # Fails.
    # variant_types = endpoints['variant_type']    # Fails.
    variant_types = tf.zeros(shape=(batch_size,), dtype=tf.int64)

    if mode == tf_estimator.ModeKeys.PREDICT:
      predictions = {
          'probabilities': endpoints['Predictions'],
          'variant': encoded_variants,
          'alt_allele_indices': features['alt_allele_indices'],
      }
      endpoints.update(predictions)

    if mode == tf_estimator.ModeKeys.EVAL:
      eval_metrics = (
          eval_metric_fn,
          [labels, endpoints['Predictions'], variant_types],
      )
    else:
      eval_metrics = None

    loss = tf.constant(0.0)
    train_op = None

    spec = tpu_estimator.TPUEstimatorSpec(
        mode=mode,
        loss=loss,
        train_op=train_op,
        eval_metrics=eval_metrics,
        predictions=endpoints,
    )
    if self.use_tpu:
      return spec
    else:
      return spec.as_estimator_spec()


class DeepVariantSmallModel(DeepVariantSlimModel):
  """A smaller of version of the DeepVariant model.

  Uses only the first layers of Inception net.
  """

  def __init__(self, representation_layer='Mixed_5d'):
    """Creates an DeepVariant CNN network based on a tf.slim model.

    Args:
      representation_layer: string. The name of the layer from the Inception net
        which will be used as an endpoint.

    Raises:
      ValueError: If any of the arguments are invalid.
    """
    super(DeepVariantSmallModel, self).__init__(
        name='small_inception',
        pretrained_model_path=(
            '/namespace/vale-project/models/classification/'
            'imagenet/inception_v3/model.ckpt-9591376'
        ),
        n_classes_model_variable='InceptionV3/Logits/Conv2d_1c_1x1/weights',
        excluded_scopes_for_incompatible_classes=[
            'InceptionV3/Logits',
            'InceptionV3/Conv2d_1a_3x3',
        ],
        excluded_scopes_for_incompatible_channels=['InceptionV3/Conv2d_1a_3x3'],
    )

    self.representation_layer = representation_layer

  def model_fn(self, features, labels, mode, params):
    """A model_fn for slim (really inception_v3), satisfying the Estimator API.

    Args:
      features: a single Tensor or dict of same (from input_fn).
      labels: a single Tensor or dict of same (from input_fn).
      mode: tf.estimator.ModeKeys.
      params: dict.

    Returns:
      EstimatorSpec or TPUEstimatorSpec depending on self.use_tpu.

    Raises:
      ValueError: If representation_layer was not found in the Inception
      architecture
    """
    # NB. The basic structure of this started from
    # //third_party/cloud_tpu/models/inception/inception_v3.py

    # TODO: get this from the model.
    num_classes = dv_constants.NUM_CLASSES

    images = features['image']
    images = self.preprocess_images(images)

    endpoints = self.create(
        images=images,
        num_classes=num_classes,
        is_training=mode == tf_estimator.ModeKeys.TRAIN,
    )

    if self.representation_layer not in endpoints.keys():
      raise ValueError(
          'Layer {} is not found Inception endpoints.'
          'Available Inception net endpoints: {}'.format(
              self.representation_layer, endpoints.keys()
          )
      )

    mid_layer = endpoints[self.representation_layer]
    # Perform 1x1 convolution similarly to the Inception architecture
    # (see 'Predictions' end points in inception_v3 architecture)

    tower = tf.nn.conv2d(
        mid_layer, 1, [1, 1], stride=1, activation_fn=tf.nn.relu
    )

    batch_size = tower.get_shape()[0].value
    tower = tf.reshape(tower, [batch_size, -1])

    with tf.compat.v1.variable_scope('denselayers'):
      with slim.arg_scope([slim.fully_connected], activation_fn=tf.nn.relu):
        logits = slim.fully_connected(tower, num_classes, scope='Dense')

    predictions = endpoints
    predictions.update({
        'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor'),
        'Logits': logits,
        'Predictions': slim.softmax(logits),
    })

    if mode == tf_estimator.ModeKeys.PREDICT:
      return self._model_fn_predict(mode, features, logits)

    # Compute loss.
    one_hot_labels = tf.one_hot(labels, num_classes, dtype=tf.int32)
    tf.compat.v1.losses.softmax_cross_entropy(
        onehot_labels=one_hot_labels,
        logits=logits,
        weights=1.0,
        label_smoothing=FLAGS.label_smoothing,
    )
    total_loss = tf.compat.v1.losses.get_total_loss(
        add_regularization_losses=True
    )

    return self.make_ops_and_estimator(
        features,
        endpoints,
        labels,
        logits,
        predictions,
        total_loss,
        mode,
        params,
    )

  def _create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(inception_v3.inception_v3_arg_scope()):
      _, endpoints = inception_v3.inception_v3(
          images, num_classes, create_aux_logits=False, is_training=is_training
      )
      return endpoints


# Our list of pre-defined model classes.
_MODEL_CLASSES = [
    DeepVariantSmallModel,
    DeepVariantInceptionV3,
    DeepVariantConstantModel,
    DeepVariantInceptionV3Embedding,
]


def all_models():
  """Gets a list of the all of the known model classes."""
  return list(_MODEL_CLASSES)


def production_models():
  """Gets a list of the models that we test extensively."""
  return [get_model('inception_v3'), get_model('inception_v3_embedding')]


def get_model(model_name, **kwargs):
  """Looks up a DeepVariantModel by name.

  Args:
    model_name: String. Looks for a pre-defined DeepVariantModel with a name
      equal to this model_name string.
    **kwargs: arguments to pass to model constructor.

  Returns:
    A DeepVariantModel instance.

  Raises:
    ValueError: If no model exists with model_name.
  """
  for model_class in _MODEL_CLASSES:
    model = model_class()
    if model_name == model.name:
      return model
  raise ValueError(
      'Unknown model_name {}, options are {}'.format(
          model_name, [model_class().name for model_class in _MODEL_CLASSES]
      )
  )
