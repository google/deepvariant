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
"""Provides an abstraction around deep learning models in DeepVariant.

This class allows us to encapsulate all of the model management, loading,
saving, and data processing in a single place so those details don't spill over
into the more general deepvariant codebase. The key thing we are aiming for here
is to make sure we can easily play with other model architectures without
modifying the surrounding training and evaluation code.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math



from absl import flags
from absl import logging

import tensorflow as tf

from deepvariant import tf_utils
from tensorflow.contrib.tpu.python.tpu import tpu_config
from tensorflow.contrib.tpu.python.tpu import tpu_estimator
from tensorflow.contrib.tpu.python.tpu import tpu_optimizer
from nets import inception
from nets import mobilenet_v1

from deepvariant import dv_constants

flags.DEFINE_float(
    'label_smoothing', 1e-6,
    'Amount of label smoothing to use. By default this is 0.0001% '
    'meaning that we expect a label error at a rate of 1 / 1,000,000')

# Training parameters.
flags.DEFINE_float('learning_rate', 0.001, 'Initial learning rate.')

flags.DEFINE_float('rmsprop_momentum', 0.9, 'Momentum.')

flags.DEFINE_float('rmsprop_decay', 0.9, 'Decay term for RMSProp.')

flags.DEFINE_float('rmsprop_epsilon', 1.0, 'Epsilon term for RMSProp.')

flags.DEFINE_float('learning_rate_decay_factor', 0.94,
                   'Learning rate decay factor.')

flags.DEFINE_float('num_epochs_per_decay', 2.0,
                   'Number of epochs after which learning rate decays.')

flags.DEFINE_float('moving_average_decay', 0.9999,
                   'The decay to use for the moving average.')

flags.DEFINE_integer(
    'save_summary_steps', 100, 'Number of steps which must have run before '
    'showing summaries.')

flags.DEFINE_integer(
    'save_interval_secs', 1000, 'Interval (in seconds) at which the model data '
    'should be checkpointed. Set to 0 to disable, -1 to ignore. '
    'Exclusive with save_interval_steps.')

flags.DEFINE_integer(
    'save_interval_steps', -1, 'Interval (in steps) at which the model data '
    'should be checkpointed. Set to 0 to disable, -1 to ignore. '
    'Exclusive with save_interval_secs.')

FLAGS = flags.FLAGS

slim = tf.contrib.slim


class UnsupportedImageDimensions(Exception):
  """Exception indicating the image dimensions aren't supported by our model."""


def binarize(labels, target_class):
  """Binarize labels and predictions.

  The labels that are not equal to target_class parameter are set to zero.

  Args:
    labels: the ground-truth labels for the examples.
    target_class: index of the class that is left as non-zero.

  Returns:
    Tensor of the same shape as labels.
  """
  labels_binary = tf.where(
      tf.equal(labels, tf.constant(target_class, dtype=tf.int64)),
      tf.zeros_like(labels), labels)
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
  return tf.metrics.recall(labels_binary, predicted_class_binary)


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
  return tf.metrics.precision(labels_binary, predicted_class_binary)


# NB. This includes only a subset of our usual metrics.
# We'll add the rest back in a subsequent change.
def eval_metric_fn(labels, predictions, unused_encoded_variants):
  """Calculate eval metrics from Tensors, on CPU host.

  Args:
    labels: the ground-truth labels for the examples.
    predictions: the predicted labels for the examples.
    unused_encoded_variants: the encoded variants.

  Returns:
    a dictionary of string name to metric.
  """

  predicted_class = tf.argmax(input=predictions, axis=1)

  metrics = {
      'Accuracy/All': tf.metrics.accuracy(labels, predicted_class),
      'Precision/All': tf.metrics.precision(labels, predicted_class),
      'Recall/All': tf.metrics.recall(labels, predicted_class),
      'FPs/All': tf.metrics.false_positives(labels, predicted_class),
      'FNs/All': tf.metrics.false_negatives(labels, predicted_class),
      'TPs/All': tf.metrics.true_positives(labels, predicted_class),
      'TNs/All': tf.metrics.true_negatives(labels, predicted_class),
      'Recall/Class1': get_class_recall(labels, predicted_class, 1),
      'Recall/Class2': get_class_recall(labels, predicted_class, 2),
      'Precision/Class1': get_class_precision(labels, predicted_class, 1),
      'Precision/Class2': get_class_precision(labels, predicted_class, 2),
  }

  return metrics

# The following two classes support loading exponential moving averages into
# their corresponding variables when a checkpoint is loaded. They're called
# as hooks by the Estimators. Note for future work: this is the documented
# way, but someone on the mailing list suggested that using the scaffold_fn
# mechanism might be better.


class LoadEMAHook(tf.train.SessionRunHook):
  """Hook to load EMA into their corresponding variables.

  This looks for the latest checkpoint in the model dir.
  """

  def __init__(self, model_dir):
    super(LoadEMAHook, self).__init__()
    self._model_dir = model_dir

  def begin(self):
    ema = tf.train.ExponentialMovingAverage(FLAGS.moving_average_decay)
    variables_to_restore = ema.variables_to_restore()
    self._load_ema = tf.contrib.framework.assign_from_checkpoint_fn(
        tf.train.latest_checkpoint(self._model_dir), variables_to_restore)

  def after_create_session(self, sess, coord):
    tf.logging.info('Reloading EMA...')
    self._load_ema(sess)


class PredictEMAHook(tf.train.SessionRunHook):
  """Hook to load EMA into their corresponding variables.

  This reads the specified checkpoint.
  """

  def __init__(self, checkpoint_path):
    super(PredictEMAHook, self).__init__()
    self._checkpoint_path = checkpoint_path

  def begin(self):
    ema = tf.train.ExponentialMovingAverage(FLAGS.moving_average_decay)
    variables_to_restore = ema.variables_to_restore()
    self._load_ema = tf.contrib.framework.assign_from_checkpoint_fn(
        self._checkpoint_path, variables_to_restore)

  def after_create_session(self, sess, coord):
    tf.logging.info('Reloading EMA...')
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

  def _create_warm_start_settings(self, start_from_checkpoint):
    """Create a proper WarmStartSettings based on start_from_checkpoint."""
    # If the special value "model_default" was passed, ask the model for
    # its default.
    if start_from_checkpoint == 'model_default':
      start_from_checkpoint = self.pretrained_model_path

    # If the path is non-False, use it.
    if start_from_checkpoint:
      logging.info('Initializing model from checkpoint at %s',
                   start_from_checkpoint)
      reader = tf.train.NewCheckpointReader(start_from_checkpoint)
      var_to_shape_map = reader.get_variable_to_shape_map()
      if tf_utils.model_num_classes(
          start_from_checkpoint,
          self.n_classes_model_variable) == dv_constants.NUM_CLASSES:
        logging.info('The model checkpoint to warm start from has the same '
                     'number of classes. If this is in training, we will '
                     'clear excluded_scopes_for_incompatible_shapes so we '
                     'include everything for '
                     'warm starting....')
        vars_to_include = var_to_shape_map.keys()
      else:
        logging.info(
            'The model checkpoint to warm start from has different '
            'number of classes. If this is in training, we will '
            'use excluded_scopes_for_incompatible_shapes=%s',
            self.excluded_scopes_for_incompatible_shapes)
        vars_to_include = [
            v for v in var_to_shape_map.keys() if not v.startswith(
                tuple(self.excluded_scopes_for_incompatible_shapes))
        ]
      return tf.estimator.WarmStartSettings(
          ckpt_to_initialize_from=start_from_checkpoint,
          vars_to_warm_start='|'.join(vars_to_include))
    else:
      # If warm_start_from is an empty string, specifically set it to None.
      logging.info('Initializing model with random parameters')
      return None

  def make_estimator(self,
                     batch_size,
                     model_dir=None,
                     max_checkpoints_to_keep=100000,
                     iterations_per_loop=100,
                     params=None,
                     unused_device_fn=None,
                     master='',
                     use_tpu=False,
                     start_from_checkpoint=None):
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
      path. According to the current implementation of Estimator, this will only
      be used in training. The inference checkpoint is loaded in a different
      place.

    Returns:
      an object implementing the tf.estimator.Estimator interface (will be a
      TPUEstimator if self.use_tpu is True).
    """
    if use_tpu is not None:
      self.use_tpu = use_tpu

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
      config = tpu_config.RunConfig(
          master=master,
          evaluation_master=master,
          model_dir=model_dir,
          log_step_count_steps=iterations_per_loop,
          keep_checkpoint_max=max_checkpoints_to_keep,
          save_checkpoints_secs=save_checkpoints_secs,
          save_checkpoints_steps=save_checkpoints_steps,
          save_summary_steps=FLAGS.save_summary_steps,
          tpu_config=tpu_config.TPUConfig(
              iterations_per_loop=iterations_per_loop))

      classifier = tpu_estimator.TPUEstimator(
          use_tpu=self.use_tpu,
          model_fn=self.model_fn,
          config=config,
          # redacted
          train_batch_size=batch_size,
          eval_batch_size=batch_size,
          predict_batch_size=batch_size,
          params=params,
          warm_start_from=warm_start_from,
      )
    else:
      config = tf.estimator.RunConfig(
          model_dir=model_dir,
          log_step_count_steps=iterations_per_loop,
          keep_checkpoint_max=max_checkpoints_to_keep,
          # device_fn=device_fn,  # Not in tf1.8?
          save_checkpoints_secs=save_checkpoints_secs,
          save_checkpoints_steps=save_checkpoints_steps,
          save_summary_steps=FLAGS.save_summary_steps,
      )
      # The TPUEstimator interface implicitly adds batch_size to the params
      # dict. Do so explicitly here, so that we can use the same model_fn.
      params_with_batch_size = {'batch_size': batch_size}
      params_with_batch_size.update(params)

      classifier = tf.estimator.Estimator(
          model_fn=self.model_fn,
          config=config,
          params=params_with_batch_size,
          warm_start_from=warm_start_from)
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
    except (ValueError, tf.errors.OpError), e:
      if self._is_bad_image_dimension_exception(e):
        _, height, width, _ = images.get_shape().as_list()
        message = (
            'Unsupported image dimensions detected: model {} was given images '
            'of w={} x h={} but a TensorFlow exception occurred while building '
            'the model, which typically indicates those dimensions are not '
            'supported by the model. The supported dimensions for {} are {}'
        ).format(self.name, width, height, self.name,
                 self.supported_dimensions_message)
        raise UnsupportedImageDimensions(message)
      else:
        raise

  def _is_bad_image_dimension_exception(self, exception):
    return any(
        x in str(exception) for x in ['Negative dimension', 'SpatialSqueeze'])

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

  # redacted

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

  def __init__(self, name, pretrained_model_path, n_classes_model_variable,
               excluded_scopes_for_incompatible_shapes):
    """Creates an DeepVariant CNN network based on a tf.slim model.

    Args:
      name: see baseclass.
      pretrained_model_path: see baseclass.
      n_classes_model_variable: str. A fully-qualitified TF variable name in the
        model that we can use to determine the shape of the output
        classification layer of the model. For example, in inception-v3 from
        slim this is 'InceptionV3/Logits/Conv2d_1c_1x1/weights'.
      excluded_scopes_for_incompatible_shapes: list of str. A list of scopes
        that will be excluded when restoring from a checkpoint to avoid loading
        incompatible shapes.

    Raises:
      ValueError: If any of the arguments are invalid.
    """
    super(DeepVariantSlimModel, self).__init__(
        name=name, pretrained_model_path=pretrained_model_path)
    if not excluded_scopes_for_incompatible_shapes:
      raise ValueError(
          'Got an empty value for '
          'excluded_scopes_for_incompatible_shapes',
          excluded_scopes_for_incompatible_shapes)
    self.n_classes_model_variable = n_classes_model_variable
    self.excluded_scopes_for_incompatible_shapes = (
        excluded_scopes_for_incompatible_shapes)

  def preprocess_images(self, images):
    """Applies preprocessing operations for Inception images.

    Because this will run in model_fn, on the accelerator, we use operations
    that efficiently execute there.

    Args:
      images: An Tensor of shape [batch_size height, width, channel] with
             uint8 values.

    Returns:
      A tensor of images of shape [batch_size height, width, channel]
      containing floating point values, with all points rescaled between
      -1 and 1 and possibly resized.
    """
    images = tf.to_float(images)
    images = tf.subtract(images, 128.0)
    images = tf.div(images, 128.0)
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
    # and //third_party/cloud_tpu/models/mobilenet/mobilenet.py

    # redacted
    num_classes = dv_constants.NUM_CLASSES

    images = features['image']
    images = self.preprocess_images(images)

    endpoints = self.create(
        images=images,
        num_classes=num_classes,
        is_training=mode == tf.estimator.ModeKeys.TRAIN)

    logits = endpoints['Logits']

    predictions = endpoints
    predictions.update({
        'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor')
    })

    if mode == tf.estimator.ModeKeys.PREDICT:
      return self._model_fn_predict(mode, features, logits)

    # Compute loss.
    one_hot_labels = tf.one_hot(labels, num_classes, dtype=tf.int32)
    tf.losses.softmax_cross_entropy(
        onehot_labels=one_hot_labels,
        logits=logits,
        weights=1.0,
        label_smoothing=FLAGS.label_smoothing)
    total_loss = tf.losses.get_total_loss(add_regularization_losses=True)

    return self.make_ops_and_estimator(features, endpoints, labels, logits,
                                       predictions, total_loss, mode, params)

  def make_ops_and_estimator(self, features, endpoints, labels, logits,
                             predictions, total_loss, mode, params):
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
    train_op = self._model_fn_train(
        mode=mode,
        total_loss=total_loss,
        batches_per_epoch=params['batches_per_epoch'],
        num_epochs_per_decay=FLAGS.num_epochs_per_decay,
        initial_learning_rate=FLAGS.learning_rate,
        learning_rate_decay_factor=FLAGS.learning_rate_decay_factor,
        rmsprop_decay=FLAGS.rmsprop_decay,
        rmsprop_momentum=FLAGS.rmsprop_momentum,
        rmsprop_epsilon=FLAGS.rmsprop_epsilon,
        moving_average_decay=FLAGS.moving_average_decay)

    eval_metrics = self._model_fn_eval(
        mode=mode,
        features=features,
        labels=labels,
        endpoints=endpoints,
        logits=logits,
        use_logits=False)

    spec = tpu_estimator.TPUEstimatorSpec(
        mode=mode,
        loss=total_loss,
        train_op=train_op,
        eval_metrics=eval_metrics,
        predictions=predictions)
    if self.use_tpu:
      return spec
    else:
      return spec.as_estimator_spec()

  def _model_fn_predict(self, mode, features, logits):
    """This is the PREDICT part of model_fn."""
    assert mode == tf.estimator.ModeKeys.PREDICT
    predictions = {
        # We don't actually use classes downstream right now.
        # 'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor'),
        # DV2 call_variants wants these passed through.
        'variant': features['variant'],
        'alt_allele_indices': features['alt_allele_indices'],
    }
    if 'label' in features:
      predictions['label'] = features['label']
    if 'locus' in features:
      predictions['locus'] = features['locus']
    if self.use_tpu:
      return tpu_estimator.TPUEstimatorSpec(mode=mode, predictions=predictions)
    else:
      return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)

  def _model_fn_eval(self, mode, features, labels, endpoints, logits,
                     use_logits):
    """This is the EVAL part of model_fn."""
    if mode != tf.estimator.ModeKeys.EVAL:
      return None
    if use_logits:
      eval_predictions = logits
    else:
      eval_predictions = endpoints['Predictions']
    encoded_variants = features['variant']
    eval_metrics = (eval_metric_fn,
                    [labels, eval_predictions, encoded_variants])
    if not self.use_tpu:
      for name, value in eval_metrics[0](*eval_metrics[1]).iteritems():
        tf.contrib.slim.summaries.add_scalar_summary(
            value, name, print_summary=True)
    return eval_metrics

  def _model_fn_train(self, mode, total_loss, batches_per_epoch,
                      num_epochs_per_decay, initial_learning_rate,
                      learning_rate_decay_factor, rmsprop_decay,
                      rmsprop_momentum, rmsprop_epsilon, moving_average_decay):
    """This is the TRAIN part of model_fn."""
    if mode != tf.estimator.ModeKeys.TRAIN:
      return None
    # Configure the learning rate using an exponetial decay.
    global_step = tf.train.get_or_create_global_step()
    decay_steps = int(1.0 * batches_per_epoch * num_epochs_per_decay)
    learning_rate = tf.train.exponential_decay(
        learning_rate=initial_learning_rate,
        global_step=global_step,
        decay_steps=decay_steps,
        decay_rate=learning_rate_decay_factor,
        staircase=True)
    # Set a minimum boundary for the learning rate.
    learning_rate = tf.maximum(
        learning_rate, 0.0001 * initial_learning_rate, name='learning_rate')
    optimizer = tf.train.RMSPropOptimizer(
        learning_rate,
        rmsprop_decay,
        momentum=rmsprop_momentum,
        epsilon=rmsprop_epsilon)
    if self.use_tpu:
      optimizer = tpu_optimizer.CrossShardOptimizer(optimizer)
    update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
    with tf.control_dependencies(update_ops):
      train_op = optimizer.minimize(total_loss, global_step=global_step)

    # NB. In the inception code this was "tf.trainable_variables()
    # + tf.moving_average_variables()", but we've settled on just
    # tf.model_variables() in the existing production DV2.
    variables_to_average = tf.model_variables()
    variable_averages = tf.train.ExponentialMovingAverage(
        decay=moving_average_decay, num_updates=global_step)
    with tf.control_dependencies([train_op]), tf.name_scope('moving_average'):
      train_op = variable_averages.apply(variables_to_average)
    tf.add_to_collection(tf.GraphKeys.UPDATE_OPS, train_op)
    return train_op

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
        excluded_scopes_for_incompatible_shapes=[
            'InceptionV3/Logits', 'InceptionV3/Conv2d_1a_3x3'
        ],
        pretrained_model_path=('/namespace/vale-project/models/classification/'
                               'imagenet/inception_v3/model.ckpt-9591376'))
    self.supported_dimensions_message = (
        'odd widths between 75-361 and any heights between 75-362')

  def _create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(inception.inception_v3_arg_scope()):
      _, endpoints = inception.inception_v3(
          images, num_classes, create_aux_logits=False, is_training=is_training)
      return endpoints


class DeepVariantInceptionV2(DeepVariantSlimModel):
  """Original DeepVariant inception_v2 network."""

  def __init__(self):
    """Creates an inception-v2 network for DeepVariant."""
    super(DeepVariantInceptionV2, self).__init__(
        name='inception_v2',
        n_classes_model_variable='InceptionV2/Logits/Conv2d_1c_1x1/weights',
        excluded_scopes_for_incompatible_shapes=[
            'InceptionV2/Logits', 'InceptionV2/Conv2d_1a_7x7'
        ],
        pretrained_model_path=('/namespace/vale-project/models/classification/'
                               'imagenet/inception_v2/model.ckpt-14284043'))

  def _create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(inception.inception_v2_arg_scope()):
      _, endpoints = inception.inception_v2(
          inputs=images, num_classes=num_classes, is_training=is_training)
      return endpoints


class DeepVariantMobileNetV1(DeepVariantSlimModel):
  """MobileNet v1 model.

  More information about MobileNets can be found in their paper:
  https://arxiv.org/pdf/1704.04861.pdf
  """

  def __init__(self):
    super(DeepVariantMobileNetV1, self).__init__(
        name='mobilenet_v1',
        n_classes_model_variable='MobilenetV1/Logits/Conv2d_1c_1x1/weights',
        excluded_scopes_for_incompatible_shapes=[
            'MobilenetV1/Logits', 'MobilenetV1/Conv2d_0'
        ],
        pretrained_model_path=('/cns/ok-d/home/howarda/slim/'
                               'mobilenet_asynch_100_224_ds_s5_cr_2.5_50/train/'
                               'model.ckpt-19527265'))

  def _create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(mobilenet_v1.mobilenet_v1_arg_scope()):
      _, endpoints = mobilenet_v1.mobilenet_v1(
          inputs=images, num_classes=num_classes, is_training=is_training)
      return endpoints


class DeepVariantDummyModel(DeepVariantModel):
  """BaseClass for dummy models that are useful for testing and benchmarking."""

  def __init__(self, name):
    """Creates a Dummy model."""
    # Note the pretrained model path isn't used but we must return a valid
    # string so here we just return "UNUSED".
    super(DeepVariantDummyModel, self).__init__(
        name=name, pretrained_model_path='UNUSED')

  def preprocess_images(self, images):
    """Preprocess images for dummy model."""
    # Note these calculations aren't necessary, but they are included here to
    # mimic the data processing pipeline used by inception. We may consider
    # removing them in a future CL, or making them optional, to reduce CPU cost
    # of this model.
    images = tf.to_float(images)
    images = tf.subtract(images, 128.0)
    images = tf.div(images, 128.0)
    return images

  @property
  def is_trainable(self):
    """A dummy model cannot be trained."""
    return False


class DeepVariantRandomGuessModel(DeepVariantDummyModel):
  """Assigns a random probability to each class.

  This model is mostly useful for testing of DeepVariant, as the evaluation of
  this model is essentially free.
  """

  def __init__(self, seed=1268458594):
    """Creates a RandomGuessing model.

    Args:
      seed: int. The random number seed to use for our tf.random_uniform op.
    """
    super(DeepVariantRandomGuessModel, self).__init__(name='random_guess')
    self.seed = seed

  def _create(self, images, num_classes, is_training):
    """The Random model emits a random uniform probability for each class."""
    batch_size = tf.shape(images)[0]
    rand_probs = tf.random_uniform(
        shape=(batch_size, num_classes), seed=self.seed)
    return {'Predictions': tf.nn.softmax(rand_probs)}

  def model_fn(self, features, labels, mode, params):
    """A model_fn for the random model."""
    # redacted
    num_classes = dv_constants.NUM_CLASSES

    # In predict-mode the last batch may be smaller; so use the measured
    # batch size.
    encoded_variants = features['variant']

    tf.set_random_seed(self.seed)
    rand_probs = tf.map_fn(
        fn=lambda _: tf.random_uniform([num_classes]),
        elems=features['image'],
        dtype=tf.float32,
    )

    # In a normal model we'd call create to get the endpoints,
    # but it's inconvenient to pass batch size in this case.
    endpoints = {'Predictions': tf.nn.softmax(rand_probs)}

    predictions = endpoints

    if mode == tf.estimator.ModeKeys.PREDICT:
      predictions = {
          'probabilities': endpoints['Predictions'],
          'variant': encoded_variants,
          'alt_allele_indices': features['alt_allele_indices'],
      }

    if mode == tf.estimator.ModeKeys.EVAL:
      eval_metrics = (eval_metric_fn,
                      [labels, endpoints['Predictions'], encoded_variants])
    else:
      eval_metrics = None

    loss = tf.constant(0.0)
    train_op = None
    spec = tpu_estimator.TPUEstimatorSpec(
        mode=mode,
        loss=loss,
        train_op=train_op,
        eval_metrics=eval_metrics,
        predictions=predictions)
    if self.use_tpu:
      return spec
    else:
      return spec.as_estimator_spec()


class DeepVariantConstantModel(DeepVariantDummyModel):
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
        'Predictions':
            tf.reshape(
                tf.tile(pred_const, [batch_size]),
                shape=(batch_size, tf.shape(pred_const)[0]))
    }

  def _create(self, images, num_classes, is_training):
    assert num_classes == len(self.predictions)
    batch_size = tf.shape(images)[0]
    pred_const = tf.constant(self.predictions)
    return self._predictions(pred_const, batch_size)

  def model_fn(self, features, labels, mode, params):
    """A model_fn for the constant model."""
    if mode == tf.estimator.ModeKeys.PREDICT:
      batch_size = tf.shape(features['image'])[0]
      logging.info('actual_batch_size %s', batch_size)
    else:
      batch_size = params['batch_size']
      logging.info('batch_size %s', batch_size)
    pred_const = tf.constant(self.predictions)
    endpoints = self._predictions(pred_const, batch_size)
    encoded_variants = features['variant']

    if mode == tf.estimator.ModeKeys.PREDICT:
      predictions = {
          'probabilities': endpoints['Predictions'],
          'variant': encoded_variants,
          'alt_allele_indices': features['alt_allele_indices'],
      }
      endpoints.update(predictions)

    if mode == tf.estimator.ModeKeys.EVAL:
      eval_metrics = (eval_metric_fn,
                      [labels, endpoints['Predictions'], encoded_variants])
    else:
      eval_metrics = None

    loss = tf.constant(0.0)
    train_op = None

    spec = tpu_estimator.TPUEstimatorSpec(
        mode=mode,
        loss=loss,
        train_op=train_op,
        eval_metrics=eval_metrics,
        predictions=endpoints)
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
        pretrained_model_path=('/namespace/vale-project/models/classification/'
                               'imagenet/inception_v3/model.ckpt-9591376'),
        n_classes_model_variable='InceptionV3/Logits/Conv2d_1c_1x1/weights',
        excluded_scopes_for_incompatible_shapes=[
            'InceptionV3/Logits', 'InceptionV3/Conv2d_1a_3x3'
        ])

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
    # and //third_party/cloud_tpu/models/mobilenet/mobilenet.py

    # redacted
    num_classes = dv_constants.NUM_CLASSES

    images = features['image']
    images = self.preprocess_images(images)

    endpoints = self.create(
        images=images,
        num_classes=num_classes,
        is_training=mode == tf.estimator.ModeKeys.TRAIN)

    if self.representation_layer not in endpoints.keys():
      raise ValueError('Layer {} is not found Inception endpoints.'
                       'Available Inception net endpoints: {}'.format(
                           self.representation_layer, endpoints.keys()))

    mid_layer = endpoints[self.representation_layer]
    # Perform 1x1 convolution similarly to the Inception architecture
    # (see 'Predictions' end points in inception_v3 architecture)

    tower = tf.contrib.layers.conv2d(
        mid_layer, 1, [1, 1], stride=1, activation_fn=tf.nn.relu)

    batch_size = tower.get_shape()[0].value
    tower = tf.reshape(tower, [batch_size, -1])

    with tf.variable_scope('denselayers'):
      with slim.arg_scope([slim.fully_connected], activation_fn=tf.nn.relu):
        logits = slim.fully_connected(tower, num_classes, scope='Dense')

    predictions = endpoints
    predictions.update({
        'classes': tf.argmax(input=logits, axis=1, output_type=tf.int32),
        'probabilities': tf.nn.softmax(logits, name='softmax_tensor'),
        'Logits': logits,
        'Predictions': slim.softmax(logits)
    })

    if mode == tf.estimator.ModeKeys.PREDICT:
      return self._model_fn_predict(mode, features, logits)

    # Compute loss.
    one_hot_labels = tf.one_hot(labels, num_classes, dtype=tf.int32)
    tf.losses.softmax_cross_entropy(
        onehot_labels=one_hot_labels,
        logits=logits,
        weights=1.0,
        label_smoothing=FLAGS.label_smoothing)
    total_loss = tf.losses.get_total_loss(add_regularization_losses=True)

    return self.make_ops_and_estimator(features, endpoints, labels, logits,
                                       predictions, total_loss, mode, params)

  def _create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(inception.inception_v3_arg_scope()):
      _, endpoints = inception.inception_v3(
          images, num_classes, create_aux_logits=False, is_training=is_training)
      return endpoints

# Our list of pre-defined models.
_MODELS = [
    DeepVariantSmallModel(),
    DeepVariantInceptionV3(),
    DeepVariantInceptionV2(),
    DeepVariantMobileNetV1(),
    DeepVariantRandomGuessModel(),
    DeepVariantConstantModel(),
]


def all_models():
  """Gets a list of the all of the known models."""
  return list(_MODELS)


def production_models():
  """Gets a list of the models that we test extensively."""
  return [get_model('inception_v3'), get_model('mobilenet_v1')]


def get_model(model_name):
  """Looks up a DeepVariantModel by name.

  Args:
    model_name: String. Looks for a pre-defined DeepVariantModel with a name
      equal to this model_name string.

  Returns:
    A DeepVariantModel instance.

  Raises:
    ValueError: If no model exists with model_name.
  """
  for model in _MODELS:
    if model_name == model.name:
      return model
  raise ValueError('Unknown model_name {}, options are {}'.format(
      model_name, [model.name for model in _MODELS]))
