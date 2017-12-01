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



import tensorflow as tf
from absl import logging
from nets import inception
from nets import mobilenet_v1
from nets import resnet_v2

from deepvariant import tf_utils

# The decay factor to use for the moving average.
MOVING_AVERAGE_DECAY = 0.9999

GENOTYPING_ERROR_RATE = 0.0001  # 0.01%

SKIP_MODEL_INITIALIZATION_IN_TEST = '__SKIP_MODEL_INITIALIZATION_IN_TEST__'

slim = tf.contrib.slim


class DeepVariantModel(object):
  """Base class for models that compute genotype likelihoods from an image.

  This class is intended for use anywhere in DeepVariant where we want to train
  or evaluate a model that computes genotype likelihoods from a pileup image. A
  bit of encapsulation helps us to try new models (beyond inception_v3) and unit
  test our code.

  The base class cannot be used directly; concrete subclasses actually implement
  specific models and all of the associated machinery to create/load/save
  models.
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
    raise NotImplementedError

  def preprocess_image(self, image):
    """Add preprocessing steps needed for this model to process image.

    Args:
      image: A (height, width, channels) 3-D Tensor of type uint8.

    Returns:
      A new image, potentially with different dimensions, based on image but
      transformed as necessary to use with this model.
    """
    raise NotImplementedError

  def initialize_from_checkpoint(self, checkpoint_path, num_classes,
                                 is_training):
    """Creates init_fn that loads a model from checkpoint in checkpoint_path.

    Load a pretrained model checkpoint from checkpoint_path. If the number of
    classes defined in that checked model is different from the number
    in our dataset, subclasses may partially load this model, but this is only
    allowed if is_training=False, otherwise an error is raised.

    Args:
      checkpoint_path: String. Path to a checkpoint.
      num_classes: The number of classes we are going to predict with this
        model. Must be integer > 1.
      is_training: boolean. Is this model being loaded for training purposes or
        for inference?

    Returns:
      An init_fn for use with slim.learning.train init_fn argument, which is
      defined as "a callable to be executed after `init_op` is called. The
      callable must accept one argument, the session being initialized."

    Raises:
      ValueError: If the checkpoint_path could not be loaded.
      ValueError: If the model endpoints in the latest checkpoint in directory
        aren't perfectly compatible with this model and the num_classes we are
        prediction and is_training=False.
    """
    raise NotImplementedError

  @property
  def is_trainable(self):
    """Returns True if this model can be trained."""
    return True

  def loss(self, endpoints, one_hot_labels):
    """Creates a loss function to train this model against one_hot_labels.

    Args:
      endpoints: Dictionary. The endpoints of this model, as returned by the
        create function call.
      one_hot_labels: One-hot encoded truth labels that we want to train this
        model to predict.

    Returns:
      A `Tensor` whose value represents the total loss.
    """
    raise NotImplementedError

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
               excluded_scopes):
    """Creates an DeepVariant CNN network based on a tf.slim model.

    Args:
      name: see baseclass.
      pretrained_model_path: see baseclass.
      n_classes_model_variable: str. A fully-qualitified TF variable name in the
        model that we can use to determine the shape of the output
        classification layer of the model. For example, in inception-v3 from
        slim this is 'InceptionV3/Logits/Conv2d_1c_1x1/weights'.
      excluded_scopes: list of str. A list of scopes that will be excluded when
        restoring from a checkpoint.

    Raises:
      ValueError: If any of the arguments are invalid.
    """
    super(DeepVariantSlimModel, self).__init__(
        name=name, pretrained_model_path=pretrained_model_path)
    if not excluded_scopes:
      raise ValueError('Got an empty value for excluded_scopes',
                       excluded_scopes)
    self.n_classes_model_variable = n_classes_model_variable
    self.excluded_scopes = excluded_scopes

  def preprocess_image(self, image):
    """Applies tf-slim preprocessing operations to image.

    Args:
      image: An single image Tensor of shape [height, width, channel] with
             uint8 values.

    Returns:
      A single image of shape [height, width, channel] containing floating point
      values, with all points rescaled between -1 and 1 and rescaled via
      resize_image_with_crop_or_pad.
    """
    image = tf.to_float(image)
    image = tf.subtract(image, 128.0)
    image = tf.div(image, 128.0)
    # redacted
    # Specifically, our current image height is 100, which is too small so that
    # in inception_v3, this final pool:
    # net = slim.avg_pool2d(net, kernel_size, padding='VALID',
    #                       scope='AvgPool_1a_{}x{}'.format(*kernel_size))
    # actually end up having a kernel_size dimension as 1, but the default
    # stride of that function is (2, 2). So, for inception_v3, we had to resize
    # the height to at least 107 if the height is smaller than that.
    if self.name == 'inception_v3':
      h, w, _ = image.get_shape().as_list()
      image = tf.image.resize_image_with_crop_or_pad(image, max(h, 107), w)
    return image

  def initialize_from_checkpoint(self, checkpoint_path, num_classes,
                                 is_training):
    """See baseclass."""
    # Note that this class supports loading from checkpoint with a different
    # number of classes than requested here, in which case the softmax endpoints
    # for the model are not loaded.
    if checkpoint_path is None:
      raise ValueError('Checkpoint cannot be None')

    if checkpoint_path == SKIP_MODEL_INITIALIZATION_IN_TEST:
      # If the checkpoint path is our magic SKIP_MODEL_INITIALIZATION_IN_TEST we
      # don't in fact do any initialization.
      return lambda sess: sess

    # Figure out how many classes this inception model was trained to predict.
    shapes = tf_utils.model_shapes(checkpoint_path,
                                   [self.n_classes_model_variable])
    # The last dimension of this tensor is the num of classes.
    checkpoint_num_classes = shapes[self.n_classes_model_variable][-1]

    exclude_scopes = []
    if checkpoint_num_classes != num_classes:
      # We have a mismatch between the number of classes the model was trained
      # with and how many we want to predict. This is fine as long as we are
      # training, as we simply don't restore the logits variables from the
      # pre-trained model and start from random weights for those values.
      if not is_training:
        raise ValueError('Checkpoint has {} classes but we want to use {} and '
                         'is_training=False'.format(checkpoint_num_classes,
                                                    num_classes))
      logging.info('Checkpoint was trained against %s '
                   'classes while our dataset is using %s, '
                   'enabling fine-tuning', checkpoint_num_classes, num_classes)
      exclude_scopes = self.excluded_scopes

    variables_to_restore = self.variables_to_restore_from_model(exclude_scopes)
    if not is_training:
      # Apply the moving averages to the variables we want to restore during
      # inference.
      variable_averages = tf.train.ExponentialMovingAverage(
          MOVING_AVERAGE_DECAY)
      for var in variables_to_restore:
        tf.add_to_collection(tf.GraphKeys.MOVING_AVERAGE_VARIABLES, var)
      variables_to_restore = variable_averages.variables_to_restore()

    # If it's training mode, we ignore any missing vars because it's ok to be
    # less strict about a starting checkpoint.
    # However, in inference, we want to make sure we still fail on missing vars
    # unless we explicity exclude them.
    return slim.assign_from_checkpoint_fn(
        checkpoint_path, variables_to_restore, ignore_missing_vars=is_training)

  def loss(self, endpoints, one_hot_labels):
    """See baseclass."""
    slim.losses.softmax_cross_entropy(
        endpoints['Logits'],
        one_hot_labels,
        label_smoothing=GENOTYPING_ERROR_RATE,
        weights=1.0)
    return slim.losses.get_total_loss()


class DeepVariantInceptionV3(DeepVariantSlimModel):
  """DeepVariant inception_v3 network."""

  def __init__(self):
    """Creates an inception-v3 network for DeepVariant."""
    super(DeepVariantInceptionV3, self).__init__(
        name='inception_v3',
        n_classes_model_variable='InceptionV3/Logits/Conv2d_1c_1x1/weights',
        excluded_scopes=['InceptionV3/Logits', 'InceptionV3/Conv2d_1a_3x3'],
        pretrained_model_path=('/namespace/vale-project/models/classification/'
                               'imagenet/inception_v3/model.ckpt-9591376'))

  def create(self, images, num_classes, is_training):
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
        excluded_scopes=['InceptionV2/Logits', 'InceptionV2/Conv2d_1a_7x7'],
        pretrained_model_path=('/namespace/vale-project/models/classification/'
                               'imagenet/inception_v2/model.ckpt-14284043'))

  def create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(inception.inception_v2_arg_scope()):
      _, endpoints = inception.inception_v2(
          inputs=images, num_classes=num_classes, is_training=is_training)
      return endpoints


class DeepVariantResnet50(DeepVariantSlimModel):
  """Resnet v2 50 model.

  References:
    See slim resnet_v2.py.
  """

  def __init__(self):
    super(DeepVariantResnet50, self).__init__(
        name='resnet_v2_50',
        n_classes_model_variable='resnet_v2_50/logits/weights',
        excluded_scopes=['resnet_v2_50/logits', 'resnet_v2_50/conv1'],
        pretrained_model_path=('/namespace/vale-project/models/classification/'
                               'imagenet/resnet_v2_50_inception_preprocessed/'
                               'model.ckpt-5136169'))

  def create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(resnet_v2.resnet_arg_scope()):
      _, endpoints = resnet_v2.resnet_v2_50(
          images, num_classes, is_training=is_training, spatial_squeeze=False)
      # Resnet's "predictions" and endpoint is (n, 1, 1, m) but we really
      # want to have an (n, m) "Predictions" endpoint.  We add a squeeze
      # op here to make that happen.
      endpoints['Predictions'] = tf.squeeze(
          endpoints['predictions'], [1, 2], name='SqueezePredictions')
      # Likewise, the endpoint "resnet_v2_50/logits" should be squeezed to
      # "Logits"
      endpoints['Logits'] = tf.squeeze(
          endpoints['resnet_v2_50/logits'], [1, 2], name='SqueezeLogits')
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
        excluded_scopes=['MobilenetV1/Logits', 'MobilenetV1/Conv2d_0'],
        pretrained_model_path=('/cns/ok-d/home/howarda/slim/'
                               'mobilenet_asynch_100_224_ds_s5_cr_2.5_50/train/'
                               'model.ckpt-19527265'))

  def create(self, images, num_classes, is_training):
    """See baseclass."""
    with slim.arg_scope(mobilenet_v1.mobilenet_v1_arg_scope()):
      _, endpoints = mobilenet_v1.mobilenet_v1(
          inputs=images, num_classes=num_classes, is_training=is_training)
      return endpoints


class DeepVariantRandomGuessModel(DeepVariantModel):
  """Assigns a random probability to each class.

  This model is mostly useful for testing of DeepVariant, as the evaluation of
  this model is essentially free.
  """

  def __init__(self, seed=1268458594):
    """Creates a RandomGuessing model.

    Args:
      seed: int. The random number seed to use for our tf.random_uniform op.
    """
    # Note the pretrained model path isn't used but we must return a valid
    # string so here we just return "UNUSED".
    super(DeepVariantRandomGuessModel, self).__init__(
        name='random_guess', pretrained_model_path='UNUSED')
    self.seed = seed

  def create(self, images, num_classes, is_training):
    """The Random model emits a random uniform probability for each class."""
    batch_size = tf.shape(images)[0]
    rand_probs = tf.random_uniform(
        shape=(batch_size, num_classes), seed=self.seed)
    return {'Predictions': tf.nn.softmax(rand_probs)}

  def preprocess_image(self, image):
    # Note this calculations aren't necessary, but they are included here to
    # mimic the data processing pipeline used by inception. We may consider
    # removing them in a future CL, or making them optional, to reduce CPU cost
    # of this model.
    image = tf.to_float(image)
    image = tf.subtract(image, 128.0)
    image = tf.div(image, 128.0)
    image = tf.reshape(image, (100, 221, 7))

    return image

  def initialize_from_checkpoint(self, checkpoint_path, num_classes,
                                 is_training):
    # No initialization is needed, so return a noop.
    return lambda sess: sess

  @property
  def is_trainable(self):
    """The RandomGuess model cannot be trained."""
    return False


# Our list of pre-defined models.
_MODELS = [
    DeepVariantInceptionV3(),
    DeepVariantInceptionV2(),
    DeepVariantMobileNetV1(),
    DeepVariantRandomGuessModel(),
    DeepVariantResnet50(),
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
