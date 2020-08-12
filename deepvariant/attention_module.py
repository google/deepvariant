# Copyright 2020 Google LLC.
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
"""Contains the attention function and attach_attention_module.

The implementation of squeeze-and-excitation (SE) block is based on the paper
https://arxiv.org/abs/1709.01507.

A use case of attach_attention_module and se_block can be found in:
attention_inception_v3.py. The modified architecture is only used in
DeepVariantAttentionInceptionV3 in modeling.py.
"""

from typing import Tuple

import tensorflow.compat.v1 as tf
import tf_slim as slim


def attach_attention_module(net: tf.Tensor,
                            attention_module: str,
                            end_point: str = None) -> Tuple[str, tf.Tensor]:
  """Attaches attention module to InceptionV3.

  Args:
    net: tf.Tensor corresponding to end_point.
    attention_module: str defining type of attention_module.
    end_point: str defining intermediate layer of model.

  Returns:
    end_point: str defining intermediate layer of model after attaching
    attention_module.
    net: tf.Tensor corresponding to new end_point.

  Raises:
    Exception: The attention_module is not defined.
  """

  if attention_module == 'se_block':
    se_block_scope = 'se_block' if end_point is None else end_point + '_SE'
    end_point = se_block_scope
    net = se_block(net, se_block_scope)
  else:
    raise Exception(
        "'{}' is not a supported attention module".format(attention_module))

  return end_point, net


def se_block(input_feature: tf.Tensor, name: str, ratio: int = 8) -> tf.Tensor:
  """Implementation of Squeeze-and-Excitation (SE) block as described in https://arxiv.org/abs/1709.01507.

  Args:
    input_feature: tf.Tensor to SE block.
    name: str defining name of SE block.
    ratio: int defining size of the bottleneck layer.

  Returns:
    output: tf.Tensor after feature recalibation using SE block.
  """

  kernel_initializer = tf.variance_scaling_initializer()
  bias_initializer = tf.constant_initializer(value=0.0)

  with tf.variable_scope(name):
    channel = input_feature.get_shape()[-1]

    # Spatial Squeeze
    squeeze = tf.reduce_mean(input_feature, axis=[1, 2], keepdims=False)
    assert squeeze.get_shape()[1:] == (channel)

    # Excitation
    excitation = slim.fully_connected(
        inputs=squeeze,
        num_outputs=int(channel // ratio),
        activation_fn=tf.nn.relu,
        weights_initializer=kernel_initializer,
        biases_initializer=bias_initializer,
        scope='bottleneck_fc')
    assert excitation.get_shape()[1:] == (channel // ratio)
    excitation = slim.fully_connected(
        inputs=excitation,
        num_outputs=int(channel),
        activation_fn=tf.nn.sigmoid,
        weights_initializer=kernel_initializer,
        biases_initializer=bias_initializer,
        scope='recover_fc')
    assert excitation.get_shape()[1:] == (channel)

    excitation = tf.expand_dims(excitation, axis=1)
    excitation = tf.expand_dims(excitation, axis=1)
    assert excitation.get_shape()[1:] == (1, 1, channel)

    output = input_feature * excitation
  return output
