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
"""Saves model ckpt files to .pb file."""

import json


from absl import flags
from absl import logging
import tensorflow as tf
import tf_slim as slim

from deepvariant import modeling
# pylint: disable=g-direct-tensorflow-import
from tensorflow.python.tools import optimize_for_inference_lib
# pylint: enable=g-direct-tensorflow-import


flags.DEFINE_string('checkpoint', None, 'Required. Path to input model.ckpt.')
flags.DEFINE_string(
    'example_info_json', None, 'Path to one *example_info.json file containing '
    'the information of the channels for the examples.')
flags.DEFINE_string('output', None, 'Required. Path to output .pb file.')
flags.DEFINE_integer('channels', 6, 'Number of channels in input tensor.')
flags.DEFINE_integer('width', 221, 'Width of the input tensor.')
flags.DEFINE_integer('height', 100, 'Height of the input tensor.')

FLAGS = flags.FLAGS


def freeze_graph(model,
                 checkpoint_path,
                 tensor_shape,
                 openvino_model_pb,
                 moving_average_decay=0.9999):
  """Converts model ckpts."""
  out_node = 'InceptionV3/Predictions/Reshape_1'
  in_node = 'input'

  inp = tf.compat.v1.placeholder(
      shape=[1] + tensor_shape, dtype=tf.float32, name=in_node)
  _ = model.create(inp, num_classes=3, is_training=False)

  ema = tf.train.ExponentialMovingAverage(moving_average_decay)
  variables_to_restore = ema.variables_to_restore()

  load_ema = slim.assign_from_checkpoint_fn(
      checkpoint_path, variables_to_restore, ignore_missing_vars=True)

  with tf.compat.v1.Session() as sess:
    sess.run(tf.compat.v1.global_variables_initializer())
    load_ema(sess)

    graph_def = sess.graph.as_graph_def()
    graph_def = tf.compat.v1.graph_util.convert_variables_to_constants(
        sess, graph_def, [out_node])
    graph_def = optimize_for_inference_lib.optimize_for_inference(
        graph_def, [in_node], [out_node], tf.float32.as_datatype_enum)

    with tf.io.gfile.GFile(openvino_model_pb, 'wb') as f:
      f.write(graph_def.SerializeToString())


def main(_):
  if FLAGS.example_info_json:
    logging.info('--example_info_json is set. '
                 'Use shape information from --example_info_json')
    example_info = json.load(tf.io.gfile.GFile(FLAGS.example_info_json, 'r'))
    tensor_shape = example_info['shape']
  else:
    tensor_shape = [FLAGS.height, FLAGS.width, FLAGS.channels]

  logging.info('Processing ckpt=%s, tensor_shape=%s.', FLAGS.checkpoint,
               tensor_shape)
  freeze_graph(
      modeling.get_model('inception_v3'), FLAGS.checkpoint, tensor_shape,
      FLAGS.output)
  logging.info('Output written to %s.', FLAGS.output)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'checkpoint',
      'output',
  ])
  tf.compat.v1.app.run()
