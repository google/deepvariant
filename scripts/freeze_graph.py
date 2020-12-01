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

import argparse
from deepvariant.modeling import get_model
import tensorflow as tf
from tensorflow.python.tools import optimize_for_inference_lib
import tf_slim as slim
print(tf.__version__)

parser = argparse.ArgumentParser(
    description='Script to create a frozen graph for DeepVariant models')
parser.add_argument('--checkpoint', required=True, help='Path to model.ckpt')
parser.add_argument('--output', required=True, help='Path to output .pb file')
parser.add_argument(
    '--moving_average_decay',
    default=0.9999,
    help='The decay to use for the moving average')
parser.add_argument(
    '--channels',
    default=6,
    type=int,
    help='Number of channels in input tensor')
parser.add_argument(
    '--width', default=221, type=int, help='Width of the input tensor')
parser.add_argument(
    '--height', default=100, type=int, help='Height of the input tensor')
args = parser.parse_args()

model = get_model('inception_v3')

out_node = 'InceptionV3/Predictions/Reshape_1'
in_node = 'input'

inp = tf.compat.v1.placeholder(
    shape=[1, args.height, args.width, args.channels],
    dtype=tf.float32,
    name=in_node)
b = model.create(inp, num_classes=3, is_training=False)

ema = tf.train.ExponentialMovingAverage(args.moving_average_decay)
variables_to_restore = ema.variables_to_restore()

load_ema = slim.assign_from_checkpoint_fn(
    args.checkpoint, variables_to_restore, ignore_missing_vars=True)

with tf.compat.v1.Session() as sess:
  sess.run(tf.compat.v1.global_variables_initializer())
  load_ema(sess)

  graph_def = sess.graph.as_graph_def()
  graph_def = tf.compat.v1.graph_util.convert_variables_to_constants(
      sess, graph_def, [out_node])
  graph_def = optimize_for_inference_lib.optimize_for_inference(
      graph_def, [in_node], [out_node], tf.float32.as_datatype_enum)

  with tf.io.gfile.GFile(args.output, 'wb') as f:
    f.write(graph_def.SerializeToString())
