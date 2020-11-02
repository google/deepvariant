import argparse

import tf_slim as slim
import tensorflow as tf
print(tf.__version__)

from tensorflow.python.tools import optimize_for_inference_lib
from deepvariant.modeling import get_model

parser = argparse.ArgumentParser(description='Script to create a frozen graph for DeepVariant models')
parser.add_argument('--checkpoint', required=True, help='Path to model.ckpt')
parser.add_argument('--output', required=True, help='Path to output .pb file')
parser.add_argument('--moving_average_decay', default=0.9999, help='The decay to use for the moving average')
parser.add_argument('--channels', default=6, type=int, help='Number of channels in input tensor')
args = parser.parse_args()

model = get_model('inception_v3')

out_node = 'InceptionV3/Predictions/Reshape_1'
in_node = 'input'

# https://github.com/google/deepvariant/blob/r1.0/deepvariant/dv_constants.py
inp = tf.compat.v1.placeholder(shape=[1, 100, 221, args.channels], dtype=tf.float32, name=in_node)
b = model.create(inp, num_classes=3, is_training=False)

ema = tf.train.ExponentialMovingAverage(args.moving_average_decay)
variables_to_restore = ema.variables_to_restore()

load_ema = slim.assign_from_checkpoint_fn(
    args.checkpoint,
    variables_to_restore,
    ignore_missing_vars=True)

with tf.compat.v1.Session() as sess:
    sess.run(tf.compat.v1.global_variables_initializer())
    load_ema(sess)

    graph_def = sess.graph.as_graph_def()
    graph_def = tf.compat.v1.graph_util.convert_variables_to_constants(sess, graph_def, [out_node])
    graph_def = optimize_for_inference_lib.optimize_for_inference(graph_def, [in_node], [out_node], tf.float32.as_datatype_enum)

    with tf.io.gfile.GFile(args.output, 'wb') as f:
        f.write(graph_def.SerializeToString())
