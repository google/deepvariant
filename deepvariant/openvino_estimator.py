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
"""An estimator that uses OpenVINO."""
import os
import subprocess
from absl import logging
import tensorflow as tf
import tf_slim as slim

try:
  from openvino.runtime import Core, AsyncInferQueue, Type  # pylint: disable=g-import-not-at-top,g-multiple-import
  from openvino.preprocess import PrePostProcessor  # pylint: disable=g-import-not-at-top
  # TODO:
  # For now, including this in the try/except block as well.
  # Fix this with a correct dep and also add unit test
  from tensorflow.python.tools import optimize_for_inference_lib  # pylint: disable=g-import-not-at-top
except ImportError:
  pass

_LOG_EVERY_N = 15000


def freeze_graph(model,
                 checkpoint_path,
                 tensor_shape,
                 openvino_model_pb,
                 moving_average_decay=0.9999):
  """Converts model ckpts."""
  logging.info('Processing ckpt=%s, tensor_shape=%s', checkpoint_path,
               tensor_shape)
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


class OpenVINOEstimator(object):
  """An estimator for OpenVINO."""

  def __init__(self, checkpoint_path, input_fn, model, openvino_model_dir=''):
    tensor_shape = input_fn.tensor_shape
    openvino_model_pb = os.path.join(openvino_model_dir, 'model.pb')
    freeze_graph(model, checkpoint_path, tensor_shape, openvino_model_pb)
    mo_tf_args = [
        'mo', '--input_model=' + openvino_model_pb, '--scale=128',
        '--mean_values', '[{}]'.format(','.join(['128'] * tensor_shape[-1]))
    ]
    if openvino_model_dir:
      mo_tf_args.append('--output_dir=' + openvino_model_dir)
    subprocess.run(mo_tf_args, check=True)
    os.remove(openvino_model_pb)

    self.ie = Core()
    self.ie.set_property('CPU', {
        'CPU_THROUGHPUT_STREAMS': 'CPU_THROUGHPUT_AUTO',
        'CPU_BIND_THREAD': 'YES'
    })
    net = self.ie.read_model(
        os.path.join(openvino_model_dir, 'model.xml'),
        os.path.join(openvino_model_dir, 'model.bin'))
    self.input_name = net.inputs[0].get_any_name()

    p = PrePostProcessor(net)
    p.input().tensor().set_element_type(Type.u8)
    net = p.build()

    compiled_model = self.ie.compile_model(net, 'CPU')
    num_requests = compiled_model.get_property(
        'OPTIMAL_NUMBER_OF_INFER_REQUESTS')
    logging.info('OpenVINO uses %d inference requests.', num_requests)
    self.infer_queue = AsyncInferQueue(compiled_model, num_requests)
    self.tf_sess = tf.compat.v1.Session()
    self.input_fn = input_fn
    self.model = model

  def __iter__(self):
    """Read input data."""
    features = tf.compat.v1.data.make_one_shot_iterator(
        self.input_fn({'batch_size': 64})).get_next()
    self.images = features['image']
    self.variant = features['variant']
    self.alt_allele_indices = features['alt_allele_indices']
    self.iter_id = 0
    self.outputs = []
    logging.info('Using OpenVINO in call_variants.')
    try:

      def completion_callback(request, inp_id):
        logging.log_every_n(
            logging.INFO,
            ('Processed %s examples in call_variants (using OpenVINO)'),
            _LOG_EVERY_N, inp_id)

        output = next(iter(request.results.values()))
        self.outputs[inp_id - 1]['probabilities'] = output.reshape(-1)

      self.infer_queue.set_callback(completion_callback)

      inp_id = 0
      while True:
        # Get next input
        inp, variant, alt_allele_indices = self.tf_sess.run(
            [self.images, self.variant, self.alt_allele_indices])

        for i in range(inp.shape[0]):
          # Start this request on new data
          inp_id += 1
          self.outputs.append({
              'probabilities': None,
              'variant': variant[i],
              'alt_allele_indices': alt_allele_indices[i]
          })
          self.infer_queue.start_async({self.input_name: inp[i:i + 1]}, inp_id)

    except (StopIteration, tf.errors.OutOfRangeError):
      self.infer_queue.wait_all()

    return self

  def __next__(self):
    if self.iter_id < len(self.outputs):
      res = self.outputs[self.iter_id]
      self.iter_id += 1
      return res
    else:
      raise StopIteration
