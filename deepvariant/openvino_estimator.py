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
import sys
from absl import logging
import tensorflow as tf
import tf_slim as slim

try:
  from openvino.inference_engine import IECore  # pylint: disable=g-import-not-at-top
  from openvino.inference_engine import StatusCode  # pylint: disable=g-import-not-at-top
  import mo_tf  # pylint: disable=g-import-not-at-top
  # redacted
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
        sys.executable, mo_tf.__file__, '--input_model=' + openvino_model_pb,
        '--scale=128', '--mean_values',
        '[{}]'.format(','.join(['128'] * tensor_shape[-1]))
    ]
    if openvino_model_dir:
      mo_tf_args.append('--output_dir=' + openvino_model_dir)
    subprocess.run(mo_tf_args, check=True)
    os.remove(openvino_model_pb)

    self.ie = IECore()
    net = self.ie.read_network(
        os.path.join(openvino_model_dir, 'model.xml'),
        os.path.join(openvino_model_dir, 'model.bin'))
    net.input_info['input'].precision = 'U8'
    self.exec_net = self.ie.load_network(
        net,
        'CPU',
        config={'CPU_THROUGHPUT_STREAMS': 'CPU_THROUGHPUT_AUTO'},
        num_requests=0)
    self.tf_sess = tf.compat.v1.Session()
    self.input_fn = input_fn
    self.model = model
    self.outputs = []

  def __iter__(self):
    """Read input data."""
    features = tf.compat.v1.data.make_one_shot_iterator(
        self.input_fn({'batch_size': 64})).get_next()
    self.images = features['image']
    self.variant = features['variant']
    self.alt_allele_indices = features['alt_allele_indices']
    self.iter_id = 0
    logging.info('Using OpenVINO in call_variants.')
    try:
      # List that maps infer requests to index of processed input.
      # -1 means that request has not been started yet.
      infer_request_input_id = [-1] * len(self.exec_net.requests)

      inp_id = 0
      while True:
        # Get next input
        inp, variant, alt_allele_indices = self.tf_sess.run(
            [self.images, self.variant, self.alt_allele_indices])

        for i in range(inp.shape[0]):
          # Get idle infer request
          infer_request_id = self.exec_net.get_idle_request_id()
          if infer_request_id < 0:
            status = self.exec_net.wait(num_requests=1)
            if status != StatusCode.OK:
              raise Exception('Wait for idle request failed!')
            infer_request_id = self.exec_net.get_idle_request_id()
            if infer_request_id < 0:
              raise Exception('Invalid request id!')

          out_id = infer_request_input_id[infer_request_id]
          request = self.exec_net.requests[infer_request_id]

          # Copy output prediction
          if out_id != -1:
            self.outputs[out_id]['probabilities'] = request.output_blobs[
                'InceptionV3/Predictions/Softmax'].buffer.reshape(-1)

          # Start this request on new data
          infer_request_input_id[infer_request_id] = inp_id
          inp_id += 1
          logging.log_every_n(
              logging.INFO,
              ('Processed %s examples in call_variants (using OpenVINO)'),
              _LOG_EVERY_N, inp_id)
          self.outputs.append({
              'probabilities': None,
              'variant': variant[i],
              'alt_allele_indices': alt_allele_indices[i]
          })

          request.async_infer({'input': inp[i:i + 1].transpose(0, 3, 1, 2)})

    except (StopIteration, tf.errors.OutOfRangeError):
      # Copy rest of outputs
      status = self.exec_net.wait()
      if status != StatusCode.OK:
        raise Exception('Wait for idle request failed!')
      for infer_request_id, out_id in enumerate(infer_request_input_id):
        if not self.outputs[out_id]['probabilities']:
          request = self.exec_net.requests[infer_request_id]
          self.outputs[out_id]['probabilities'] = request.output_blobs[
              'InceptionV3/Predictions/Softmax'].buffer.reshape(-1)
    return self

  def __next__(self):
    if self.iter_id < len(self.outputs):
      res = self.outputs[self.iter_id]
      self.iter_id += 1
      return res
    else:
      raise StopIteration
