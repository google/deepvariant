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
import queue
import threading
import tensorflow as tf

try:
  from openvino.inference_engine import IECore  # pylint: disable=g-import-not-at-top
  from openvino.inference_engine import StatusCode  # pylint: disable=g-import-not-at-top
except ImportError:
  pass


class OpenVINOEstimator(object):
  """An estimator for OpenVINO."""

  def __init__(self, model_xml, model_bin, *, input_fn, model):
    self.ie = IECore()
    net = self.ie.read_network(model_xml, model_bin)
    net.input_info['input'].precision = 'U8'
    self.exec_net = self.ie.load_network(
        net,
        'CPU',
        config={'CPU_THROUGHPUT_STREAMS': 'CPU_THROUGHPUT_AUTO'},
        num_requests=0)
    self.tf_sess = tf.compat.v1.Session()
    self.input_fn = input_fn
    self.model = model
    self.outputs = {}
    self.results = queue.Queue()
    self.process_thread = threading.Thread(target=self._process)
    self.features = tf.compat.v1.data.make_one_shot_iterator(
        self.input_fn({'batch_size': 64})).get_next()

  def _process(self):
    """Read input data."""
    self.images = self.features['image']
    self.variant = self.features['variant']
    self.alt_allele_indices = self.features['alt_allele_indices']

    try:
      # List that maps infer requests to index of processed input.
      # -1 means that request has not been started yet.
      infer_request_input_id = [-1] * len(self.exec_net.requests)

      inp_id = 0
      iter_id = 0
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
          self.outputs[inp_id] = {
              'probabilities': None,
              'variant': variant[i],
              'alt_allele_indices': alt_allele_indices[i]
          }
          inp_id += 1
          request.async_infer({'input': inp[i:i + 1].transpose(0, 3, 1, 2)})

          while self.outputs:
            if self.outputs[iter_id]['probabilities'] is not None:
              self.results.put(self.outputs.pop(iter_id))
              iter_id += 1
            else:
              break

    except (StopIteration, tf.errors.OutOfRangeError):
      # Copy rest of outputs
      status = self.exec_net.wait()
      if status != StatusCode.OK:
        raise Exception('Wait for idle request failed!')
      for infer_request_id, out_id in enumerate(infer_request_input_id):
        if not self.outputs[out_id]['probabilities']:
          request = self.exec_net.requests[infer_request_id]
          res = self.outputs[out_id]
          res['probabilities'] = request.output_blobs[
              'InceptionV3/Predictions/Softmax'].buffer.reshape(-1)
          self.results.put(res)

  def __iter__(self):
    self.process_thread.start()
    return self

  def __next__(self):
    if self.process_thread.isAlive() or not self.results.empty():
      return self.results.get()
    else:
      raise StopIteration
