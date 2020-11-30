import os
import subprocess
import tensorflow as tf
import tf_slim as slim
from tensorflow.python.tools import optimize_for_inference_lib
from threading import Thread
from queue import Queue

try:
  from openvino.inference_engine import IECore, StatusCode
  import mo_tf
  openvino_available = True
except:
  openvino_available = False
  pass


def freeze_graph(model, checkpoint_path, num_channels, moving_average_decay=0.9999):
  out_node = 'InceptionV3/Predictions/Reshape_1'
  in_node = 'input'

  # https://github.com/google/deepvariant/blob/r1.0/deepvariant/dv_constants.py
  inp = tf.compat.v1.placeholder(shape=[1, 100, 221, num_channels], dtype=tf.float32, name=in_node)
  b = model.create(inp, num_classes=3, is_training=False)

  ema = tf.train.ExponentialMovingAverage(moving_average_decay)
  variables_to_restore = ema.variables_to_restore()

  load_ema = slim.assign_from_checkpoint_fn(
      checkpoint_path,
      variables_to_restore,
      ignore_missing_vars=True)

  with tf.compat.v1.Session() as sess:
      sess.run(tf.compat.v1.global_variables_initializer())
      load_ema(sess)

      graph_def = sess.graph.as_graph_def()
      graph_def = tf.compat.v1.graph_util.convert_variables_to_constants(sess, graph_def, [out_node])
      graph_def = optimize_for_inference_lib.optimize_for_inference(graph_def, [in_node], [out_node], tf.float32.as_datatype_enum)

      with tf.io.gfile.GFile('model.pb', 'wb') as f:
          f.write(graph_def.SerializeToString())


class OpenVINOEstimator(object):
  @staticmethod
  def is_available():
    return openvino_available


  def __init__(self, checkpoint_path, num_channels, *, input_fn, model):
    freeze_graph(model, checkpoint_path, num_channels)
    subprocess.run([mo_tf.__file__, '--input_model=model.pb', '--scale=128',
                    '--mean_values', '[{}]'.format(','.join(['128'] * num_channels))])
    os.remove('model.pb')

    self.ie = IECore()
    net = self.ie.read_network('model.xml', 'model.bin')
    net.input_info['input'].precision = 'U8'
    self.exec_net = self.ie.load_network(net, 'CPU',
                                         config={'CPU_THROUGHPUT_STREAMS': 'CPU_THROUGHPUT_AUTO'},
                                         num_requests=0)
    self.tf_sess = tf.compat.v1.Session()
    self.input_fn = input_fn
    self.model = model
    self.outputs = {}
    self.results = Queue()
    self.process_thread = Thread(target=self._process)
    self.features = tf.compat.v1.data.make_one_shot_iterator(self.input_fn({'batch_size': 64})).get_next()


  def _process(self):
    # Read input data
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
        inp, variant, alt_allele_indices = self.tf_sess.run([self.images, self.variant, self.alt_allele_indices])

        for i in range(inp.shape[0]):
          # Get idle infer request
          infer_request_id = self.exec_net.get_idle_request_id()
          if infer_request_id < 0:
            status = self.exec_net.wait(num_requests=1)
            if status != StatusCode.OK:
                raise Exception("Wait for idle request failed!")
            infer_request_id = self.exec_net.get_idle_request_id()
            if infer_request_id < 0:
                raise Exception("Invalid request id!")

          out_id = infer_request_input_id[infer_request_id]
          request = self.exec_net.requests[infer_request_id]

          # Copy output prediction
          if out_id != -1:
            self.outputs[out_id]['probabilities'] = request.output_blobs['InceptionV3/Predictions/Softmax'].buffer.reshape(-1)

          # Start this request on new data
          infer_request_input_id[infer_request_id] = inp_id
          self.outputs[inp_id] = {
            'probabilities': None,
            'variant': variant[i],
            'alt_allele_indices': alt_allele_indices[i]
          }
          inp_id += 1
          request.async_infer({'input': inp[i:i+1].transpose(0, 3, 1, 2)})

          while self.outputs:
            if not self.outputs[iter_id]['probabilities'] is None:
              self.results.put(self.outputs.pop(iter_id))
              iter_id += 1
            else:
              break


    except (StopIteration, tf.errors.OutOfRangeError):
      # Copy rest of outputs
      status = self.exec_net.wait()
      if status != StatusCode.OK:
          raise Exception("Wait for idle request failed!")
      for infer_request_id, out_id in enumerate(infer_request_input_id):
        if not self.outputs[out_id]['probabilities']:
          request = self.exec_net.requests[infer_request_id]
          res = self.outputs[out_id]
          res['probabilities'] = request.output_blobs['InceptionV3/Predictions/Softmax'].buffer.reshape(-1)
          self.results.put(res)


  def __iter__(self):
    self.process_thread.start()
    return self


  def __next__(self):
    if self.process_thread.isAlive() or not self.results.empty():
      return self.results.get()
    else:
      raise StopIteration
