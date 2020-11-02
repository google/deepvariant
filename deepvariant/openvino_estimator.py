import tensorflow as tf

try:
  from openvino.inference_engine import IECore, StatusCode
except:
  pass

class OpenVINOEstimator(object):
  def __init__(self, model_xml, model_bin, *, input_fn, model):
    self.ie = IECore()
    net = self.ie.read_network(model_xml, model_bin)
    self.exec_net = self.ie.load_network(net, 'CPU',
                                         config={'CPU_THROUGHPUT_STREAMS': 'CPU_THROUGHPUT_AUTO'},
                                         num_requests=0)
    self.tf_sess = tf.compat.v1.Session()
    self.input_fn = input_fn
    self.model = model
    self.outputs = []


  def __iter__(self):
    # Read input data
    features = tf.compat.v1.data.make_one_shot_iterator(
        self.input_fn({'batch_size': 1})).get_next()
    images = features['image']
    self.images = self.model.preprocess_images(images)
    self.variant = features['variant']
    self.alt_allele_indices = features['alt_allele_indices']
    self.iter_id = 0

    try:
      # List that maps infer requests to index of processed input.
      # -1 means that request has not been started yet.
      infer_request_input_id = [-1] * len(self.exec_net.requests)

      inp_id = 0
      while True:
        # Get next input
        inp, variant, alt_allele_indices = self.tf_sess.run([self.images, self.variant, self.alt_allele_indices])

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
        inp_id += 1
        self.outputs.append({
          'probabilities': None,
          'variant': variant[0],
          'alt_allele_indices': alt_allele_indices[0]
        })
        request.async_infer({'input': inp.transpose(0, 3, 1, 2)})

    except (StopIteration, tf.errors.OutOfRangeError):
      # Copy rest of outputs
      status = self.exec_net.wait()
      if status != StatusCode.OK:
          raise Exception("Wait for idle request failed!")
      for infer_request_id, out_id in enumerate(infer_request_input_id):
        if not self.outputs[out_id]['probabilities']:
          request = self.exec_net.requests[infer_request_id]
          self.outputs[out_id]['probabilities'] = request.output_blobs['InceptionV3/Predictions/Softmax'].buffer.reshape(-1)

    return self


  def __next__(self):
    if self.iter_id < len(self.outputs):
      res = self.outputs[self.iter_id]
      self.iter_id += 1
      return res
    else:
      raise StopIteration
