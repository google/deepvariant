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
"""Code for calling variants with a trained DeepVariant model."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time



from tensorflow import flags
import numpy as np
import tensorflow as tf

from absl import logging

from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import io_utils
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import variant_utils
from deepvariant import logging_level
from deepvariant import modeling
from deepvariant import tf_utils
from deepvariant.protos import deepvariant_pb2

_ALLOW_EXECUTION_HARDWARE = [
    'auto',  # Default, no validation.
    'cpu',  # Don't use accelerators, even if available.
    'accelerator',  # Must be hardware acceleration or an error will be raised.
]

# The number of digits past the decimal point that genotype likelihoods are
# rounded to, for numerical stability.
_GL_PRECISION = 10

FLAGS = flags.FLAGS

flags.DEFINE_string(
    'examples', None,
    'Required. tf.Example protos containing DeepVariant candidate variants in '
    'TFRecord format, as emitted by make_examples.')
flags.DEFINE_string(
    'outfile', None,
    'Required. Destination path where we will write output candidate variants '
    'with additional likelihood information in TFRecord format of '
    'CallVariantsOutput protos.')
flags.DEFINE_string(
    'checkpoint', None,
    'Required. Path to the TensorFlow model checkpoint to use to evaluate '
    'candidate variant calls.')
flags.DEFINE_integer(
    'batch_size', 512,
    'Number of candidate variant tensors to batch together during inference. '
    'Larger batches use more memory but are more computational efficient.')
flags.DEFINE_integer('max_batches', None,
                     'Max. batches to evaluate. Defaults to all.')
flags.DEFINE_integer('num_readers', 8,
                     'Number of parallel readers to create for examples.')
flags.DEFINE_string('model_name', 'inception_v3',
                    'The name of the model architecture of --checkpoint.')
flags.DEFINE_boolean('include_debug_info', False,
                     'If true, include extra debug info in the output.')
flags.DEFINE_string(
    'execution_hardware', 'auto',
    'When in cpu mode, call_variants will not place any ops on the GPU, even '
    'if one is available. In accelerator mode call_variants validates that at '
    'least some hardware accelerator (GPU/TPU) was available for us. This '
    'option is primarily for QA purposes to allow users to validate their '
    'accelerator environment is correctly configured. In auto mode, the '
    'default, op placement is entirely left up to TensorFlow')


class ExecutionHardwareError(Exception):
  pass


def prepare_inputs(source_path, model, batch_size, num_readers=None):
  """Prepares image and encoded_variant ops.

  Reads image / encoded_variant tuples from source_path, extracting the image
  and encoded_variant tensors from source_path. The image is decoded from its
  png encoding and preprocessed with model.preprocess_image as well. Every
  example in source_path is read once (num_epoch=1).

  Args:
    source_path: Path to a TFRecord file containing deepvariant tf.Example
      protos.
    model: A DeepVariantModel whose preprocess_image function will be used on
      image.
    batch_size: int > 0. Size of batches to use during inference.
    num_readers: int > 0 or None. Number of parallel readers to use to read
      examples from source_path. If None, uses FLAGS.num_readers instead.

  Returns:
    A tuple of (image, encoded_variant, encoded_alt_allele_indices) TF ops.
    Image is a [height, width, channel] tensor.
    encoded_variants is a tf.string tensor containing a serialized Variant proto
    describing the variant call associated with image.
    encoded_alt_allele_indices is a tf.string tensor containing a serialized
    CallVariantsOutput.AltAlleleIndices proto containing the
    alternate alleles indices used as "alt" when constructing the image.
  """
  if not num_readers:
    num_readers = FLAGS.num_readers

  tensor_shape = tf_utils.get_shape_from_examples_path(source_path)

  def _parse_single_example(serialized_example):
    """Parses serialized example into a dictionary of de-serialized features."""
    features = tf.parse_single_example(
        serialized_example,
        features={
            'image/encoded': tf.FixedLenFeature([], tf.string),
            'variant/encoded': tf.FixedLenFeature([], tf.string),
            # deepvariant_pb2.CallVariantsOutput.AltAlleleIndices
            'alt_allele_indices/encoded': tf.FixedLenFeature([], tf.string),
        })
    return features

  with tf.name_scope('input'):

    def _preprocess_image(features):
      """Preprocess images (decode, reshape, and apply model-specific steps)."""
      image = features['image/encoded']
      # Bypassing the reshaping and preprocessing if there is no tensor_shape.
      # Currently that could happen when the input file is empty.
      if tensor_shape:
        image = tf.reshape(tf.decode_raw(image, tf.uint8), tensor_shape)
        image = model.preprocess_image(image)
      features['image/encoded'] = image
      return features

    files = tf.gfile.Glob(io_utils.NormalizeToShardedFilePattern(source_path))
    reader_options = io_utils.make_tfrecord_options(files)
    if reader_options.compression_type == (
        tf.python_io.TFRecordCompressionType.GZIP):
      compression_type = 'GZIP'
    else:
      compression_type = None
    dataset = tf.data.TFRecordDataset(files, compression_type=compression_type)
    dataset = dataset.map(
        _parse_single_example, num_parallel_calls=FLAGS.num_readers)
    dataset = dataset.map(
        _preprocess_image, num_parallel_calls=FLAGS.num_readers)
    dataset = dataset.prefetch(10 * batch_size)
    dataset = dataset.batch(batch_size)
    iterator = dataset.make_one_shot_iterator()
    features = iterator.get_next()
    return (features['image/encoded'], features['variant/encoded'],
            features['alt_allele_indices/encoded'])


def round_gls(gls, precision=None):
  """Returns genotype likelihoods rounded to the desired precision level.

  Args:
    gls: A list of floats. The input genotype likelihoods at any precision.
    precision: Positive int. The number of places past the decimal point to
      round to. If None, no rounding is performed.

  Returns:
    A list of floats rounded to the desired precision.

  Raises:
    ValueError: The input gls do not sum to nearly 1.
  """
  if abs(sum(gls) - 1) > 1e-6:
    raise ValueError(
        'Invalid genotype likelihoods do not sum to one: sum({}) = {}'.format(
            gls, sum(gls)))
  if precision is None:
    return gls

  min_ix = 0
  min_gl = gls[0]
  for ix, gl in enumerate(gls):
    if gl < min_gl:
      min_gl = gl
      min_ix = ix

  rounded_gls = [round(gl, precision) for gl in gls]
  rounded_gls[min_ix] = max(
      0.0,
      round(1 - sum(rounded_gls[:min_ix] + rounded_gls[min_ix + 1:]),
            precision))
  return rounded_gls


def call_batch(sess, writer, encoded_variants, encoded_alt_allele_indices,
               predictions):
  """Calls variants by computing the genotype likelihoods predictions.

  This function runs TF to get the values for the predictions for each
  encoded_variants. The resulting genotype likelihoods incorporated into the
  decoded versions of encoded_variants, populating the QUAL, genotype, PL, and
  GQ fields.

  NOTE: This function isn't yet fully operational or complete. We don't
  actually update Variant, we don't handle multi-allelic variants, etc. The
  API will need to change as we incorporate these improvements.

  Args:
    sess: A TensorFlow session where we can evaluate encoded_variants and
      predictions.
    writer: A object with a write() function that will be called for each
      encoded_variant and genotype likelihoods.
    encoded_variants: A [batch_size, 1] tensor of strings. Each string should be
      a third_party.nucleus.protos.Variant encoded protobuf.
    encoded_alt_allele_indices: a tf.string tensor containing a serialized
      CallVariantsOutput.AltAlleleIndices proto containing the
      alternate alleles indices used as "alt" when constructing the image.
    predictions: A [batch_size, 3] tensor of floats. These are the predicted
      genotype likelihoods (p00, p0x, pxx) for some alt allele x, in the same
      order as encoded_variants.

  Returns:
    The number of variants called and written out.
  """
  # After we run the input encoded_variants and predictions they are now
  # concrete objects, not promises, so we can start worthing their values.
  encoded_variants, encoded_alt_allele_indices, predictions = sess.run(
      [encoded_variants, encoded_alt_allele_indices, predictions])

  # Walk over the variants / prediction pairs, creating our calls.
  for encoded_variant, one_encoded_alt_allele_indices, gls in zip(
      encoded_variants, encoded_alt_allele_indices, predictions):
    # Round the genotype probabilities to a stable precision level.
    rounded_gls = round_gls(gls, precision=_GL_PRECISION)
    cvo = _create_cvo_proto(encoded_variant, rounded_gls,
                            one_encoded_alt_allele_indices)
    writer.write(cvo.SerializeToString())

  return len(encoded_variants)


def _create_cvo_proto(encoded_variant, gls, encoded_alt_allele_indices):
  """Returns a CallVariantsOutput proto from the relevant input information."""
  variant = variants_pb2.Variant.FromString(encoded_variant)
  alt_allele_indices = (
      deepvariant_pb2.CallVariantsOutput.AltAlleleIndices.FromString(
          encoded_alt_allele_indices))
  debug_info = None
  if FLAGS.include_debug_info:
    debug_info = deepvariant_pb2.CallVariantsOutput.DebugInfo(
        has_insertion=variant_utils.has_insertion(variant),
        has_deletion=variant_utils.has_deletion(variant),
        is_snp=variant_utils.is_snp(variant),
        predicted_label=np.argmax(gls))
  call_variants_output = deepvariant_pb2.CallVariantsOutput(
      variant=variant,
      alt_allele_indices=alt_allele_indices,
      genotype_probabilities=gls,
      debug_info=debug_info)
  return call_variants_output


def call_variants(examples_filename,
                  checkpoint_path,
                  model,
                  output_file,
                  execution_hardware='auto',
                  batch_size=16,
                  max_batches=None):
  """Main driver of call_variants."""
  # Read a single TFExample to make sure we're not loading an older version.
  example_format = tf_utils.get_format_from_examples_path(examples_filename)
  if example_format is None:
    logging.warning('Unable to read any records from %s. Output will contain '
                    'zero records.', examples_filename)
    io_utils.write_tfrecords([], output_file)
    return
  elif example_format != 'raw':
    raise ValueError('The TF examples in {} has image/format \'{}\' '
                     '(expected \'raw\') which means you might need to rerun '
                     'make_examples to generate the examples again.'.format(
                         examples_filename, example_format))

  if execution_hardware not in _ALLOW_EXECUTION_HARDWARE:
    raise ValueError(
        'Unexpected execution_hardware={} value. Allowed values are {}'.format(
            execution_hardware, ','.join(_ALLOW_EXECUTION_HARDWARE)))

  with tf.Graph().as_default():
    images, encoded_variants, encoded_alt_allele_indices = prepare_inputs(
        examples_filename, model, batch_size, FLAGS.num_readers)

    # Create our model and extract the predictions from the model endpoints.
    predictions = model.create(images, 3, is_training=False)['Predictions']

    # The op for initializing the variables.
    init_op = tf.group(tf.global_variables_initializer(),
                       tf.local_variables_initializer())

    device_count = {'GPU': 0, 'TPU': 0} if execution_hardware == 'cpu' else {}
    config = tf.ConfigProto(device_count=device_count)
    with tf.Session(config=config) as sess:
      sess.run(init_op)

      # Initial the model from the provided checkpoint using our session.
      logging.info('Initializing model from %s', checkpoint_path)
      model.initialize_from_checkpoint(checkpoint_path, 3, False)(sess)

      if execution_hardware == 'accelerator':
        if not any(dev.device_type != 'CPU' for dev in sess.list_devices()):
          raise ExecutionHardwareError(
              'execution_hardware is set to accelerator, but no accelerator '
              'was found')

      logging.info('Writing calls to %s', output_file)
      writer, _ = io_utils.make_proto_writer(output_file)
      with writer:
        start_time = time.time()
        try:
          n_batches = 0
          n_examples = 0
          while max_batches is None or n_batches < max_batches:
            n_called = call_batch(sess, writer, encoded_variants,
                                  encoded_alt_allele_indices, predictions)

            duration = time.time() - start_time
            n_batches += 1
            n_examples += n_called
            logging.info(
                ('Processed %s examples in %s batches [%.2f sec per 100]'),
                n_examples, n_batches, (100 * duration) / n_examples)
        except tf.errors.OutOfRangeError:
          logging.info('Done evaluating variants')


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: call_variants does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)), errors.CommandLineError)
    del argv  # Unused.
    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()

    model = modeling.get_model(FLAGS.model_name)
    call_variants(
        examples_filename=FLAGS.examples,
        checkpoint_path=FLAGS.checkpoint,
        model=model,
        execution_hardware=FLAGS.execution_hardware,
        output_file=FLAGS.outfile,
        max_batches=FLAGS.max_batches,
        batch_size=FLAGS.batch_size)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'outfile',
      'checkpoint',
  ])
  tf.app.run()
