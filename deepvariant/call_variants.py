# Copyright 2017 Google LLC.
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

import os
import time



from absl import flags
from absl import logging
import numpy as np
import tensorflow as tf
from tensorflow import estimator as tf_estimator

from deepvariant import data_providers
from deepvariant import dv_utils
from deepvariant import logging_level
from deepvariant import modeling
from deepvariant.protos import deepvariant_pb2
from google.protobuf import text_format
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import variant_utils


tf.compat.v1.disable_eager_execution()

_ALLOW_EXECUTION_HARDWARE = [
    'auto',  # Default, no validation.
    'cpu',  # Don't use accelerators, even if available.
    'accelerator',  # Must be hardware acceleration or an error will be raised.
]

# The number of digits past the decimal point that genotype likelihoods are
# rounded to, for numerical stability.
_GL_PRECISION = 10

# This number is estimated by the following logic:
# For a sample with 10,000,000 examples, if we log every 50,000 examples,
# there will be 200 lines per sample.
_LOG_EVERY_N = 50000

FLAGS = flags.FLAGS

flags.DEFINE_string(
    'examples',
    None,
    (
        'Required. tf.Example protos containing DeepVariant candidate variants'
        ' in TFRecord format, as emitted by make_examples. Can be a'
        ' comma-separated list of files, and the file names can contain'
        ' wildcard characters.'
    ),
)
flags.DEFINE_string(
    'outfile',
    None,
    (
        'Required. Destination path where we will write output candidate'
        ' variants with additional likelihood information in TFRecord format of'
        ' CallVariantsOutput protos.'
    ),
)
flags.DEFINE_string(
    'checkpoint',
    None,
    (
        'Required. Path to the TensorFlow model checkpoint to use to evaluate '
        'candidate variant calls.'
    ),
)
flags.DEFINE_integer(
    'batch_size',
    512,
    (
        'Number of candidate variant tensors to batch together during'
        ' inference. Larger batches use more memory but are more computational'
        ' efficient.'
    ),
)
flags.DEFINE_integer(
    'max_batches', None, 'Max. batches to evaluate. Defaults to all.'
)
flags.DEFINE_integer(
    'num_readers', 8, 'Number of parallel readers to create for examples.'
)
flags.DEFINE_string(
    'model_name',
    'inception_v3',
    'The name of the model architecture of --checkpoint.',
)
flags.DEFINE_boolean(
    'include_debug_info',
    False,
    'If true, include extra debug info in the output.',
)
flags.DEFINE_boolean(
    'debugging_true_label_mode',
    False,
    (
        'If true, read the true labels from examples and add to '
        'output. Note that the program will crash if the input '
        'examples do not have the label field. '
        'When true, this will also fill everything when '
        '--include_debug_info is set to true.'
    ),
)
flags.DEFINE_string(
    'execution_hardware',
    'auto',
    (
        'When in cpu mode, call_variants will not place any ops on the GPU,'
        ' even if one is available. In accelerator mode call_variants validates'
        ' that at least some hardware accelerator (GPU/TPU) was available for'
        ' us. This option is primarily for QA purposes to allow users to'
        ' validate their accelerator environment is correctly configured. In'
        ' auto mode, the default, op placement is entirely left up to'
        ' TensorFlow.  In tpu mode, use and require TPU.'
    ),
)
flags.DEFINE_string(
    'config_string',
    None,
    (
        'String representation of a tf.ConfigProto message, with'
        ' comma-separated key: value pairs, such as "allow_soft_placement:'
        ' True". The value can itself be another message, such as "gpu_options:'
        ' {per_process_gpu_memory_fraction: 0.5}".'
    ),
)

# Cloud TPU Cluster Resolvers
flags.DEFINE_string(
    'gcp_project',
    None,
    (
        'Project name for the Cloud TPU-enabled project. If not specified, we '
        'will attempt to automatically detect the GCE project from metadata.'
    ),
)
flags.DEFINE_string(
    'tpu_zone',
    None,
    (
        'GCE zone where the Cloud TPU is located in. If not specified, we '
        'will attempt to automatically detect the GCE project from metadata.'
    ),
)
flags.DEFINE_string(
    'tpu_name',
    None,
    (  # pylint: disable=line-too-long
        'Name of the Cloud TPU for Cluster Resolvers. You must specify either '
        'this flag or --master. An empty value corresponds to no Cloud TPU.'
        ' See '
        'https://www.tensorflow.org/api_docs/python/tf/distribute/cluster_resolver/TPUClusterResolver'
    ),
)

flags.DEFINE_string(
    'master',
    None,
    (
        'GRPC URL of the master (e.g. grpc://ip.address.of.tpu:8470). You '
        'must specify either this flag or --tpu_name.'
    ),
)

flags.DEFINE_boolean('use_tpu', False, 'Use tpu if available.')

flags.DEFINE_boolean('use_openvino', False, 'Use Intel OpenVINO as backend.')

flags.DEFINE_string(
    'openvino_model_dir',
    '',
    'If set, use this directory to save the temporary model file for OpenVINO.',
)

flags.DEFINE_string(
    'kmp_blocktime',
    '0',
    (
        'Value to set the KMP_BLOCKTIME environment variable to for efficient'
        ' MKL inference. See'
        ' https://www.tensorflow.org/performance/performance_guide for more'
        ' information. The default value is 0, which provides the best'
        ' performance in our tests. Set this flag to "" to not set the'
        ' variable.'
    ),
)


class ExecutionHardwareError(Exception):
  pass


def prepare_inputs(source_path, use_tpu=False, num_readers=None):
  """Return a tf.data input_fn from the source_path.

  Args:
    source_path: Path to a TFRecord file containing deepvariant tf.Example
      protos.
    use_tpu: boolean.  Use the tpu code path.
    num_readers: int > 0 or None. Number of parallel readers to use to read
      examples from source_path. If None, uses FLAGS.num_readers instead.

  Returns:
    A tf input_fn yielding batches of image, encoded_variant,
    encoded_alt_allele_indices.

    The image is a [batch_size, height, width, channel] tensor. The
    encoded_variants is a tf.string or tpu-encoded tensor containing a
    serialized Variant proto describing the variant call associated with
    image. The encoded_alt_allele_indices is a tf.string or tpu-encoded
    tensor containing a serialized CallVariantsOutput.AltAlleleIndices proto
    containing the alternate alleles indices used as "alt" when constructing
    the image.
  """
  if not num_readers:
    num_readers = FLAGS.num_readers

  return data_providers.get_input_fn_from_filespec(
      input_file_spec=source_path,
      mode=tf_estimator.ModeKeys.PREDICT,
      use_tpu=use_tpu,
      input_read_threads=num_readers,
      debugging_true_label_mode=FLAGS.debugging_true_label_mode,
  )


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
            gls, sum(gls)
        )
    )
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
      round(
          1 - sum(rounded_gls[:min_ix] + rounded_gls[min_ix + 1 :]), precision
      ),
  )
  return rounded_gls


def write_variant_call(writer, prediction, use_tpu):
  """Write the variant call based on prediction.

  Args:
    writer: A object with a write() function that will be called for each
      encoded_variant and genotype likelihoods.
    prediction: A [3] tensor of floats. These are the predicted genotype
      likelihoods (p00, p0x, pxx) for some alt allele x, in the same order as
      encoded_variants.
      use_tpu: bool.  Decode the tpu specific encoding of prediction.

  Returns:
    The return status from writer.
  """
  encoded_variant = prediction['variant']
  if use_tpu:
    encoded_variant = dv_utils.int_tensor_to_string(encoded_variant)

  encoded_alt_allele_indices = prediction['alt_allele_indices']
  if use_tpu:
    encoded_alt_allele_indices = dv_utils.int_tensor_to_string(
        encoded_alt_allele_indices
    )

  rounded_gls = round_gls(prediction['probabilities'], precision=_GL_PRECISION)

  # Write it out.
  true_labels = prediction['label'] if FLAGS.debugging_true_label_mode else None
  cvo = _create_cvo_proto(
      encoded_variant,
      rounded_gls,
      encoded_alt_allele_indices,
      true_labels,
      logits=prediction.get('logits'),
      prelogits=prediction.get('prelogits'),
  )
  return writer.write(cvo)


def _create_cvo_proto(
    encoded_variant,
    gls,
    encoded_alt_allele_indices,
    true_labels=None,
    logits=None,
    prelogits=None,
):
  """Returns a CallVariantsOutput proto from the relevant input information."""
  variant = variants_pb2.Variant.FromString(encoded_variant)
  alt_allele_indices = (
      deepvariant_pb2.CallVariantsOutput.AltAlleleIndices.FromString(
          encoded_alt_allele_indices
      )
  )
  debug_info = None
  if FLAGS.include_debug_info or FLAGS.debugging_true_label_mode:
    if prelogits is not None:
      assert prelogits.shape == (1, 1, 2048)
      prelogits = prelogits[0][0]
    debug_info = deepvariant_pb2.CallVariantsOutput.DebugInfo(
        has_insertion=variant_utils.has_insertion(variant),
        has_deletion=variant_utils.has_deletion(variant),
        is_snp=variant_utils.is_snp(variant),
        predicted_label=np.argmax(gls),
        true_label=true_labels,
        logits=logits,
        prelogits=prelogits,
    )
  call_variants_output = deepvariant_pb2.CallVariantsOutput(
      variant=variant,
      alt_allele_indices=alt_allele_indices,
      genotype_probabilities=gls,
      debug_info=debug_info,
  )
  return call_variants_output


def call_variants(
    examples_filename,
    checkpoint_path,
    model,
    output_file,
    execution_hardware='auto',
    batch_size=16,
    max_batches=None,
    use_tpu=False,
    master='',
):
  """Main driver of call_variants."""
  if FLAGS.kmp_blocktime:
    os.environ['KMP_BLOCKTIME'] = FLAGS.kmp_blocktime
    logging.vlog(
        3, 'Set KMP_BLOCKTIME to {}'.format(os.environ['KMP_BLOCKTIME'])
    )

  # Read a single TFExample to make sure we're not loading an older version.
  first_example = dv_utils.get_one_example_from_examples_path(examples_filename)
  if first_example is None:
    logging.warning(
        'Unable to read any records from %s. Output will contain zero records.',
        examples_filename,
    )
    tfrecord.write_tfrecords([], output_file)
    return

  example_info_json = dv_utils.get_example_info_json_filename(
      examples_filename, 0
  )
  example_shape, example_channels_enum = (
      dv_utils.get_shape_and_channels_from_json(example_info_json)
  )

  # Check if the checkpoint_path has the same shape.
  if checkpoint_path is not None and example_shape is not None:
    reader = tf.compat.v1.train.NewCheckpointReader(checkpoint_path)
    shape_map_for_layers = reader.get_variable_to_shape_map()
    first_layer = 'InceptionV3/Conv2d_1a_3x3/weights'
    # For a shape map of [3, 3, 6, 32] for the Conv2d_1a_3x3 layer, the 6
    # is the number of channels.
    num_channels_in_checkpoint_model = shape_map_for_layers[first_layer][2]
    if num_channels_in_checkpoint_model != example_shape[2]:
      raise ValueError(
          'The number of channels in examples and checkpoint '
          'should match, but the checkpoint has {} channels while '
          'the examples have {}.'.format(
              num_channels_in_checkpoint_model, example_shape[2]
          )
      )
    input_info_file = os.path.join(
        os.path.dirname(checkpoint_path), 'model.ckpt.example_info.json'
    )
    ckpt_shape, ckpt_channels_enum = dv_utils.get_shape_and_channels_from_json(
        input_info_file
    )

    if ckpt_shape is not None and ckpt_channels_enum is not None:
      if example_shape != ckpt_shape:
        raise ValueError(
            f'Shape mismatch in {example_info_json} and {input_info_file}.'
        )
      if example_channels_enum != ckpt_channels_enum:
        raise ValueError(
            f'Channels mismatch in {example_info_json} and {input_info_file}.'
        )
    else:
      # We can consider more strictly enforcing this.
      logging.warning(
          'Starting from v1.4.0, we recommend having a '
          'model.ckpt.example_info.json file with your model.'
      )

  # Check accelerator status.
  if execution_hardware not in _ALLOW_EXECUTION_HARDWARE:
    raise ValueError(
        'Unexpected execution_hardware={} value. Allowed values are {}'.format(
            execution_hardware, ','.join(_ALLOW_EXECUTION_HARDWARE)
        )
    )
  init_op = tf.group(
      tf.compat.v1.global_variables_initializer(),
      tf.compat.v1.local_variables_initializer(),
  )

  config = tf.compat.v1.ConfigProto()
  if FLAGS.config_string is not None:
    text_format.Parse(FLAGS.config_string, config)
  if execution_hardware == 'cpu':
    # Don't overwrite entire dictionary.
    config.device_count['GPU'] = 0
    config.device_count['TPU'] = 0

  # Perform sanity check.
  with tf.compat.v1.Session(config=config) as sess:
    sess.run(init_op)
    if execution_hardware == 'accelerator':
      if not any(dev.device_type != 'CPU' for dev in sess.list_devices()):
        raise ExecutionHardwareError(
            'execution_hardware is set to accelerator, but no accelerator '
            'was found'
        )
    # TODO. Sort out auto-detection of TPU. Just calling
    # sess.list_devices here doesn't return the correct answer. That can only
    # work later, after the device (on the other VM) has been initialized,
    # which is generally not yet.

  # Prepare input stream and estimator.
  tf_dataset = prepare_inputs(source_path=examples_filename, use_tpu=use_tpu)
  estimator = model.make_estimator(
      batch_size=batch_size,
      master=master,
      use_tpu=use_tpu,
      session_config=config,
      include_debug_info=FLAGS.include_debug_info,
  )

  # Instantiate the prediction "stream", and select the EMA values from
  # the model.
  if checkpoint_path is None:
    # Unit tests use this branch.
    predict_hooks = []
  else:
    predict_hooks = [h(checkpoint_path) for h in model.session_predict_hooks()]
  predictions = iter(
      estimator.predict(
          input_fn=tf_dataset,
          checkpoint_path=checkpoint_path,
          hooks=predict_hooks,
      )
  )

  # The following code is introduced to by in sync with call_variants_keras
  # where we use multiple writers to write outpts.
  # If output file is already sharded then don't dynamically shard.
  if sharded_file_utils.is_sharded_filename(output_file):
    logging.info('Output is already sharded, so dynamic sharding is disabled.')
  else:
    # For call_variants, we always use one writer process.
    total_writer_process = 1
    # Convert output filename to sharded output filename.
    filename_pattern = output_file.replace(
        '.tfrecord.gz', '@' + str(total_writer_process) + '.tfrecord.gz'
    )
    output_file = sharded_file_utils.maybe_generate_sharded_filenames(
        filename_pattern
    )[0]
  # Consume predictions one at a time and write them to output_file.
  logging.info('Writing calls to %s', output_file)
  writer = tfrecord.Writer(output_file)
  with writer:
    start_time = time.time()
    n_examples, n_batches = 0, 0
    while max_batches is None or n_batches <= max_batches:
      try:
        prediction = next(predictions)
      except (StopIteration, tf.errors.OutOfRangeError):
        break
      write_variant_call(writer, prediction, use_tpu)
      n_examples += 1
      n_batches = n_examples // batch_size + 1
      duration = time.time() - start_time

      logging.log_every_n(
          logging.INFO,
          'Processed %s examples in %s batches [%.3f sec per 100]',
          _LOG_EVERY_N,
          n_examples,
          n_batches,
          (100 * duration) / n_examples,
      )
    # One last log to capture the extra examples.
    logging.info(
        'Processed %s examples in %s batches [%.3f sec per 100]',
        n_examples,
        n_batches,
        (100 * duration) / n_examples,
    )

    logging.info(
        'Done calling variants from a total of %d examples.', n_examples
    )


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: call_variants does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)),
          errors.CommandLineError,
      )
    del argv  # Unused.
    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()

    if FLAGS.use_tpu:
      master = dv_utils.resolve_master(
          FLAGS.master, FLAGS.tpu_name, FLAGS.tpu_zone, FLAGS.gcp_project
      )
    else:
      master = ''

    model = modeling.get_model(FLAGS.model_name)
    call_variants(
        examples_filename=FLAGS.examples,
        checkpoint_path=FLAGS.checkpoint,
        model=model,
        execution_hardware=FLAGS.execution_hardware,
        output_file=FLAGS.outfile,
        max_batches=FLAGS.max_batches,
        batch_size=FLAGS.batch_size,
        master=master,
        use_tpu=FLAGS.use_tpu,
    )


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'outfile',
      'checkpoint',
  ])
  tf.compat.v1.app.run()
