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
"""Calling variants with a trained DeepVariant TF2/Keras model."""

import multiprocessing
import os
import time
from typing import Any, Sequence



from absl import flags
from absl import logging
import numpy as np
import tensorflow as tf

from deepvariant import dv_utils
from deepvariant import keras_modeling as modeling
from deepvariant import logging_level
from deepvariant.protos import deepvariant_pb2
from absl import app
import multiprocessing
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import vis


_ALLOW_EXECUTION_HARDWARE = [
    'auto',  # Default, no validation.
    'cpu',  # Don't use accelerators, even if available.
    'accelerator',  # Must be hardware acceleration or an error will be raised.
]

# The number of digits past the decimal point that genotype likelihoods are
# rounded to, for numerical stability.
_GL_PRECISION = 10

_MAX_WRITER_THREADS = 16

_LOG_EVERY_N_BATCHES = 50

_DEFAULT_INPUT_READ_THREADS = 32
_DEFAULT_PREFETCH_BUFFER_BYTES = 16 * 1000 * 1000

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
    1024,
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
flags.DEFINE_boolean(
    'include_debug_info',
    False,
    'If true, include extra debug info in the output, including the original '
    'image_encoded.',
)
flags.DEFINE_list(
    'activation_layers',
    [],
    'A list of activation layer names which we add to the debug info output.'
    ' Needs include_debug_info flag to be True.',
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
flags.DEFINE_integer(
    'writer_threads',
    0,
    (
        'Number of threads to use for writing. Set 0 to autodetect.'
        ' In autodetect mode, 1 thread is used in CPU inference and'
        ' all cpus are used when GPU is available. If set to a specific value'
        ' other than 0 then autodetect is disabled. Maximum 16 processes'
        ' can be used for writing. Default: 0'
    ),
)
_LIMIT = flags.DEFINE_integer(
    'limit', 0, 'If set to > 0, limit processing to <= limit examples.'
)


class ExecutionHardwareError(Exception):
  pass


class CustomCallback(tf.keras.callbacks.Callback):
  """Custom callbacks for `predict`."""

  def on_predict_batch_begin(self, batch, logs=None):
    logging.log_every_n(
        logging.INFO,
        'Begin `predict` on batch %d.',
        _LOG_EVERY_N_BATCHES,
        batch,
    )

  def on_predict_batch_end(self, batch, logs=None):
    logging.log_every_n(
        logging.INFO, 'End `predict` on batch %d.', _LOG_EVERY_N_BATCHES, batch
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


def write_variant_call(
    writer, prediction, include_debug_info, debugging_true_label_mode
):
  """Write the variant call based on prediction.

  Args:
    writer: A object with a write() function that will be called for each
      encoded_variant and genotype likelihoods.
    prediction: A [3] tensor of floats. These are the predicted genotype
      likelihoods (p00, p0x, pxx) for some alt allele x, in the same order as
      encoded_variants.
    include_debug_info: If true, include debug information.
    debugging_true_label_mode: If true, include true label from the example.

  Returns:
    The return status from writer.
  """
  encoded_variant = prediction['variant']
  encoded_alt_allele_indices = prediction['alt_allele_indices']
  rounded_gls = round_gls(prediction['probabilities'], precision=_GL_PRECISION)

  # Write it out.
  true_labels = None
  if debugging_true_label_mode:
    if 'label' not in prediction:
      logging.warning('No label is available for debugging_true_label_mode.')
    else:
      true_labels = prediction['label']
  image_encoded = prediction['image_encoded'] if include_debug_info else None
  layer_outputs_encoded = (
      prediction['layer_outputs_encoded'] if include_debug_info else None
  )
  cvo = _create_cvo_proto(
      encoded_variant,
      rounded_gls,
      encoded_alt_allele_indices,
      true_labels,
      logits=prediction.get('logits'),
      prelogits=prediction.get('prelogits'),
      image_encoded=image_encoded,
      include_debug_info=include_debug_info,
      debugging_true_label_mode=debugging_true_label_mode,
      layer_outputs_encoded=layer_outputs_encoded,
      pileup_curation=prediction.get('pileup_curation'),
  )
  return writer.write(cvo)


def add_pileup_curation_to_debug_info(
    debug_info: deepvariant_pb2.CallVariantsOutput.DebugInfo,
    pileup_curation: vis.PileupCuration,
) -> deepvariant_pb2.CallVariantsOutput.DebugInfo:
  """Adds pileup curation to debug_info."""
  if pileup_curation is None:
    return debug_info
  debug_info.pileup_curation.diff_category = pileup_curation.diff_category.value
  debug_info.pileup_curation.base_quality = pileup_curation.base_quality.value
  debug_info.pileup_curation.mapping_quality = (
      pileup_curation.mapping_quality.value
  )
  debug_info.pileup_curation.strand_bias = pileup_curation.strand_bias.value
  debug_info.pileup_curation.read_support = pileup_curation.read_support.value
  return debug_info


def _create_cvo_proto(
    encoded_variant,
    gls,
    encoded_alt_allele_indices,
    true_labels=None,
    logits=None,
    prelogits=None,
    image_encoded=None,
    include_debug_info=False,
    debugging_true_label_mode=False,
    layer_outputs_encoded=None,
    pileup_curation=None,
):
  """Returns a CallVariantsOutput proto from the relevant input information."""
  variant = variants_pb2.Variant.FromString(encoded_variant)
  alt_allele_indices = (
      deepvariant_pb2.CallVariantsOutput.AltAlleleIndices.FromString(
          encoded_alt_allele_indices
      )
  )
  debug_info = None
  if include_debug_info or debugging_true_label_mode:
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
        image_encoded=image_encoded,
        layer_output_encoded=layer_outputs_encoded,
    )
    debug_info = add_pileup_curation_to_debug_info(debug_info, pileup_curation)

  call_variants_output = deepvariant_pb2.CallVariantsOutput(
      variant=variant,
      alt_allele_indices=alt_allele_indices,
      genotype_probabilities=gls,
      debug_info=debug_info,
  )
  return call_variants_output


# TODO: Consider creating one data loading function to re-use simliar
#                code with training in train_inceptionv3.py.
def get_dataset(
    path,
    example_shape,
    batch_size,
    include_debug_info,
    debugging_true_label_mode,
):
  """Parse TFRecords, do image preprocessing, and return the image dataset for inference and the variant/alt-allele dataset for writing the variant calls."""

  proto_features = {
      'image/encoded': tf.io.FixedLenFeature((), tf.string),
      'variant/encoded': tf.io.FixedLenFeature((), tf.string),
      'alt_allele_indices/encoded': tf.io.FixedLenFeature((), tf.string),
  }
  if debugging_true_label_mode:
    proto_features['label'] = tf.io.FixedLenFeature((1), tf.int64)

  def _parse_example(example):
    """Parses a serialized tf.Example."""
    parsed_features = tf.io.parse_single_example(
        serialized=example, features=proto_features
    )
    image_encoded = (
        parsed_features['image/encoded'] if include_debug_info else b''
    )
    image = tf.io.decode_raw(parsed_features['image/encoded'], tf.uint8)
    image = tf.reshape(image, example_shape)
    image = tf.cast(image, tf.float32)
    image = dv_utils.preprocess_images(image)
    variant = parsed_features['variant/encoded']
    alt_allele_indices = parsed_features['alt_allele_indices/encoded']
    optional_label = None
    if debugging_true_label_mode:
      optional_label = parsed_features['label']
    return image_encoded, image, variant, alt_allele_indices, optional_label

  ds = tf.data.TFRecordDataset.list_files(
      sharded_file_utils.normalize_to_sharded_file_pattern(path), shuffle=False
  )

  def load_dataset(filename):
    dataset = tf.data.TFRecordDataset(
        filename,
        buffer_size=_DEFAULT_PREFETCH_BUFFER_BYTES,
        compression_type='GZIP',
    )
    return dataset

  ds = ds.interleave(
      load_dataset,
      cycle_length=_DEFAULT_INPUT_READ_THREADS,
      num_parallel_calls=tf.data.AUTOTUNE,
  )
  if _LIMIT.value > 0:
    ds = ds.take(_LIMIT.value)

  enc_image_variant_alt_allele_ds = ds.map(
      map_func=_parse_example, num_parallel_calls=tf.data.AUTOTUNE
  )

  enc_image_variant_alt_allele_ds = enc_image_variant_alt_allele_ds.batch(
      batch_size=batch_size
  ).prefetch(tf.data.AUTOTUNE)
  return enc_image_variant_alt_allele_ds


def post_processing(
    output_file: str,
    output_queue: Any,
    include_debug_info: bool,
    debugging_true_label_mode: bool,
) -> None:
  """Post processing of called variants.

  Args:
    output_file: Path to output file where outputs will be written.
    output_queue: Multiprocessing queue to fetch predictions from.
    include_debug_info: If true, include debug information.
    debugging_true_label_mode: If true, include true label from the example.
  """
  writer = tfrecord.Writer(output_file, compression_type='GZIP')
  n_examples = 0
  n_batches = 0
  while True:
    try:
      # If no records found in 3 minutes, then break.
      item = output_queue.get(timeout=180)
      if item is None:
        break
      (
          predictions,
          image_encodes,
          variants,
          alt_allele_indices_list,
          optional_label_list,
          layer_outputs_encoded,
          optional_pileup_curation_list,
      ) = item
      if optional_label_list is not None:
        optional_label_list = optional_label_list.numpy()
      for i, probabilities in enumerate(predictions):
        pred = {
            'probabilities': probabilities,
            'variant': variants[i].numpy(),
            'image_encoded': image_encodes[i].numpy(),
            'alt_allele_indices': alt_allele_indices_list[i].numpy(),
            'layer_outputs_encoded': layer_outputs_encoded,
        }
        if include_debug_info:
          if optional_pileup_curation_list is not None:
            pred['pileup_curation'] = optional_pileup_curation_list[i]
          if debugging_true_label_mode and optional_label_list is not None:
            pred['label'] = optional_label_list[i][0]
        write_variant_call(
            writer, pred, include_debug_info, debugging_true_label_mode
        )
        n_examples += 1
      n_batches += 1
    except TimeoutError as error:
      logging.warning('Writer timeout occurred %s.', str(error))
      break
  writer.close()


def call_variants(
    examples_filename: str,
    checkpoint_path: str,
    output_file: str,
    writer_threads: int,
    kmp_blocktime: str,
    batch_size: int,
    include_debug_info: bool,
    debugging_true_label_mode: bool,
    activation_layers: Sequence[str],
):
  """Main driver of call_variants."""
  first_example = dv_utils.get_one_example_from_examples_path(examples_filename)
  if first_example is None:
    logging.warning(
        'Unable to read any records from %s. Output shards will contain zero'
        ' records.',
        examples_filename,
    )
    # Write empty shards
    total_writer_process = 1
    output_file = output_file.replace(
        '.tfrecord.gz', '@' + str(total_writer_process) + '.tfrecord.gz'
    )
    paths = sharded_file_utils.maybe_generate_sharded_filenames(output_file)
    for path in paths:
      tfrecord.write_tfrecords([], path, compression_type='GZIP')
    return

  # See if GPU is available
  is_gpu_available = True if tf.config.list_physical_devices('GPU') else False

  # This means we are in autodetect mode.
  if writer_threads == 0:
    # If GPU is available then use all CPUs for writing.
    total_writer_process = (
        multiprocessing.cpu_count() if is_gpu_available else 1
    )
  else:
    total_writer_process = writer_threads

  # If output file is already sharded then don't dynamically shard.
  if sharded_file_utils.is_sharded_filename(output_file):
    logging.info('Output is already sharded, so dynamic sharding is disabled.')
    total_writer_process = 1
    paths = [output_file]
  else:
    # Use maximum _MAX_WRITER_THREADS threads, this is to lower the overhead of
    # spinning up and shutting down many processes.
    total_writer_process = min(total_writer_process, _MAX_WRITER_THREADS)
    # Convert output filename to sharded output filename.
    output_file = output_file.replace(
        '.tfrecord.gz', '@' + str(total_writer_process) + '.tfrecord.gz'
    )
    paths = sharded_file_utils.maybe_generate_sharded_filenames(output_file)

  output_queue = multiprocessing.Queue()
  all_processes = []
  for process_id in range(0, total_writer_process):
    post_processing_process = multiprocessing.get_context().Process(
        target=post_processing,
        args=(
            paths[process_id],
            output_queue,
            include_debug_info,
            debugging_true_label_mode,
        ),
    )
    all_processes.append(post_processing_process)
    post_processing_process.start()

  logging.info('Total %d writing processes started.', len(all_processes))

  if kmp_blocktime:
    os.environ['KMP_BLOCKTIME'] = kmp_blocktime
    logging.vlog(
        3, 'Set KMP_BLOCKTIME to {}'.format(os.environ['KMP_BLOCKTIME'])
    )

  example_info_json = dv_utils.get_example_info_json_filename(
      examples_filename, 0
  )
  example_shape = dv_utils.get_shape_and_channels_from_json(example_info_json)[
      0
  ]

  if example_shape is None:
    logging.info(
        (
            'Unable to read shape information from %s. Directly read from '
            'examples instead.'
        ),
        example_info_json,
    )
    example_shape = dv_utils.example_image_shape(first_example)

  logging.info('Shape of input examples: %s', str(example_shape))
  use_saved_model = tf.io.gfile.exists(checkpoint_path) and tf.io.gfile.exists(
      f'{checkpoint_path}/saved_model.pb'
  )
  logging.info('Use saved model: %s', str(use_saved_model))

  if checkpoint_path is not None:
    if use_saved_model:
      model = tf.saved_model.load(checkpoint_path)
      model_example_info_json = f'{checkpoint_path}/example_info.json'
      model_example_shape = dv_utils.get_shape_and_channels_from_json(
          model_example_info_json
      )
      # This usually happens in multi-sample cases where the sample_name
      # is added to the prefix, so we remove it.
      if not tf.io.gfile.exists(example_info_json):
        # Find the sample name by finding file prefix and then the last bit
        # would be the sample name.
        filename_prefix = os.path.basename(example_info_json).split('.')[0]
        sample_name = filename_prefix.split('_')[-1]
        expected_filename = example_info_json.replace('_' + sample_name, '')
        if not tf.io.gfile.exists(expected_filename):
          raise ValueError(
              f'File {example_info_json} or {expected_filename} not found.'
              'Please check make_examples output.'
          )
        example_info_json = expected_filename

      input_example_shape = dv_utils.get_shape_and_channels_from_json(
          example_info_json
      )
      # These checks make sure we are using the right model with right input.
      if model_example_shape[0] != input_example_shape[0]:
        # The following has been changed from ValueError to Warning because of
        # internal
        logging.warning(
            'Input shape %s and model shape %s does not match.',
            str(input_example_shape[0]),
            str(model_example_shape[0]),
        )
      if model_example_shape[1] != input_example_shape[1]:
        # The following has been changed from ValueError to Warning because of
        # internal
        logging.warning(
            'Input channels %s and model channels %s do not match.',
            str(input_example_shape[1]),
            str(model_example_shape[1]),
        )
    else:
      model = modeling.inceptionv3(
          example_shape, init_backbone_with_imagenet=False
      )
      model.load_weights(checkpoint_path).expect_partial()

    enc_image_variant_alt_allele_ds = get_dataset(
        examples_filename,
        example_shape,
        batch_size,
        include_debug_info,
        debugging_true_label_mode,
    )

    activation_models = {}
    if include_debug_info and activation_layers and not use_saved_model:
      activation_models = {
          layer_name: modeling.get_activations_model(model, layer_name)
          for layer_name in activation_layers
      }

    batch_no = 0
    n_examples = 0
    n_batches = 0
    start_time = time.time()
    for batch in enc_image_variant_alt_allele_ds:
      # Here are the elements in `batch`. To see where the reading logic came
      # from, please see the `get_dataset` function, specifically the
      # `_parse_example` function within it.
      # TODO: consider using these as variable names to make this
      # function more readable.
      # batch[0]: image_encodes
      # batch[1]: images_in_batch
      # batch[2]: variants
      # batch[3]: alt_allele_indices_list
      # batch[4]: optional_label_list
      if use_saved_model:
        predictions = model.signatures['serving_default'](batch[1])
        predictions = predictions['classification'].numpy()
      else:
        if not is_gpu_available:
          # This is faster on CPU but slower on GPU.
          predictions = model.predict_on_batch(batch[1])
        else:
          # This is faster on GPU but slower on CPU.
          predictions = model(batch[1], training=False).numpy()

      layer_outputs_encoded = None
      if include_debug_info and activation_layers:
        if not use_saved_model:
          if not is_gpu_available:
            # This is faster on CPU but slower on GPU.
            layer_outputs_encoded = {
                layer: activation_model.predict_on_batch(batch[1]).tobytes()
                for layer, activation_model in activation_models.items()
            }
          else:
            # This is faster on GPU but slower on CPU.
            layer_outputs_encoded = {
                layer: (
                    activation_model(batch[1], training=False).numpy().tobytes()
                )
                for layer, activation_model in activation_models.items()
            }
        else:
          logging.warning(
              'Activation layers: %s are not saved in CVO. Change to ckpt'
              'model file to get activation layers outputs.',
              ','.join(activation_layers),
          )

      pileup_curation_in_batch = None
      if include_debug_info:
        pileup_curation_in_batch = [
            vis.curate_pileup(
                vis.split_3d_array_into_channels(
                    dv_utils.unpreprocess_images(image.numpy())
                )
            )
            for image in batch[1]
        ]
      batch_no += 1
      duration = time.time() - start_time
      n_examples += len(predictions)
      n_batches += 1
      logging.log_every_n(
          logging.INFO,
          'Predicted %s examples in %s batches [%.3f sec per 100].',
          _LOG_EVERY_N_BATCHES,
          n_examples,
          n_batches,
          (100 * duration) / n_examples,
      )
      output_queue.put((
          predictions,
          batch[0],  # image_encodes
          batch[2],  # variants
          batch[3],  # alt_allele_indices_list
          batch[4],  # optional_label_list
          layer_outputs_encoded,
          pileup_curation_in_batch,
      ))

    # Put none values in the queue so the running processes can terminate.
    for _ in range(0, total_writer_process):
      output_queue.put(None)
    for post_processing_process in all_processes:
      post_processing_process.join()


def main(argv=()):
  env = os.environ.copy()
  logging.info('call_variants: env = %s', env)
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
    # Make sure output filename is consistent and can be used for multi-writing.
    if not sharded_file_utils.is_sharded_filename(
        FLAGS.outfile
    ) and not FLAGS.outfile.endswith('.tfrecord.gz'):
      raise ValueError('Output filename must end with .tfrecord.gz')

    if FLAGS.activation_layers:
      if not FLAGS.include_debug_info:
        logging.warning(
            '"--include_debug_info" need to be set True to have the input'
            ' activation_layers: %s included in CVO.',
            ','.join(FLAGS.activation_layers),
        )

    call_variants(
        examples_filename=FLAGS.examples,
        checkpoint_path=FLAGS.checkpoint,
        output_file=FLAGS.outfile,
        writer_threads=FLAGS.writer_threads,
        kmp_blocktime=FLAGS.kmp_blocktime,
        batch_size=FLAGS.batch_size,
        include_debug_info=FLAGS.include_debug_info,
        debugging_true_label_mode=FLAGS.debugging_true_label_mode,
        activation_layers=FLAGS.activation_layers,
    )
    logging.info('Complete: call_variants.')


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'outfile',
      'checkpoint',
  ])
  logging.use_python_logging()
  app.run(main)
