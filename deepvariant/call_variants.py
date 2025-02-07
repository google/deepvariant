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

import copy
import itertools
import json
import multiprocessing
import os
import time
from typing import Any, Sequence

from absl import flags
from absl import logging
import numpy as np
import tensorflow as tf


# pylint: disable=C0301
from bazel_tools.tools.python.runfiles import runfiles
r = runfiles.Create()
filename = os.path.join(r.Rlocation('com_google_deepvariant/deepvariant/examples_from_stream.so'))
dv_stream_dataset = tf.load_op_library(filename)
from deepvariant import dv_utils
from deepvariant import keras_modeling as modeling
from deepvariant import logging_level
from deepvariant.protos import deepvariant_pb2
from tensorflow.python.platform import gfile
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

_EXAMPLES = flags.DEFINE_string(
    'examples',
    None,
    (
        'Required. tf.Example protos containing DeepVariant candidate variants'
        ' in TFRecord format, as emitted by make_examples. Can be a'
        ' comma-separated list of files, and the file names can contain'
        ' wildcard characters.'
    ),
)
_OUTFILE = flags.DEFINE_string(
    'outfile',
    None,
    (
        'Required. Destination path where we will write output candidate'
        ' variants with additional likelihood information in TFRecord format of'
        ' CallVariantsOutput protos.'
    ),
)
_CHECKPOINT = flags.DEFINE_string(
    'checkpoint',
    None,
    (
        'Required. Path to the TensorFlow model checkpoint to use to evaluate '
        'candidate variant calls.'
    ),
)
_BATCH_SIZE = flags.DEFINE_integer(
    'batch_size',
    1024,
    (
        'Number of candidate variant tensors to batch together during'
        ' inference. Larger batches use more memory but are more computational'
        ' efficient.'
    ),
)
_MAX_BATCHES = flags.DEFINE_integer(
    'max_batches', None, 'Max. batches to evaluate. Defaults to all.'
)
_NUM_READERS = flags.DEFINE_integer(
    'num_readers', 8, 'Number of parallel readers to create for examples.'
)
_INCLUDE_DEBUG_INFO = flags.DEFINE_boolean(
    'include_debug_info',
    False,
    'If true, include extra debug info in the output, including the original '
    'image_encoded.',
)
_ACTIVATION_LAYERS = flags.DEFINE_list(
    'activation_layers',
    [],
    'A list of activation layer names which we add to the debug info output.'
    ' Needs include_debug_info flag to be True.',
)
_DEBUGGING_TRUE_LABEL_MODE = flags.DEFINE_boolean(
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
_EXECUTION_HARDWARE = flags.DEFINE_string(
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
_CONFIG_STRING = flags.DEFINE_string(
    'config_string',
    None,
    (
        'String representation of a tf.ConfigProto message, with'
        ' comma-separated key: value pairs, such as "allow_soft_placement:'
        ' True". The value can itself be another message, such as "gpu_options:'
        ' {per_process_gpu_memory_fraction: 0.5}".'
    ),
)
_KMP_BLOCKTIME = flags.DEFINE_string(
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
_WRITER_THREADS = flags.DEFINE_integer(
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
_NUM_INPUT_SHARDS = flags.DEFINE_integer(
    'num_input_shards',
    None,
    'Number of input shards. This flag is required when examples are'
    ' streammed.',
)
_SHM_PREFIX = flags.DEFINE_string(
    'shm_prefix',
    'deepvariant',
    'Shared memory files prefix.',
)
_STREAM_EXAMPLES = flags.DEFINE_boolean(
    'stream_examples',
    False,
    'If true, examples are streammed from make_example processes.',
)
_ALLOW_EMPTY_EXAMPLES = flags.DEFINE_boolean(
    'allow_empty_examples',
    True,
    'If true, call_variants will not crash if the examples are empty. This may'
    ' reasonably happen when the small model is used, or when processing a'
    ' small region.',
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


# This class implements the stream dataset.
class FromStreamDataset(tf.data.Dataset):
  """Dataset implementation for streaming examples from shared memory."""

  def __init__(self, shm_prefix, num_shards):
    resource = dv_stream_dataset.stream_examples_init(shm_prefix, num_shards)
    self._resource = resource
    self.num_empty = tf.Variable(0)
    self.num_examples = tf.Variable(0, dtype=tf.int64)
    self.num_shards = tf.Variable(num_shards)

    # Generates consecutive numbers.
    dataset = tf.data.Dataset.counter()
    # Stream is stopped once all make_examples shards are processed.
    dataset = dataset.take_while(
        lambda v: tf.less(self.num_empty, self.num_shards)
    )
    # Retrieving examples from the stream.
    dataset = dataset.map(
        lambda i: dv_stream_dataset.stream_examples_next(self._resource, i),
        num_parallel_calls=tf.data.AUTOTUNE,
    )
    # Filter operation keeps the counter of completed shards.
    dataset = dataset.filter(self.filter_fn)
    # Unbatching is needed because examples are streammed in a variable size
    # batches.
    dataset = dataset.unbatch()

    self._dataset = dataset
    super().__init__(self._dataset._variant_tensor)

  def filter_fn(self, x):
    if tf.equal(tf.shape(x.image)[0], 0):
      self.num_empty.assign_add(1)
      return False
    return True

  def _inputs(self):
    return []

  @property
  def element_spec(self):
    return self._dataset.element_spec


# TODO: Consider creating one data loading function to re-use simliar
#                code with training in train_inceptionv3.py.
def get_dataset(
    path,
    example_shape,
    channel_indices,
    batch_size,
    include_debug_info,
    debugging_true_label_mode,
    use_data_from_stream=False,
    shm_prefix=None,
    num_shards=None,
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
    image = dv_utils.preprocess_images(image, channel_indices)
    variant = parsed_features['variant/encoded']
    alt_allele_indices = parsed_features['alt_allele_indices/encoded']
    optional_label = None
    if debugging_true_label_mode:
      optional_label = parsed_features['label']
    return image_encoded, image, variant, alt_allele_indices, optional_label

  def _parse_example_from_stream(blob):
    """Parses a data from shared memory buffer."""
    image = tf.io.decode_raw(blob.image, tf.uint8)
    image = tf.reshape(image, example_shape)
    image = dv_utils.preprocess_images(image, channel_indices)
    variant = blob.variant
    alt_allele_indices = blob.alt_allele_idx
    return b'', image, variant, alt_allele_indices, None

  # Dataset from stream.
  if use_data_from_stream:
    ds = FromStreamDataset(shm_prefix, num_shards)
    enc_image_variant_alt_allele_ds = ds.map(
        map_func=_parse_example_from_stream, num_parallel_calls=tf.data.AUTOTUNE
    )
    enc_image_variant_alt_allele_ds = enc_image_variant_alt_allele_ds.batch(
        batch_size=batch_size
    ).prefetch(tf.data.AUTOTUNE)
    return enc_image_variant_alt_allele_ds
  # Dataset from TFRecord files.
  else:
    ds = tf.data.TFRecordDataset.list_files(
        sharded_file_utils.normalize_to_sharded_file_pattern(path),
        shuffle=False,
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
      # If no records found in 5 minutes, then break.
      item = output_queue.get(timeout=300)
      if item is None:
        break
      (
          predictions,
          image_encodes,
          variants,
          alt_allele_indices_list,
          optional_label_list,
          layer_outputs,
          optional_pileup_curation_list,
      ) = item
      for i, probabilities in enumerate(predictions):
        pred = {
            'probabilities': probabilities,
            'variant': variants[i],
            'image_encoded': image_encodes[i],
            'alt_allele_indices': alt_allele_indices_list[i],
            'layer_outputs_encoded': {
                layer_name: output[i].tobytes()
                for layer_name, output in layer_outputs.items()
            },
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
  logging.info(
      'Wrote %d examples over %d batches. Closing writer.',
      n_examples,
      n_batches,
  )
  writer.close()


def write_empty_output_file(examples_filename: str, output_file: str) -> None:
  """Writes an empty output file."""
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


def load_model_and_check_shape(
    checkpoint_path: str,
    examples_filename: str,
    first_example: Any,
    use_saved_model: bool,
    use_examples_from_stream: bool,
) -> tuple[Any, Any]:
  """Ensure that the example shape matches the model shape."""

  example_shape = []
  example_info_json: str = ''
  if not use_examples_from_stream:
    example_info_json = dv_utils.get_example_info_json_filename(
        examples_filename, 0
    )
    example_shape, _, _ = dv_utils.get_shape_and_channels_from_json(
        example_info_json
    )

    if example_shape is None:
      logging.info(
          (
              'Unable to read shape information from %s. Directly read from '
              'examples instead.'
          ),
          example_info_json,
      )
      example_shape = dv_utils.example_image_shape(first_example)
  # end of (if not use_examples_from_stream)

  if use_saved_model:
    model = tf.saved_model.load(checkpoint_path)
    model_example_info_json = f'{checkpoint_path}/example_info.json'
    model_example_shape, _, channel_indices = (
        dv_utils.get_shape_and_channels_from_json(model_example_info_json)
    )
    if channel_indices:
      raise ValueError('Channel ablation is not supported for saved model.')
    # This usually happens in multi-sample cases where the sample_name
    # is added to the prefix, so we remove it.
    if use_examples_from_stream:
      example_shape = model_example_shape
    else:
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

    # These checks make sure we are using the right model with right input.
    if not use_examples_from_stream:
      if model_example_shape[0] != example_shape[0]:
        # The following has been changed from ValueError to Warning because of
        # internal
        logging.warning(
            'Input shape %s and model shape %s does not match.',
            str(example_shape[0]),
            str(model_example_shape[0]),
        )
      if model_example_shape[1] != example_shape[1]:
        # The following has been changed from ValueError to Warning because of
        # internal
        logging.warning(
            'Input channels %s and model channels %s do not match.',
            str(example_shape[1]),
            str(model_example_shape[1]),
        )
  else:
    model_example_info_json = None
    # If example_shape could not be inferred from the examples, then we try to
    # infer it from the model directory. This is the case when examples from
    # stream is used.
    model_dir = os.path.dirname(checkpoint_path)
    for dirname, subdir, fnames in gfile.Walk(model_dir):
      if subdir:
        continue
      for fname in fnames:
        if fname.endswith('example_info.json'):
          model_example_info_json = f'{dirname}/{fname}'
          break
    if model_example_info_json:
      example_shape, _, channel_indices = (
          dv_utils.get_shape_and_channels_from_json(
              f'{model_example_info_json}'
          )
      )
    else:
      raise ValueError(
          'Could not infer example shape from examples or model directory.'
          'example_info.json was not found in the model directory.'
      )

    # If we are using channel ablation, the example shape of incoming
    # data remains the same, but the model example shape needs to be
    # adjusted to reflect the channels that are being ablated.
    # Update shape to reflect ablation.
    if channel_indices:
      model_example_shape = copy.copy(example_shape)
      model_example_shape[2] = len(channel_indices)
      logging.info(
          'Channel Ablation updates model input shape %s.',
          str(model_example_shape),
      )
    else:
      model_example_shape = example_shape

    logging.info('model_example_shape: %s', model_example_shape)
    model = modeling.inceptionv3(
        model_example_shape, init_backbone_with_imagenet=False
    )
    model.load_weights(checkpoint_path).expect_partial()
  return example_shape, model


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
    use_dataset_from_stream: bool,
    shm_prefix: str,
    num_shards: int,
    allow_empty_examples: bool,
):
  """Main driver of call_variants."""
  first_example = None
  if not use_dataset_from_stream:
    first_example = dv_utils.get_one_example_from_examples_path(
        examples_filename
    )
    if first_example is None:
      write_empty_output_file(examples_filename, output_file)
      if allow_empty_examples:
        logging.info('Allowing empty examples. Exiting.')
        return
      else:
        raise ValueError(
            'Examples file is empty. If you want to allow empty examples,'
            ' please set --allow_empty_examples to True.'
        )

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

  writer_queues = []
  for _ in range(total_writer_process):
    writer_queues.append(multiprocessing.Queue())
  writer_queues_iterator = itertools.cycle(writer_queues)

  all_processes = []
  for process_id in range(0, total_writer_process):
    post_processing_process = multiprocessing.get_context().Process(
        target=post_processing,
        args=(
            paths[process_id],
            writer_queues[process_id],
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

  use_saved_model = tf.io.gfile.exists(checkpoint_path) and tf.io.gfile.exists(
      f'{checkpoint_path}/saved_model.pb'
  )
  logging.info('Use saved model: %s', str(use_saved_model))

  if checkpoint_path is None:
    raise ValueError('Checkpoint filename must be specified.')

  channel_indices = []
  if not use_saved_model:
    if not os.path.isdir(checkpoint_path):
      example_info_json_path = os.path.join(
          os.path.dirname(checkpoint_path), 'example_info.json'
      )
    else:
      example_info_json_path = os.path.join(
          checkpoint_path, 'example_info.json'
      )
    example_info = json.loads(tf.io.gfile.GFile(example_info_json_path).read())

    for idx, channel_enum in enumerate(example_info['channels']):
      if channel_enum not in example_info.get('ablation_channels', []):
        channel_indices.append(idx)

  example_shape, model = load_model_and_check_shape(
      checkpoint_path,
      examples_filename,
      first_example,
      use_saved_model,
      _STREAM_EXAMPLES.value,
  )

  if not example_shape:
    raise ValueError(
        'Could not infer example shape from examples or model directory.'
    )

  logging.info('example_shape: %s', example_shape)
  enc_image_variant_alt_allele_ds = get_dataset(
      examples_filename,
      example_shape,
      channel_indices,
      batch_size,
      include_debug_info,
      debugging_true_label_mode,
      use_dataset_from_stream,
      shm_prefix,
      num_shards,
  )

  activation_model = model
  if include_debug_info and activation_layers:
    if not use_saved_model:
      activation_model = modeling.get_activations_model(
          model, activation_layers
      )

  batch_no = 0
  n_examples = 0
  n_batches = 0
  start_time = time.time()
  for (
      image_encodes,
      images_in_batch,
      variants,
      alt_allele_indices_list,
      optional_label_list,
  ) in enc_image_variant_alt_allele_ds:
    # These elements per iteration are read from the `get_dataset` function,
    # specifically the `_parse_example` function within it.
    if use_saved_model:
      predictions = model.signatures['serving_default'](images_in_batch)
      predictions = predictions['classification'].numpy()
    else:
      if not is_gpu_available:
        # This is faster on CPU but slower on GPU.
        predictions = model.predict_on_batch(images_in_batch)
      else:
        # This is faster on GPU but slower on CPU.
        predictions = model(images_in_batch, training=False).numpy()

    layer_outputs = {}
    if include_debug_info and activation_layers:
      if not use_saved_model:
        if not is_gpu_available:
          # This is faster on CPU but slower on GPU.
          layer_outputs = activation_model.predict_on_batch(images_in_batch)
        else:
          # This is faster on GPU but slower on CPU.
          layer_outputs = activation_model(images_in_batch, training=False)
          layer_outputs = {
              layer_name: output.numpy()
              for layer_name, output in layer_outputs.items()
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
          for image in images_in_batch
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
    if optional_label_list is not None:
      optional_label_list = optional_label_list.numpy()
    next(writer_queues_iterator).put((
        predictions,
        image_encodes.numpy(),
        variants.numpy(),
        alt_allele_indices_list.numpy(),
        optional_label_list,
        layer_outputs,
        pileup_curation_in_batch,
    ))

  # Put none values in each queue so the running processes can terminate.
  for queue in writer_queues:
    queue.put(None)
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

    if not _CHECKPOINT.value:
      raise ValueError('Checkpoint filename must be specified.')

    logging_level.set_from_flag()
    # Make sure output filename is consistent and can be used for multi-writing.
    outfile = _OUTFILE.value
    if not outfile:
      raise ValueError('Output filename must be specified.')

    if not sharded_file_utils.is_sharded_filename(
        outfile
    ) and not outfile.endswith('.tfrecord.gz'):
      raise ValueError('Output filename must end with .tfrecord.gz')

    if _ACTIVATION_LAYERS.value:
      if not _INCLUDE_DEBUG_INFO.value:
        logging.warning(
            '"--include_debug_info" need to be set True to have the input'
            ' activation_layers: %s included in CVO.',
            ','.join(_ACTIVATION_LAYERS.value),
        )

    call_variants(
        examples_filename=_EXAMPLES.value,
        checkpoint_path=_CHECKPOINT.value,
        output_file=_OUTFILE.value,
        writer_threads=_WRITER_THREADS.value,
        kmp_blocktime=_KMP_BLOCKTIME.value,
        batch_size=_BATCH_SIZE.value,
        include_debug_info=_INCLUDE_DEBUG_INFO.value,
        debugging_true_label_mode=_DEBUGGING_TRUE_LABEL_MODE.value,
        activation_layers=_ACTIVATION_LAYERS.value,
        use_dataset_from_stream=_STREAM_EXAMPLES.value,
        shm_prefix=_SHM_PREFIX.value,
        num_shards=_NUM_INPUT_SHARDS.value,
        allow_empty_examples=_ALLOW_EMPTY_EXAMPLES.value,
    )
    logging.info('Complete: call_variants.')


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'outfile',
      'checkpoint',
  ])
  logging.use_python_logging()
  app.run(main)
