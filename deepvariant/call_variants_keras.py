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
"""Experimental code for calling variants with a trained DeepVariant TF2/Keras model.

Added in v1.4.0 but not officially supported.

TODO: Write a unit test suite like call_variants_test.py.
"""

import multiprocessing
import os
import time
from typing import Any


from absl import flags
from absl import logging
import numpy as np
import tensorflow as tf

from deepvariant import dv_utils
from deepvariant import keras_modeling as modeling
from deepvariant import logging_level
from deepvariant.protos import deepvariant_pb2
from absl import app

from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import variant_utils


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

_LOG_EVERY_N_BATCHES = 100

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


def write_variant_call(writer, prediction):
  """Write the variant call based on prediction.

  Args:
    writer: A object with a write() function that will be called for each
      encoded_variant and genotype likelihoods.
    prediction: A [3] tensor of floats. These are the predicted genotype
      likelihoods (p00, p0x, pxx) for some alt allele x, in the same order as
      encoded_variants.

  Returns:
    The return status from writer.
  """
  encoded_variant = prediction['variant']
  encoded_alt_allele_indices = prediction['alt_allele_indices']
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


# TODO: Consider creating one data loading function to re-use simliar
#                code with training in train_inceptionv3.py.
def get_dataset(path, example_shape):
  """Parse TFRecords, do image preprocessing, and return the image dataset for inference and the variant/alt-allele dataset for writing the variant calls."""

  proto_features = {
      'image/encoded': tf.io.FixedLenFeature((), tf.string),
      'variant/encoded': tf.io.FixedLenFeature((), tf.string),
      'alt_allele_indices/encoded': tf.io.FixedLenFeature((), tf.string),
  }

  def _parse_example(example):
    """Parses a serialized tf.Example."""
    parsed_features = tf.io.parse_single_example(
        serialized=example, features=proto_features
    )
    image = tf.io.decode_raw(parsed_features['image/encoded'], tf.uint8)
    image = tf.reshape(image, example_shape)
    image = tf.cast(image, tf.float32)
    image = dv_utils.preprocess_images(image)
    variant = parsed_features['variant/encoded']
    alt_allele_indices = parsed_features['alt_allele_indices/encoded']
    return image, variant, alt_allele_indices

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

  image_variant_alt_allele_ds = ds.map(
      map_func=_parse_example, num_parallel_calls=tf.data.AUTOTUNE
  )

  image_variant_alt_allele_ds = image_variant_alt_allele_ds.batch(
      batch_size=FLAGS.batch_size
  ).prefetch(tf.data.AUTOTUNE)
  return image_variant_alt_allele_ds


def post_processing(output_file: str, output_queue: Any) -> None:
  """Post processing of called variants.

  Args:
    output_file: Path to output file where outputs will be written.
    output_queue: Multiprocessing queue to fetch predictions from.
  """
  writer = tfrecord.Writer(output_file)
  n_examples = 0
  n_batches = 0
  start_time = time.time()
  while True:
    item = output_queue.get()
    if item is None:
      break
    predictions, variants, alt_allele_indices_list = item
    for probabilities, variant, alt_allele_indices in zip(
        predictions, variants, alt_allele_indices_list
    ):
      pred = {
          'probabilities': probabilities,
          'variant': variant.numpy(),
          'alt_allele_indices': alt_allele_indices.numpy(),
      }
      write_variant_call(writer, pred)
      n_examples += 1
      duration = time.time() - start_time
      logging.log_every_n(
          logging.INFO,
          'Processed %s examples in %s batches [%.3f sec per 100]',
          _LOG_EVERY_N,
          n_examples,
          n_batches,
          (100 * duration) / n_examples,
      )
    n_batches += 1
  writer.close()
  duration = time.time() - start_time
  logging.info(
      'Processed %s examples in %s batches [%.3f sec per 100]',
      n_examples,
      n_batches,
      (100 * duration) / n_examples,
  )
  logging.info('Done calling variants from a total of %d examples.', n_examples)


def call_variants(
    examples_filename: str, checkpoint_path: str, output_file: str
):
  """Main driver of call_variants."""
  output_queue = multiprocessing.Queue()
  post_processing_process = multiprocessing.get_context().Process(
      target=post_processing,
      args=(
          output_file,
          output_queue,
      ),
  )
  post_processing_process.start()
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

  if checkpoint_path is not None:
    model = modeling.inceptionv3(
        example_shape, init_backbone_with_imagenet=False
    )
    model.load_weights(checkpoint_path).expect_partial()

    image_variant_alt_allele_ds = get_dataset(examples_filename, example_shape)

    batch_no = 0
    for batch in image_variant_alt_allele_ds:
      predictions = model.predict_on_batch(batch[0])
      batch_no += 1
      output_queue.put((predictions, batch[1], batch[2]))
    output_queue.put(None)
    post_processing_process.join()


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

    call_variants(
        examples_filename=FLAGS.examples,
        checkpoint_path=FLAGS.checkpoint,
        output_file=FLAGS.outfile,
    )


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'outfile',
      'checkpoint',
  ])
  logging.use_python_logging()
  app.run(main)
