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
"""Tests for deepvariant .call_variants_slim."""

import collections
import errno
import sys
from unittest import mock



from absl import flags
from absl import logging
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
import numpy as np
import tensorflow as tf
from tensorflow import estimator as tf_estimator

from deepvariant import call_variants_slim
from deepvariant import dv_utils
from deepvariant import modeling
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io import tfrecord
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import variant_utils

FLAGS = flags.FLAGS

# NB. This entire collection of tests will be invoked with '--use_tpu=' 'true'
# and 'false' by the BUILD file, and a tpu device will be allocated when
# necessary.


def setUpModule():
  testdata.init()


# For tests that don't actually want to read a real checkpoint,
# return a fake one.  The estimator understands None to mean
# that all the variables should be left uninitialized.
_LEAVE_MODEL_UNINITIALIZED = None


# Return the stream of batched images from a dataset.
def _get_infer_batches(tf_dataset, model, batch_size):
  """Provides batches of pileup images from this dataset.

  This instantiates an iterator on the dataset, and returns the
  image, variant, alt_allele_indices, features in batches. It calls
  model.preprocess_images on the images (but note that we will be moving
  that step into model_fn for the Estimator api).

  Args:
    tf_dataset: DeepVariantInput.
    model: DeepVariantModel.
    batch_size: int.  The batch size.

  Returns:
    (image, variant, alt_allele_indices)

  Raises:
    ValueError: if the dataset has the wrong mode.
  """
  if tf_dataset.mode != tf_estimator.ModeKeys.PREDICT:
    raise ValueError(
        'tf_dataset.mode is {} but must be PREDICT.'.format(tf_dataset.mode)
    )

  params = dict(batch_size=batch_size)
  features = tf.compat.v1.data.make_one_shot_iterator(
      tf_dataset(params)
  ).get_next()

  images = features['image']
  if tf_dataset.tensor_shape:
    # tensor_shape will be None if the input was an empty file.
    images = model.preprocess_images(images)
  variant = features['variant']
  alt_allele_indices = features['alt_allele_indices']

  return images, variant, alt_allele_indices


class CallVariantsEndToEndTests(
    tf.compat.v1.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
):

  def setUp(self):
    super().setUp()
    self.checkpoint_dir = tf.compat.v1.test.get_temp_dir()

  def assertCallVariantsEmitsNRecordsForInceptionV3(
      self, filename, num_examples
  ):
    outfile = test_utils.test_tmpfile('inception_v3.call_variants.tfrecord')
    model = modeling.get_model('inception_v3')
    checkpoint_path = _LEAVE_MODEL_UNINITIALIZED

    call_variants_slim.call_variants(
        examples_filename=filename,
        checkpoint_path=checkpoint_path,
        model=model,
        output_file=outfile,
        batch_size=4,
        max_batches=None,
    )
    call_variants_outputs = list(
        tfrecord.read_tfrecords(outfile, deepvariant_pb2.CallVariantsOutput)
    )
    # Check that we have the right number of output protos.
    self.assertEqual(len(call_variants_outputs), num_examples)

  def assertCallVariantsEmitsNRecordsForConstantModel(
      self, filename, num_examples
  ):
    checkpoint_path = _LEAVE_MODEL_UNINITIALIZED
    outfile = test_utils.test_tmpfile('call_variants.tfrecord')
    model = modeling.get_model('constant')
    call_variants_slim.call_variants(
        examples_filename=filename,
        checkpoint_path=checkpoint_path,
        model=model,
        output_file=outfile,
        batch_size=4,
        max_batches=None,
        primary='',
        use_tpu=FLAGS.use_tpu,
    )
    call_variants_outputs = list(
        tfrecord.read_tfrecords(outfile, deepvariant_pb2.CallVariantsOutput)
    )
    # Check that we have the right number of output protos.
    self.assertEqual(len(call_variants_outputs), num_examples)

  def test_call_end2end_with_empty_shards(self):
    # Get only up to 10 examples.
    examples = list(
        tfrecord.read_tfrecords(
            testdata.GOLDEN_CALLING_EXAMPLES, max_records=10
        )
    )
    # Write to 15 shards, which means there will be multiple empty shards.
    source_path = test_utils.test_tmpfile('sharded@{}'.format(15))
    tfrecord.write_tfrecords(examples, source_path)
    self.assertCallVariantsEmitsNRecordsForConstantModel(
        source_path, len(examples)
    )

  def test_call_end2end_empty_first_shard(self):
    # Get only up to 10 examples.
    examples = list(
        tfrecord.read_tfrecords(
            testdata.GOLDEN_CALLING_EXAMPLES, max_records=10
        )
    )
    empty_first_file = test_utils.test_tmpfile('empty_1st_shard-00000-of-00002')
    tfrecord.write_tfrecords([], empty_first_file)
    second_file = test_utils.test_tmpfile('empty_1st_shard-00001-of-00002')
    tfrecord.write_tfrecords(examples, second_file)
    self.assertCallVariantsEmitsNRecordsForConstantModel(
        test_utils.test_tmpfile('empty_1st_shard@2'), len(examples)
    )

  def test_call_end2end_zero_record_file_for_inception_v3(self):
    zero_record_file = test_utils.test_tmpfile('zero_record_file')
    tfrecord.write_tfrecords([], zero_record_file)
    self.assertCallVariantsEmitsNRecordsForInceptionV3(
        test_utils.test_tmpfile('zero_record_file'), 0
    )

  def _call_end2end_helper(self, examples_path, model, shard_inputs):
    examples = list(tfrecord.read_tfrecords(examples_path))

    if shard_inputs:
      # Create a sharded version of our golden examples.
      source_path = test_utils.test_tmpfile('sharded@{}'.format(3))
      tfrecord.write_tfrecords(examples, source_path)
    else:
      source_path = examples_path

    # If we point the test at a headless server, it will often be 2x2,
    # which has 8 replicas.  Otherwise a smaller batch size is fine.
    if FLAGS.use_tpu:
      batch_size = 8
    else:
      batch_size = 4

    if model.name == 'constant':
      # For the constant model we can run everything.
      max_batches = None
    else:
      # For all other models we only run a single batch for inference.
      max_batches = 1

    outfile = test_utils.test_tmpfile('call_variants.tfrecord')
    call_variants_slim.call_variants(
        examples_filename=source_path,
        checkpoint_path=_LEAVE_MODEL_UNINITIALIZED,
        model=model,
        output_file=outfile,
        batch_size=batch_size,
        max_batches=max_batches,
        primary='',
        use_tpu=FLAGS.use_tpu,
    )

    call_variants_outputs = list(
        tfrecord.read_tfrecords(outfile, deepvariant_pb2.CallVariantsOutput)
    )

    return call_variants_outputs, examples, batch_size, max_batches

  @parameterized.parameters(model for model in modeling.production_models())
  @flagsaver.flagsaver
  def test_call_end2end_with_labels(self, model):
    FLAGS.debugging_true_label_mode = True
    (call_variants_outputs, examples, batch_size, max_batches) = (
        self._call_end2end_helper(
            testdata.GOLDEN_TRAINING_EXAMPLES, model, False
        )
    )
    # Check that we have the right number of output protos.
    self.assertEqual(
        len(call_variants_outputs),
        batch_size * max_batches if max_batches else len(examples),
    )

    # Checks that at least some of the `true_label`s are filled.
    self.assertTrue(
        any(cvo.debug_info.true_label > 0 for cvo in call_variants_outputs)
    )

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(
      (model, shard_inputs, include_debug_info)
      for shard_inputs in [False, True]
      for model in modeling.production_models()
      for include_debug_info in [False, True]
  )  # pylint: enable=g-complex-comprehension
  @flagsaver.flagsaver
  def test_call_end2end(self, model, shard_inputs, include_debug_info):
    FLAGS.include_debug_info = include_debug_info
    (call_variants_outputs, examples, batch_size, max_batches) = (
        self._call_end2end_helper(
            testdata.GOLDEN_CALLING_EXAMPLES, model, shard_inputs
        )
    )
    # Check that we have the right number of output protos.
    self.assertEqual(
        len(call_variants_outputs),
        batch_size * max_batches if max_batches else len(examples),
    )

    # Check that our CallVariantsOutput (CVO) have the following critical
    # properties:
    # - we have one CVO for each example we processed.
    # - the variant in the CVO is exactly what was in the example.
    # - the alt_allele_indices of the CVO match those of its corresponding
    #   example.
    # - there are 3 genotype probabilities and these are between 0.0 and 1.0.
    # We can only do this test when processing all of the variants (max_batches
    # is None), since we processed all of the examples with that model.
    if max_batches is None:
      self.assertCountEqual(
          [cvo.variant for cvo in call_variants_outputs],
          [dv_utils.example_variant(ex) for ex in examples],
      )

    # Check the CVO debug_info: not filled if include_debug_info is False;
    # else, filled by logic based on CVO.
    if not include_debug_info:
      for cvo in call_variants_outputs:
        self.assertEqual(
            cvo.debug_info, deepvariant_pb2.CallVariantsOutput.DebugInfo()
        )
    else:
      for cvo in call_variants_outputs:
        self.assertEqual(
            cvo.debug_info.has_insertion,
            variant_utils.has_insertion(cvo.variant),
        )
        self.assertEqual(
            cvo.debug_info.has_deletion, variant_utils.has_deletion(cvo.variant)
        )
        self.assertEqual(
            cvo.debug_info.is_snp, variant_utils.is_snp(cvo.variant)
        )
        self.assertEqual(
            cvo.debug_info.predicted_label,
            np.argmax(cvo.genotype_probabilities),
        )
        self.assertEqual(len(cvo.debug_info.logits), 3)
        self.assertEqual(len(cvo.debug_info.prelogits), 2048)

    def example_matches_call_variants_output(example, call_variants_output):
      return (
          dv_utils.example_variant(example) == call_variants_output.variant
          and dv_utils.example_alt_alleles_indices(example)
          == call_variants_output.alt_allele_indices.indices
      )

    for call_variants_output in call_variants_outputs:
      # Find all matching examples.
      matches = [
          ex
          for ex in examples
          if example_matches_call_variants_output(ex, call_variants_output)
      ]
      # We should have exactly one match.
      self.assertEqual(len(matches), 1)
      example = matches[0]
      # Check that we've faithfully copied in the alt alleles (though currently
      # as implemented we find our example using this information so it cannot
      # fail). Included here in case that changes in the future.
      self.assertEqual(
          list(dv_utils.example_alt_alleles_indices(example)),
          list(call_variants_output.alt_allele_indices.indices),
      )
      # We should have exactly three genotype probabilities (assuming our
      # ploidy == 2).
      self.assertEqual(len(call_variants_output.genotype_probabilities), 3)
      # These are probabilities so they should be between 0 and 1.
      self.assertTrue(
          0 <= gp <= 1 for gp in call_variants_output.genotype_probabilities
      )

  @parameterized.parameters(model for model in modeling.production_models())
  def test_call_variants_with_no_shape(self, model):
    # Read one good record from a valid file.
    example = next(tfrecord.read_tfrecords(testdata.GOLDEN_CALLING_EXAMPLES))
    # Remove image/shape.
    del example.features.feature['image/shape']
    source_path = test_utils.test_tmpfile('make_examples_out_noshape.tfrecord')
    tfrecord.write_tfrecords([example], source_path)
    with self.assertRaisesRegex(
        ValueError,
        (
            'Invalid image/shape: we expect to find an image/shape '
            'field with length 3.'
        ),
    ):
      ds = call_variants_slim.prepare_inputs(source_path)
      _ = list(_get_infer_batches(ds, model=model, batch_size=1))

  def test_call_variants_with_empty_input(self):
    source_path = test_utils.test_tmpfile('empty.tfrecord')
    tfrecord.write_tfrecords([], source_path)
    # Make sure that prepare_inputs don't crash on empty input.
    ds = call_variants_slim.prepare_inputs(source_path)
    m = modeling.get_model('constant')

    # The API specifies that OutOfRangeError is thrown in this case.
    batches = list(_get_infer_batches(ds, model=m, batch_size=1))
    with self.test_session() as sess:
      sess.run(tf.compat.v1.local_variables_initializer())
      sess.run(tf.compat.v1.global_variables_initializer())
      try:
        _ = sess.run(batches)
      except tf.errors.OutOfRangeError:
        pass


class CallVariantsUnitTests(
    tf.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
):

  @classmethod
  def setUpClass(cls):
    super().setUpClass()
    cls.examples = list(
        tfrecord.read_tfrecords(testdata.GOLDEN_CALLING_EXAMPLES)
    )
    cls.variants = [dv_utils.example_variant(ex) for ex in cls.examples]
    cls.model = modeling.get_model('constant')

  @parameterized.parameters(
      ('not_sharded', 'not_sharded'),
      ('sharded@3', 'sharded@3'),
      ('sharded@3', 'sharded-?????-of-00003'),
      ('asterisks@2', 'asterisks-*-of-00002'),
  )
  def test_prepare_inputs(self, filename_to_write, file_string_input):
    source_path = test_utils.test_tmpfile(filename_to_write)
    tfrecord.write_tfrecords(self.examples, source_path)
    # file_string_input could be a comma-separated list. Add the prefix to all
    # of them, and join it back to a string.
    file_string_input = ','.join(
        [test_utils.test_tmpfile(f) for f in file_string_input.split(',')]
    )

    with self.test_session() as sess:
      sess.run(tf.compat.v1.local_variables_initializer())
      sess.run(tf.compat.v1.global_variables_initializer())

      ds = call_variants_slim.prepare_inputs(file_string_input)
      _, variants, _ = _get_infer_batches(ds, model=self.model, batch_size=1)

      seen_variants = []
      try:
        while True:
          seen_variants.extend(sess.run(variants))
      except tf.errors.OutOfRangeError:
        pass

      self.assertCountEqual(
          self.variants, variant_utils.decode_variants(seen_variants)
      )

  @parameterized.parameters(
      (None, [3.592555731302127e-5, 0.99992620944976807, 3.78809563699178e-5]),
      (2, [0.0, 1.0, 0.0]),
      (10, [3.59096e-5, 0.9999262094, 3.7881e-5]),
  )
  def test_round_gls(self, precision, expected):
    test_data = [3.592555731302127e-5, 0.99992620944976807, 3.78809563699178e-5]
    actual = call_variants_slim.round_gls(test_data, precision)
    self.assertEqual(actual, expected)

  @parameterized.parameters('auto', 'cpu')
  def test_call_variants_non_accelerated_execution_runs(
      self, execution_hardware
  ):
    if FLAGS.use_tpu:
      # predict batch size must be divisible by number of replicas.
      batch_size = 2
    else:
      batch_size = 1
    outfile = test_utils.test_tmpfile('call_variants_cpu_only.tfrecord')
    call_variants_slim.call_variants(
        examples_filename=testdata.GOLDEN_CALLING_EXAMPLES,
        checkpoint_path=_LEAVE_MODEL_UNINITIALIZED,
        model=self.model,
        execution_hardware=execution_hardware,
        max_batches=1,
        batch_size=batch_size,
        output_file=outfile,
        use_tpu=FLAGS.use_tpu,
    )

  @parameterized.parameters(
      dict(hardware_env='auto', devices=['cpu'], expect_exception=False),
      dict(hardware_env='auto', devices=['gpu'], expect_exception=False),
      dict(hardware_env='auto', devices=['tpu'], expect_exception=False),
      dict(hardware_env='cpu', devices=['cpu'], expect_exception=False),
      dict(hardware_env='cpu', devices=['gpu'], expect_exception=False),
      dict(hardware_env='cpu', devices=['tpu'], expect_exception=False),
      dict(hardware_env='accelerator', devices=['cpu'], expect_exception=True),
      dict(hardware_env='accelerator', devices=['gpu'], expect_exception=False),
      dict(hardware_env='accelerator', devices=['tpu'], expect_exception=False),
      dict(
          hardware_env='accelerator',
          devices=['cpu', 'gpu'],
          expect_exception=False,
      ),
      dict(
          hardware_env='accelerator',
          devices=['cpu', 'tpu'],
          expect_exception=False,
      ),
      dict(
          hardware_env='accelerator',
          devices=['cpu', 'gpu', 'tpu'],
          expect_exception=False,
      ),
  )
  def test_call_variants_execution_hardware(
      self, hardware_env, devices, expect_exception
  ):
    # We cannot access the full _DeviceAttribute as it's not exported. So use a
    # namedtuple with the same field names instead.
    device = collections.namedtuple('_DeviceAttribute', ['name', 'device_type'])

    # Mocking the list_devices call means the framework attempts to use a bogus
    # TPU device, which fails, so don't do that.  Handle the TPU case elsewhere.
    if 'tpu' in devices or FLAGS.use_tpu:
      return

    with mock.patch.object(
        call_variants_slim.tf.compat.v1.Session, 'list_devices'
    ) as mock_ld:
      mock_ld.return_value = [
          device(name=dt + '/' + str(i), device_type=dt.upper())
          for i, dt in enumerate(devices)
      ]

      # Only run the tpu cases when we have an actual tpu device, supplied
      # by the flags from the BUILD rule.
      def _run():
        call_variants_slim.call_variants(
            use_tpu=FLAGS.use_tpu,
            examples_filename=testdata.GOLDEN_CALLING_EXAMPLES,
            checkpoint_path=_LEAVE_MODEL_UNINITIALIZED,
            model=self.model,
            execution_hardware=hardware_env,
            max_batches=1,
            batch_size=1,
            output_file=test_utils.test_tmpfile('zzz.tfrecord'),
        )

      if expect_exception:
        with self.assertRaises(call_variants_slim.ExecutionHardwareError):
          _run()
      else:
        _run()

  def test_catches_bad_argv(self):
    with (
        mock.patch.object(logging, 'error') as mock_logging,
        mock.patch.object(sys, 'exit') as mock_exit,
    ):
      call_variants_slim.main(['call_variants_slim.py', 'extra_arg'])
    mock_logging.assert_called_once_with(
        'Command line parsing failure: call_variants does not accept '
        'positional arguments but some are present on the command line: '
        "\"['call_variants_slim.py', 'extra_arg']\"."
    )
    mock_exit.assert_called_once_with(errno.ENOENT)


if __name__ == '__main__':
  absltest.main()
