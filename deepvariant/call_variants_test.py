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
"""Tests for deepvariant .deepvariant."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import errno
import sys



from absl.testing import parameterized
import mock
import numpy as np
import six
import tensorflow as tf

from absl import logging
from deepvariant import call_variants
from deepvariant import modeling
from deepvariant import test_utils
from deepvariant import tf_utils
from deepvariant.core import io_utils
from deepvariant.core import variantutils
from deepvariant.protos import deepvariant_pb2
from deepvariant.testing import flagsaver

FLAGS = tf.flags.FLAGS


def setUpModule():
  test_utils.init()


class CallVariantsEndToEndTests(
    six.with_metaclass(parameterized.TestGeneratorMetaclass, tf.test.TestCase)):

  def assertCallVariantsEmitsNRecordsForRandomGuess(self, filename,
                                                    num_examples):
    outfile = test_utils.test_tmpfile('call_variants.tfrecord')
    model = modeling.get_model('random_guess')
    call_variants.call_variants(
        examples_filename=filename,
        checkpoint_path=modeling.SKIP_MODEL_INITIALIZATION_IN_TEST,
        model=model,
        output_file=outfile,
        batch_size=4,
        max_batches=None)
    call_variants_outputs = list(
        io_utils.read_tfrecords(outfile, deepvariant_pb2.CallVariantsOutput))
    # Check that we have the right number of output protos.
    self.assertEqual(len(call_variants_outputs), num_examples)

  def test_call_end2end_with_empty_shards(self):
    # Get only up to 10 examples.
    examples = list(
        io_utils.read_tfrecords(
            test_utils.GOLDEN_CALLING_EXAMPLES, max_records=10))
    # Write to 15 shards, which means there will be multiple empty shards.
    source_path = test_utils.test_tmpfile('sharded@{}'.format(15))
    io_utils.write_tfrecords(examples, source_path)
    self.assertCallVariantsEmitsNRecordsForRandomGuess(source_path,
                                                       len(examples))

  def test_call_end2end_empty_first_shard(self):
    # Get only up to 10 examples.
    examples = list(
        io_utils.read_tfrecords(
            test_utils.GOLDEN_CALLING_EXAMPLES, max_records=10))
    empty_first_file = test_utils.test_tmpfile('empty_1st_shard-00000-of-00002')
    io_utils.write_tfrecords([], empty_first_file)
    second_file = test_utils.test_tmpfile('empty_1st_shard-00001-of-00002')
    io_utils.write_tfrecords(examples, second_file)
    self.assertCallVariantsEmitsNRecordsForRandomGuess(
        test_utils.test_tmpfile('empty_1st_shard@2'), len(examples))

  @parameterized.parameters((model, shard_inputs, include_debug_info)
                            for shard_inputs in [False, True]
                            for model in modeling.production_models()
                            for include_debug_info in [False, True])
  @flagsaver.FlagSaver
  def test_call_end2end(self, model, shard_inputs, include_debug_info):
    FLAGS.include_debug_info = include_debug_info
    examples = list(io_utils.read_tfrecords(test_utils.GOLDEN_CALLING_EXAMPLES))

    if shard_inputs:
      # Create a sharded version of our golden examples.
      source_path = test_utils.test_tmpfile('sharded@{}'.format(3))
      io_utils.write_tfrecords(examples, source_path)
    else:
      source_path = test_utils.GOLDEN_CALLING_EXAMPLES

    batch_size = 4
    if model.name == 'random_guess':
      # For the random guess model we can run everything.
      max_batches = None
    else:
      # For all other models we only run a single batch for inference.
      max_batches = 1

    outfile = test_utils.test_tmpfile('call_variants.tfrecord')
    call_variants.call_variants(
        examples_filename=source_path,
        checkpoint_path=modeling.SKIP_MODEL_INITIALIZATION_IN_TEST,
        model=model,
        output_file=outfile,
        batch_size=batch_size,
        max_batches=max_batches)

    call_variants_outputs = list(
        io_utils.read_tfrecords(outfile, deepvariant_pb2.CallVariantsOutput))

    # Check that we have the right number of output protos.
    self.assertEqual(
        len(call_variants_outputs), batch_size * max_batches
        if max_batches else len(examples))

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
      self.assertItemsEqual([cvo.variant for cvo in call_variants_outputs],
                            [tf_utils.example_variant(ex) for ex in examples])

    # Check the CVO debug_info: not filled if include_debug_info is False;
    # else, filled by logic based on CVO.
    if not include_debug_info:
      for cvo in call_variants_outputs:
        self.assertEqual(cvo.debug_info,
                         deepvariant_pb2.CallVariantsOutput.DebugInfo())
    else:
      for cvo in call_variants_outputs:
        self.assertEqual(cvo.debug_info.has_insertion,
                         variantutils.has_insertion(cvo.variant))
        self.assertEqual(cvo.debug_info.has_deletion,
                         variantutils.has_deletion(cvo.variant))
        self.assertEqual(cvo.debug_info.is_snp, variantutils.is_snp(
            cvo.variant))
        self.assertEqual(cvo.debug_info.predicted_label,
                         np.argmax(cvo.genotype_probabilities))

    def example_matches_call_variants_output(example, call_variants_output):
      return (tf_utils.example_variant(example) == call_variants_output.variant
              and tf_utils.example_alt_alleles_indices(
                  example) == call_variants_output.alt_allele_indices.indices)

    for call_variants_output in call_variants_outputs:
      # Find all matching examples.
      matches = [
          ex for ex in examples
          if example_matches_call_variants_output(ex, call_variants_output)
      ]
      # We should have exactly one match.
      self.assertEqual(len(matches), 1)
      example = matches[0]
      # Check that we've faithfully copied in the alt alleles (though currently
      # as implemented we find our example using this information so it cannot
      # fail). Included here in case that changes in the future.
      self.assertEqual(
          list(tf_utils.example_alt_alleles_indices(example)),
          list(call_variants_output.alt_allele_indices.indices))
      # We should have exactly three genotype probabilities (assuming our
      # ploidy == 2).
      self.assertEqual(len(call_variants_output.genotype_probabilities), 3)
      # These are probabilities so they should be between 0 and 1.
      self.assertTrue(
          0 <= gp <= 1 for gp in call_variants_output.genotype_probabilities)

  @parameterized.parameters((model, bad_format)
                            for model in modeling.production_models()
                            for bad_format in ['', 'png'])
  def test_call_variants_with_invalid_format(self, model, bad_format):
    # Read one good record from a valid file.
    example = next(io_utils.read_tfrecords(test_utils.GOLDEN_CALLING_EXAMPLES))
    # Overwrite the image/format field to be an invalid value
    # (anything but 'raw').
    example.features.feature['image/format'].bytes_list.value[0] = bad_format
    source_path = test_utils.test_tmpfile('make_examples_output.tfrecord')
    io_utils.write_tfrecords([example], source_path)
    outfile = test_utils.test_tmpfile('call_variants_invalid_format.tfrecord')

    with self.assertRaises(ValueError):
      call_variants.call_variants(
          examples_filename=source_path,
          checkpoint_path=modeling.SKIP_MODEL_INITIALIZATION_IN_TEST,
          model=model,
          output_file=outfile,
          batch_size=1,
          max_batches=1)

  @parameterized.parameters(model for model in modeling.production_models())
  def test_call_variants_with_no_shape(self, model):
    # Read one good record from a valid file.
    example = next(io_utils.read_tfrecords(test_utils.GOLDEN_CALLING_EXAMPLES))
    # Remove image/shape.
    del example.features.feature['image/shape']
    source_path = test_utils.test_tmpfile('make_examples_out_noshape.tfrecord')
    io_utils.write_tfrecords([example], source_path)
    with self.assertRaisesRegexp(
        ValueError, 'Invalid image/shape: we expect to find an image/shape '
        'field with length 3.'):
      call_variants.prepare_inputs(source_path, model, batch_size=1)

  def test_call_variants_with_empty_input(self):
    source_path = test_utils.test_tmpfile('empty.tfrecord')
    io_utils.write_tfrecords([], source_path)
    # Make sure that prepare_inputs don't crash on empty input.
    call_variants.prepare_inputs(
        source_path, modeling.get_model('random_guess'), batch_size=1)


class CallVariantsUnitTests(
    six.with_metaclass(parameterized.TestGeneratorMetaclass, tf.test.TestCase)):

  @classmethod
  def setUpClass(cls):
    cls.examples = list(
        io_utils.read_tfrecords(test_utils.GOLDEN_CALLING_EXAMPLES))
    cls.variants = [tf_utils.example_variant(ex) for ex in cls.examples]
    cls.model = modeling.get_model('random_guess')

  @parameterized.parameters(
      ('not_sharded', False),
      ('sharded@3', False),
      ('sharded@3', True),
  )
  def test_prepare_inputs(self, filename, expand_to_file_pattern):
    source_path = test_utils.test_tmpfile(filename)
    io_utils.write_tfrecords(self.examples, source_path)
    if expand_to_file_pattern:
      # Transform foo@3 to foo-?????-of-00003.
      source_path = io_utils.NormalizeToShardedFilePattern(source_path)

    with self.test_session() as sess:
      _, variants, _ = call_variants.prepare_inputs(
          source_path, self.model, batch_size=1)
      sess.run(tf.local_variables_initializer())
      sess.run(tf.global_variables_initializer())

      seen_variants = []
      try:
        while True:
          seen_variants.extend(sess.run(variants))
      except tf.errors.OutOfRangeError:
        pass

      self.assertItemsEqual(self.variants,
                            variantutils.decode_variants(seen_variants))

  @parameterized.parameters(
      (None, [3.592555731302127e-5, 0.99992620944976807, 3.78809563699178e-5]),
      (2, [0.0, 1.0, 0.0]),
      (10, [3.59096e-5, 0.9999262094, 3.7881e-5]),
  )
  def test_round_gls(self, precision, expected):
    test_data = [3.592555731302127e-5, 0.99992620944976807, 3.78809563699178e-5]
    actual = call_variants.round_gls(test_data, precision)
    self.assertEqual(actual, expected)

  @parameterized.parameters('auto', 'cpu')
  def test_call_variants_non_accelerated_execution_runs(self,
                                                        execution_hardware):
    # This doesn't mock out the list_devices call so it's worth keeping
    # despite being very similar to the parameterized test below.
    outfile = test_utils.test_tmpfile('call_variants_cpu_only.tfrecord')
    call_variants.call_variants(
        examples_filename=test_utils.GOLDEN_CALLING_EXAMPLES,
        checkpoint_path=modeling.SKIP_MODEL_INITIALIZATION_IN_TEST,
        model=self.model,
        execution_hardware=execution_hardware,
        max_batches=1,
        batch_size=1,
        output_file=outfile)

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
          expect_exception=False),
      dict(
          hardware_env='accelerator',
          devices=['cpu', 'tpu'],
          expect_exception=False),
      dict(
          hardware_env='accelerator',
          devices=['cpu', 'gpu', 'tpu'],
          expect_exception=False),
  )
  def test_call_variants_execution_hardware(self, hardware_env, devices,
                                            expect_exception):
    # We cannot access the full _DeviceAttribute as it's not exported. So use a
    # namedtuple with the same field names instead.
    device = collections.namedtuple('_DeviceAttribute', ['name', 'device_type'])

    with mock.patch.object(call_variants.tf.Session, 'list_devices') as mock_ld:
      mock_ld.return_value = [
          device(name=dt + '/' + str(i), device_type=dt.upper())
          for i, dt in enumerate(devices)
      ]

      def _run():
        call_variants.call_variants(
            examples_filename=test_utils.GOLDEN_CALLING_EXAMPLES,
            checkpoint_path=modeling.SKIP_MODEL_INITIALIZATION_IN_TEST,
            model=self.model,
            execution_hardware=hardware_env,
            max_batches=1,
            batch_size=1,
            output_file=test_utils.test_tmpfile('zzz.tfrecord'))

      if expect_exception:
        with self.assertRaises(call_variants.ExecutionHardwareError):
          _run()
      else:
        _run()

  def test_catches_bad_argv(self):
    with mock.patch.object(logging, 'error') as mock_logging,\
        mock.patch.object(sys, 'exit') as mock_exit:
      call_variants.main(['call_variants.py', 'extra_arg'])
    mock_logging.assert_called_once_with(
        'Command line parsing failure: call_variants does not accept '
        'positional arguments but some are present on the command line: '
        '"[\'call_variants.py\', \'extra_arg\']".')
    mock_exit.assert_called_once_with(errno.ENOENT)


if __name__ == '__main__':
  tf.test.main()
