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
"""Tests for genomics.deepvariant.model_eval."""

import os
from unittest import mock


from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
import six
import tensorflow as tf
from tensorflow import estimator as tf_estimator

from deepvariant import data_providers_test
from deepvariant import model_eval
from deepvariant import testdata
from deepvariant.testing import tf_test_utils

FLAGS = flags.FLAGS

# Note that this test suite is invoked twice, with --use_tpu set both ways.


def setUpModule():
  testdata.init()


class ModelEvalTest(
    six.with_metaclass(parameterized.TestGeneratorMetaclass, tf.test.TestCase)):

  def setUp(self):
    self.checkpoint_dir = tf.compat.v1.test.get_temp_dir()
    # Use this to generate a random name.  The framework
    # will create the directory under self.checkpoint_dir.
    self.eval_name = os.path.basename(tf.compat.v1.test.get_temp_dir())

  @parameterized.parameters(['inception_v3'])
  @flagsaver.flagsaver
  @mock.patch('deepvariant.data_providers.'
              'get_input_fn_from_dataset')
  def test_end2end(self, model_name, mock_get_input_fn_from_dataset):
    """End-to-end test of model_eval."""
    tf_test_utils.write_fake_checkpoint('inception_v3', self.test_session(),
                                        self.checkpoint_dir,
                                        FLAGS.moving_average_decay)

    # Start up eval, loading that checkpoint.
    FLAGS.batch_size = 2
    FLAGS.checkpoint_dir = self.checkpoint_dir
    FLAGS.eval_name = self.eval_name
    FLAGS.max_ckpt_to_evaluate = 0
    FLAGS.max_examples = 2
    FLAGS.best_checkpoint_metric = 'F1/All'
    FLAGS.model_name = model_name
    FLAGS.dataset_config_pbtxt = '/path/to/mock.pbtxt'
    FLAGS.master = ''
    # Always try to read in compressed inputs to stress that case. Uncompressed
    # inputs are certain to work. This test is expensive to run, so we want to
    # minimize the number of times we need to run this.
    mock_get_input_fn_from_dataset.return_value = (
        data_providers_test.make_golden_dataset(
            compressed_inputs=True, use_tpu=FLAGS.use_tpu))
    model_eval.main(0)
    mock_get_input_fn_from_dataset.assert_called_once_with(
        dataset_config_filename=FLAGS.dataset_config_pbtxt,
        mode=tf_estimator.ModeKeys.EVAL,
        use_tpu=FLAGS.use_tpu)
    self.assertTrue(
        tf_test_utils.check_file_exists(
            'best_checkpoint.txt', eval_name=self.eval_name))
    self.assertTrue(
        tf_test_utils.check_file_exists(
            'best_checkpoint.metrics', eval_name=self.eval_name))
    self.assertEqual(
        tf.train.load_checkpoint(self.checkpoint_dir).get_tensor('global_step'),
        0)

  # Using a constant model, check that running an eval returns the expected
  # metrics.
  @flagsaver.flagsaver
  @mock.patch(
      'deepvariant.model_eval.checkpoints_iterator')
  @mock.patch('deepvariant.data_providers.'
              'get_input_fn_from_dataset')
  def test_fixed_eval_sees_the_same_evals(self, mock_get_input_fn_from_dataset,
                                          mock_checkpoints_iterator):
    dataset = data_providers_test.make_golden_dataset(use_tpu=FLAGS.use_tpu)
    n_checkpoints = 3
    checkpoints = [
        tf_test_utils.write_fake_checkpoint(
            'constant',
            self.test_session(),
            self.checkpoint_dir,
            FLAGS.moving_average_decay,
            global_step=i,
            name='model' + str(i)) for i in range(n_checkpoints)
    ]

    # Setup our mocks.
    mock_checkpoints_iterator.return_value = checkpoints
    mock_get_input_fn_from_dataset.return_value = dataset

    # Start up eval, loading that checkpoint.
    FLAGS.batch_size = 2
    FLAGS.checkpoint_dir = self.checkpoint_dir
    FLAGS.eval_name = self.eval_name
    FLAGS.max_ckpt_to_evaluate = n_checkpoints - 1
    FLAGS.model_name = 'constant'
    FLAGS.dataset_config_pbtxt = '/path/to/mock.pbtxt'
    FLAGS.master = ''
    model_eval.main(0)
    self.assertEqual(
        tf.train.load_checkpoint(self.checkpoint_dir).get_tensor('global_step'),
        n_checkpoints - 1)
    self.assertEqual(mock_get_input_fn_from_dataset.call_args_list, [
        mock.call(
            use_tpu=FLAGS.use_tpu,
            dataset_config_filename=FLAGS.dataset_config_pbtxt,
            mode=tf_estimator.ModeKeys.EVAL)
    ])

    metrics = [
        model_eval.read_metrics(checkpoint, eval_name=FLAGS.eval_name)
        for checkpoint in checkpoints
    ]

    # Check that our metrics are what we expect them to be.
    # See internal for details on how to compute these counts:
    # Counts of labels in our golden dataset:
    #  1 0
    # 12 1
    # 35 2
    expected_values_for_all_exact = {
        # We have 12 correct calls [there are 12 variants with a label of 1] and
        # 1 label 0 + 35 with a label of 2, so we have an accuracy of 12 / 48,
        # which is 0.25.
        'Accuracy/All': 0.25,
        # We don't have any FNs because we call everything het.
        'FNs/All': 0,
        # Two of our labels are 0, which we call het, giving us 2 FP.
        'FPs/All': 1.0,
        # We call everything as het, so the recall has to be 1.
        'Recall/All': 1.0,
        # TODO: include this metric when we upgrade to TF 1.5.
        # # We don't call anything but hets, so TNs has to be 0.
        # 'TNs/All': 0,
        # We find 47 positives, so this has to be 47.
        'TPs/All': 47,
    }
    for key, expected_value in expected_values_for_all_exact.items():
      print(str(key) + '=' + str(metrics[0][key]))

    for key, expected_value in expected_values_for_all_exact.items():
      self.assertEqual(metrics[0][key], expected_value)

    expected_values_for_all_close = {
        # We called 47 / 48 correctly ~ 0.979167
        'Precision/All': 0.979167,
        # We called (2 * 47 / 48) / (1 + 47 / 48) correctly ~ 0.989474
        'F1/All': 0.989474,
    }
    for key, expected_value in expected_values_for_all_close.items():
      self.assertAlmostEqual(metrics[0][key], expected_value, places=6)

    for m1, m2 in zip(metrics, metrics[1:]):
      # Remove global_step from comparison first.
      m1.pop('global_step', None)
      m2.pop('global_step', None)
      self.assertEqual(m1, m2)


if __name__ == '__main__':
  absltest.main()
