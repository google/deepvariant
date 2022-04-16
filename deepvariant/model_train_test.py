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
"""Tests for deepvariant .model_train."""

import json
from unittest import mock
import uuid



from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
import six
import tensorflow as tf
from tensorflow import estimator as tf_estimator

from deepvariant import data_providers_test
from deepvariant import dv_constants
from deepvariant import model_train
from deepvariant import modeling
from deepvariant import testdata
from deepvariant.testing import tf_test_utils


FLAGS = flags.FLAGS
MOCK_SENTINEL_RETURN_VALUE = 'mocked_return_value'

# Note that this test suite is invoked twice, with --use_tpu set both ways.


def setUpModule():
  testdata.init()


class ModelTrainTest(parameterized.TestCase, tf.test.TestCase):

  def _run_tiny_training(self, model_name, dataset, warm_start_from=''):
    """Runs one training step. This function always starts a new train_dir."""
    with mock.patch(
        'deepvariant.data_providers.'
        'get_input_fn_from_dataset') as mock_get_input_fn_from_dataset:
      with mock.patch('deepvariant.model_train.'
                      'copy_over_example_info_json') as _:
        mock_get_input_fn_from_dataset.return_value = dataset
        FLAGS.train_dir = tf_test_utils.test_tmpdir(uuid.uuid4().hex)
        FLAGS.batch_size = 2
        FLAGS.model_name = model_name
        FLAGS.save_interval_secs = -1
        FLAGS.save_interval_steps = 1
        FLAGS.number_of_steps = 1
        FLAGS.dataset_config_pbtxt = '/path/to/mock.pbtxt'
        FLAGS.start_from_checkpoint = warm_start_from
        FLAGS.master = ''
        model_train.parse_and_run()
        # We have a checkpoint after training.
        mock_get_input_fn_from_dataset.assert_called_once_with(
            dataset_config_filename=FLAGS.dataset_config_pbtxt,
            mode=tf_estimator.ModeKeys.TRAIN,
            use_tpu=mock.ANY,
            max_examples=None,
        )
        self.assertIsNotNone(tf.train.latest_checkpoint(FLAGS.train_dir))

  @mock.patch('deepvariant'
              '.modeling.tf.compat.v1.losses.softmax_cross_entropy')
  @mock.patch('deepvariant'
              '.modeling.tf.compat.v1.losses.get_total_loss')
  def test_loss(self, mock_total_loss, mock_cross):
    labels = [[0, 1, 0], [1, 0, 0]]
    logits = 'Logits'
    smoothing = 0.01
    actual = model_train.loss(logits, labels, smoothing)
    mock_total_loss.assert_called_once_with()
    self.assertEqual(actual, mock_total_loss.return_value)
    mock_cross.assert_called_once_with(
        logits, labels, label_smoothing=smoothing, weights=1.0)

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters((model.name, compressed_inputs)
                            for model in modeling.production_models()
                            if model.is_trainable
                            for compressed_inputs in [True, False])
  # pylint: enable=g-complex-comprehension
  @flagsaver.flagsaver
  def test_end2end(self, model_name, compressed_inputs):
    """End-to-end test of model_train script."""
    self._run_tiny_training(
        model_name=model_name,
        dataset=data_providers_test.make_golden_dataset(
            compressed_inputs=compressed_inputs, use_tpu=FLAGS.use_tpu))

  @flagsaver.flagsaver
  def test_end2end_inception_v3_warm_up_from(self):
    """End-to-end test of model_train script."""
    checkpoint_dir = tf_test_utils.test_tmpdir('inception_v3_warm_up_from')
    tf_test_utils.write_fake_checkpoint('inception_v3', self.test_session(),
                                        checkpoint_dir)
    self._run_tiny_training(
        model_name='inception_v3',
        dataset=data_providers_test.make_golden_dataset(use_tpu=FLAGS.use_tpu),
        warm_start_from=checkpoint_dir + '/model')

  @flagsaver.flagsaver
  def test_end2end_inception_v3_warm_up_allow_different_num_channels(self):
    """End-to-end test of model_train script."""
    FLAGS.allow_warmstart_from_different_num_channels = True
    checkpoint_dir = tf_test_utils.test_tmpdir(
        'inception_v3_warm_up_allow_different_num_channels')
    tf_test_utils.write_fake_checkpoint(
        'inception_v3',
        self.test_session(),
        checkpoint_dir,
        num_channels=dv_constants.PILEUP_NUM_CHANNELS + 1)
    self._run_tiny_training(
        model_name='inception_v3',
        dataset=data_providers_test.make_golden_dataset(use_tpu=FLAGS.use_tpu),
        warm_start_from=checkpoint_dir + '/model')

  @flagsaver.flagsaver
  def test_end2end_inception_v3_warm_up_by_default_fail_diff_num_channels(self):
    """End-to-end test of model_train script."""
    checkpoint_dir = tf_test_utils.test_tmpdir(
        'test_end2end_inception_v3_warm_up_by_default_fail_diff_num_channels')
    tf_test_utils.write_fake_checkpoint(
        'inception_v3',
        self.test_session(),
        checkpoint_dir,
        num_channels=dv_constants.PILEUP_NUM_CHANNELS + 1)
    with self.assertRaisesRegex(
        ValueError,
        r'Shape of variable InceptionV3/Conv2d_1a_3x3/weights:0 \(\(.*\)\) '
        r'doesn\'t match with shape of tensor '
        r'InceptionV3/Conv2d_1a_3x3/weights \(\[.*\]\) from checkpoint reader.'
    ):
      self._run_tiny_training(
          model_name='inception_v3',
          dataset=data_providers_test.make_golden_dataset(
              use_tpu=FLAGS.use_tpu),
          warm_start_from=checkpoint_dir + '/model')

  @flagsaver.flagsaver
  def test_end2end_inception_v3_failed_warm_up_from(self):
    """End-to-end test of model_train script with a non-existent path."""
    # Internal TF raises tf.errors.OpError, public TF raises ValueError.
    with self.assertRaises((tf.errors.OpError, ValueError)):
      self._run_tiny_training(
          model_name='inception_v3',
          dataset=data_providers_test.make_golden_dataset(
              use_tpu=FLAGS.use_tpu),
          warm_start_from='this/path/does/not/exist')

  @flagsaver.flagsaver
  def test_end2end_inception_v3_embedding_invalid_embedding_size(self):
    """End-to-end test of model_train script with an invalid embedding size."""
    with six.assertRaisesRegex(
        self, ValueError, 'Expected seq_type_embedding_size '
        'to be a positive number but saw -100 '
        'instead.'):
      FLAGS.seq_type_embedding_size = -100
      self._run_tiny_training(
          model_name='inception_v3_embedding',
          dataset=data_providers_test.make_golden_dataset(
              use_tpu=FLAGS.use_tpu))

  @parameterized.parameters((False), (True))
  @flagsaver.flagsaver
  @mock.patch('deepvariant.model_train.'
              'tf.compat.v1.train.replica_device_setter')
  @mock.patch('deepvariant.model_train.run')
  def test_main_internal(self, use_tpu, mock_run, mock_device_setter):
    FLAGS.master = 'some_master'
    FLAGS.use_tpu = use_tpu
    FLAGS.ps_tasks = 10
    FLAGS.task = 5

    model_train.parse_and_run()

    mock_device_setter.assert_called_once_with(10)
    mock_run.assert_called_once_with(
        'some_master' if use_tpu else '',
        False,
        device_fn=mock.ANY,
        use_tpu=mock.ANY)

  @mock.patch('deepvariant.model_train.os.environ')
  @mock.patch('deepvariant.model_train.'
              'tf.compat.v1.train.replica_device_setter')
  @mock.patch('deepvariant.model_train.run')
  def test_main_tfconfig_local(self, mock_run, mock_device_setter,
                               mock_environ):
    mock_environ.get.return_value = '{}'
    model_train.parse_and_run()

    mock_device_setter.assert_called_once_with(0)
    mock_run.assert_called_once_with(
        '', True, device_fn=mock.ANY, use_tpu=mock.ANY)

  @parameterized.named_parameters(
      ('master', 'master', 0, True, '/job:master/task:0'),
      ('worker', 'worker', 10, False, '/job:worker/task:10'),
  )
  @mock.patch(
      'deepvariant.model_train.tf.distribute.Server')
  @mock.patch('deepvariant.model_train.os.environ')
  @mock.patch('deepvariant.model_train.'
              'tf.compat.v1.train.replica_device_setter')
  @mock.patch('deepvariant.model_train.run')
  def test_main_tfconfig_dist(self, job_name, task_index, expected_is_chief,
                              expected_worker, mock_run, mock_device_setter,
                              mock_environ, mock_server):
    tf_config = {
        'cluster': {
            'ps': ['ps1:800', 'ps2:800']
        },
        'task': {
            'type': job_name,
            'index': task_index,
        },
    }

    class FakeServer(object):
      target = 'some-target'

    mock_environ.get.return_value = json.dumps(tf_config)
    mock_server.return_value = FakeServer()

    model_train.parse_and_run()

    mock_device_setter.assert_called_once_with(
        2, worker_device=expected_worker, cluster=mock.ANY)
    mock_run.assert_called_once_with(
        'some-target', expected_is_chief, device_fn=mock.ANY, use_tpu=mock.ANY)

  @parameterized.parameters(
      ('master', 'some-master'),
      ('task', 10),
      ('ps_tasks', 5),
  )
  @flagsaver.flagsaver
  @mock.patch('deepvariant.model_train.os.environ')
  def test_main_invalid_args(self, flag_name, flag_value, mock_environ):
    # Ensure an exception is raised if flags and TF_CONFIG are set.
    tf_config = {
        'cluster': {
            'ps': ['ps1:800', 'ps2:800']
        },
        'task': {
            'type': 'master',
            'index': 0,
        },
    }

    mock_environ.get.return_value = json.dumps(tf_config)
    setattr(FLAGS, flag_name, flag_value)
    self.assertRaises(ValueError, model_train.parse_and_run)


if __name__ == '__main__':
  absltest.main()
