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
"""Tests for deepvariant .model_train."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json



from absl.testing import absltest
from absl.testing import parameterized
import mock
import tensorflow as tf

from deepvariant import data_providers_test
from deepvariant import model_train
from deepvariant import modeling
from deepvariant import test_utils
from deepvariant.testing import flagsaver

FLAGS = tf.flags.FLAGS
MOCK_SENTINEL_RETURN_VALUE = 'mocked_return_value'


def setUpModule():
  test_utils.init()


class ModelTrainTest(parameterized.TestCase):

  @flagsaver.FlagSaver
  def test_training_works_with_compressed_inputs(self):
    """End-to-end test of model_train script."""
    self._run_tiny_training(
        model_name='mobilenet_v1',
        dataset=data_providers_test.make_golden_dataset(compressed_inputs=True))

  def _run_tiny_training(self, model_name, dataset):
    with mock.patch(
        'deepvariant.data_providers.get_dataset'
    ) as mock_get_dataset:
      mock_get_dataset.return_value = dataset
      FLAGS.train_dir = test_utils.test_tmpfile(model_name)
      FLAGS.batch_size = 2
      FLAGS.model_name = model_name
      FLAGS.save_interval_secs = 0
      FLAGS.number_of_steps = 1
      FLAGS.dataset_config_pbtxt = '/path/to/mock.pbtxt'
      FLAGS.start_from_checkpoint = ''
      model_train.parse_and_run()
      # We have a checkpoint after training.
      mock_get_dataset.assert_called_once_with(FLAGS.dataset_config_pbtxt)
      self.assertIsNotNone(tf.train.latest_checkpoint(FLAGS.train_dir))

  @parameterized.parameters(
      model.name for model in modeling.production_models() if model.is_trainable
  )
  @flagsaver.FlagSaver
  def test_end2end(self, model_name):
    """End-to-end test of model_train script."""
    self._run_tiny_training(
        model_name=model_name,
        dataset=data_providers_test.make_golden_dataset())

  @parameterized.parameters(
      (None, None),
      ('', None),
      ('/path/to/file', MOCK_SENTINEL_RETURN_VALUE),
      ('USE_FLAG_VALUE', MOCK_SENTINEL_RETURN_VALUE),
  )
  def test_model_init_function(self, path, expected):
    model = mock.Mock(spec=modeling.DeepVariantModel)
    model.initialize_from_checkpoint.return_value = MOCK_SENTINEL_RETURN_VALUE
    self.assertEqual(expected, model_train.model_init_function(model, 3, path))
    if expected:
      model.initialize_from_checkpoint.assert_called_once_with(
          path, 3, is_training=True)
    else:
      test_utils.assert_not_called_workaround(model.initialize_from_checkpoint)

  @flagsaver.FlagSaver
  @mock.patch('deepvariant.model_train.'
              'tf.train.replica_device_setter')
  @mock.patch('deepvariant.model_train.run')
  def test_main_internal(self, mock_run, mock_device_setter):
    FLAGS.master = 'some_master'
    FLAGS.ps_tasks = 10
    FLAGS.task = 5

    model_train.parse_and_run()

    mock_device_setter.assert_called_once_with(10)
    mock_run.assert_called_once_with('some_master', False, device_fn=mock.ANY)

  @mock.patch('deepvariant.model_train.os.environ')
  @mock.patch('deepvariant.model_train.'
              'tf.train.replica_device_setter')
  @mock.patch('deepvariant.model_train.run')
  def test_main_tfconfig_local(self, mock_run, mock_device_setter,
                               mock_environ):
    mock_environ.get.return_value = '{}'
    model_train.parse_and_run()

    mock_device_setter.assert_called_once_with(0)
    mock_run.assert_called_once_with('', True, device_fn=mock.ANY)

  @parameterized.named_parameters(
      ('master', 'master', 0, True, '/job:master/task:0'),
      ('worker', 'worker', 10, False, '/job:worker/task:10'),
  )
  @mock.patch(
      'deepvariant.model_train.tf.train.Server')
  @mock.patch('deepvariant.model_train.os.environ')
  @mock.patch('deepvariant.model_train.'
              'tf.train.replica_device_setter')
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
        'some-target', expected_is_chief, device_fn=mock.ANY)

  @parameterized.parameters(
      ('master', 'some-master'),
      ('task', 10),
      ('ps_tasks', 5),
  )
  @flagsaver.FlagSaver
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
