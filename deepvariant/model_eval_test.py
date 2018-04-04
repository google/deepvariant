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
"""Tests for genomics.deepvariant.model_eval."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os



from tensorflow import flags
from absl.testing import absltest
from absl.testing import parameterized
import mock
import numpy.testing as npt
import six
import tensorflow as tf

from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import variant_utils
from deepvariant import data_providers_test
from deepvariant import model_eval
from deepvariant import modeling
from deepvariant import pileup_image
from deepvariant import testdata
from deepvariant.testing import flagsaver

FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


class ModelEvalTest(
    six.with_metaclass(parameterized.TestGeneratorMetaclass, tf.test.TestCase)):

  def test_select_variants_weights(self):
    variants = [
        test_utils.make_variant(start=10, alleles=['C', 'T']),
        test_utils.make_variant(start=11, alleles=['C', 'TA']),
        test_utils.make_variant(start=12, alleles=['C', 'A']),
        test_utils.make_variant(start=13, alleles=['CA', 'T']),
    ]
    encoded = tf.constant([v.SerializeToString() for v in variants])

    with self.test_session() as sess:
      sess.run(tf.global_variables_initializer())
      op = model_eval.select_variants_weights(
          variant_utils.is_snp, encoded, name='tf_is_snp')
      self.assertTrue(op.name.startswith('tf_is_snp'))
      npt.assert_array_equal(op.eval(), [1.0, 0.0, 1.0, 0.0])

  def test_calling_metrics(self):

    def make_mock_metric(name):
      # pylint: disable=unused-argument
      def _side_effect(predictions, labels, weights):
        if weights:
          return name + ':' + ','.join(str(int(w)) for w in weights)
        else:
          return name + ':None'

      return mock.MagicMock(side_effect=_side_effect)

    predictions = tf.constant([0, 1, 2, 0])
    labels = tf.constant([0, 2, 1, 1])
    metrics = {
        'm1': make_mock_metric('mock_metric1'),
        'm2': make_mock_metric('mock_metric2')
    }
    selectors = {'s1': [1, 1, 1, 1], 's2': [0, 0, 0, 0], 's3': None}

    # The returned dictionary has the expected keys and values.
    self.assertEqual({
        'm1/s1': 'mock_metric1:1,1,1,1',
        'm1/s2': 'mock_metric1:0,0,0,0',
        'm1/s3': 'mock_metric1:None',
        'm2/s1': 'mock_metric2:1,1,1,1',
        'm2/s2': 'mock_metric2:0,0,0,0',
        'm2/s3': 'mock_metric2:None',
    },
                     model_eval.calling_metrics(
                         metrics_map=metrics,
                         selectors_map=selectors,
                         predictions=predictions,
                         labels=labels))

    # Check that our mocked metrics have all of the calls we.
    for mocked in metrics.values():
      self.assertEqual([
          mock.call(
              predictions=predictions, labels=labels, weights=selectors[x])
          for x in selectors
      ], mocked.call_args_list)

  def setUp(self):
    self.checkpoint_dir = tf.test.get_temp_dir()

  def _write_fake_checkpoint(self,
                             model_name,
                             checkpoint_dir=None,
                             name='model'):
    if checkpoint_dir is None:
      checkpoint_dir = self.checkpoint_dir

    path = os.path.join(checkpoint_dir, name)
    with self.test_session() as sess:
      model = modeling.get_model(model_name)
      # Needed to protect ourselves for models without an input image shape.
      h, w = getattr(model, 'input_image_shape', (100, 221))
      images = tf.placeholder(
          tf.float32, shape=(4, h, w, pileup_image.DEFAULT_NUM_CHANNEL))
      model.create(images, num_classes=3, is_training=True)
      # This is gross, but necessary as model_eval assumes the model was trained
      # with model_train which uses exp moving averages. Unfortunately we cannot
      # just call into model_train as it uses FLAGS which conflict with the
      # flags in use by model_eval. So we inline the creation of the EMA here.
      variable_averages = tf.train.ExponentialMovingAverage(
          FLAGS.moving_average_decay, tf.train.get_or_create_global_step())
      tf.add_to_collection(tf.GraphKeys.UPDATE_OPS,
                           variable_averages.apply(
                               tf.contrib.framework.get_model_variables()))
      sess.run(tf.global_variables_initializer())
      save = tf.train.Saver(tf.contrib.framework.get_variables())
      save.save(sess, path)
    return path

  @parameterized.parameters(['inception_v3'])
  @flagsaver.FlagSaver
  @mock.patch(
      'deepvariant.data_providers.get_dataset')
  def test_end2end(self, model_name, mock_get_dataset):
    """End-to-end test of model_eval."""
    self._write_fake_checkpoint(model_name)

    # Start up eval, loading that checkpoint.
    FLAGS.batch_size = 2
    FLAGS.checkpoint_dir = self.checkpoint_dir
    FLAGS.eval_dir = tf.test.get_temp_dir()
    FLAGS.max_evaluations = 1
    FLAGS.max_examples = 2
    FLAGS.model_name = model_name
    FLAGS.dataset_config_pbtxt = '/path/to/mock.pbtxt'
    # Always try to read in compressed inputs to stress that case. Uncompressed
    # inputs are certain to work. This test is expensive to run, so we want to
    # minimize the number of times we need to run this.
    mock_get_dataset.return_value = data_providers_test.make_golden_dataset(
        compressed_inputs=True)
    model_eval.main(0)
    mock_get_dataset.assert_called_once_with(FLAGS.dataset_config_pbtxt)

  @flagsaver.FlagSaver
  @mock.patch(
      'deepvariant.model_eval.checkpoints_iterator')
  @mock.patch(
      'deepvariant.data_providers.get_dataset')
  def test_fixed_eval_sees_the_same_evals(self, mock_get_dataset,
                                          mock_checkpoints_iterator):
    dataset = data_providers_test.make_golden_dataset()
    n_checkpoints = 3
    checkpoints = [
        self._write_fake_checkpoint('constant', name='model' + str(i))
        for i in range(n_checkpoints)
    ]

    # Setup our mocks.
    mock_checkpoints_iterator.return_value = checkpoints
    mock_get_dataset.return_value = dataset

    # Start up eval, loading that checkpoint.
    FLAGS.batch_size = 2
    FLAGS.checkpoint_dir = self.checkpoint_dir
    FLAGS.eval_dir = tf.test.get_temp_dir()
    FLAGS.max_evaluations = n_checkpoints
    FLAGS.model_name = 'constant'
    FLAGS.dataset_config_pbtxt = '/path/to/mock.pbtxt'
    model_eval.main(0)

    self.assertEqual(mock_get_dataset.call_args_list,
                     [mock.call(FLAGS.dataset_config_pbtxt)] * n_checkpoints)

    metrics = [
        model_eval.read_metrics(checkpoint, eval_dir=FLAGS.eval_dir)
        for checkpoint in checkpoints
    ]

    # Check that our metrics are what we expect them to be.
    # See b/62864044 for details on how to compute these counts:
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
        # One of our labels is 0, which we call het, giving us 1 FP.
        'FPs/All': 1.0,
        # We call everything as het, so the recall has to be 1.
        'Recall/All': 1.0,
        # redacted
        # # We don't call anything but hets, so TNs has to be 0.
        # 'TNs/All': 0,
        # We find all positives, so this has to be 47.
        'TPs/All': 47,
    }
    for key, expected_value in expected_values_for_all_exact.iteritems():
      self.assertEqual(metrics[0][key], expected_value)

    expected_values_for_all_close = {
        # We called 47 / 48 correctly.
        'Precision/All': 47. / 48,
    }
    for key, expected_value in expected_values_for_all_close.iteritems():
      self.assertAlmostEqual(metrics[0][key], expected_value, places=6)

    for m1, m2 in zip(metrics, metrics[1:]):
      self.assertEqual(m1, m2)


if __name__ == '__main__':
  absltest.main()
