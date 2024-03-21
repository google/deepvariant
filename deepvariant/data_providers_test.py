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
"""Tests for learning.genomics.deepvariant.data_provider."""




from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant import data_providers
from deepvariant import dv_config
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import testdata
from third_party.nucleus.io import tfrecord


def setUpModule():
  testdata.init()


# Return a DeepVariantInput attached to the golden training data.
def get_golden_dataset(
    config=dv_config.get_config('exome+test'),
    mode='train',
):
  return iter(
      data_providers.input_fn(
          path=testdata.GOLDEN_TRAINING_EXAMPLES,
          config=config,
          mode=mode,
      )
  )


class ParseExampleTest(parameterized.TestCase):

  @parameterized.parameters(
      dict(mode='train', expected_keys={'image', 'label', 'sample_weight'}),
      dict(mode='tune', expected_keys={'image', 'label', 'sample_weight'}),
      dict(
          mode='predict',
          expected_keys={
              'image',
              'sample_weight',
              'variant/encoded',
              'alt_allele_indices/encoded',
          },
      ),
      dict(
          mode='debug',
          expected_keys={
              'image',
              'label',
              'locus',
              'variant_type',
              'sample_weight',
              'variant/encoded',
              'alt_allele_indices/encoded',
              'sequencing_type',
          },
      ),
  )
  def test_parse_example(self, mode, expected_keys):
    path = testdata.GOLDEN_TRAINING_EXAMPLES
    ds = tf.data.TFRecordDataset(path, compression_type='GZIP')
    item = ds.take(1).get_single_element()
    input_shape = dv_utils.get_shape_from_examples_path(path)
    config = dv_config.get_config('exome')
    parse_example = data_providers.create_parse_example_fn(
        config,
        mode,
    )
    output = parse_example(item, input_shape)
    self.assertIsInstance(output, dict)
    self.assertSetEqual(set(output.keys()), expected_keys)


class ClassWeightsTest(absltest.TestCase):

  def test_parse_example(self):
    config = dv_config.get_config('exome+test')
    config.batch_size = 128
    config.class_weights = '1,1,10'
    batch = next(get_golden_dataset(config=config, mode='train'))
    _, _, sample_weights = batch
    # Test for sample weights greater than 1
    self.assertGreater(
        tf.reduce_sum(tf.cast(sample_weights > 1, tf.float32)), 0
    )


class CreateExamplesTest(absltest.TestCase):

  def setUp(self):
    super().setUp()
    self.config = dv_config.get_config('exome+test')

  def test_invalid_mode(self):
    with self.assertRaisesRegex(ValueError, 'Mode must be set to'):
      _ = data_providers.input_fn(
          path=testdata.GOLDEN_TRAINING_EXAMPLES,
          config=self.config,
          mode='invalid_mode',
      )

  def test_create_input_dataset(self):
    ds = data_providers.input_fn(
        path=testdata.GOLDEN_TRAINING_EXAMPLES,
        config=self.config,
        mode='train',
    )
    item = ds.take(1).get_single_element()
    self.assertIsInstance(item, tuple)


class DataProviderTest(parameterized.TestCase):

  def test_expected_number_of_loci(self):
    batch = get_golden_dataset(
        config=dv_config.get_config('exome'),
        mode='debug',
    ).get_next()
    loci = set(batch['locus'].numpy())

    expected_loci = {
        example.features.feature['locus'].bytes_list.value[0]
        for example in tfrecord.read_tfrecords(
            testdata.GOLDEN_TRAINING_EXAMPLES, compression_type='GZIP'
        )
    }
    self.assertSetEqual(loci, expected_loci)

  @parameterized.parameters(['train', 'tune', 'predict', 'debug'])
  def test_mode_and_dims(self, mode):
    config = dv_config.get_config('exome+test')
    ds = get_golden_dataset(mode=mode)

    # Check that the batch size matches our config
    if mode in ['train', 'tune']:
      images, labels, sample_weights = ds.get_next()
    else:
      examples = ds.get_next()
      images, labels, sample_weights = (
          examples['image'],
          examples.get('label', None),
          examples['sample_weight'],
      )
    self.assertEqual(config.batch_size, images.shape[0])
    self.assertEqual(config.batch_size, sample_weights.shape[0])

    if labels is not None:
      self.assertEqual(config.batch_size, labels.shape[0])
      # Check that our labels are the correct shape
      self.assertTrue(dv_constants.NUM_CLASSES, labels.shape[1])
      # Check that our labels are one-hot encoded
      self.assertTrue(tf.reduce_all(tf.reduce_sum(labels, axis=1) == 1))


if __name__ == '__main__':
  absltest.main()
