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
"""Tests for deepvariant.tf_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized
import mock
import tensorflow as tf

from tensorflow.core.example import example_pb2
from deepvariant import test_utils
from deepvariant import tf_utils
from deepvariant.core import io_utils
from deepvariant.core.genomics import variants_pb2


class TFUtilsTest(parameterized.TestCase):

  def setUp(self):
    self.alts = ['A']
    self.variant = variants_pb2.Variant(
        reference_name='1',
        start=10,
        end=11,
        reference_bases='C',
        alternate_bases=self.alts)
    self.encoded_image = 'encoded_image_data'
    self.default_shape = [5, 5, 7]
    self.default_format = 'raw'

  def testModelShapes(self):
    # Builds a graph.
    v0 = tf.Variable([[1, 2, 3], [4, 5, 6]], dtype=tf.float32, name='v0')
    v1 = tf.Variable(
        [[[1], [2]], [[3], [4]], [[5], [6]]], dtype=tf.float32, name='v1')
    init_all_op = tf.initialize_all_variables()
    save = tf.train.Saver({'v0': v0, 'v1': v1})
    save_path = test_utils.test_tmpfile('ckpt_for_debug_string')
    with tf.Session() as sess:
      sess.run(init_all_op)
      # Saves a checkpoint.
      save.save(sess, save_path)

      # Model shapes without any variable requests gives you all variables.
      self.assertEqual({
          'v0': (2, 3),
          'v1': (3, 2, 1)
      }, tf_utils.model_shapes(save_path))
      # Asking for v0 gives you only v0's shape.
      self.assertEqual({'v0': (2, 3)}, tf_utils.model_shapes(save_path, ['v0']))
      # Asking for v1 gives you only v1's shape.
      self.assertEqual({
          'v1': (3, 2, 1)
      }, tf_utils.model_shapes(save_path, ['v1']))

      # Verifies model_shapes() fails for non-existent tensors.
      with self.assertRaisesRegexp(KeyError, 'v3'):
        tf_utils.model_shapes(save_path, ['v3'])

  def testMakeExample(self):
    example = tf_utils.make_example(self.variant, self.alts, self.encoded_image,
                                    self.default_shape, self.default_format)

    self.assertEqual(self.encoded_image,
                     tf_utils.example_encoded_image(example))
    self.assertEqual(
        'raw', example.features.feature['image/format'].bytes_list.value[0])
    self.assertEqual(self.variant, tf_utils.example_variant(example))
    self.assertEqual('1:11-11', tf_utils.example_locus(example))
    self.assertEqual([0], tf_utils.example_alt_alleles_indices(example))
    self.assertEqual('1:11:C->A', tf_utils.example_key(example))

  def testMakeExampleMultiAllelic(self):
    alts = ['AA', 'CC', 'GG']
    self.variant.alternate_bases[:] = alts
    # Providing GG, AA checks that we're sorting the indices.
    example = tf_utils.make_example(self.variant, ['GG', 'AA'], 'foo',
                                    self.default_shape, self.default_format)
    self.assertEqual([0, 2], tf_utils.example_alt_alleles_indices(example))
    self.assertEqual(['AA', 'GG'], tf_utils.example_alt_alleles(example))
    self.assertEqual('1:11:C->AA/GG', tf_utils.example_key(example))

  def testAltAllelesWithVariant(self):
    alts = list(self.variant.alternate_bases)
    example = tf_utils.make_example(self.variant, alts, 'foo',
                                    self.default_shape, self.default_format)
    self.assertEqual([0], tf_utils.example_alt_alleles_indices(example))
    with mock.patch(
        'deepvariant.tf_utils.example_variant'
    ) as mock_ex_variant:
      # Providing variant directly avoids the call to example_variant().
      self.assertEqual(alts,
                       tf_utils.example_alt_alleles(
                           example, variant=self.variant))
      mock_ex_variant.assert_not_called()

      # Checks that we load the variant if needed and that our mock is working.
      mock_ex_variant.return_value = self.variant
      self.assertEqual(alts, tf_utils.example_alt_alleles(example))
      mock_ex_variant.assert_called_once_with(example)

  def assertIsNotAFeature(self, label, example):
    self.assertNotIn(label, example.features.feature)

  def testExampleSetLabel(self):
    example = tf_utils.make_example(self.variant, self.alts, self.encoded_image,
                                    self.default_shape, self.default_format)

    self.assertIsNotAFeature('label', example)
    for label in [0, 1, 2]:
      tf_utils.example_set_label(example, label)
      self.assertEqual(label, tf_utils.example_label(example))

  def testExampleSetTruthVariant(self):
    example = tf_utils.make_example(self.variant, self.alts, self.encoded_image,
                                    self.default_shape, self.default_format)
    full_tvariant = variants_pb2.Variant(
        variant_set_id='variant_set_id',
        id='id',
        names=['name1'],
        created=1234,
        reference_name='1',
        start=10,
        end=11,
        reference_bases='C',
        alternate_bases=['A'],
        filter=['PASS'],
        quality=1234.5,
        calls=[
            variants_pb2.VariantCall(
                call_set_id='call_set_id',
                call_set_name='call_set_name',
                genotype=[0, 1],
                phaseset='phaseset',
                genotype_likelihood=[0.1, 0.2, 0.3])
        ])
    test_utils.set_list_values(full_tvariant.info['key'], [1])
    test_utils.set_list_values(full_tvariant.calls[0].info['key'], [2])

    simple_tvariant = variants_pb2.Variant(
        reference_name='1',
        start=10,
        end=11,
        reference_bases='C',
        alternate_bases=['A'],
        filter=['PASS'],
        quality=1234.5,
        calls=[
            variants_pb2.VariantCall(
                call_set_name='call_set_name', genotype=[0, 1])
        ])
    test_utils.set_list_values(simple_tvariant.calls[0].info['key'], [2])

    self.assertIsNotAFeature('truth_variant/encoded', example)
    tf_utils.example_set_truth_variant(example, full_tvariant, simplify=False)
    self.assertEqual(full_tvariant, tf_utils.example_truth_variant(example))

    # Check that reencoding with simplify=True produces the simplified version.
    tf_utils.example_set_truth_variant(example, full_tvariant, simplify=True)
    self.assertEqual(simple_tvariant, tf_utils.example_truth_variant(example))

  def testExampleImageShape(self):
    example = tf_utils.make_example(self.variant, self.alts, self.encoded_image,
                                    self.default_shape, self.default_format)
    self.assertEqual(self.default_shape, tf_utils.example_image_shape(example))

  def testFailedExampleImageShape(self):
    # Create an empty example that doesn't have the required image/shape field.
    example = example_pb2.Example()
    with self.assertRaisesRegexp(ValueError,
                                 'Invalid image/shape: we expect to find an '
                                 'image/shape field with length 3.'):
      tf_utils.example_image_shape(example)

  @parameterized.parameters(
      ('test_shape.gz', 'test_shape.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape@1.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-?????-of-00001.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-*.gz'), ('output', 'output'),
      ('test_shape-00000-of-00001', 'test_shape@1'),
      ('test_shape-00000-of-00001', 'test_shape-?????-of-00001'),
      ('test_shape-00000-of-00001', 'test_shape-*'))
  def testGetShapeFromExamplesPath(self, file_name_to_write,
                                   tfrecord_path_to_match):
    example = example_pb2.Example()
    valid_shape = [1, 2, 3]
    example.features.feature['image/shape'].int64_list.value.extend(valid_shape)
    output_file = test_utils.test_tmpfile(file_name_to_write)
    io_utils.write_tfrecords([example], output_file)
    tf_utils.get_shape_from_examples_path(
        test_utils.test_tmpfile(tfrecord_path_to_match))

  @parameterized.parameters(
      ('test_shape.gz', 'test_shape.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape@1.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-?????-of-00001.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-*.gz'), ('output', 'output'),
      ('test_shape-00000-of-00001', 'test_shape@1'),
      ('test_shape-00000-of-00001', 'test_shape-?????-of-00001'),
      ('test_shape-00000-of-00001', 'test_shape-*'))
  def testGetNoneShapeFromEmptyExamplesPath(self, file_name_to_write,
                                            tfrecord_path_to_match):
    output_file = test_utils.test_tmpfile(file_name_to_write)
    io_utils.write_tfrecords([], output_file)
    self.assertIsNone(
        tf_utils.get_shape_from_examples_path(
            test_utils.test_tmpfile(tfrecord_path_to_match)))

  def testGetShapeFromExamplesPathInvalidPath(self):
    # This calls tf.gfile.Glob, which will raise errors.OpError,
    # at least on a Posix filesystem.  Other filesystems might
    # not fail like that, and will return an empty list, which
    # is turned into a different exception.
    with self.assertRaisesRegexp(Exception, '/this/path/does/not'):
      tf_utils.get_shape_from_examples_path('/this/path/does/not/exist')


if __name__ == '__main__':
  absltest.main()
