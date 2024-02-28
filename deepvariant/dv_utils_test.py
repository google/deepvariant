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
"""Tests for deepvariant.dv_utils."""



from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tensorflow as tf

from deepvariant import dv_utils
from deepvariant.protos import deepvariant_pb2
from tensorflow.python.platform import gfile
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from tensorflow.core.example import example_pb2


class DVUtilsTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.alts = ['A']
    self.variant = variants_pb2.Variant(
        reference_name='1',
        start=10,
        end=11,
        reference_bases='C',
        alternate_bases=self.alts,
    )
    self.encoded_image = b'encoded_image_data'
    self.default_shape = [5, 5, 7]

  def testFailedExampleImageShape(self):
    # Create an empty example that doesn't have the required image/shape field.
    example = example_pb2.Example()
    with self.assertRaisesRegex(
        ValueError,
        (
            'Invalid image/shape: we expect to find an '
            'image/shape field with length 3.'
        ),
    ):
      dv_utils.example_image_shape(example)

  @parameterized.parameters(
      ('test_shape.gz', 'test_shape.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape@1.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-?????-of-00001.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-*.gz'),
      ('output', 'output'),
      ('test_shape-00000-of-00001', 'test_shape@1'),
      ('test_shape-00000-of-00001', 'test_shape-?????-of-00001'),
      ('test_shape-00000-of-00001', 'test_shape-*'),
  )
  def testGetShapeFromExamplesPath(
      self, file_name_to_write, tfrecord_path_to_match
  ):
    example = example_pb2.Example()
    valid_shape = [1, 2, 3]
    example.features.feature['image/shape'].int64_list.value.extend(valid_shape)
    output_file = test_utils.test_tmpfile(file_name_to_write)
    tfrecord.write_tfrecords([example], output_file, compression_type='GZIP')
    self.assertEqual(
        valid_shape,
        dv_utils.get_shape_from_examples_path(
            test_utils.test_tmpfile(tfrecord_path_to_match)
        ),
    )
    # clean up
    gfile.Remove(output_file)

  @parameterized.parameters(
      ('test_shape.gz', 'test_shape.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape@1.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-?????-of-00001.gz'),
      ('test_shape-00000-of-00001.gz', 'test_shape-*.gz'),
      ('output', 'output'),
      ('test_shape-00000-of-00001', 'test_shape@1'),
      ('test_shape-00000-of-00001', 'test_shape-?????-of-00001'),
      ('test_shape-00000-of-00001', 'test_shape-*'),
  )
  def testGetNoneShapeFromEmptyExamplesPath(
      self, file_name_to_write, tfrecord_path_to_match
  ):
    output_file = test_utils.test_tmpfile(file_name_to_write)
    tfrecord.write_tfrecords([], output_file, compression_type='GZIP')
    self.assertIsNone(
        dv_utils.get_shape_from_examples_path(
            test_utils.test_tmpfile(tfrecord_path_to_match)
        )
    )
    # Clean up
    gfile.Remove(output_file)

  @parameterized.parameters(
      ('/this/path/does/not/exist', '/this/path/does/not'),
      ('/bad/pathA/a,/bad/pathB/b', '/bad/pathA'),
  )
  def testGetShapeFromExamplesPathInvalidPath(
      self, source_paths, expected_partial_message
  ):
    # This calls tf.io.gfile.Glob, which will raise errors.OpError,
    # at least on a Posix filesystem.  Other filesystems might
    # not fail like that, and will return an empty list, which
    # is turned into a different exception.
    with self.assertRaisesRegex(Exception, expected_partial_message):
      dv_utils.get_shape_from_examples_path(source_paths)

  def testStringToIntTensor(self):
    with tf.compat.v1.Session() as sess:
      s = '\001\002\003\004\005\006\007'
      it = dv_utils.string_to_int_tensor(s)
      x = sess.run(it)
      a = x[0]
      self.assertLen(s, a)
      b = list(x[1 : a + 1])
      self.assertEqual(b, [1, 2, 3, 4, 5, 6, 7])

  def testIntTensorToString(self):
    with tf.compat.v1.Session() as sess:
      s = b'\001\002\003\004\005\006\007'
      it = dv_utils.string_to_int_tensor(s)
      x = sess.run(it)
      t = dv_utils.int_tensor_to_string(x)
      self.assertEqual(t, s)

  def test_tfexample_conversion(self):
    # This tests the partial conversion of a CallVariantsOutput proto.
    input_callvariant = (
        b'\nm2\x01G:\x01AR\x19\n\tBAM_FNAME\x12\x0c\n\n\x1a\x08test.bamZ;\x12'
        b'\x12\n\x03VAF\x12\x0b\n\t\x11\x00\x00\x00\x00\x00\x00\xf0?\x12\n\n'
        b'\x02DP\x12\x04\n\x028\x13\x12\x0e\n\x02AD\x12\x08\n\x028\x00\n\x028'
        b'\x13:\x02\x01\x01J\x05HG003h\xae\xea0r\x04chr1\x80\x01\xad\xea0\x12'
        b'\x03\n\x01\x00\x1a\x18\xcc\xe7~\xc5\x98?\xf8>\x00\x00\x00\xc0\x0cZ"?'
        b'\x00\x00\x00\xe0\xa9\xfe\xef?"\x0e\x08\x02 \x01(\x02B\x06d\xfe\xfeF'
        b'\x982'
    )
    expected_example = (
        b'\n\x84\x01\n\x17\n\x0bimage/shape\x12\x08\x1a\x06\n\x04d\xdd\x01\x07'
        b'\n\x0e\n\x05label\x12\x05\x1a\x03\n\x01\x02\n\x1f\n\x05locus\x12\x16'
        b'\n\x14\n\x12chr1:800045-800046\n\x1b\n\rimage/encoded\x12\n\n\x08\n'
        b'\x06d\xfe\xfeF\x982\n\x1b\n\x12alt_allele_indices\x12\x05\x1a\x03\n'
        b'\x01\x00'
    )
    a_cvo = deepvariant_pb2.CallVariantsOutput.FromString(input_callvariant)

    output_example = dv_utils.call_variant_to_tfexample(a_cvo, [100, 221, 7])
    expected = example_pb2.Example.FromString(expected_example)
    self.assertEqual(output_example, expected)

  def test_preprocess_images(self):
    # Create a test input tensor.
    test_input = tf.constant([0, 128, 255], dtype=tf.uint8)
    test_input = tf.reshape(
        test_input, [1, 3]
    )  # Reshaping to simulate an image tensor.

    # Created the corresponding expected output.
    expected_output = tf.constant([-1, 0, 0.9921875], dtype=tf.float32)
    expected_output = tf.reshape(expected_output, [1, 3])

    # Call the preprocess_images function.
    output = dv_utils.preprocess_images(test_input)

    # Check if the output is as expected.
    self.assertTrue(tf.reduce_all(tf.equal(output, expected_output)))

    # Check if the output is in the correct range [-1, 1].
    self.assertTrue(
        tf.reduce_all(output >= -1.0) and tf.reduce_all(output <= 1.0)
    )

    # Check if the output dtype is float32.
    self.assertEqual(output.dtype, tf.float32)

  def test_unpreprocess_images(self):
    # Create a test input array.
    test_input = np.array([-1, 0, 0.9921875], dtype=np.float32)
    # Reshaping to simulate an image array.
    test_input = np.reshape(test_input, [1, 3])

    # Created the corresponding expected output.
    expected_output = np.array([0, 128, 255], dtype=np.float32)
    expected_output = np.reshape(expected_output, [1, 3])

    # Call the unpreprocess_images function.
    output = dv_utils.unpreprocess_images(test_input)

    # Check if the output is as expected.
    np.testing.assert_array_equal(output, expected_output)


if __name__ == '__main__':
  tf.config.run_functions_eagerly(True)
  absltest.main()
