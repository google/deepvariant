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
import tensorflow as tf

from deepvariant import dv_utils
from tensorflow.python.platform import gfile
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from tensorflow.core.example import example_pb2


class TFUtilsTest(parameterized.TestCase):

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

  def testModelShapes(self):
    # Builds a graph.
    v0 = tf.Variable([[1, 2, 3], [4, 5, 6]], dtype=tf.float32, name='v0')
    v1 = tf.Variable(
        [[[1], [2]], [[3], [4]], [[5], [6]]], dtype=tf.float32, name='v1'
    )
    init_all_op = tf.compat.v1.initialize_all_variables()
    save = tf.compat.v1.train.Saver({'v0': v0, 'v1': v1})
    save_path = test_utils.test_tmpfile('ckpt_for_debug_string')
    with tf.compat.v1.Session() as sess:
      sess.run(init_all_op)
      # Saves a checkpoint.
      save.save(sess, save_path)

      # Model shapes without any variable requests gives you all variables.
      self.assertEqual(
          {'v0': (2, 3), 'v1': (3, 2, 1)}, dv_utils.model_shapes(save_path)
      )
      # Asking for v0 gives you only v0's shape.
      self.assertEqual({'v0': (2, 3)}, dv_utils.model_shapes(save_path, ['v0']))
      # Asking for v1 gives you only v1's shape.
      self.assertEqual(
          {'v1': (3, 2, 1)}, dv_utils.model_shapes(save_path, ['v1'])
      )

      # Verifies model_shapes() fails for non-existent tensors.
      with self.assertRaisesRegex(KeyError, 'v3'):
        dv_utils.model_shapes(save_path, ['v3'])

  def testModelNumClasses(self):
    # Builds a graph.
    class_variable_name = 'class_variable_name'
    v0 = tf.Variable([[1, 2, 3]], dtype=tf.int32, name='class_variable_name')
    v1 = tf.Variable(
        [[[1], [2]], [[3], [4]], [[5], [6]]], dtype=tf.float32, name='v1'
    )
    init_all_op = tf.compat.v1.initialize_all_variables()
    save = tf.compat.v1.train.Saver({class_variable_name: v0, 'v1': v1})
    save_path = test_utils.test_tmpfile('ckpt_for_debug_classes')
    with tf.compat.v1.Session() as sess:
      sess.run(init_all_op)
      # Saves a checkpoint.
      save.save(sess, save_path)

      # If you pass in the correct class_variable_name, you'll find the number
      # of classes.
      self.assertEqual(
          3, dv_utils.model_num_classes(save_path, class_variable_name)
      )
      # If the class variable name doesn't existin the checkpoint, return None.
      self.assertEqual(
          None, dv_utils.model_num_classes(save_path, 'non-existent-var')
      )
      # If the checkpoint doesn't exist, return none.
      self.assertIsNone(dv_utils.model_num_classes(None, class_variable_name))

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
    tfrecord.write_tfrecords([example], output_file)
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
    tfrecord.write_tfrecords([], output_file)
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

  def testCompressionTypeOfFiles(self):
    self.assertEqual(
        'GZIP', dv_utils.compression_type_of_files(['/tmp/foo.tfrecord.gz'])
    )
    self.assertIsNone(dv_utils.compression_type_of_files(['/tmp/foo.tfrecord']))


if __name__ == '__main__':
  tf.compat.v1.disable_eager_execution()
  absltest.main()
