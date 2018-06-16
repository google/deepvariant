# Copyright 2018 Google Inc.  All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License
"""Tests for shuffle_tfrecords_beam.

To run:
1) Install Beam using the instructions at
   https://beam.apache.org/get-started/quickstart-py/

2) Install TensorFlow using the instructions at
   https://www.tensorflow.org/install/

3) Run
$ python path/to/shuffle_tfrecords_beam_test.py
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# We must import tensorflow here; if we import it after the others it
# tries to use the google 3.third_party version, which we do not want.
import tensorflow as tf  # pylint: disable=g-bad-import-order

import tempfile
import unittest

import shuffle_tfrecords_beam


def _write_tfrecords_file(filename, records):
  writer = tf.python_io.TFRecordWriter(filename)
  for r in records:
    writer.write(r)
  writer.close()


def _write_random_tfrecords_file(num_records):
  with tempfile.NamedTemporaryFile(delete=False) as f:
    f.close()
    _write_tfrecords_file(f.name, map(str, range(num_records)))
    return f.name


def _read_tfrecords_file(filename):
  reader = tf.python_io.tf_record_iterator(path=filename)
  return list(reader)


class ShuffleTFRecordsTest(unittest.TestCase):
  """Unittest for shuffle_tfrecords_beam.main."""

  def test_main(self):
    input_tfrecords = _write_random_tfrecords_file(1000)

    fout1 = input_tfrecords + '.out1'
    shuffle_tfrecords_beam.main([
        '--input_pattern_list={}'.format(input_tfrecords),
        '--output_pattern={}'.format(fout1), '--runner=DirectRunner'
    ])
    out1 = _read_tfrecords_file(fout1 + '-00000-of-00001')

    self.assertEqual(1000, len(out1))

    fout2 = input_tfrecords + '.out2'
    shuffle_tfrecords_beam.main([
        '--input_pattern_list={}'.format(input_tfrecords),
        '--output_pattern={}'.format(fout2), '--runner=DirectRunner'
    ])
    out2 = _read_tfrecords_file(fout2 + '-00000-of-00001')

    self.assertEqual(out1, out2)


if __name__ == '__main__':
  unittest.main()
