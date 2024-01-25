# Copyright 2023 Google LLC.
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
import os

from absl import flags
from absl.testing import flagsaver
from absl.testing import parameterized
import ml_collections
import tensorflow as tf

from deepvariant import call_variants
from deepvariant import dv_utils
from deepvariant import keras_modeling
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2
from absl import absltest
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord


FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


class CallVariantsTest(parameterized.TestCase):
  """Test cases for call variants."""

  @parameterized.named_parameters(
      dict(
          testcase_name="end2end_default",
          model="inception_v3",
          shard_input=False,
      ),
      dict(
          testcase_name="end2end_shard_input",
          model="inception_v3",
          shard_input=True,
      ),
      dict(
          testcase_name="end2end_no_debug_info_no_shard",
          model="inception_v3",
          shard_input=False,
      ),
      dict(
          testcase_name="end2end_shard_input_with_debug",
          model="inception_v3",
          shard_input=True,
          include_debug_info=True,
      ),
  )
  @flagsaver.flagsaver
  def test_call_variants_end2end(
      self,
      model,
      shard_input,
      include_debug_info=False,
  ):
    # Load in test data and get input shape
    calling_testdata_path = testdata.GOLDEN_CALLING_EXAMPLES
    example_info_json_path = dv_utils.get_example_info_json_filename(
        calling_testdata_path, None
    )
    if shard_input:
      # Input a sharded version of our golden examples
      calling_testdata_path = testdata.GOLDEN_CALLING_EXAMPLES_SHARDED
      example_info_json_path = dv_utils.get_example_info_json_filename(
          calling_testdata_path, 0
      )

    # Load and save a model with random weights
    input_checkpoint_dir = os.path.join(
        self.create_tempdir("input"), "saved_model"
    )
    tf.io.gfile.makedirs(input_checkpoint_dir)
    tf.io.gfile.copy(
        example_info_json_path,
        os.path.join(input_checkpoint_dir, "example_info.json"),
        overwrite=True,
    )
    input_shape = dv_utils.get_shape_from_examples_path(calling_testdata_path)

    config = ml_collections.ConfigDict()
    with config.unlocked() as config:
      config.model_type = model
    model = keras_modeling.get_model(config)(
        input_shape, weights=None, init_backbone_with_imagenet=False
    )
    model.save(input_checkpoint_dir)

    # set up output directory
    output_dir = self.create_tempdir()
    output_tfrecord = os.path.join(output_dir, "output.tfrecord.gz")

    # Run end to end variant calling
    FLAGS.batch_size = 4
    FLAGS.include_debug_info = include_debug_info
    FLAGS.outfile = output_tfrecord
    FLAGS.examples = calling_testdata_path
    FLAGS.checkpoint = input_checkpoint_dir
    call_variants.main()

    # Assert
    sharded_output_files = tf.io.gfile.listdir(output_dir)
    self.assertLen(sharded_output_files, 1)

    only_sharded_output_filepath = os.path.join(
        output_dir, sharded_output_files[0]
    )
    call_variants_outputs = list(
        tfrecord.read_tfrecords(
            only_sharded_output_filepath,
            proto=deepvariant_pb2.CallVariantsOutput,
        )
    )

    # Check that we have the right number of output protos
    if shard_input:
      num_examples = len(
          list(
              tf.data.Dataset.list_files(
                  sharded_file_utils.normalize_to_sharded_file_pattern(
                      calling_testdata_path
                  ),
                  shuffle=False,
              ).interleave(
                  lambda file: tf.data.TFRecordDataset(
                      file, compression_type="GZIP"
                  ),
                  num_parallel_calls=tf.data.AUTOTUNE,
              )
          )
      )
    else:
      num_examples = len(list(tfrecord.read_tfrecords(calling_testdata_path)))
    self.assertLen(call_variants_outputs, num_examples)
    if include_debug_info:
      self.assertNotEmpty(call_variants_outputs[0].debug_info.image_encoded)

  @parameterized.named_parameters(
      dict(
          testcase_name="round_gls (precision=None)",
          precision=None,
          expected_value=[0.2102311329, 0.099768768, 0.6899999991],
      ),
      dict(
          testcase_name="round_gls (precision=1)",
          precision=1,
          expected_value=[0.2, 0.1, 0.7],
      ),
      dict(
          testcase_name="round_gl (precision=2)",
          precision=2,
          expected_value=[0.21, 0.10, 0.69],
      ),
  )
  def test_round_gls_with_precision(self, precision, expected_value):
    input_value = [0.2102311329, 0.099768768, 0.6899999991]

    actual_value = call_variants.round_gls(input_value, precision)

    self.assertEqual(actual_value, expected_value)


if __name__ == "__main__":
  absltest.main()
