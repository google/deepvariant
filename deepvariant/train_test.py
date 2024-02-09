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

from absl.testing import flagsaver
from absl.testing import parameterized
from google.protobuf import text_format
import tensorflow as tf

from deepvariant import dv_config
from deepvariant import testdata
from deepvariant import train
from deepvariant.protos import deepvariant_pb2


def setUpModule():
  testdata.init()


class TrainTest(parameterized.TestCase):
  """Test cases training deepvariant."""

  @parameterized.named_parameters(
      dict(
          testcase_name="deepvariant_base_train_test",
          config_name="base:test",
          testdata_tfrecord_filename="golden.training_examples.tfrecord.gz",
          num_examples=49,
      ),
      dict(
          testcase_name="deepvariant_exome_train_test",
          config_name="exome:test",
          testdata_tfrecord_filename="golden.training_examples.tfrecord.gz",
          num_examples=49,
      ),
  )
  def test_train_and_evaluate(
      self, config_name, testdata_tfrecord_filename, num_examples
  ):
    config = dv_config.get_config(config_name)

    # Create deepvariant dataset config for training
    # Filepath names need to be resolved at runtime during tests o ensure they
    # can be opened and used for tests. Creating a temporary directory with the
    # dataset config defined below at runtime allows the test logic to be
    # contained here and allows us to use testdata training data.
    dataset_config = deepvariant_pb2.DeepVariantDatasetConfig(
        name=config_name,
        num_examples=num_examples,
        tfrecord_path=testdata.deepvariant_testdata(testdata_tfrecord_filename),
    )
    dataset_config_dir = os.path.join(
        self.create_tempdir("input"), "testdata.dataset_config.pbtxt"
    )
    with tf.io.gfile.GFile(dataset_config_dir, "w") as f:
      f.write(text_format.MessageToString(dataset_config))
    with config.unlocked() as config:
      config.train_dataset_pbtxt = dataset_config_dir
      config.tune_dataset_pbtxt = dataset_config_dir

    # Create temporary output path to contain output files from traiing
    # including model checkpoints.
    experiment_dir = self.create_tempdir().full_path

    # Run the model using specified config and flags
    with flagsaver.as_parsed(xm_runlocal="True", experiment_dir=experiment_dir):
      train.train(config)

    # Asset checkpoints and example_info are present as expected in the output
    # directory.
    test_output_directory = os.path.join(experiment_dir, "wu_-1")
    self.assertNotEmpty(
        tf.io.gfile.listdir(os.path.join(test_output_directory, "checkpoints"))
    )
    self.assertTrue(
        tf.io.gfile.exists(
            os.path.join(
                test_output_directory, "checkpoints", "example_info.json"
            )
        )
    )


if __name__ == "__main__":
  tf.test.main()
