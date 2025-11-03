# Copyright 2024 Google LLC.
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

from absl.testing import absltest
import tensorflow as tf

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.small_model import make_small_model_examples
from deepvariant.small_model import small_model_config
from deepvariant.small_model import train_utils
from third_party.nucleus.protos import variants_pb2

# Fake data for testing
MAIN_SAMPLE = 'main_sample'
FAKE_VARIANT_HET = variants_pb2.Variant(
    reference_name='chr9',
    start=5000,
    end=5001,
    reference_bases='A',
    alternate_bases=['C'],
    calls=[],
)
FAKE_VARIANT_CALL_HET = deepvariant_pb2.DeepVariantCall(
    variant=FAKE_VARIANT_HET,
    ref_support_ext=deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
        read_infos=[
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name='read_1',
                mapping_quality=60,
                average_base_quality=30,
                is_reverse_strand=False,
                sample_name=MAIN_SAMPLE,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name='read_2',
                mapping_quality=20,
                average_base_quality=35,
                is_reverse_strand=True,
                sample_name=MAIN_SAMPLE,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name='read_3',
                mapping_quality=40,
                average_base_quality=25,
                is_reverse_strand=True,
                sample_name=MAIN_SAMPLE,
            ),
        ]
    ),
    allele_support_ext={
        'C': deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
            read_infos=[
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name='read_4',
                    mapping_quality=60,
                    average_base_quality=50,
                    is_reverse_strand=False,
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name='read_5',
                    mapping_quality=30,
                    average_base_quality=60,
                    is_reverse_strand=False,
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name='read_6',
                    mapping_quality=60,
                    average_base_quality=40,
                    is_reverse_strand=False,
                    sample_name=MAIN_SAMPLE,
                ),
            ]
        )
    },
    allele_frequency_at_position={
        4997: 0,
        4998: 0,
        4999: 50,
        5000: 87,
        5001: 30,
        5002: 40,
        5003: 0,
    },
)
FAKE_VARIANT_CALL_HET_LABEL = variant_labeler.VariantLabel(
    is_confident=True,
    variant=FAKE_VARIANT_HET,
    genotype=(0, 1),
)


def _get_test_config():
  config = small_model_config.get_config('wgs')
  config.batch_size = 1
  config.epochs = 1
  config.read_ahead_buffer = '1M'
  config.tfrecord_buffer_size = 1
  config.tfrecord_num_parallel_calls = 1
  config.interleave_cycle_length = 1
  config.shuffle_buffer_elements = 1
  config.map_parallel_calls = 1
  config.prefetch_buffer_size = 1
  config.logging_frequency = 1
  config.num_train_samples = 1
  config.num_tune_samples = 1
  config.k_folds = 0
  config.model_params.hidden_layer_sizes = [16]
  config.model_params.vaf_context_window_size = 3
  return config


class TrainUtilsTest(absltest.TestCase):

  def test_e2e_training(self):
    config = _get_test_config()
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=config.model_params.vaf_context_window_size,
            sample_names=[MAIN_SAMPLE],
            expand_by_haplotype=config.model_params.expand_by_haplotype,
        )
    )
    training_examples = small_model_example_factory.encode_training_examples(
        [
            (FAKE_VARIANT_CALL_HET, FAKE_VARIANT_CALL_HET_LABEL),
        ],
        {},
        sample_order=[0],
    )
    temp_dir = self.create_tempdir().full_path
    tfrecord_path = os.path.join(temp_dir, 'examples.tfrecord.gz')
    with tf.io.TFRecordWriter(tfrecord_path, options='GZIP') as writer:
      for example in training_examples:
        writer.write(example.SerializeToString())

    config.train_tfrecord_directory = tfrecord_path
    config.tune_tfrecord_directory = tfrecord_path

    model_runner = train_utils.MirroredModelRunner(
        config,
        working_dir=temp_dir,
        tensorboard_dir=os.path.join(temp_dir, 'tensorboard'),
        warm_start_checkpoint=None,
        is_test=True,
    )
    model_runner.train()


if __name__ == '__main__':
  absltest.main()
