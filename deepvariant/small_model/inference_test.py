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
from unittest import mock
from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf
from deepvariant.protos import deepvariant_pb2
from deepvariant.small_model import inference
from deepvariant.small_model import make_small_model_examples
from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2


FAKE_VARIANT_CALL_1 = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_name="chr9",
        reference_bases="A",
        alternate_bases=["C"],
        calls=[variants_pb2.VariantCall(genotype=[-1])],
    )
)
FAKE_VARIANT_CALL_2 = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_name="chr10",
        reference_bases="T",
        alternate_bases=["C"],
    )
)
FAKE_VARIANT_CALL_INSERTION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_name="chr11",
        reference_bases="T",
        alternate_bases=["CCT"],
        calls=[variants_pb2.VariantCall(genotype=[-1])],
    )
)
FAKE_VARIANT_CALL_DELETION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_name="chr11",
        reference_bases="TCC",
        alternate_bases=["T"],
    )
)
FAKE_EXAMPLE = tf.train.Example(
    features=tf.train.Features(
        feature={
            make_small_model_examples.FEATURES_ENCODED: tf.train.Feature(
                int64_list=tf.train.Int64List(value=[])
            ),
        }
    )
)


class SmallModelVariantCallerTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.mock_classifier = mock.MagicMock()
    self.mock_classifier.predict_on_batch.return_value = [
        (0.0, 0.999, 0.0),
        (0.0, 0.0, 0.1),
        (0.0, 0.0, 0.999),
        (0.1, 0.0, 0.0),
    ]
    self.variant_caller = inference.SmallModelVariantCaller(
        classifier=self.mock_classifier,
        snp_gq_threshold=10,
        indel_gq_threshold=10,
        batch_size=4,
    )

  def test_call_variants(self):
    call_variant_outputs, filtered_candidates = (
        self.variant_caller.call_variants(
            candidates_with_alt_allele_indices=[
                (
                    FAKE_VARIANT_CALL_1,
                    make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
                ),
                (
                    FAKE_VARIANT_CALL_2,
                    make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
                ),
                (
                    FAKE_VARIANT_CALL_INSERTION,
                    make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
                ),
                (
                    FAKE_VARIANT_CALL_DELETION,
                    make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
                ),
            ],
            examples=[
                tf.train.Example(),
                tf.train.Example(),
                tf.train.Example(),
                tf.train.Example(),
            ],
        )
    )
    self.assertEqual(
        call_variant_outputs,
        [
            deepvariant_pb2.CallVariantsOutput(
                variant=variants_pb2.Variant(
                    reference_name=FAKE_VARIANT_CALL_1.variant.reference_name,
                    reference_bases="A",
                    alternate_bases=["C"],
                    calls=[
                        variants_pb2.VariantCall(
                            genotype=[-1],
                            info={
                                "MID": struct_pb2.ListValue(
                                    values=[
                                        struct_pb2.Value(
                                            string_value=inference.SMALL_MODEL_ID
                                        )
                                    ]
                                )
                            },
                        )
                    ],
                ),
                genotype_probabilities=[0, 0.999, 0],
                debug_info=None,
                alt_allele_indices=deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
                    indices=[0]
                ),
            ),
            deepvariant_pb2.CallVariantsOutput(
                variant=variants_pb2.Variant(
                    reference_name=FAKE_VARIANT_CALL_INSERTION.variant.reference_name,
                    reference_bases="T",
                    alternate_bases=["CCT"],
                    calls=[
                        variants_pb2.VariantCall(
                            genotype=[-1],
                            info={
                                "MID": struct_pb2.ListValue(
                                    values=[
                                        struct_pb2.Value(
                                            string_value=inference.SMALL_MODEL_ID
                                        )
                                    ]
                                )
                            },
                        )
                    ],
                ),
                genotype_probabilities=[0, 0, 0.999],
                debug_info=None,
                alt_allele_indices=deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
                    indices=[0]
                ),
            ),
        ],
    )
    self.assertEqual(
        filtered_candidates, [FAKE_VARIANT_CALL_2, FAKE_VARIANT_CALL_DELETION]
    )

  def test_call_variants_small_batch(self):
    with mock.patch.object(
        type(self.mock_classifier),
        "__call__",
        lambda self, x, training: [
            (0.0, 0.999, 0.0),
            (0.0, 0.0, 0.1),
        ],
    ):
      call_variant_outputs, filtered_candidates = (
          self.variant_caller.call_variants(
              candidates_with_alt_allele_indices=[
                  (
                      FAKE_VARIANT_CALL_1,
                      make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
                  ),
                  (
                      FAKE_VARIANT_CALL_2,
                      make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
                  ),
              ],
              examples=[
                  tf.train.Example(),
                  tf.train.Example(),
              ],
          )
      )
    self.assertLen(call_variant_outputs, 1)
    self.assertLen(filtered_candidates, 1)


if __name__ == "__main__":
  absltest.main()
