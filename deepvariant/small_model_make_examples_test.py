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
"""Tests for make_candidate_summaries."""

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import small_model_make_examples
from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2

FAKE_VARIANT_HET = variants_pb2.Variant(
    reference_name="chr9",
    start=5000,
    end=5001,
    reference_bases="A",
    alternate_bases=["C"],
    calls=[
        variants_pb2.VariantCall(
            info={
                "AD": struct_pb2.ListValue(
                    values=[
                        struct_pb2.Value(int_value=14),
                        struct_pb2.Value(int_value=16),
                    ]
                ),
                "DP": struct_pb2.ListValue(
                    values=[struct_pb2.Value(int_value=30)]
                ),
                "VAF": struct_pb2.ListValue(
                    values=[struct_pb2.Value(number_value=0.875)]
                ),
            }
        )
    ],
)
FAKE_VARIANT_CALL_HET = deepvariant_pb2.DeepVariantCall(
    variant=FAKE_VARIANT_HET
)
FAKE_VARIANT_CALL_HET_LABEL = variant_labeler.VariantLabel(
    is_confident=True,
    variant=FAKE_VARIANT_HET,
    genotype=(0, 1),
)


class SmallModelMakeExamplesTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name="contig_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.CONTIG,
          expected_value="chr9",
      ),
      dict(
          testcase_name="start_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.START,
          expected_value=5000,
      ),
      dict(
          testcase_name="end_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.END,
          expected_value=5001,
      ),
      dict(
          testcase_name="ref_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.REF,
          expected_value="A",
      ),
      dict(
          testcase_name="alt_1_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.ALT_1,
          expected_value="C",
      ),
      dict(
          testcase_name="num_reads_supports_ref_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.NUM_READS_SUPPORTS_REF,
          expected_value=14,
      ),
      dict(
          testcase_name="num_reads_supports_alt_1_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.NUM_READS_SUPPORTS_ALT_1,
          expected_value=16,
      ),
      dict(
          testcase_name="total_depth_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.TOTAL_DEPTH,
          expected_value=30,
      ),
      dict(
          testcase_name="variant_allele_frequency_1_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.VARIANT_ALLELE_FREQUENCY_1,
          expected_value=87,
      ),
      dict(
          testcase_name="genotype_feature",
          candidate=FAKE_VARIANT_CALL_HET,
          label=FAKE_VARIANT_CALL_HET_LABEL,
          feature=small_model_make_examples.SmallModelFeature.GENOTYPE,
          expected_value=small_model_make_examples.GenotypeEncoding.HET.value,
      ),
  )
  def test_get_feature_from_candidate_or_label(
      self, candidate, label, feature, expected_value
  ):
    self.assertEqual(
        small_model_make_examples.get_feature_from_candidate_or_label(
            feature,
            candidate,
            label,
        ),
        expected_value,
    )

  def test_get_example_feature_columns(self):
    # Note: for now, testing only the identifying columns and the first four
    # features since the features will evolve.
    num_test_columns = len(small_model_make_examples.IDENTIFYING_FEATURES) + 4
    expected_columns = [
        "contig",
        "start",
        "end",
        "ref",
        "alt_1",
        "num_reads_supports_ref",
        "num_reads_supports_alt_1",
        "total_depth",
        "variant_allele_frequency_1",
    ]

    example_feature_columns = (
        small_model_make_examples.get_example_feature_columns()
    )

    self.assertEqual(
        example_feature_columns[:num_test_columns],
        expected_columns,
    )
    # Last column should always be the truth feature.
    self.assertEqual(example_feature_columns[-1], "genotype")

  def test_generate_training_examples(self):
    num_test_columns = len(small_model_make_examples.IDENTIFYING_FEATURES) + 4

    small_model_examples = (
        small_model_make_examples.generate_training_examples([
            (FAKE_VARIANT_CALL_HET, FAKE_VARIANT_CALL_HET_LABEL),
        ])
    )
    self.assertEqual(
        small_model_examples[0][:num_test_columns],
        [
            "chr9",
            5000,
            5001,
            "A",
            "C",
            14,
            16,
            30,
            87,
        ],
    )
    # Last column should be truth value
    self.assertEqual(
        small_model_examples[0][-1],
        small_model_make_examples.GenotypeEncoding.HET.value,
    )


if __name__ == "__main__":
  absltest.main()
