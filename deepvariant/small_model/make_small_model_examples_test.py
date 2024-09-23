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
from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.small_model import make_small_model_examples
from third_party.nucleus.protos import variants_pb2

FAKE_VARIANT_HET = variants_pb2.Variant(
    reference_name="chr9",
    start=5000,
    end=5001,
    reference_bases="A",
    alternate_bases=["C"],
    calls=[],
)
FAKE_VARIANT_CALL_HET = deepvariant_pb2.DeepVariantCall(
    variant=FAKE_VARIANT_HET,
    ref_support_ext=deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
        read_infos=[
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_1",
                mapping_quality=60,
                average_base_quality=30,
                is_reverse_strand=False,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_2",
                mapping_quality=20,
                average_base_quality=35,
                is_reverse_strand=True,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_3",
                mapping_quality=40,
                average_base_quality=25,
                is_reverse_strand=True,
            ),
        ]
    ),
    allele_support_ext={
        "C": deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
            read_infos=[
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_4",
                    mapping_quality=60,
                    average_base_quality=50,
                    is_reverse_strand=False,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_5",
                    mapping_quality=30,
                    average_base_quality=60,
                    is_reverse_strand=False,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_6",
                    mapping_quality=60,
                    average_base_quality=40,
                    is_reverse_strand=False,
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

FAKE_VARIANT_INSERTION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_bases="A",
        alternate_bases=["AAC"],
    )
)
FAKE_VARIANT_DELETION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_bases="AACC",
        alternate_bases=["C"],
    )
)
READ_PHASES = {
    "read_1": 1,
    "read_2": 2,
    "read_3": 0,
    "read_4": 1,
    "read_5": 1,
    "read_6": 2,
}


class SmallModelMakeExamplesTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name="contig_feature",
          feature=make_small_model_examples.SmallModelFeature.CONTIG.value,
          expected_value="chr9",
      ),
      dict(
          testcase_name="start_feature",
          feature=make_small_model_examples.SmallModelFeature.START.value,
          expected_value=5000,
      ),
      dict(
          testcase_name="end_feature",
          feature=make_small_model_examples.SmallModelFeature.END.value,
          expected_value=5001,
      ),
      dict(
          testcase_name="ref_feature",
          feature=make_small_model_examples.SmallModelFeature.REF.value,
          expected_value="A",
      ),
      dict(
          testcase_name="alt_1_feature",
          feature=make_small_model_examples.SmallModelFeature.ALT_1.value,
          expected_value="C",
      ),
      dict(
          testcase_name="num_reads_supports_ref_feature",
          feature=make_small_model_examples.SmallModelFeature.NUM_READS_SUPPORTS_REF.value,
          expected_value=3,
      ),
      dict(
          testcase_name="num_reads_supports_alt_1_feature",
          feature=make_small_model_examples.SmallModelFeature.NUM_READS_SUPPORTS_ALT_1.value,
          expected_value=3,
      ),
      dict(
          testcase_name="total_depth_feature",
          feature=make_small_model_examples.SmallModelFeature.TOTAL_DEPTH.value,
          expected_value=6,
      ),
      dict(
          testcase_name="variant_allele_frequency_1_feature",
          feature=make_small_model_examples.SmallModelFeature.VARIANT_ALLELE_FREQUENCY_1.value,
          expected_value=50,
      ),
      dict(
          testcase_name="genotype_feature",
          feature=make_small_model_examples.SmallModelFeature.GENOTYPE.value,
          expected_value=make_small_model_examples.GenotypeEncoding.HET.value,
      ),
      dict(
          testcase_name="ref_mapping_quality_feature",
          feature=make_small_model_examples.SmallModelFeature.REF_MAPPING_QUALITY.value,
          expected_value=40,
      ),
      dict(
          testcase_name="alt_1_mapping_quality_feature",
          feature=make_small_model_examples.SmallModelFeature.ALT_1_MAPPING_QUALITY.value,
          expected_value=50,
      ),
      dict(
          testcase_name="ref_base_quality_feature",
          feature=make_small_model_examples.SmallModelFeature.REF_BASE_QUALITY.value,
          expected_value=30,
      ),
      dict(
          testcase_name="alt_1_base_quality_feature",
          feature=make_small_model_examples.SmallModelFeature.ALT_1_BASE_QUALITY.value,
          expected_value=50,
      ),
      dict(
          testcase_name="variant_allele_frequency_at_minus_3_feature",
          feature="variant_allele_frequency_at_minus_3",
          expected_value=0,
      ),
      dict(
          testcase_name="variant_allele_frequency_at_minus_1_feature",
          feature="variant_allele_frequency_at_minus_1",
          expected_value=50,
      ),
      dict(
          testcase_name="variant_allele_frequency_at_plus_0_feature",
          feature="variant_allele_frequency_at_plus_0",
          expected_value=87,
      ),
      dict(
          testcase_name="variant_allele_frequency_at_plus_1_feature",
          feature="variant_allele_frequency_at_plus_1",
          expected_value=30,
      ),
      dict(
          testcase_name="variant_allele_frequency_at_plus_2_feature",
          feature="variant_allele_frequency_at_plus_2",
          expected_value=40,
      ),
      dict(
          testcase_name="is_snp_feature",
          feature=make_small_model_examples.SmallModelFeature.IS_SNP.value,
          expected_value=1,
      ),
      dict(
          testcase_name="is_insertion",
          feature=make_small_model_examples.SmallModelFeature.IS_INSERTION.value,
          expected_value=0,
      ),
      dict(
          testcase_name="is_deletion_feature",
          feature=make_small_model_examples.SmallModelFeature.IS_DELETION.value,
          expected_value=0,
      ),
      dict(
          testcase_name="insertion_length_feature",
          feature=make_small_model_examples.SmallModelFeature.INSERTION_LENGTH.value,
          expected_value=0,
      ),
      dict(
          testcase_name="deletion_length_feature",
          feature=make_small_model_examples.SmallModelFeature.DELETION_LENGTH.value,
          expected_value=0,
      ),
      dict(
          testcase_name="ref_reverse_strand_ratio_feature",
          feature=make_small_model_examples.SmallModelFeature.REF_REVERSE_STRAND_RATIO.value,
          expected_value=66,
      ),
      dict(
          testcase_name="alt_1_reverse_strand_ratio_feature",
          feature=make_small_model_examples.SmallModelFeature.ALT_1_REVERSE_STRAND_RATIO.value,
          expected_value=0,
      ),
  )
  def test_feature_encoder_encode_feature(self, feature, expected_value):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            FAKE_VARIANT_CALL_HET, FAKE_VARIANT_CALL_HET_LABEL
        ).encode_feature(
            feature,
        ),
        expected_value,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name="is_snp_feature_insertion",
          feature=make_small_model_examples.SmallModelFeature.IS_SNP.value,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_snp_feature_deletion",
          feature=make_small_model_examples.SmallModelFeature.IS_SNP.value,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_insertion_feature_insertion",
          feature=make_small_model_examples.SmallModelFeature.IS_INSERTION.value,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=1,
      ),
      dict(
          testcase_name="is_insertion_feature_deletion",
          feature=make_small_model_examples.SmallModelFeature.IS_INSERTION.value,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_deletion_feature_insertion",
          feature=make_small_model_examples.SmallModelFeature.IS_DELETION.value,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_deletion_feature_deletion",
          feature=make_small_model_examples.SmallModelFeature.IS_DELETION.value,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=1,
      ),
      dict(
          testcase_name="insertion_length_feature_insertion",
          feature=make_small_model_examples.SmallModelFeature.INSERTION_LENGTH.value,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=2,
      ),
      dict(
          testcase_name="insertion_length_feature_deletion",
          feature=make_small_model_examples.SmallModelFeature.INSERTION_LENGTH.value,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="deletion_length_feature_insertion",
          feature=make_small_model_examples.SmallModelFeature.DELETION_LENGTH.value,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=0,
      ),
      dict(
          testcase_name="deletion_length_feature_deletion",
          feature=make_small_model_examples.SmallModelFeature.DELETION_LENGTH.value,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=3,
      ),
  )
  def test_encode_indel_feature(self, feature, candidate, expected_value):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            candidate, label=None
        ).encode_feature(
            feature,
        ),
        expected_value,
    )

  @parameterized.parameters(
      dict(
          window_size=0,
          expand_by_haplotype=False,
          expected_columns=[
              "contig",
              "start",
              "end",
              "ref",
              "alt_1",
              "genotype",
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
          ],
      ),
      dict(
          window_size=3,
          expand_by_haplotype=False,
          expected_columns=[
              "contig",
              "start",
              "end",
              "ref",
              "alt_1",
              "genotype",
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
          ],
      ),
      dict(
          window_size=11,
          expand_by_haplotype=False,
          expected_columns=[
              "contig",
              "start",
              "end",
              "ref",
              "alt_1",
              "genotype",
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "variant_allele_frequency_at_minus_5",
              "variant_allele_frequency_at_minus_4",
              "variant_allele_frequency_at_minus_3",
              "variant_allele_frequency_at_minus_2",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
              "variant_allele_frequency_at_plus_2",
              "variant_allele_frequency_at_plus_3",
              "variant_allele_frequency_at_plus_4",
              "variant_allele_frequency_at_plus_5",
          ],
      ),
      dict(
          window_size=3,
          expand_by_haplotype=True,
          expected_columns=[
              "contig",
              "start",
              "end",
              "ref",
              "alt_1",
              "genotype",
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
              "num_reads_supports_ref_hp_0",
              "num_reads_supports_alt_1_hp_0",
              "total_depth_hp_0",
              "variant_allele_frequency_1_hp_0",
              "ref_mapping_quality_hp_0",
              "alt_1_mapping_quality_hp_0",
              "ref_base_quality_hp_0",
              "alt_1_base_quality_hp_0",
              "ref_reverse_strand_ratio_hp_0",
              "alt_1_reverse_strand_ratio_hp_0",
              "num_reads_supports_ref_hp_1",
              "num_reads_supports_alt_1_hp_1",
              "total_depth_hp_1",
              "variant_allele_frequency_1_hp_1",
              "ref_mapping_quality_hp_1",
              "alt_1_mapping_quality_hp_1",
              "ref_base_quality_hp_1",
              "alt_1_base_quality_hp_1",
              "ref_reverse_strand_ratio_hp_1",
              "alt_1_reverse_strand_ratio_hp_1",
              "num_reads_supports_ref_hp_2",
              "num_reads_supports_alt_1_hp_2",
              "total_depth_hp_2",
              "variant_allele_frequency_1_hp_2",
              "ref_mapping_quality_hp_2",
              "alt_1_mapping_quality_hp_2",
              "ref_base_quality_hp_2",
              "alt_1_base_quality_hp_2",
              "ref_reverse_strand_ratio_hp_2",
              "alt_1_reverse_strand_ratio_hp_2",
          ],
      ),
  )
  def test_small_model_example_factory_training_features(
      self, window_size, expand_by_haplotype, expected_columns
  ):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=window_size,
            expand_by_haplotype=expand_by_haplotype,
        )
    )
    self.assertEqual(
        small_model_example_factory.training_features,
        expected_columns,
    )

  @parameterized.parameters(
      dict(
          window_size=0,
          expand_by_haplotype=False,
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
          ],
      ),
      dict(
          window_size=3,
          expand_by_haplotype=False,
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
          ],
      ),
      dict(
          window_size=11,
          expand_by_haplotype=False,
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "variant_allele_frequency_at_minus_5",
              "variant_allele_frequency_at_minus_4",
              "variant_allele_frequency_at_minus_3",
              "variant_allele_frequency_at_minus_2",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
              "variant_allele_frequency_at_plus_2",
              "variant_allele_frequency_at_plus_3",
              "variant_allele_frequency_at_plus_4",
              "variant_allele_frequency_at_plus_5",
          ],
      ),
      dict(
          window_size=3,
          expand_by_haplotype=True,
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt_1",
              "total_depth",
              "variant_allele_frequency_1",
              "ref_mapping_quality",
              "alt_1_mapping_quality",
              "ref_base_quality",
              "alt_1_base_quality",
              "ref_reverse_strand_ratio",
              "alt_1_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
              "num_reads_supports_ref_hp_0",
              "num_reads_supports_alt_1_hp_0",
              "total_depth_hp_0",
              "variant_allele_frequency_1_hp_0",
              "ref_mapping_quality_hp_0",
              "alt_1_mapping_quality_hp_0",
              "ref_base_quality_hp_0",
              "alt_1_base_quality_hp_0",
              "ref_reverse_strand_ratio_hp_0",
              "alt_1_reverse_strand_ratio_hp_0",
              "num_reads_supports_ref_hp_1",
              "num_reads_supports_alt_1_hp_1",
              "total_depth_hp_1",
              "variant_allele_frequency_1_hp_1",
              "ref_mapping_quality_hp_1",
              "alt_1_mapping_quality_hp_1",
              "ref_base_quality_hp_1",
              "alt_1_base_quality_hp_1",
              "ref_reverse_strand_ratio_hp_1",
              "alt_1_reverse_strand_ratio_hp_1",
              "num_reads_supports_ref_hp_2",
              "num_reads_supports_alt_1_hp_2",
              "total_depth_hp_2",
              "variant_allele_frequency_1_hp_2",
              "ref_mapping_quality_hp_2",
              "alt_1_mapping_quality_hp_2",
              "ref_base_quality_hp_2",
              "alt_1_base_quality_hp_2",
              "ref_reverse_strand_ratio_hp_2",
              "alt_1_reverse_strand_ratio_hp_2",
          ],
      ),
  )
  def test_small_model_example_factory_inference_features(
      self, window_size, expand_by_haplotype, expected_columns
  ):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=window_size,
            expand_by_haplotype=expand_by_haplotype,
        )
    )
    self.assertEqual(
        small_model_example_factory.inference_features,
        expected_columns,
    )

  def test_encode_training_examples(self):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=3
        )
    )
    small_model_examples = small_model_example_factory.encode_training_examples(
        [
            (FAKE_VARIANT_CALL_HET, FAKE_VARIANT_CALL_HET_LABEL),
        ],
        {},
    )
    self.assertEqual(
        small_model_examples[0],
        [
            "chr9",
            5000,
            5001,
            "A",
            "C",
            1,
            3,
            3,
            6,
            50,
            40,
            50,
            30,
            50,
            66,
            0,
            1,
            0,
            0,
            0,
            0,
            50,
            87,
            30,
        ],
    )

  def test_encode_inference_examples(self):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=5
        )
    )
    skipped, kept, examples = (
        small_model_example_factory.encode_inference_examples(
            [FAKE_VARIANT_CALL_HET], {}
        )
    )
    self.assertEqual(skipped, [])
    self.assertEqual(kept, [FAKE_VARIANT_CALL_HET])
    self.assertEqual(
        examples[0],
        [
            3,
            3,
            6,
            50,
            40,
            50,
            30,
            50,
            66,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            50,
            87,
            30,
            40,
        ],
    )

  def test_encode_inference_examples_expand_by_haplotype(self):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=5,
            expand_by_haplotype=True,
        )
    )
    skipped, kept, examples = (
        small_model_example_factory.encode_inference_examples(
            [FAKE_VARIANT_CALL_HET],
            READ_PHASES,
        )
    )
    self.assertEqual(skipped, [])
    self.assertEqual(kept, [FAKE_VARIANT_CALL_HET])
    self.assertEqual(
        examples[0],
        [
            3,
            3,
            6,
            50,
            40,
            50,
            30,
            50,
            66,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            50,
            87,
            30,
            40,
            1,
            0,
            1,
            0,
            40,
            0,
            25,
            0,
            100,
            0,
            1,
            2,
            3,
            66,
            60,
            45,
            30,
            55,
            0,
            0,
            1,
            1,
            2,
            50,
            20,
            60,
            35,
            40,
            100,
            0,
        ],
    )


if __name__ == "__main__":
  absltest.main()
