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
import tensorflow as tf

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.small_model import make_small_model_examples
from third_party.nucleus.protos import variants_pb2

MAIN_SAMPLE = "main_sample"
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
                sample_name=MAIN_SAMPLE,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_2",
                mapping_quality=20,
                average_base_quality=35,
                is_reverse_strand=True,
                sample_name=MAIN_SAMPLE,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_3",
                mapping_quality=40,
                average_base_quality=25,
                is_reverse_strand=True,
                sample_name=MAIN_SAMPLE,
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
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_5",
                    mapping_quality=30,
                    average_base_quality=60,
                    is_reverse_strand=False,
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_6",
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

FAKE_VARIANT_INSERTION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_bases="A",
        alternate_bases=["AAC"],
    ),
    allele_support_ext={
        "AAC": deepvariant_pb2.DeepVariantCall.SupportingReadsExt()
    },
)
FAKE_VARIANT_DELETION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_bases="AACC",
        alternate_bases=["C"],
    ),
    allele_support_ext={
        "C": deepvariant_pb2.DeepVariantCall.SupportingReadsExt()
    },
)
FAKE_VARIANT_MULTIALLELIC_INSERTION = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_bases="A",
        alternate_bases=["AC", "ACC", "ACCC"],
    ),
    ref_support_ext=deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
        read_infos=[
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_1",
                sample_name=MAIN_SAMPLE,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_2",
                sample_name=MAIN_SAMPLE,
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_3",
                sample_name=MAIN_SAMPLE,
            ),
        ]
    ),
    allele_support_ext={
        "AC": deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
            read_infos=[
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_4",
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_5",
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_6",
                    sample_name=MAIN_SAMPLE,
                ),
            ]
        ),
        "ACC": deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
            read_infos=[
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_7",
                    sample_name=MAIN_SAMPLE,
                ),
            ]
        ),
        "ACCC": deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
            read_infos=[
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_8",
                    sample_name=MAIN_SAMPLE,
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_9",
                    sample_name=MAIN_SAMPLE,
                ),
            ]
        ),
    },
)
FAKE_VARIANT_MULTIALLELIC_INSERTION_LABEL_HETALT = variant_labeler.VariantLabel(
    is_confident=True,
    variant=FAKE_VARIANT_MULTIALLELIC_INSERTION.variant,
    genotype=(1, 2),
)
READ_PHASES = {
    "read_1": 1,
    "read_2": 2,
    "read_3": 0,
    "read_4": 1,
    "read_5": 1,
    "read_6": 2,
}
FAKE_VARIANT_MULTI_SAMPLE_SNP = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(
        reference_bases="A",
        alternate_bases=["T"],
    ),
    ref_support_ext=deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
        read_infos=[
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_1",
                sample_name="sample_1",
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_2",
                sample_name="sample_1",
            ),
            deepvariant_pb2.DeepVariantCall.ReadSupport(
                read_name="read_3",
                sample_name="sample_2",
            ),
        ]
    ),
    allele_support_ext={
        "T": deepvariant_pb2.DeepVariantCall.SupportingReadsExt(
            read_infos=[
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_4",
                    sample_name="sample_1",
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_5",
                    sample_name="sample_2",
                ),
                deepvariant_pb2.DeepVariantCall.ReadSupport(
                    read_name="read_6",
                    sample_name="sample_3",
                ),
            ]
        ),
    },
)
FAKE_VARIANT_MULTI_SAMPLE_SNP_LABEL_HET = variant_labeler.VariantLabel(
    is_confident=True,
    variant=FAKE_VARIANT_MULTI_SAMPLE_SNP.variant,
    genotype=(0, 1),
)
BASE_FEATURE_VALUES_TOGETHER = [3, 3, 6, 6, 50, 50, 0, 0, 0, 0, 0, 0]
BASE_FEATURE_VALUES_SAMPLE_1 = [2, 1, 3, 6, 16, 33, 0, 0, 0, 0, 0, 0]
BASE_FEATURE_VALUES_SAMPLE_2 = [1, 1, 2, 6, 16, 50, 0, 0, 0, 0, 0, 0]


class SmallModelMakeExamplesTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name="contig_feature",
          feature=make_small_model_examples.IdentifyingFeature.CONTIG,
          expected_value="chr9",
      ),
      dict(
          testcase_name="start_feature",
          feature=make_small_model_examples.IdentifyingFeature.START,
          expected_value="5000",
      ),
      dict(
          testcase_name="end_feature",
          feature=make_small_model_examples.IdentifyingFeature.END,
          expected_value="5001",
      ),
      dict(
          testcase_name="ref_feature",
          feature=make_small_model_examples.IdentifyingFeature.REF,
          expected_value="A",
      ),
      dict(
          testcase_name="alt_feature",
          feature=make_small_model_examples.IdentifyingFeature.ALT,
          expected_value="C",
      ),
  )
  def test_feature_encoder_encode_identifying_feature(
      self, feature, expected_value
  ):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            FAKE_VARIANT_CALL_HET,
            make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
        ).encode_identifying_feature(
            feature,
        ),
        expected_value,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name="num_reads_supports_ref_feature",
          feature=make_small_model_examples.BaseFeature.NUM_READS_SUPPORTS_REF,
          expected_value=3,
      ),
      dict(
          testcase_name="num_reads_supports_alt_feature",
          feature=make_small_model_examples.BaseFeature.NUM_READS_SUPPORTS_ALT,
          expected_value=3,
      ),
      dict(
          testcase_name="alt_indices_depth_feature",
          feature=make_small_model_examples.BaseFeature.ALT_INDICES_DEPTH,
          expected_value=6,
      ),
      dict(
          testcase_name="alt_indices_depth_feature_multiallelic",
          feature=make_small_model_examples.BaseFeature.ALT_INDICES_DEPTH,
          expected_value=7,
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          alt_allele_indices=(0, 1),
      ),
      dict(
          testcase_name="total_depth_feature",
          feature=make_small_model_examples.BaseFeature.TOTAL_DEPTH,
          expected_value=6,
      ),
      dict(
          testcase_name="total_depth_feature_multiallelic",
          feature=make_small_model_examples.BaseFeature.TOTAL_DEPTH,
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          expected_value=9,
      ),
      dict(
          testcase_name="variant_allele_frequency_feature",
          feature=make_small_model_examples.BaseFeature.VARIANT_ALLELE_FREQUENCY,
          expected_value=50,
      ),
      dict(
          testcase_name="alt_indices_variant_allele_frequency_feature",
          feature=make_small_model_examples.BaseFeature.ALT_INDICES_VARIANT_ALLELE_FREQUENCY,
          expected_value=50,
      ),
      dict(
          testcase_name=(
              "alt_indices_variant_allele_frequency_feature_multiallelic"
          ),
          feature=make_small_model_examples.BaseFeature.ALT_INDICES_VARIANT_ALLELE_FREQUENCY,
          expected_value=57,
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          alt_allele_indices=(0, 1),
      ),
      dict(
          testcase_name="ref_mapping_quality_feature",
          feature=make_small_model_examples.BaseFeature.REF_MAPPING_QUALITY,
          expected_value=40,
      ),
      dict(
          testcase_name="alt_mapping_quality_feature",
          feature=make_small_model_examples.BaseFeature.ALT_MAPPING_QUALITY,
          expected_value=50,
      ),
      dict(
          testcase_name="ref_base_quality_feature",
          feature=make_small_model_examples.BaseFeature.REF_BASE_QUALITY,
          expected_value=30,
      ),
      dict(
          testcase_name="alt_base_quality_feature",
          feature=make_small_model_examples.BaseFeature.ALT_BASE_QUALITY,
          expected_value=50,
      ),
      dict(
          testcase_name="ref_reverse_strand_ratio_feature",
          feature=make_small_model_examples.BaseFeature.REF_REVERSE_STRAND_RATIO,
          expected_value=66,
      ),
      dict(
          testcase_name="alt_reverse_strand_ratio_feature",
          feature=make_small_model_examples.BaseFeature.ALT_REVERSE_STRAND_RATIO,
          expected_value=0,
      ),
  )
  def test_feature_encoder_encode_base_feature(
      self,
      feature,
      expected_value,
      candidate=FAKE_VARIANT_CALL_HET,
      alt_allele_indices=make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
  ):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            candidate,
            alt_allele_indices,
        ).encode_base_feature(
            feature,
        ),
        expected_value,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name="het_label",
          candidate=FAKE_VARIANT_CALL_HET,
          alt_allele_indices=(0,),
          label=FAKE_VARIANT_CALL_HET_LABEL,
          expected_value=[0, 1, 0],
      ),
      dict(
          testcase_name="hetalt_label_with_alt_allele_indices_0",
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          alt_allele_indices=(0,),
          label=FAKE_VARIANT_MULTIALLELIC_INSERTION_LABEL_HETALT,
          expected_value=[0, 1, 0],
      ),
      dict(
          testcase_name="hetalt_label_with_alt_allele_indices_1",
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          alt_allele_indices=(1,),
          label=FAKE_VARIANT_MULTIALLELIC_INSERTION_LABEL_HETALT,
          expected_value=[0, 1, 0],
      ),
      dict(
          testcase_name="hetalt_label_with_alt_allele_indices_2",
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          alt_allele_indices=(2,),
          label=FAKE_VARIANT_MULTIALLELIC_INSERTION_LABEL_HETALT,
          expected_value=[1, 0, 0],
      ),
      dict(
          testcase_name="hetalt_label_with_alt_allele_indices_1,2",
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          alt_allele_indices=(1, 2),
          label=FAKE_VARIANT_MULTIALLELIC_INSERTION_LABEL_HETALT,
          expected_value=[0, 1, 0],
      ),
  )
  def test_feature_encoder_one_hot_encode_label(
      self, candidate, alt_allele_indices, label, expected_value
  ):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            candidate,
            alt_allele_indices,
        ).encode_label(label),
        expected_value,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name="is_snp_feature",
          feature=make_small_model_examples.VariantFeature.IS_SNP,
          candidate=FAKE_VARIANT_CALL_HET,
          expected_value=1,
      ),
      dict(
          testcase_name="is_insertion",
          feature=make_small_model_examples.VariantFeature.IS_INSERTION,
          candidate=FAKE_VARIANT_CALL_HET,
          expected_value=0,
      ),
      dict(
          testcase_name="is_deletion_feature",
          feature=make_small_model_examples.VariantFeature.IS_DELETION,
          candidate=FAKE_VARIANT_CALL_HET,
          expected_value=0,
      ),
      dict(
          testcase_name="insertion_length_feature",
          feature=make_small_model_examples.VariantFeature.INSERTION_LENGTH,
          candidate=FAKE_VARIANT_CALL_HET,
          expected_value=0,
      ),
      dict(
          testcase_name="deletion_length_feature",
          feature=make_small_model_examples.VariantFeature.DELETION_LENGTH,
          candidate=FAKE_VARIANT_CALL_HET,
          expected_value=0,
      ),
      dict(
          testcase_name="is_snp_feature_insertion",
          feature=make_small_model_examples.VariantFeature.IS_SNP,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_snp_feature_deletion",
          feature=make_small_model_examples.VariantFeature.IS_SNP,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_insertion_feature_insertion",
          feature=make_small_model_examples.VariantFeature.IS_INSERTION,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=1,
      ),
      dict(
          testcase_name="is_insertion_feature_deletion",
          feature=make_small_model_examples.VariantFeature.IS_INSERTION,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_deletion_feature_insertion",
          feature=make_small_model_examples.VariantFeature.IS_DELETION,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_deletion_feature_deletion",
          feature=make_small_model_examples.VariantFeature.IS_DELETION,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=1,
      ),
      dict(
          testcase_name="insertion_length_feature_insertion",
          feature=make_small_model_examples.VariantFeature.INSERTION_LENGTH,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=2,
      ),
      dict(
          testcase_name="insertion_length_feature_deletion",
          feature=make_small_model_examples.VariantFeature.INSERTION_LENGTH,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="deletion_length_feature_insertion",
          feature=make_small_model_examples.VariantFeature.DELETION_LENGTH,
          candidate=FAKE_VARIANT_INSERTION,
          expected_value=0,
      ),
      dict(
          testcase_name="deletion_length_feature_deletion",
          feature=make_small_model_examples.VariantFeature.DELETION_LENGTH,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=3,
      ),
      dict(
          testcase_name="is_multiallelic_feature_biallelic",
          feature=make_small_model_examples.VariantFeature.IS_MULTIALLELIC,
          candidate=FAKE_VARIANT_DELETION,
          expected_value=0,
      ),
      dict(
          testcase_name="is_multiallelic_feature_multiallelic",
          feature=make_small_model_examples.VariantFeature.IS_MULTIALLELIC,
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          expected_value=1,
      ),
      dict(
          testcase_name="is_multiple_alt_alleles_feature_alt_allele_indices_0",
          feature=make_small_model_examples.VariantFeature.IS_MULTIPLE_ALT_ALLELES,
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          expected_value=0,
          alt_allele_indices=(0,),
      ),
      dict(
          testcase_name=(
              "is_multiple_alt_alleles_feature_alt_allele_indices_1,2"
          ),
          feature=make_small_model_examples.VariantFeature.IS_MULTIPLE_ALT_ALLELES,
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          expected_value=1,
          alt_allele_indices=(1, 2),
      ),
  )
  def test_encode_variant_feature(
      self,
      feature,
      candidate,
      expected_value,
      alt_allele_indices=make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
  ):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            candidate,
            alt_allele_indices,
        ).encode_variant_feature(
            feature,
        ),
        expected_value,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name="vaf_context_window_size_3",
          vaf_context_window_size=3,
          expected_value={
              "variant_allele_frequency_at_minus_1": 50,
              "variant_allele_frequency_at_plus_0": 87,
              "variant_allele_frequency_at_plus_1": 30,
          },
      ),
      dict(
          testcase_name="vaf_context_window_size_5",
          vaf_context_window_size=5,
          expected_value={
              "variant_allele_frequency_at_minus_2": 0,
              "variant_allele_frequency_at_minus_1": 50,
              "variant_allele_frequency_at_plus_0": 87,
              "variant_allele_frequency_at_plus_1": 30,
              "variant_allele_frequency_at_plus_2": 40,
          },
      ),
      dict(
          testcase_name="vaf_context_window_size_7",
          vaf_context_window_size=7,
          expected_value={
              "variant_allele_frequency_at_minus_3": 0,
              "variant_allele_frequency_at_minus_2": 0,
              "variant_allele_frequency_at_minus_1": 50,
              "variant_allele_frequency_at_plus_0": 87,
              "variant_allele_frequency_at_plus_1": 30,
              "variant_allele_frequency_at_plus_2": 40,
              "variant_allele_frequency_at_plus_3": 0,
          },
      ),
  )
  def test_feature_encoder_encode_variant_allele_frequency_at_position(
      self, vaf_context_window_size, expected_value
  ):
    self.assertEqual(
        make_small_model_examples.FeatureEncoder(
            FAKE_VARIANT_CALL_HET,
            make_small_model_examples.DEFAULT_ALT_ALLELE_INDICES,
        ).encode_variant_allele_frequency_at_position(
            vaf_context_window_size,
        ),
        expected_value,
    )

  @parameterized.named_parameters(
      dict(
          testcase_name="biallelic_candidate",
          candidate=FAKE_VARIANT_CALL_HET,
          expected_value=[(0,)],
      ),
      dict(
          testcase_name="multiallelic_candidate",
          candidate=FAKE_VARIANT_MULTIALLELIC_INSERTION,
          expected_value=[(0,), (1,), (2,), (0, 1), (0, 2), (1, 2)],
      ),
  )
  def test_get_set_of_allele_indices(self, candidate, expected_value):
    self.assertEqual(
        make_small_model_examples.get_set_of_allele_indices(candidate),
        expected_value,
    )

  @parameterized.parameters(
      dict(
          window_size=0,
          expand_by_haplotype=False,
          sample_names=[MAIN_SAMPLE],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
          ],
      ),
      dict(
          window_size=3,
          expand_by_haplotype=False,
          sample_names=[MAIN_SAMPLE],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
          ],
      ),
      dict(
          window_size=11,
          expand_by_haplotype=False,
          sample_names=[MAIN_SAMPLE],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
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
          sample_names=[MAIN_SAMPLE],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
              "variant_allele_frequency_at_minus_1",
              "variant_allele_frequency_at_plus_0",
              "variant_allele_frequency_at_plus_1",
              "main_sample_num_reads_supports_ref_hp_0",
              "main_sample_num_reads_supports_alt_hp_0",
              "main_sample_alt_indices_depth_hp_0",
              "main_sample_total_depth_hp_0",
              "main_sample_variant_allele_frequency_hp_0",
              "main_sample_alt_indices_variant_allele_frequency_hp_0",
              "main_sample_ref_mapping_quality_hp_0",
              "main_sample_alt_mapping_quality_hp_0",
              "main_sample_ref_base_quality_hp_0",
              "main_sample_alt_base_quality_hp_0",
              "main_sample_ref_reverse_strand_ratio_hp_0",
              "main_sample_alt_reverse_strand_ratio_hp_0",
              "main_sample_num_reads_supports_ref_hp_1",
              "main_sample_num_reads_supports_alt_hp_1",
              "main_sample_alt_indices_depth_hp_1",
              "main_sample_total_depth_hp_1",
              "main_sample_variant_allele_frequency_hp_1",
              "main_sample_alt_indices_variant_allele_frequency_hp_1",
              "main_sample_ref_mapping_quality_hp_1",
              "main_sample_alt_mapping_quality_hp_1",
              "main_sample_ref_base_quality_hp_1",
              "main_sample_alt_base_quality_hp_1",
              "main_sample_ref_reverse_strand_ratio_hp_1",
              "main_sample_alt_reverse_strand_ratio_hp_1",
              "main_sample_num_reads_supports_ref_hp_2",
              "main_sample_num_reads_supports_alt_hp_2",
              "main_sample_alt_indices_depth_hp_2",
              "main_sample_total_depth_hp_2",
              "main_sample_variant_allele_frequency_hp_2",
              "main_sample_alt_indices_variant_allele_frequency_hp_2",
              "main_sample_ref_mapping_quality_hp_2",
              "main_sample_alt_mapping_quality_hp_2",
              "main_sample_ref_base_quality_hp_2",
              "main_sample_alt_base_quality_hp_2",
              "main_sample_ref_reverse_strand_ratio_hp_2",
              "main_sample_alt_reverse_strand_ratio_hp_2",
          ],
      ),
      dict(
          window_size=0,
          expand_by_haplotype=False,
          sample_names=["sample_1", "sample_2"],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "sample_1_num_reads_supports_ref",
              "sample_1_num_reads_supports_alt",
              "sample_1_alt_indices_depth",
              "sample_1_total_depth",
              "sample_1_variant_allele_frequency",
              "sample_1_alt_indices_variant_allele_frequency",
              "sample_1_ref_mapping_quality",
              "sample_1_alt_mapping_quality",
              "sample_1_ref_base_quality",
              "sample_1_alt_base_quality",
              "sample_1_ref_reverse_strand_ratio",
              "sample_1_alt_reverse_strand_ratio",
              "sample_2_num_reads_supports_ref",
              "sample_2_num_reads_supports_alt",
              "sample_2_alt_indices_depth",
              "sample_2_total_depth",
              "sample_2_variant_allele_frequency",
              "sample_2_alt_indices_variant_allele_frequency",
              "sample_2_ref_mapping_quality",
              "sample_2_alt_mapping_quality",
              "sample_2_ref_base_quality",
              "sample_2_alt_base_quality",
              "sample_2_ref_reverse_strand_ratio",
              "sample_2_alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
          ],
      ),
      dict(
          window_size=0,
          expand_by_haplotype=False,
          sample_names=["sample_2", "sample_1"],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "sample_2_num_reads_supports_ref",
              "sample_2_num_reads_supports_alt",
              "sample_2_alt_indices_depth",
              "sample_2_total_depth",
              "sample_2_variant_allele_frequency",
              "sample_2_alt_indices_variant_allele_frequency",
              "sample_2_ref_mapping_quality",
              "sample_2_alt_mapping_quality",
              "sample_2_ref_base_quality",
              "sample_2_alt_base_quality",
              "sample_2_ref_reverse_strand_ratio",
              "sample_2_alt_reverse_strand_ratio",
              "sample_1_num_reads_supports_ref",
              "sample_1_num_reads_supports_alt",
              "sample_1_alt_indices_depth",
              "sample_1_total_depth",
              "sample_1_variant_allele_frequency",
              "sample_1_alt_indices_variant_allele_frequency",
              "sample_1_ref_mapping_quality",
              "sample_1_alt_mapping_quality",
              "sample_1_ref_base_quality",
              "sample_1_alt_base_quality",
              "sample_1_ref_reverse_strand_ratio",
              "sample_1_alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
          ],
      ),
      dict(
          window_size=0,
          expand_by_haplotype=True,
          sample_names=["sample_1", "sample_2"],
          expected_columns=[
              "num_reads_supports_ref",
              "num_reads_supports_alt",
              "alt_indices_depth",
              "total_depth",
              "variant_allele_frequency",
              "alt_indices_variant_allele_frequency",
              "ref_mapping_quality",
              "alt_mapping_quality",
              "ref_base_quality",
              "alt_base_quality",
              "ref_reverse_strand_ratio",
              "alt_reverse_strand_ratio",
              "sample_1_num_reads_supports_ref",
              "sample_1_num_reads_supports_alt",
              "sample_1_alt_indices_depth",
              "sample_1_total_depth",
              "sample_1_variant_allele_frequency",
              "sample_1_alt_indices_variant_allele_frequency",
              "sample_1_ref_mapping_quality",
              "sample_1_alt_mapping_quality",
              "sample_1_ref_base_quality",
              "sample_1_alt_base_quality",
              "sample_1_ref_reverse_strand_ratio",
              "sample_1_alt_reverse_strand_ratio",
              "sample_2_num_reads_supports_ref",
              "sample_2_num_reads_supports_alt",
              "sample_2_alt_indices_depth",
              "sample_2_total_depth",
              "sample_2_variant_allele_frequency",
              "sample_2_alt_indices_variant_allele_frequency",
              "sample_2_ref_mapping_quality",
              "sample_2_alt_mapping_quality",
              "sample_2_ref_base_quality",
              "sample_2_alt_base_quality",
              "sample_2_ref_reverse_strand_ratio",
              "sample_2_alt_reverse_strand_ratio",
              "is_snp",
              "is_insertion",
              "is_deletion",
              "insertion_length",
              "deletion_length",
              "is_multiallelic",
              "is_multiple_alt_alleles",
              "sample_1_num_reads_supports_ref_hp_0",
              "sample_1_num_reads_supports_alt_hp_0",
              "sample_1_alt_indices_depth_hp_0",
              "sample_1_total_depth_hp_0",
              "sample_1_variant_allele_frequency_hp_0",
              "sample_1_alt_indices_variant_allele_frequency_hp_0",
              "sample_1_ref_mapping_quality_hp_0",
              "sample_1_alt_mapping_quality_hp_0",
              "sample_1_ref_base_quality_hp_0",
              "sample_1_alt_base_quality_hp_0",
              "sample_1_ref_reverse_strand_ratio_hp_0",
              "sample_1_alt_reverse_strand_ratio_hp_0",
              "sample_1_num_reads_supports_ref_hp_1",
              "sample_1_num_reads_supports_alt_hp_1",
              "sample_1_alt_indices_depth_hp_1",
              "sample_1_total_depth_hp_1",
              "sample_1_variant_allele_frequency_hp_1",
              "sample_1_alt_indices_variant_allele_frequency_hp_1",
              "sample_1_ref_mapping_quality_hp_1",
              "sample_1_alt_mapping_quality_hp_1",
              "sample_1_ref_base_quality_hp_1",
              "sample_1_alt_base_quality_hp_1",
              "sample_1_ref_reverse_strand_ratio_hp_1",
              "sample_1_alt_reverse_strand_ratio_hp_1",
              "sample_1_num_reads_supports_ref_hp_2",
              "sample_1_num_reads_supports_alt_hp_2",
              "sample_1_alt_indices_depth_hp_2",
              "sample_1_total_depth_hp_2",
              "sample_1_variant_allele_frequency_hp_2",
              "sample_1_alt_indices_variant_allele_frequency_hp_2",
              "sample_1_ref_mapping_quality_hp_2",
              "sample_1_alt_mapping_quality_hp_2",
              "sample_1_ref_base_quality_hp_2",
              "sample_1_alt_base_quality_hp_2",
              "sample_1_ref_reverse_strand_ratio_hp_2",
              "sample_1_alt_reverse_strand_ratio_hp_2",
              "sample_2_num_reads_supports_ref_hp_0",
              "sample_2_num_reads_supports_alt_hp_0",
              "sample_2_alt_indices_depth_hp_0",
              "sample_2_total_depth_hp_0",
              "sample_2_variant_allele_frequency_hp_0",
              "sample_2_alt_indices_variant_allele_frequency_hp_0",
              "sample_2_ref_mapping_quality_hp_0",
              "sample_2_alt_mapping_quality_hp_0",
              "sample_2_ref_base_quality_hp_0",
              "sample_2_alt_base_quality_hp_0",
              "sample_2_ref_reverse_strand_ratio_hp_0",
              "sample_2_alt_reverse_strand_ratio_hp_0",
              "sample_2_num_reads_supports_ref_hp_1",
              "sample_2_num_reads_supports_alt_hp_1",
              "sample_2_alt_indices_depth_hp_1",
              "sample_2_total_depth_hp_1",
              "sample_2_variant_allele_frequency_hp_1",
              "sample_2_alt_indices_variant_allele_frequency_hp_1",
              "sample_2_ref_mapping_quality_hp_1",
              "sample_2_alt_mapping_quality_hp_1",
              "sample_2_ref_base_quality_hp_1",
              "sample_2_alt_base_quality_hp_1",
              "sample_2_ref_reverse_strand_ratio_hp_1",
              "sample_2_alt_reverse_strand_ratio_hp_1",
              "sample_2_num_reads_supports_ref_hp_2",
              "sample_2_num_reads_supports_alt_hp_2",
              "sample_2_alt_indices_depth_hp_2",
              "sample_2_total_depth_hp_2",
              "sample_2_variant_allele_frequency_hp_2",
              "sample_2_alt_indices_variant_allele_frequency_hp_2",
              "sample_2_ref_mapping_quality_hp_2",
              "sample_2_alt_mapping_quality_hp_2",
              "sample_2_ref_base_quality_hp_2",
              "sample_2_alt_base_quality_hp_2",
              "sample_2_ref_reverse_strand_ratio_hp_2",
              "sample_2_alt_reverse_strand_ratio_hp_2",
          ],
      ),
  )
  def test_small_model_example_factory_model_features(
      self, window_size, expand_by_haplotype, sample_names, expected_columns
  ):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=window_size,
            expand_by_haplotype=expand_by_haplotype,
            sample_names=sample_names,
        )
    )
    self.assertEqual(
        small_model_example_factory.model_features,
        expected_columns,
    )

  def test_encode_training_examples(self):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=3,
            sample_names=[MAIN_SAMPLE],
        )
    )
    training_examples = small_model_example_factory.encode_training_examples(
        [
            (FAKE_VARIANT_CALL_HET, FAKE_VARIANT_CALL_HET_LABEL),
        ],
        {},
        sample_order=[0],
    )
    self.assertEqual(
        training_examples[0],
        tf.train.Example(
            features=tf.train.Features(
                feature={
                    make_small_model_examples.FEATURES_ENCODED: (
                        tf.train.Feature(
                            int64_list=tf.train.Int64List(
                                value=[
                                    3,
                                    3,
                                    6,
                                    6,
                                    50,
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
                                    0,
                                    50,
                                    87,
                                    30,
                                ]
                            )
                        )
                    ),
                    make_small_model_examples.IDS_ENCODED: tf.train.Feature(
                        bytes_list=tf.train.BytesList(
                            value=[
                                b"chr9",
                                b"5000",
                                b"5001",
                                b"A",
                                b"C",
                            ]
                        )
                    ),
                    make_small_model_examples.LABEL_ENCODED: tf.train.Feature(
                        int64_list=tf.train.Int64List(value=[0, 1, 0])
                    ),
                    make_small_model_examples.GENOTYPE_ENCODED: (
                        tf.train.Feature(
                            int64_list=tf.train.Int64List(value=[0, 1])
                        )
                    ),
                }
            )
        ),
    )

  def test_encode_inference_examples(self):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=5,
            sample_names=[MAIN_SAMPLE],
        )
    )
    example_set = small_model_example_factory.encode_inference_examples(
        [FAKE_VARIANT_CALL_HET, FAKE_VARIANT_MULTIALLELIC_INSERTION], {}, [0]
    )
    self.assertEqual(example_set.skipped_candidates, [])
    self.assertEqual(
        example_set.candidates_with_alt_allele_indices,
        [
            (FAKE_VARIANT_CALL_HET, (0,)),
            (FAKE_VARIANT_MULTIALLELIC_INSERTION, (0,)),
            (FAKE_VARIANT_MULTIALLELIC_INSERTION, (1,)),
            (FAKE_VARIANT_MULTIALLELIC_INSERTION, (2,)),
            (FAKE_VARIANT_MULTIALLELIC_INSERTION, (0, 1)),
            (FAKE_VARIANT_MULTIALLELIC_INSERTION, (0, 2)),
            (FAKE_VARIANT_MULTIALLELIC_INSERTION, (1, 2)),
        ],
    )
    self.assertEqual(
        example_set.inference_examples[0],
        [
            3,
            3,
            6,
            6,
            50,
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
            sample_names=[MAIN_SAMPLE],
            expand_by_haplotype=True,
        )
    )
    example_set = small_model_example_factory.encode_inference_examples(
        [FAKE_VARIANT_CALL_HET], READ_PHASES, [0]
    )
    self.assertEqual(example_set.skipped_candidates, [])
    self.assertEqual(
        example_set.candidates_with_alt_allele_indices,
        [
            (FAKE_VARIANT_CALL_HET, (0,)),
        ],
    )
    self.assertEqual(
        example_set.inference_examples[0],
        [
            3,
            3,
            6,
            6,
            50,
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
            0,
            0,
            50,
            87,
            30,
            40,
            1,
            0,
            1,
            6,
            0,
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
            6,
            33,
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
            6,
            16,
            50,
            20,
            60,
            35,
            40,
            100,
            0,
        ],
    )

  @parameterized.parameters(
      dict(
          sample_order=[0, 1],
          expected_values=[
              BASE_FEATURE_VALUES_TOGETHER,
              BASE_FEATURE_VALUES_SAMPLE_1,
              BASE_FEATURE_VALUES_SAMPLE_2,
          ],
      ),
      dict(
          sample_order=[1, 0],
          expected_values=[
              BASE_FEATURE_VALUES_TOGETHER,
              BASE_FEATURE_VALUES_SAMPLE_2,
              BASE_FEATURE_VALUES_SAMPLE_1,
          ],
      ),
  )
  def test_encode_inference_examples_multi_sample(
      self, sample_order, expected_values
  ):
    small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            vaf_context_window_size=0,
            sample_names=["sample_1", "sample_2"],
        )
    )
    example_set = small_model_example_factory.encode_inference_examples(
        [FAKE_VARIANT_MULTI_SAMPLE_SNP], {}, sample_order
    )
    self.assertEqual(example_set.skipped_candidates, [])
    self.assertEqual(
        example_set.inference_examples[0],
        expected_values[0]
        + expected_values[1]
        + expected_values[2]
        + [1, 0, 0, 0, 0, 0, 0],
    )


if __name__ == "__main__":
  absltest.main()
