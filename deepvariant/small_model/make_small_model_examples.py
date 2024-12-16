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
"""Module for generating small model examples."""

import enum
from typing import Dict, List, Sequence, Tuple

import tensorflow as tf

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import variant_utils

# Example encoding keys
FEATURES_ENCODED = 'features/encoded'
IDS_ENCODED = 'ids/encoded'
LABEL_ENCODED = 'label/encoded'


class GenotypeEncoding(enum.Enum):
  REF = 0
  HET = 1
  HOM_ALT = 2


class Haplotype(enum.Enum):
  HP_0 = 0
  HP_1 = 1
  HP_2 = 2


class SmallModelFeature(enum.Enum):
  """Ordered list of features used in the small model."""


class IdentifyingFeature(SmallModelFeature):
  """Features used to identify each candidate."""

  CONTIG = 'contig'
  START = 'start'
  END = 'end'
  REF = 'ref'
  ALT_1 = 'alt_1'


class TruthFeature(SmallModelFeature):
  """A single truth label for the candidate."""

  GENOTYPE = 'genotype'


class BaseFeature(SmallModelFeature):
  """Set of basic features that can be computed over a subset of reads."""

  NUM_READS_SUPPORTS_REF = 'num_reads_supports_ref'
  NUM_READS_SUPPORTS_ALT_1 = 'num_reads_supports_alt_1'
  TOTAL_DEPTH = 'total_depth'
  VARIANT_ALLELE_FREQUENCY_1 = 'variant_allele_frequency_1'
  REF_MAPPING_QUALITY = 'ref_mapping_quality'
  ALT_1_MAPPING_QUALITY = 'alt_1_mapping_quality'
  REF_BASE_QUALITY = 'ref_base_quality'
  ALT_1_BASE_QUALITY = 'alt_1_base_quality'
  REF_REVERSE_STRAND_RATIO = 'ref_reverse_strand_ratio'
  ALT_1_REVERSE_STRAND_RATIO = 'alt_1_reverse_strand_ratio'


class VariantFeature(SmallModelFeature):
  """Set of features that describe the variant candidate."""

  IS_SNP = 'is_snp'
  IS_INSERTION = 'is_insertion'
  IS_DELETION = 'is_deletion'
  INSERTION_LENGTH = 'insertion_length'
  DELETION_LENGTH = 'deletion_length'


ENCODING_BY_GENOTYPE = {
    (0, 0): GenotypeEncoding.REF.value,
    (0, 1): GenotypeEncoding.HET.value,
    (1, 0): GenotypeEncoding.HET.value,
    (1, 1): GenotypeEncoding.HOM_ALT.value,
}
FAKE_CANDIDATE = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(alternate_bases=[''])
)


def _mean_for_attribute(
    read_infos: Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport],
    attribute_name: str,
    multiplier: int = 1,
) -> int:
  """Calculates the mean of an attribute for a list of reads."""
  if not read_infos:
    return 0
  return (
      multiplier
      * sum(getattr(r, attribute_name) for r in read_infos)
      // len(read_infos)
  )


def _filter_by_haplotype(
    read_infos: Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport],
    read_phases: Dict[str, int],
    haplotype: int | None = None,
) -> Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport]:
  """Filters the reads by haplotype."""
  return [r for r in read_infos if read_phases.get(r.read_name, 0) == haplotype]


def _get_context_allele_frequency_offsets(
    vaf_context_window_size: int,
) -> Sequence[int]:
  """Returns the list of offsets for the given context window size."""
  if not vaf_context_window_size:
    return []
  half_window_size = vaf_context_window_size // 2
  return list(range(-half_window_size, half_window_size + 1))


class FeatureEncoder:
  """Class for encoding small model features for a single candidate."""

  def __init__(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      haplotype: int | None = None,
      read_phases: Dict[str, int] | None = None,
  ):
    """Initializes the feature encoder.

    Args:
      candidate: the candidate proto from which to extract the feature.
      haplotype: (optional) the haplotype tag to use for filtering reads.
      read_phases: (optional)a dictionary of read names to haplotype phases.
    """
    self.candidate = candidate
    self.haplotype = haplotype
    self.ref_read_infos = self.candidate.ref_support_ext.read_infos
    if haplotype is not None and read_phases is not None:
      self.ref_read_infos = _filter_by_haplotype(
          self.ref_read_infos, read_phases, haplotype
      )
    self.ref_read_infos_count = len(self.ref_read_infos)

    alt_base = candidate.variant.alternate_bases[0]
    self.alt_1_read_infos = []
    if alt_base in candidate.allele_support_ext:
      self.alt_1_read_infos = candidate.allele_support_ext[alt_base].read_infos
    if haplotype is not None and read_phases is not None:
      self.alt_1_read_infos = _filter_by_haplotype(
          self.alt_1_read_infos, read_phases, haplotype
      )
    self.alt_1_read_infos_count = len(self.alt_1_read_infos)

  def _get_contig(self) -> str:
    """Returns the contig of the candidate."""
    return self.candidate.variant.reference_name

  def _get_start(self) -> str:
    """Returns the start position of the candidate."""
    return str(self.candidate.variant.start)

  def _get_end(self) -> str:
    """Returns the end position of the candidate."""
    return str(self.candidate.variant.end)

  def _get_ref(self) -> str:
    """Returns the reference base of the candidate."""
    return self.candidate.variant.reference_bases

  def _get_alt_1(self) -> str:
    """Returns the first alternate base of the candidate."""
    return self.candidate.variant.alternate_bases[0]

  def _get_num_reads_supports_ref(self) -> int:
    """Returns the number of reads supporting the reference base."""
    return self.ref_read_infos_count

  def _get_num_reads_supports_alt_1(self) -> int:
    """Returns the number of reads supporting the alternate base."""
    return self.alt_1_read_infos_count

  def _get_total_depth(self) -> int:
    """Returns the total depth of the candidate."""
    return self.ref_read_infos_count + self.alt_1_read_infos_count

  def _get_variant_allele_frequency_1(self) -> int:
    """Returns the variant allele frequency of the candidate."""
    dp = self._get_total_depth()
    if dp == 0:
      return 0
    return 100 * self.alt_1_read_infos_count // dp

  def _get_ref_mapping_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.ref_read_infos, 'mapping_quality')

  def _get_alt_1_mapping_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.alt_1_read_infos, 'mapping_quality')

  def _get_ref_base_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.ref_read_infos, 'average_base_quality')

  def _get_alt_1_base_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.alt_1_read_infos, 'average_base_quality')

  def _get_ref_reverse_strand_ratio(self) -> int:
    """Returns the reverse strand ratio of the candidate."""
    return _mean_for_attribute(
        self.ref_read_infos, 'is_reverse_strand', multiplier=100
    )

  def _get_alt_1_reverse_strand_ratio(self) -> int:
    """Returns the reverse strand ratio of the candidate."""
    return _mean_for_attribute(
        self.alt_1_read_infos, 'is_reverse_strand', multiplier=100
    )

  def _get_is_snp(self) -> int:
    """Returns 1 if the candidate is a SNP, 0 otherwise."""
    return int(variant_utils.is_snp(self.candidate.variant))

  def _get_is_insertion(self) -> int:
    """Returns 1 if the candidate is an insertion, 0 otherwise."""
    return int(
        variant_utils.is_insertion(
            self.candidate.variant.reference_bases,
            self.candidate.variant.alternate_bases[0],
        )
    )

  def _get_is_deletion(self) -> int:
    """Returns 1 if the candidate is an insertion, 0 otherwise."""
    return int(
        variant_utils.is_deletion(
            self.candidate.variant.reference_bases,
            self.candidate.variant.alternate_bases[0],
        )
    )

  def _get_insertion_length(self) -> int:
    """Returns the insertion length of the candidate."""
    insertion_length = len(self.candidate.variant.alternate_bases[0]) - len(
        self.candidate.variant.reference_bases
    )
    return max(0, insertion_length)

  def _get_deletion_length(self) -> int:
    """Returns the deletion length of the candidate."""
    deletion_length = len(self.candidate.variant.reference_bases) - len(
        self.candidate.variant.alternate_bases[0]
    )
    return max(0, deletion_length)

  def encode_base_feature(
      self,
      feature: BaseFeature,
  ) -> int:
    """Encodes the given base feature for the candidate.

    Args:
      feature: specifies which feature to extract

    Returns:
      The extracted value for that base feature for the candidate.
    """
    if feature == BaseFeature.NUM_READS_SUPPORTS_REF:
      return self._get_num_reads_supports_ref()
    elif feature == BaseFeature.NUM_READS_SUPPORTS_ALT_1:
      return self._get_num_reads_supports_alt_1()
    elif feature == BaseFeature.TOTAL_DEPTH:
      return self._get_total_depth()
    elif feature == BaseFeature.VARIANT_ALLELE_FREQUENCY_1:
      return self._get_variant_allele_frequency_1()
    elif feature == BaseFeature.REF_MAPPING_QUALITY:
      return self._get_ref_mapping_quality()
    elif feature == BaseFeature.ALT_1_MAPPING_QUALITY:
      return self._get_alt_1_mapping_quality()
    elif feature == BaseFeature.REF_BASE_QUALITY:
      return self._get_ref_base_quality()
    elif feature == BaseFeature.ALT_1_BASE_QUALITY:
      return self._get_alt_1_base_quality()
    elif feature == BaseFeature.REF_REVERSE_STRAND_RATIO:
      return self._get_ref_reverse_strand_ratio()
    elif feature == BaseFeature.ALT_1_REVERSE_STRAND_RATIO:
      return self._get_alt_1_reverse_strand_ratio()
    else:
      raise ValueError(f'{feature.value} does not map to a callable.')

  def encode_variant_feature(
      self,
      feature: VariantFeature,
  ) -> int:
    """Encodes the given variant feature for the candidate.

    Args:
      feature: specifies which feature to extract

    Returns:
      The extracted value for that variant feature for the candidate.
    """
    if feature == VariantFeature.IS_SNP:
      return self._get_is_snp()
    elif feature == VariantFeature.IS_INSERTION:
      return self._get_is_insertion()
    elif feature == VariantFeature.IS_DELETION:
      return self._get_is_deletion()
    elif feature == VariantFeature.INSERTION_LENGTH:
      return self._get_insertion_length()
    elif feature == VariantFeature.DELETION_LENGTH:
      return self._get_deletion_length()
    else:
      raise ValueError(f'{feature.value} does not map to a callable.')

  def encode_identifying_feature(
      self,
      feature: IdentifyingFeature,
  ) -> str:
    """Encodes the given identifying feature for the candidate.

    Args:
      feature: specifies which feature to extract

    Returns:
      The extracted value for that identifying feature for the candidate.
    """
    if feature == IdentifyingFeature.CONTIG:
      return self._get_contig()
    elif feature == IdentifyingFeature.START:
      return self._get_start()
    elif feature == IdentifyingFeature.END:
      return self._get_end()
    elif feature == IdentifyingFeature.REF:
      return self._get_ref()
    elif feature == IdentifyingFeature.ALT_1:
      return self._get_alt_1()
    else:
      raise ValueError(f'{feature.value} does not map to a callable.')

  def encode_variant_allele_frequency_at_position(
      self,
      vaf_context_window_size: int,
  ) -> dict[str, int]:
    """Returns the VAF for the candidate at the given position.

    Args:
      vaf_context_window_size: the context window size to use for the VAF.

    Returns:
      A dictionary mapping context VAF feature names to values.
    """
    if not vaf_context_window_size:
      return {}
    allele_frequencies_at_positions = {}
    for offset in _get_context_allele_frequency_offsets(
        vaf_context_window_size
    ):
      position = self.candidate.variant.start + offset
      direction = 'minus' if offset < 0 else 'plus'
      feature_name = f'variant_allele_frequency_at_{direction}_{abs(offset)}'
      allele_frequencies_at_positions[feature_name] = (
          self.candidate.allele_frequency_at_position.get(position, 0)
      )

    return allele_frequencies_at_positions

  def encode_label(
      self,
      label: variant_labeler.VariantLabel,
  ) -> Sequence[int]:
    """Returns a one-hot encoded genotype.

    Args:
      label: the label to encode into a genotype.

    Returns:
      A one-hot encoded genotype.
    """
    value = [0 for _ in GenotypeEncoding]
    value[ENCODING_BY_GENOTYPE[label.genotype]] = 1
    return value


class SmallModelExampleFactory:
  """Class for making small model examples."""

  def __init__(
      self,
      vaf_context_window_size: int,
      accept_snps: bool = True,
      accept_indels: bool = True,
      expand_by_haplotype: bool = False,
      model_features: Sequence[str] | None = None,
  ):
    self.accept_snps = accept_snps
    self.accept_indels = accept_indels
    self.expand_by_haplotype = expand_by_haplotype
    self.vaf_context_window_size = vaf_context_window_size
    self.model_features = self._get_and_validate_model_features(model_features)

  def _get_and_validate_model_features(
      self, model_features: Sequence[str] | None
  ) -> Sequence[str]:
    """Returns the model features for the small model."""
    all_features = list(
        self._encode_candidate_feature_dict(FAKE_CANDIDATE).keys()
    )
    if model_features:
      if not set(model_features).issubset(all_features):
        unrecognized_features = set(model_features).difference(all_features)
        raise ValueError(
            'The specified small model features are invalid:'
            f' {",".join(unrecognized_features)}'
        )
      return model_features
    else:
      return all_features

  def _pass_candidate_to_small_model(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
  ) -> bool:
    """Determines if the candidate is eligible for the small model."""
    if not variant_utils.is_biallelic(candidate.variant):
      return False
    elif variant_utils.is_snp(candidate.variant):
      return self.accept_snps
    else:
      return self.accept_indels

  def _encode_candidate_feature_dict(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      read_phases: dict[str, int] | None = None,
  ) -> dict[str, int]:
    """Encodes all model features for a given candidate into a key-value dict."""
    feature_encoder = FeatureEncoder(candidate)
    candidate_example = {}
    candidate_example.update({
        feature.value: feature_encoder.encode_base_feature(feature)
        for feature in BaseFeature
    })
    candidate_example.update({
        feature.value: feature_encoder.encode_variant_feature(feature)
        for feature in VariantFeature
    })
    candidate_example.update(
        feature_encoder.encode_variant_allele_frequency_at_position(
            self.vaf_context_window_size
        )
    )
    if self.expand_by_haplotype:
      for haplotype in Haplotype:
        haplotype_feature_encoder = FeatureEncoder(
            candidate, haplotype.value, read_phases
        )
        for feature in BaseFeature:
          candidate_example.update({
              f'{feature.value}_hp_{haplotype.value}': (
                  haplotype_feature_encoder.encode_base_feature(feature)
              )
          })
    return candidate_example

  def _encode_model_features(self, candidate, read_phases) -> Sequence[int]:
    """Encodes a candidate example into the selected model features.

    Args:
      candidate: The candidate to be encoded.
      read_phases: A dictionary mapping read names to haplotype tags.

    Returns:
      A feature vector.
    """
    encoded_candidate = self._encode_candidate_feature_dict(
        candidate, read_phases
    )
    return [encoded_candidate[feature] for feature in self.model_features]

  def _encode_candidate_example(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      label: variant_labeler.VariantLabel | None,
      read_phases: Dict[str, int],
  ) -> tf.train.Example:
    """Encodes a candidate example."""
    feature_encoder = FeatureEncoder(candidate)
    return tf.train.Example(
        features=tf.train.Features(
            feature={
                FEATURES_ENCODED: tf.train.Feature(
                    int64_list=tf.train.Int64List(
                        value=self._encode_model_features(
                            candidate, read_phases
                        )
                    )
                ),
                IDS_ENCODED: tf.train.Feature(
                    bytes_list=tf.train.BytesList(
                        value=[
                            feature_encoder.encode_identifying_feature(
                                feature
                            ).encode()
                            for feature in IdentifyingFeature
                        ]
                    )
                ),
                LABEL_ENCODED: tf.train.Feature(
                    int64_list=tf.train.Int64List(
                        value=feature_encoder.encode_label(label)
                    )
                ),
            }
        )
    )

  def encode_training_examples(
      self,
      candidates_with_label: Sequence[
          Tuple[deepvariant_pb2.DeepVariantCall, variant_labeler.VariantLabel]
      ],
      read_phases: Dict[str, int],
  ) -> Sequence[tf.train.Example]:
    """Generates examples from the given candidates for training.

    Args:
      candidates_with_label: List of candidates with labels to be processed into
        examples.
      read_phases: A dictionary mapping read names to haplotype tags.

    Returns:
      A list of encoded candidate examples.
    """
    candidate_examples = []
    for candidate, label in candidates_with_label:
      if not self._pass_candidate_to_small_model(candidate):
        continue
      candidate_example = self._encode_candidate_example(
          candidate, label, read_phases
      )
      candidate_examples.append(candidate_example)
    return candidate_examples

  def encode_inference_examples(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      read_phases: Dict[str, int],
  ) -> Tuple[
      List[deepvariant_pb2.DeepVariantCall],
      List[deepvariant_pb2.DeepVariantCall],
      List[Sequence[int]],
  ]:
    """Generates examples from the given candidates for inference.

    Args:
      candidates: List of candidates to be processed into a summary.
      read_phases: A dictionary mapping read names to haplotype tags.

    Returns:
      A tuple containing:
        A list of all candidates that were skipped.
        A list of all candidates for which examples were generated.
        A list of encoded candidate examples.
    """
    skipped_candidates = []
    kept_candidates = []
    examples = []
    for candidate in candidates:
      if not self._pass_candidate_to_small_model(candidate):
        skipped_candidates.append(candidate)
        continue
      kept_candidates.append(candidate)
      candidate_example = self._encode_model_features(candidate, read_phases)
      examples.append(candidate_example)
    return skipped_candidates, kept_candidates, examples
