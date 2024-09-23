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
import re
from typing import Dict, List, Sequence, Tuple, Union

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.util import variant_utils


class GenotypeEncoding(enum.Enum):
  REF = 0
  HET = 1
  HOM_ALT = 2


class SmallModelFeature(enum.Enum):
  """Ordered list of features used in the small model."""

  # Features used to label and track examples
  CONTIG = 'contig'
  START = 'start'
  END = 'end'
  REF = 'ref'
  ALT_1 = 'alt_1'
  # Truth feature
  GENOTYPE = 'genotype'
  # Features passed to the model
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
  IS_SNP = 'is_snp'
  IS_INSERTION = 'is_insertion'
  IS_DELETION = 'is_deletion'
  INSERTION_LENGTH = 'insertion_length'
  DELETION_LENGTH = 'deletion_length'


IDENTIFYING_FEATURES = (
    SmallModelFeature.CONTIG,
    SmallModelFeature.START,
    SmallModelFeature.END,
    SmallModelFeature.REF,
    SmallModelFeature.ALT_1,
)
TRUTH_FEATURE = SmallModelFeature.GENOTYPE
VARIANT_ALLELE_FREQUENCY_AT_PREFIX = 'variant_allele_frequency_at'

# This set of features gets calculated per haplotype tag.
HAPLOTYPE_FEATURE_SET = [
    SmallModelFeature.NUM_READS_SUPPORTS_REF,
    SmallModelFeature.NUM_READS_SUPPORTS_ALT_1,
    SmallModelFeature.TOTAL_DEPTH,
    SmallModelFeature.VARIANT_ALLELE_FREQUENCY_1,
    SmallModelFeature.REF_MAPPING_QUALITY,
    SmallModelFeature.ALT_1_MAPPING_QUALITY,
    SmallModelFeature.REF_BASE_QUALITY,
    SmallModelFeature.ALT_1_BASE_QUALITY,
    SmallModelFeature.REF_REVERSE_STRAND_RATIO,
    SmallModelFeature.ALT_1_REVERSE_STRAND_RATIO,
]


ENCODING_BY_GENOTYPE = {
    (0, 0): GenotypeEncoding.REF.value,
    (0, 1): GenotypeEncoding.HET.value,
    (1, 0): GenotypeEncoding.HET.value,
    (1, 1): GenotypeEncoding.HOM_ALT.value,
}


def filter_by_haplotype(
    read_infos: Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport],
    read_phases: Dict[str, int],
    haplotype: int | None = None,
) -> Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport]:
  """Filters read infos by haplotype."""
  return [r for r in read_infos if read_phases.get(r.read_name, 0) == haplotype]


class FeatureEncoder:
  """Class for encoding small model features for a single candidate."""

  def __init__(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      label: Union[variant_labeler.VariantLabel, None] = None,
      haplotype: int | None = None,
      read_phases: Dict[str, int] | None = None,
  ):
    """Initializes the feature encoder.

    Args:
      candidate: the candidate proto from which to extract the feature.
      label: the label proto that contains the truth information.
      haplotype: (optional) the haplotype tag to use for filtering reads.
      read_phases: (optional)a dictionary of read names to haplotype phases.
    """
    self.candidate = candidate
    self.label = label
    self.haplotype = haplotype
    self.ref_read_infos = self.candidate.ref_support_ext.read_infos
    if haplotype is not None and read_phases is not None:
      self.ref_read_infos = filter_by_haplotype(
          self.ref_read_infos, read_phases, haplotype
      )
    self.ref_read_infos_count = len(self.ref_read_infos)

    alt_base = candidate.variant.alternate_bases[0]
    self.alt_1_read_infos = []
    if alt_base in candidate.allele_support_ext:
      self.alt_1_read_infos = candidate.allele_support_ext[alt_base].read_infos
    if haplotype is not None and read_phases is not None:
      self.alt_1_read_infos = filter_by_haplotype(
          self.alt_1_read_infos, read_phases, haplotype
      )
    self.alt_1_read_infos_count = len(self.alt_1_read_infos)

  def _get_contig(self) -> str:
    """Returns the contig of the candidate."""
    return self.candidate.variant.reference_name

  def _get_start(self) -> int:
    """Returns the start position of the candidate."""
    return self.candidate.variant.start

  def _get_end(self) -> int:
    """Returns the end position of the candidate."""
    return self.candidate.variant.end

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
    if not self.ref_read_infos:
      return 0
    return (
        sum(r.mapping_quality for r in self.ref_read_infos)
        // self.ref_read_infos_count
    )

  def _get_alt_1_mapping_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    if not self.alt_1_read_infos:
      return 0
    return (
        sum(r.mapping_quality for r in self.alt_1_read_infos)
        // self.alt_1_read_infos_count
    )

  def _get_ref_base_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    if not self.ref_read_infos:
      return 0
    return (
        sum(r.average_base_quality for r in self.ref_read_infos)
        // self.ref_read_infos_count
    )

  def _get_alt_1_base_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    if not self.alt_1_read_infos:
      return 0
    return (
        sum(r.average_base_quality for r in self.alt_1_read_infos)
        // self.alt_1_read_infos_count
    )

  def _get_ref_reverse_strand_ratio(self) -> int:
    """Returns the reverse strand ratio of the candidate."""
    if not self.ref_read_infos:
      return 0
    return (
        100
        * sum(r.is_reverse_strand for r in self.ref_read_infos)
        // self.ref_read_infos_count
    )

  def _get_alt_1_reverse_strand_ratio(self) -> int:
    """Returns the reverse strand ratio of the candidate."""
    if not self.alt_1_read_infos:
      return 0
    return (
        100
        * sum(r.is_reverse_strand for r in self.alt_1_read_infos)
        // self.alt_1_read_infos_count
    )

  def _get_variant_allele_frequency_at_position(
      self,
      feature_name: str,
  ) -> int:
    """Returns the VAF for the candidate at the given position."""
    offset = int(re.findall('[0-9]+', feature_name)[0])
    plus_or_minus = -1 if 'minus' in feature_name else 1
    position = self.candidate.variant.start + (plus_or_minus * offset)
    return self.candidate.allele_frequency_at_position.get(position, 0)

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

  def encode_feature(
      self,
      feature: str,
  ) -> Union[int, str]:
    """Encodes the given feature from a candidate or label.

    Args:
      feature: specifies which feature to extract

    Returns:
      The extracted value for that feature for the candidate.
    """
    if feature == SmallModelFeature.CONTIG.value:
      return self._get_contig()
    elif feature == SmallModelFeature.START.value:
      return self._get_start()
    elif feature == SmallModelFeature.END.value:
      return self._get_end()
    elif feature == SmallModelFeature.REF.value:
      return self._get_ref()
    elif feature == SmallModelFeature.ALT_1.value:
      return self._get_alt_1()
    elif feature == SmallModelFeature.NUM_READS_SUPPORTS_REF.value:
      return self._get_num_reads_supports_ref()
    elif feature == SmallModelFeature.NUM_READS_SUPPORTS_ALT_1.value:
      return self._get_num_reads_supports_alt_1()
    elif feature == SmallModelFeature.TOTAL_DEPTH.value:
      return self._get_total_depth()
    elif feature == SmallModelFeature.VARIANT_ALLELE_FREQUENCY_1.value:
      return self._get_variant_allele_frequency_1()
    elif feature == SmallModelFeature.REF_MAPPING_QUALITY.value:
      return self._get_ref_mapping_quality()
    elif feature == SmallModelFeature.ALT_1_MAPPING_QUALITY.value:
      return self._get_alt_1_mapping_quality()
    elif feature == SmallModelFeature.REF_BASE_QUALITY.value:
      return self._get_ref_base_quality()
    elif feature == SmallModelFeature.ALT_1_BASE_QUALITY.value:
      return self._get_alt_1_base_quality()
    elif feature == SmallModelFeature.REF_REVERSE_STRAND_RATIO.value:
      return self._get_ref_reverse_strand_ratio()
    elif feature == SmallModelFeature.ALT_1_REVERSE_STRAND_RATIO.value:
      return self._get_alt_1_reverse_strand_ratio()
    elif feature.startswith(VARIANT_ALLELE_FREQUENCY_AT_PREFIX):
      return self._get_variant_allele_frequency_at_position(feature)
    elif feature == SmallModelFeature.IS_SNP.value:
      return self._get_is_snp()
    elif feature == SmallModelFeature.IS_INSERTION.value:
      return self._get_is_insertion()
    elif feature == SmallModelFeature.IS_DELETION.value:
      return self._get_is_deletion()
    elif feature == SmallModelFeature.INSERTION_LENGTH.value:
      return self._get_insertion_length()
    elif feature == SmallModelFeature.DELETION_LENGTH.value:
      return self._get_deletion_length()
    elif feature == SmallModelFeature.GENOTYPE.value and self.label:
      return ENCODING_BY_GENOTYPE[self.label.genotype]
    else:
      raise ValueError(f'{feature} does not map to a callable.')


class SmallModelExampleFactory:
  """Class for making small model examples."""
  HP_TAGS = (0, 1, 2)

  def __init__(
      self,
      vaf_context_window_size: int,
      accept_snps: bool = True,
      accept_indels: bool = True,
      expand_by_haplotype: bool = False,
  ):
    self.accept_snps = accept_snps
    self.accept_indels = accept_indels
    self.expand_by_haplotype = expand_by_haplotype
    context_vaf_features = self._get_context_vaf_features(
        vaf_context_window_size
    )
    self._haplotype_features = []
    if expand_by_haplotype:
      self._haplotype_features = self._get_haplotype_feature_set()

    self._training_features = [f.value for f in SmallModelFeature]
    self._inference_features = [
        feature.value
        for feature in SmallModelFeature
        if feature not in IDENTIFYING_FEATURES and feature != TRUTH_FEATURE
    ]
    self._inference_features.extend(context_vaf_features)
    self._training_features.extend(context_vaf_features)

  @property
  def inference_features(self) -> list[str]:
    return self._inference_features + self._haplotype_features

  @property
  def training_features(self) -> list[str]:
    return self._training_features + self._haplotype_features

  def _get_haplotype_feature_set(self):
    haplotype_features = []
    for hp_tag in self.HP_TAGS:
      haplotype_features.extend(
          [f'{feature.value}_hp_{hp_tag}' for feature in HAPLOTYPE_FEATURE_SET]
      )
    return haplotype_features

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

  def _get_context_vaf_features(
      self,
      vaf_context_window_size: int,
  ) -> List[str]:
    """Adds context VAF to model_features based on the window size."""
    if not vaf_context_window_size:
      return []
    half_window_size = vaf_context_window_size // 2
    context_vaf_features = []
    for offset in range(-half_window_size, half_window_size + 1):
      prefix = 'minus' if offset < 0 else 'plus'
      context_vaf_features.append(
          f'{VARIANT_ALLELE_FREQUENCY_AT_PREFIX}_{prefix}_{abs(offset)}'
      )
    return context_vaf_features

  def _encode_candidate_example(
      self, candidate, label, read_phases, is_training
  ):
    """Encodes a candidate example.

    Args:
      candidate: The candidate to be encoded.
      label: The label of the candidate.
      read_phases: A dictionary mapping read names to haplotype tags.
      is_training: Whether the example is for training or inference.

    Returns:
      A list of encoded features.
    """
    feature_set = (
        self._training_features if is_training else self._inference_features
    )
    feature_encoder = FeatureEncoder(candidate, label)
    candidate_example = [
        feature_encoder.encode_feature(feature) for feature in feature_set
    ]
    if self.expand_by_haplotype:
      for haplotype in self.HP_TAGS:
        feature_encoder = FeatureEncoder(
            candidate, label, haplotype, read_phases
        )
        candidate_example.extend([
            feature_encoder.encode_feature(feature.value)
            for feature in HAPLOTYPE_FEATURE_SET
        ])
    return candidate_example

  def encode_training_examples(
      self,
      candidates_with_label: Sequence[
          Tuple[deepvariant_pb2.DeepVariantCall, variant_labeler.VariantLabel]
      ],
      read_phases: Dict[str, int],
  ) -> Sequence[Sequence[Union[str, int]]]:
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
          candidate, label, read_phases, is_training=True
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
      List[Sequence[Union[str, int]]],
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
      candidate_example = self._encode_candidate_example(
          candidate, None, read_phases, is_training=False
      )
      examples.append(candidate_example)
    return skipped_candidates, kept_candidates, examples
