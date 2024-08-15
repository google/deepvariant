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
from typing import List, Sequence, Tuple, Union

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

ENCODING_BY_GENOTYPE = {
    (0, 0): GenotypeEncoding.REF.value,
    (0, 1): GenotypeEncoding.HET.value,
    (1, 0): GenotypeEncoding.HET.value,
    (1, 1): GenotypeEncoding.HOM_ALT.value,
}


def _get_contig(candidate: deepvariant_pb2.DeepVariantCall) -> str:
  """Returns the contig of the candidate."""
  return candidate.variant.reference_name


def _get_start(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns the start position of the candidate."""
  return candidate.variant.start


def _get_end(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns the end position of the candidate."""
  return candidate.variant.end


def _get_ref(candidate: deepvariant_pb2.DeepVariantCall) -> str:
  """Returns the reference base of the candidate."""
  return candidate.variant.reference_bases


def _get_alt_1(candidate: deepvariant_pb2.DeepVariantCall) -> str:
  """Returns the first alternate base of the candidate."""
  return candidate.variant.alternate_bases[0]


def _get_num_reads_supports_ref(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the number of reads supporting the reference base."""
  return candidate.variant.calls[0].info['AD'].values[0].int_value


def _get_num_reads_supports_alt_1(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the number of reads supporting the alternate base."""
  return candidate.variant.calls[0].info['AD'].values[1].int_value


def _get_total_depth(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the total depth of the candidate."""
  return candidate.variant.calls[0].info['DP'].values[0].int_value


def _get_variant_allele_frequency_1(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the variant allele frequency of the candidate."""
  return int(
      candidate.variant.calls[0].info['VAF'].values[0].number_value * 100
  )


def _get_ref_mapping_quality(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns the mapping quality of the candidate."""
  read_infos = candidate.ref_support_ext.read_infos
  if not read_infos:
    return 0
  return sum(r.mapping_quality for r in read_infos) // len(read_infos)


def _get_alt_1_mapping_quality(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the mapping quality of the candidate."""
  alt_base = candidate.variant.alternate_bases[0]
  read_infos = []
  if alt_base in candidate.allele_support_ext:
    read_infos = candidate.allele_support_ext[alt_base].read_infos
  if not read_infos:
    return 0
  return sum(r.mapping_quality for r in read_infos) // len(read_infos)


def _get_ref_base_quality(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns the mapping quality of the candidate."""
  read_infos = candidate.ref_support_ext.read_infos
  if not read_infos:
    return 0
  return sum(r.average_base_quality for r in read_infos) // len(read_infos)


def _get_alt_1_base_quality(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the mapping quality of the candidate."""
  alt_base = candidate.variant.alternate_bases[0]
  read_infos = []
  if alt_base in candidate.allele_support_ext:
    read_infos = candidate.allele_support_ext[alt_base].read_infos
  if not read_infos:
    return 0
  return sum(r.average_base_quality for r in read_infos) // len(read_infos)


def _get_variant_allele_frequency_at_position(
    feature_name: str,
    candidate: deepvariant_pb2.DeepVariantCall,
) -> int:
  """Returns the VAF for the candidate at the given position."""
  offset = int(re.findall('[0-9]+', feature_name)[0])
  plus_or_minus = -1 if 'minus' in feature_name else 1
  position = candidate.variant.start + (plus_or_minus * offset)
  return candidate.allele_frequency_at_position.get(position, 0)


def _get_is_snp(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns 1 if the candidate is a SNP, 0 otherwise."""
  return int(variant_utils.is_snp(candidate.variant))


def _get_is_insertion(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns 1 if the candidate is an insertion, 0 otherwise."""
  return int(
      variant_utils.is_insertion(
          candidate.variant.reference_bases,
          candidate.variant.alternate_bases[0],
      )
  )


def _get_is_deletion(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns 1 if the candidate is an insertion, 0 otherwise."""
  return int(
      variant_utils.is_deletion(
          candidate.variant.reference_bases,
          candidate.variant.alternate_bases[0],
      )
  )


def _get_insertion_length(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns the insertion length of the candidate."""
  insertion_length = len(candidate.variant.alternate_bases[0]) - len(
      candidate.variant.reference_bases
  )
  return max(0, insertion_length)


def _get_deletion_length(candidate: deepvariant_pb2.DeepVariantCall) -> int:
  """Returns the deletion length of the candidate."""
  deletion_length = len(candidate.variant.reference_bases) - len(
      candidate.variant.alternate_bases[0]
  )
  return max(0, deletion_length)


def encode_feature(
    feature: str,
    candidate: deepvariant_pb2.DeepVariantCall,
    label: Union[variant_labeler.VariantLabel, None],
) -> Union[int, str]:
  """Encodes the given feature from a candidate or label.

  Args:
    feature: specifies which feature to extract
    candidate: the candidate proto from which to extract the feature.
    label: the variant label proto that contains the truth information.

  Returns:
    The extracted value for that feature for the candidate.
  """
  if feature == SmallModelFeature.CONTIG.value:
    return _get_contig(candidate)
  elif feature == SmallModelFeature.START.value:
    return _get_start(candidate)
  elif feature == SmallModelFeature.END.value:
    return _get_end(candidate)
  elif feature == SmallModelFeature.REF.value:
    return _get_ref(candidate)
  elif feature == SmallModelFeature.ALT_1.value:
    return _get_alt_1(candidate)
  elif feature == SmallModelFeature.NUM_READS_SUPPORTS_REF.value:
    return _get_num_reads_supports_ref(candidate)
  elif feature == SmallModelFeature.NUM_READS_SUPPORTS_ALT_1.value:
    return _get_num_reads_supports_alt_1(candidate)
  elif feature == SmallModelFeature.TOTAL_DEPTH.value:
    return _get_total_depth(candidate)
  elif feature == SmallModelFeature.VARIANT_ALLELE_FREQUENCY_1.value:
    return _get_variant_allele_frequency_1(candidate)
  elif feature == SmallModelFeature.REF_MAPPING_QUALITY.value:
    return _get_ref_mapping_quality(candidate)
  elif feature == SmallModelFeature.ALT_1_MAPPING_QUALITY.value:
    return _get_alt_1_mapping_quality(candidate)
  elif feature == SmallModelFeature.REF_BASE_QUALITY.value:
    return _get_ref_base_quality(candidate)
  elif feature == SmallModelFeature.ALT_1_BASE_QUALITY.value:
    return _get_alt_1_base_quality(candidate)
  elif feature.startswith(VARIANT_ALLELE_FREQUENCY_AT_PREFIX):
    return _get_variant_allele_frequency_at_position(feature, candidate)
  elif feature == SmallModelFeature.IS_SNP.value:
    return _get_is_snp(candidate)
  elif feature == SmallModelFeature.IS_INSERTION.value:
    return _get_is_insertion(candidate)
  elif feature == SmallModelFeature.IS_DELETION.value:
    return _get_is_deletion(candidate)
  elif feature == SmallModelFeature.INSERTION_LENGTH.value:
    return _get_insertion_length(candidate)
  elif feature == SmallModelFeature.DELETION_LENGTH.value:
    return _get_deletion_length(candidate)
  elif feature == SmallModelFeature.GENOTYPE.value and label:
    return ENCODING_BY_GENOTYPE[label.genotype]
  else:
    raise ValueError(f'{feature} does not map to a callable.')


class SmallModelExampleFactory:
  """Class for making small model examples."""

  def __init__(
      self,
      vaf_context_window_size: int,
      accept_snps: bool = True,
      accept_indels: bool = True,
  ):
    self.accept_snps = accept_snps
    self.accept_indels = accept_indels
    context_vaf_features = self._get_context_vaf_features(
        vaf_context_window_size
    )
    self.training_features = [f.value for f in SmallModelFeature]
    self.inference_features = [
        feature.value
        for feature in SmallModelFeature
        if feature not in IDENTIFYING_FEATURES and feature != TRUTH_FEATURE
    ]
    self.inference_features.extend(context_vaf_features)
    self.training_features.extend(context_vaf_features)

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

  def encode_training_examples(
      self,
      candidates_with_label: Sequence[
          Tuple[deepvariant_pb2.DeepVariantCall, variant_labeler.VariantLabel]
      ],
  ) -> Sequence[Sequence[Union[str, int]]]:
    """Generates examples from the given candidates for training.

    Args:
      candidates_with_label: List of candidates with labels to be processed into
        examples.

    Returns:
      A list of encoded candidate examples.
    """
    candidate_examples = []
    for candidate, label in candidates_with_label:
      if not self._pass_candidate_to_small_model(candidate):
        continue
      candidate_example = [
          encode_feature(feature, candidate, label)
          for feature in self.training_features
      ]
      candidate_examples.append(candidate_example)
    return candidate_examples

  def encode_inference_examples(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
  ) -> Tuple[
      List[deepvariant_pb2.DeepVariantCall],
      List[deepvariant_pb2.DeepVariantCall],
      List[Sequence[Union[str, int]]],
  ]:
    """Generates examples from the given candidates for inference.

    Args:
      candidates: List of candidates to be processed into a summary.

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
      candidate_example = [
          encode_feature(feature, candidate, None)
          for feature in self.inference_features
      ]
      examples.append(candidate_example)
    return skipped_candidates, kept_candidates, examples
