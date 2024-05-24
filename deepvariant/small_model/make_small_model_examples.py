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
from typing import Sequence, Tuple, Union

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
  # Features passed to the model
  NUM_READS_SUPPORTS_REF = 'num_reads_supports_ref'
  NUM_READS_SUPPORTS_ALT_1 = 'num_reads_supports_alt_1'
  TOTAL_DEPTH = 'total_depth'
  VARIANT_ALLELE_FREQUENCY_1 = 'variant_allele_frequency_1'
  # Truth feature
  GENOTYPE = 'genotype'


IDENTIFYING_FEATURES = (
    SmallModelFeature.CONTIG,
    SmallModelFeature.START,
    SmallModelFeature.END,
    SmallModelFeature.REF,
    SmallModelFeature.ALT_1,
)
TRUTH_FEATURE = SmallModelFeature.GENOTYPE
MODEL_FEATURES = [
    feature
    for feature in SmallModelFeature
    if feature not in IDENTIFYING_FEATURES and feature != TRUTH_FEATURE
]


ENCODING_BY_GENOTYPE = {
    (0, 0): GenotypeEncoding.REF.value,
    (0, 1): GenotypeEncoding.HET.value,
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
  return candidate.variant.reference_bases[0]


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


def _pass_candidate_to_small_model(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> bool:
  """Determines if the candidate is eligible for the small model."""
  return (
      variant_utils.is_snp(candidate.variant)
      and len(candidate.variant.alternate_bases) == 1
  )


def encode_feature(
    feature: SmallModelFeature,
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
  if feature == SmallModelFeature.CONTIG:
    return _get_contig(candidate)
  elif feature == SmallModelFeature.START:
    return _get_start(candidate)
  elif feature == SmallModelFeature.END:
    return _get_end(candidate)
  elif feature == SmallModelFeature.REF:
    return _get_ref(candidate)
  elif feature == SmallModelFeature.ALT_1:
    return _get_alt_1(candidate)
  elif feature == SmallModelFeature.NUM_READS_SUPPORTS_REF:
    return _get_num_reads_supports_ref(candidate)
  elif feature == SmallModelFeature.NUM_READS_SUPPORTS_ALT_1:
    return _get_num_reads_supports_alt_1(candidate)
  elif feature == SmallModelFeature.TOTAL_DEPTH:
    return _get_total_depth(candidate)
  elif feature == SmallModelFeature.VARIANT_ALLELE_FREQUENCY_1:
    return _get_variant_allele_frequency_1(candidate)
  elif feature == SmallModelFeature.GENOTYPE and label:
    return ENCODING_BY_GENOTYPE[label.genotype]
  else:
    raise ValueError(f'{feature} does not map to a callable.')


def get_example_feature_columns() -> Sequence[str]:
  """Produces an ordered list of all feature as columns in TSV example.

  Returns:
    A list of feature names.
  """
  return [f.value for f in SmallModelFeature]


def generate_training_examples(
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
    if not _pass_candidate_to_small_model(candidate):
      continue
    candidate_example = [
        encode_feature(feature, candidate, label)
        for feature in SmallModelFeature
    ]
    candidate_examples.append(candidate_example)
  return candidate_examples


def generate_inference_examples(
    candidates: Sequence[deepvariant_pb2.DeepVariantCall],
) -> Tuple[
    Sequence[deepvariant_pb2.DeepVariantCall],
    Sequence[deepvariant_pb2.DeepVariantCall],
    Sequence[Sequence[Union[str, int]]],
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
    if not _pass_candidate_to_small_model(candidate):
      skipped_candidates.append(candidate)
      continue
    kept_candidates.append(candidate)
    candidate_example = [
        encode_feature(feature, candidate, None) for feature in MODEL_FEATURES
    ]
    examples.append(candidate_example)
  return skipped_candidates, kept_candidates, examples
