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

import dataclasses
import enum
import itertools
from typing import Sequence

import tensorflow as tf

from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import variant_utils

# Example encoding keys
FEATURES_ENCODED = 'features/encoded'
IDS_ENCODED = 'ids/encoded'
LABEL_ENCODED = 'label/encoded'
GENOTYPE_ENCODED = 'genotype/encoded'


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
  ALT = 'alt'


class TruthFeature(SmallModelFeature):
  """A single truth label for the candidate."""

  GENOTYPE = 'genotype'


class BaseFeature(SmallModelFeature):
  """Set of basic features that can be computed over a subset of reads."""

  NUM_READS_SUPPORTS_REF = 'num_reads_supports_ref'
  NUM_READS_SUPPORTS_ALT = 'num_reads_supports_alt'
  ALT_INDICES_DEPTH = 'alt_indices_depth'
  TOTAL_DEPTH = 'total_depth'
  VARIANT_ALLELE_FREQUENCY = 'variant_allele_frequency'
  ALT_INDICES_VARIANT_ALLELE_FREQUENCY = 'alt_indices_variant_allele_frequency'
  REF_MAPPING_QUALITY = 'ref_mapping_quality'
  ALT_MAPPING_QUALITY = 'alt_mapping_quality'
  REF_BASE_QUALITY = 'ref_base_quality'
  ALT_BASE_QUALITY = 'alt_base_quality'
  REF_REVERSE_STRAND_RATIO = 'ref_reverse_strand_ratio'
  ALT_REVERSE_STRAND_RATIO = 'alt_reverse_strand_ratio'


class VariantFeature(SmallModelFeature):
  """Set of features that describe the variant candidate."""

  IS_SNP = 'is_snp'
  IS_INSERTION = 'is_insertion'
  IS_DELETION = 'is_deletion'
  INSERTION_LENGTH = 'insertion_length'
  DELETION_LENGTH = 'deletion_length'
  IS_MULTIALLELIC = 'is_multiallelic'
  IS_MULTIPLE_ALT_ALLELES = 'is_multiple_alt_alleles'


ENCODING_BY_GENOTYPE = {
    (0, 0): GenotypeEncoding.REF.value,
    (0, 1): GenotypeEncoding.HET.value,
    (1, 0): GenotypeEncoding.HET.value,
    (1, 1): GenotypeEncoding.HOM_ALT.value,
    (0,): GenotypeEncoding.REF.value,
    (1,): GenotypeEncoding.HOM_ALT.value,
}
FAKE_CANDIDATE = deepvariant_pb2.DeepVariantCall(
    variant=variants_pb2.Variant(alternate_bases=['']),
    allele_support_ext={
        '': deepvariant_pb2.DeepVariantCall.SupportingReadsExt(read_infos=[])
    },
)
DEFAULT_ALT_ALLELE_INDICES = (0,)


def _mean_for_attribute(
    read_infos: Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport],
    attribute_name: str,
    multiplier: int = 1,
) -> int:
  """Calculates the mean of an attribute for a list of reads."""
  if not read_infos:
    return 0
  return (
      multiplier * sum(getattr(r, attribute_name) for r in read_infos)
  ) // len(read_infos)


def _filter_by_haplotype(
    read_infos: Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport],
    read_phases: dict[str, int],
    haplotype: int | None = None,
) -> Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport]:
  """Filters the reads by haplotype."""
  return [r for r in read_infos if read_phases.get(r.read_name, 0) == haplotype]


def _filter_by_sample(
    read_infos: Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport],
    sample_name: str,
) -> Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport]:
  """Filters the reads by sample."""
  return [r for r in read_infos if r.sample_name == sample_name]


def _get_context_allele_frequency_offsets(
    vaf_context_window_size: int,
) -> Sequence[int]:
  """Returns the list of offsets for the given context window size."""
  if not vaf_context_window_size:
    return []
  half_window_size = vaf_context_window_size // 2
  return list(range(-half_window_size, half_window_size + 1))


def get_set_of_allele_indices(
    candidate: deepvariant_pb2.DeepVariantCall,
) -> list[tuple[int, ...]]:
  """Returns the complete set of allele indices for the candidate."""
  num_alt_alleles = len(candidate.variant.alternate_bases)
  biallelic = [(i,) for i in range(num_alt_alleles)]
  multiallelic = list(itertools.combinations(range(num_alt_alleles), 2))
  return biallelic + multiallelic


def _get_alt_read_infos(
    candidate: deepvariant_pb2.DeepVariantCall,
    alt_allele_indices: tuple[int, ...],
) -> Sequence[deepvariant_pb2.DeepVariantCall.ReadSupport]:
  """Returns the read infos for the given allele index, where -1 is ref."""
  read_infos = []
  for allele_index in alt_allele_indices:
    alt_base = candidate.variant.alternate_bases[allele_index]
    if alt_base not in candidate.allele_support_ext:
      raise ValueError(f'Alt base "{alt_base}" not found in candidate.')
    read_infos.extend(candidate.allele_support_ext[alt_base].read_infos)
  return read_infos


def get_exclude_alleles(
    candidate: deepvariant_pb2.DeepVariantCall,
    alt_allele_indices: tuple[int, ...],
) -> list[str]:
  """Returns the list of alleles to exclude from the candidate."""
  return [
      alternate_bases
      for i, alternate_bases in enumerate(candidate.variant.alternate_bases)
      if i not in alt_allele_indices
  ]


class FeatureEncoder:
  """Class for encoding small model features for a single candidate."""

  def __init__(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      alt_allele_indices: tuple[int, ...],
      sample: str | None = None,
      haplotype: int | None = None,
      read_phases: dict[str, int] | None = None,
  ):
    """Initializes the feature encoder.

    Args:
      candidate: the candidate proto from which to extract the feature.
      alt_allele_indices: the indices of the alt alleles to consider when
        computing the feature values.
      sample: (optional) the sample name to use for filtering reads.
      haplotype: (optional) the haplotype tag to use for filtering reads.
      read_phases: (optional)a dictionary of read names to haplotype phases.
    """
    self.candidate = candidate
    self.haplotype = haplotype
    self.alt_allele_indices = alt_allele_indices
    self.ref_read_infos = self.candidate.ref_support_ext.read_infos
    self.alt_read_infos = _get_alt_read_infos(
        candidate, self.alt_allele_indices
    )
    self.all_alt_read_infos = _get_alt_read_infos(
        candidate, tuple(range(len(candidate.variant.alternate_bases)))
    )
    if sample:
      self.ref_read_infos = _filter_by_sample(self.ref_read_infos, sample)
      self.alt_read_infos = _filter_by_sample(self.alt_read_infos, sample)
      self.all_alt_read_infos = _filter_by_sample(
          self.all_alt_read_infos, sample
      )
    if haplotype is not None and read_phases is not None:
      self.ref_read_infos = _filter_by_haplotype(
          self.ref_read_infos, read_phases, haplotype
      )
      self.alt_read_infos = _filter_by_haplotype(
          self.alt_read_infos, read_phases, haplotype
      )
      self.all_alt_read_infos = _filter_by_haplotype(
          self.all_alt_read_infos, read_phases, haplotype
      )
    self.ref_read_infos_count = len(self.ref_read_infos)
    self.alt_read_infos_count = len(self.alt_read_infos)
    self.all_alt_read_infos_count = len(self.all_alt_read_infos)
    self.exclude_alleles = get_exclude_alleles(candidate, alt_allele_indices)

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

  def _get_alt(self) -> str:
    """Returns the first alternate base of the candidate."""
    return '|'.join(
        self.candidate.variant.alternate_bases[allele_index]
        for allele_index in self.alt_allele_indices
    )

  def _get_num_reads_supports_ref(self) -> int:
    """Returns the number of reads supporting the reference base."""
    return self.ref_read_infos_count

  def _get_num_reads_supports_alt(self) -> int:
    """Returns the number of reads supporting the alternate base."""
    return self.alt_read_infos_count

  def _get_alt_indices_depth(self) -> int:
    """Returns the depth of the alt indices."""
    return self.ref_read_infos_count + self.alt_read_infos_count

  def _get_total_depth(self) -> int:
    """Returns the total depth of the candidate."""
    return len(self.candidate.ref_support_ext.read_infos) + sum(
        len(r.read_infos) for r in self.candidate.allele_support_ext.values()
    )

  def _get_alt_indices_variant_allele_frequency(self) -> int:
    """Returns the ref allele frequency of the candidate."""
    dp = self._get_alt_indices_depth()
    if dp == 0:
      return 0
    return 100 * self.alt_read_infos_count // dp

  def _get_variant_allele_frequency(self) -> int:
    """Returns the variant allele frequency of the candidate."""
    dp = self._get_total_depth()
    if dp == 0:
      return 0
    return 100 * self.alt_read_infos_count // dp

  def _get_ref_mapping_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.ref_read_infos, 'mapping_quality')

  def _get_alt_mapping_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.alt_read_infos, 'mapping_quality')

  def _get_ref_base_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.ref_read_infos, 'average_base_quality')

  def _get_alt_base_quality(self) -> int:
    """Returns the mapping quality of the candidate."""
    return _mean_for_attribute(self.alt_read_infos, 'average_base_quality')

  def _get_ref_reverse_strand_ratio(self) -> int:
    """Returns the reverse strand ratio of the candidate."""
    return _mean_for_attribute(
        self.ref_read_infos, 'is_reverse_strand', multiplier=100
    )

  def _get_alt_reverse_strand_ratio(self) -> int:
    """Returns the reverse strand ratio of the candidate."""
    return _mean_for_attribute(
        self.alt_read_infos, 'is_reverse_strand', multiplier=100
    )

  def _get_is_snp(self) -> int:
    """Returns 1 if the candidate is a SNP, 0 otherwise."""
    return int(
        variant_utils.is_snp(self.candidate.variant, self.exclude_alleles)
    )

  def _get_is_insertion(self) -> int:
    """Returns 1 if the candidate is an insertion, 0 otherwise."""
    return int(
        variant_utils.variant_is_insertion(
            self.candidate.variant, self.exclude_alleles
        )
    )

  def _get_is_deletion(self) -> int:
    """Returns 1 if the candidate is an insertion, 0 otherwise."""
    return int(
        variant_utils.variant_is_deletion(
            self.candidate.variant, self.exclude_alleles
        )
    )

  def _get_is_multiallelic(self) -> int:
    """Returns 1 if the candidate is a multiallelic, 0 otherwise."""
    return int(variant_utils.is_multiallelic(self.candidate.variant))

  def _get_insertion_length(self) -> int:
    """Returns the insertion length of the candidate."""
    insertion_lengths = [
        len(self.candidate.variant.alternate_bases[i])
        - len(self.candidate.variant.reference_bases)
        for i in self.alt_allele_indices
    ]
    return max(0, max(insertion_lengths))

  def _get_deletion_length(self) -> int:
    """Returns the deletion length of the candidate."""
    deletion_lengths = [
        len(self.candidate.variant.reference_bases)
        - len(self.candidate.variant.alternate_bases[i])
        for i in self.alt_allele_indices
    ]
    return max(0, max(deletion_lengths))

  def _get_is_multiple_alt_alleles(self) -> int:
    """Returns the length of the alt allele indices."""
    return int(len(self.alt_allele_indices) > 1)

  def _genotype_label(self, label: variant_labeler.VariantLabel) -> int:
    """Returns the genotype of the candidate."""
    if variant_utils.is_biallelic(self.candidate.variant):
      return ENCODING_BY_GENOTYPE[label.genotype]
    return label.label_for_alt_alleles(self.alt_allele_indices)

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
    elif feature == BaseFeature.NUM_READS_SUPPORTS_ALT:
      return self._get_num_reads_supports_alt()
    elif feature == BaseFeature.ALT_INDICES_DEPTH:
      return self._get_alt_indices_depth()
    elif feature == BaseFeature.TOTAL_DEPTH:
      return self._get_total_depth()
    elif feature == BaseFeature.ALT_INDICES_VARIANT_ALLELE_FREQUENCY:
      return self._get_alt_indices_variant_allele_frequency()
    elif feature == BaseFeature.VARIANT_ALLELE_FREQUENCY:
      return self._get_variant_allele_frequency()
    elif feature == BaseFeature.REF_MAPPING_QUALITY:
      return self._get_ref_mapping_quality()
    elif feature == BaseFeature.ALT_MAPPING_QUALITY:
      return self._get_alt_mapping_quality()
    elif feature == BaseFeature.REF_BASE_QUALITY:
      return self._get_ref_base_quality()
    elif feature == BaseFeature.ALT_BASE_QUALITY:
      return self._get_alt_base_quality()
    elif feature == BaseFeature.REF_REVERSE_STRAND_RATIO:
      return self._get_ref_reverse_strand_ratio()
    elif feature == BaseFeature.ALT_REVERSE_STRAND_RATIO:
      return self._get_alt_reverse_strand_ratio()
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
    elif feature == VariantFeature.IS_MULTIALLELIC:
      return self._get_is_multiallelic()
    elif feature == VariantFeature.IS_MULTIPLE_ALT_ALLELES:
      return self._get_is_multiple_alt_alleles()
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
    elif feature == IdentifyingFeature.ALT:
      return self._get_alt()
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
    genotype_encoding = self._genotype_label(label)
    value = [0 for _ in GenotypeEncoding]
    value[genotype_encoding] = 1
    return value


@dataclasses.dataclass
class InferenceExampleSet:
  """A set of small model examples."""

  skipped_candidates: list[deepvariant_pb2.DeepVariantCall] = dataclasses.field(
      default_factory=list
  )
  candidates_with_alt_allele_indices: list[
      tuple[deepvariant_pb2.DeepVariantCall, tuple[int, ...]]
  ] = dataclasses.field(default_factory=list)
  inference_examples: list[Sequence[int]] = dataclasses.field(
      default_factory=list
  )


class SmallModelExampleFactory:
  """Class for making small model examples."""

  def __init__(
      self,
      vaf_context_window_size: int,
      sample_names: Sequence[str],
      accept_snps: bool = True,
      accept_indels: bool = True,
      accept_multiallelics: bool = True,
      expand_by_haplotype: bool = False,
      model_features: Sequence[str] | None = None,
  ):
    self.sample_names = sample_names
    self.accept_snps = accept_snps
    self.accept_indels = accept_indels
    self.accept_multiallelics = accept_multiallelics
    self.expand_by_haplotype = expand_by_haplotype
    self.vaf_context_window_size = vaf_context_window_size
    self.model_features = self._get_and_validate_model_features(model_features)

  def _get_and_validate_model_features(
      self,
      model_features: Sequence[str] | None,
  ) -> Sequence[str]:
    """Returns the model features for the small model."""
    all_features = list(
        self._encode_candidate_feature_dict(
            candidate=FAKE_CANDIDATE,
            alt_allele_indices=DEFAULT_ALT_ALLELE_INDICES,
            read_phases=None,
            sample_order=list(range(len(self.sample_names))),
        ).keys()
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
      return self.accept_multiallelics
    elif variant_utils.is_snp(candidate.variant):
      return self.accept_snps
    else:
      return self.accept_indels

  def _encode_candidate_feature_dict(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      alt_allele_indices: tuple[int, ...],
      read_phases: dict[str, int] | None,
      sample_order: Sequence[int],
  ) -> dict[str, int]:
    """Encodes all model features for a given candidate into a key-value dict."""
    candidate_example = {}
    feature_encoder = FeatureEncoder(candidate, alt_allele_indices, sample=None)
    candidate_example.update({
        feature.value: feature_encoder.encode_base_feature(feature)
        for feature in BaseFeature
    })
    if len(self.sample_names) > 1:
      for sample_index in sample_order:
        sample_name = self.sample_names[sample_index]
        sample_feature_encoder = FeatureEncoder(
            candidate, alt_allele_indices, sample_name
        )
        candidate_example.update({
            f'{sample_name}_{feature.value}': (
                sample_feature_encoder.encode_base_feature(feature)
            )
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
      for sample_index in sample_order:
        sample_name = self.sample_names[sample_index]
        for haplotype in Haplotype:
          haplotype_feature_encoder = FeatureEncoder(
              candidate,
              alt_allele_indices,
              sample_name,
              haplotype.value,
              read_phases,
          )
          for feature in BaseFeature:
            candidate_example.update({
                f'{sample_name}_{feature.value}_hp_{haplotype.value}': (
                    haplotype_feature_encoder.encode_base_feature(feature)
                )
            })
    return candidate_example

  def _encode_model_features(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      alt_allele_indices: tuple[int, ...],
      read_phases: dict[str, int] | None,
      sample_order: Sequence[int],
  ) -> Sequence[int]:
    """Encodes a candidate example into the selected model features.

    Args:
      candidate: The candidate to be encoded.
      alt_allele_indices: The alt-allele indices to use for the candidate.
      read_phases: A dictionary mapping read names to haplotype tags.
      sample_order: The order in which the samples are to be encoded.

    Returns:
      A feature vector.
    """
    encoded_candidate = self._encode_candidate_feature_dict(
        candidate, alt_allele_indices, read_phases, sample_order
    )
    return [
        value
        for feature, value in encoded_candidate.items()
        if feature in self.model_features
    ]

  def _encode_training_example(
      self,
      candidate: deepvariant_pb2.DeepVariantCall,
      alt_allele_indices: tuple[int, ...],
      label: variant_labeler.VariantLabel,
      read_phases: dict[str, int],
      sample_order: Sequence[int],
  ) -> tf.train.Example:
    """Encodes a candidate example."""
    feature_encoder = FeatureEncoder(candidate, alt_allele_indices)
    return tf.train.Example(
        features=tf.train.Features(
            feature={
                FEATURES_ENCODED: tf.train.Feature(
                    int64_list=tf.train.Int64List(
                        value=self._encode_model_features(
                            candidate,
                            alt_allele_indices,
                            read_phases,
                            sample_order,
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
                GENOTYPE_ENCODED: tf.train.Feature(
                    int64_list=tf.train.Int64List(value=label.genotype)
                ),
            }
        )
    )

  def encode_training_examples(
      self,
      candidates_with_label: Sequence[
          tuple[deepvariant_pb2.DeepVariantCall, variant_labeler.VariantLabel]
      ],
      read_phases: dict[str, int],
      sample_order: Sequence[int],
  ) -> Sequence[tf.train.Example]:
    """Generates examples from the given candidates for training.

    Args:
      candidates_with_label: list of candidates with labels to be processed into
        examples.
      read_phases: A dictionary mapping read names to haplotype tags.
      sample_order: The order in which the samples are to be encoded.

    Returns:
      A list of encoded candidate examples.
    """
    training_examples = []
    for candidate, label in candidates_with_label:
      if not self._pass_candidate_to_small_model(candidate):
        continue
      alt_allele_indices_set = get_set_of_allele_indices(candidate)
      for alt_allele_indices in alt_allele_indices_set:
        candidate_example = self._encode_training_example(
            candidate, alt_allele_indices, label, read_phases, sample_order
        )
        training_examples.append(candidate_example)
    return training_examples

  def encode_inference_examples(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      read_phases: dict[str, int],
      sample_order: Sequence[int],
  ) -> InferenceExampleSet:
    """Generates examples from the given candidates for inference.

    Args:
      candidates: list of candidates to be processed into a summary.
      read_phases: A dictionary mapping read names to haplotype tags.
      sample_order: The order in which the samples are to be encoded.

    Returns:
      A InferenceExampleSet containing:
        skipped_candidates: A list of all candidates that were skipped.
        candidates_with_alt_allele_indices: A list of candidates and
        alt-allele-indices pairs for which
        examples were generated.
        inference_examples: A list of encoded candidate examples.
    """
    example_set = InferenceExampleSet()
    for candidate in candidates:
      if not self._pass_candidate_to_small_model(candidate):
        example_set.skipped_candidates.append(candidate)
        continue
      alt_allele_indices_set = get_set_of_allele_indices(candidate)
      for alt_allele_indices in alt_allele_indices_set:
        example_set.candidates_with_alt_allele_indices.append(
            (candidate, alt_allele_indices)
        )
        candidate_example = self._encode_model_features(
            candidate, alt_allele_indices, read_phases, sample_order
        )
        example_set.inference_examples.append(candidate_example)
    return example_set
