# Copyright 2025 Google LLC.
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
"""A labeler that combines positional and haplotype-based labeling."""

from deepvariant.labeler import haplotype_labeler
from deepvariant.labeler import positional_labeler
from deepvariant.labeler import variant_labeler
from third_party.nucleus.util import variant_utils


class CombinedLabeler(variant_labeler.VariantLabeler):
  """A labeler that combines positional and haplotype-based labeling.

  This labeler gives precedence to the haplotype-based labeler's output.
  It falls back to the positional labeler's output only when the
  haplotype labeler either fails to produce a label or calls a variant
  as homozygous-reference (0/0). This strategy aims to leverage the higher
  accuracy of the haplotype labeler in complex cases, while using the
  positional labeler as a backup.
  """

  def __init__(
      self,
      truth_vcf_reader,
      ref_reader,
      confident_regions,
      max_group_size=haplotype_labeler._MAX_GROUP_SIZE,
      max_separation=haplotype_labeler._MAX_SEPARATION_WITHIN_VARIANT_GROUP,
      max_gt_options_product=haplotype_labeler._MAX_GT_OPTIONS_PRODUCT,
  ):
    """Initializes a new CombinedLabeler.

    Args:
      truth_vcf_reader: A VcfReader object for the truth variants.
      ref_reader: A FastaReader object for the reference genome.
      confident_regions: A RangeSet of confident regions.
      max_group_size: The maximum number of variants to group together for
        haplotype labeling.
      max_separation: The maximum separation between variants within a group for
        haplotype labeling.
      max_gt_options_product: The maximum product of genotype options for
        haplotype labeling.
    """
    super().__init__(
        truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions
    )
    self._positional_labeler = positional_labeler.PositionalVariantLabeler(
        truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions
    )
    self._haplotype_labeler = haplotype_labeler.HaplotypeLabeler(
        truth_vcf_reader=truth_vcf_reader,
        ref_reader=ref_reader,
        confident_regions=confident_regions,
        max_group_size=max_group_size,
        max_separation=max_separation,
        max_gt_options_product=max_gt_options_product,
    )

  def variant_key(self, variant):
    return f'{variant.reference_name}-{variant.start}-{variant.reference_bases}'

  def label_variants(self, variants, region):
    """Labels the given variants.

    Args:
      variants: An iterable of variants to label.
      region: The region containing the variants.

    Yields:
      variant_labeler.VariantLabel objects for each input variant.
    """
    variants = list(variants)
    positional_labels = self._positional_labeler.label_variants(
        variants, region
    )
    haplotype_labels = self._haplotype_labeler.label_variants(variants, region)

    positional_labels_by_key = {
        self.variant_key(l.variant): l for l in positional_labels
    }
    haplotype_labels_by_key = {
        self.variant_key(l.variant): l for l in haplotype_labels
    }

    for variant in variants:
      key = self.variant_key(variant)
      positional_label = positional_labels_by_key.get(key)
      haplotype_label = haplotype_labels_by_key.get(key)
      if not positional_label and not haplotype_label:
        continue

      if variant_utils.is_snp(variant):
        if haplotype_label:
          yield haplotype_label
      else:
        if haplotype_label and haplotype_label.genotype != (0, 0):
          yield haplotype_label
        elif positional_label:
          yield positional_label
