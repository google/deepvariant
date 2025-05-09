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

r"""Script to preprocess a truth VCF file.

This script pre-processes a truth VCF file to normalize the alleles and set the
genotype correctly. Mostly used for the T2T truth VCFs that have overlapping
variants with star alleles.

Please note that this script should only be used for truth VCFs that are
completely phased. Without phasing inormation, this script will not be able
to normalize the alleles correctly.

Example usage:
blaze run -c opt //learning/genomics/deepvariant/opensource_only/tools:preprocess_truth -- \
--truth_vcf /data/kishwar/program_outputs/Q100_truth/update_20250117/GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z2k.intersected.vcf.gz \
--output_vcf /data/kishwar/program_outputs/Q100_truth/update_20250117/GRCh38_HG2-T2TQ100-V1.1_smvar.star_normalized.vcf.gz \
--logtostderr
"""

import collections
from collections.abc import Sequence
from typing import Any
from absl import app
from absl import flags
from absl import logging
import pysam

_TRUTH_VCF = flags.DEFINE_string(
    'truth_vcf', None, 'Input truth VCF file containing truth phased variants.'
)
_OUTPUT_VCF = flags.DEFINE_string(
    'output_vcf', None, 'Output VCF file containing truth phased variants.'
)


def get_genotype(variant: pysam.VariantRecord) -> Any:
  """Returns the genotype of a variant."""
  for sample in variant.samples:
    return variant.samples[sample]['GT']


def normalize_alleles(
    variant: pysam.VariantRecord,
    allele_map: collections.defaultdict[int, list[str]],
) -> collections.defaultdict[int, list[str]]:
  """Normalizes a variant."""
  ref_allele = variant.alleles[0]
  genotype = list(get_genotype(variant))
  genotype_index = 1
  # Now go through the alts and add them to the map.
  for allele_index in genotype:
    alt_allele = variant.alleles[allele_index]
    if alt_allele == '*':
      genotype_index += 1
      continue
    # Check if the reference allele has common suffix with the alt allele
    common_suffix_length = 0
    if len(alt_allele) > 1 and len(ref_allele) > 1:
      i = len(ref_allele) - 1
      j = len(alt_allele) - 1
      while i > 0 and j > 0 and ref_allele[i] == alt_allele[j]:
        common_suffix_length += 1
        i -= 1
        j -= 1
    allele_stops_at = variant.stop - common_suffix_length
    alt_allele_normalized = alt_allele
    ref_allele_normalized = ref_allele
    ref_allele_stop = variant.stop
    if common_suffix_length > 0:
      alt_allele_normalized = alt_allele[:-common_suffix_length]
      ref_allele_normalized = ref_allele[:-common_suffix_length]
      ref_allele_stop = variant.stop - common_suffix_length

    for pos in range(variant.start, allele_stops_at):
      if pos not in allele_map:
        raise ValueError(
            'Alt allele at position %d is not in the allele map for variant %s.'
            % pos,
            variant,
        )
      base_index = pos - variant.start
      if base_index < len(alt_allele_normalized):
        alt_allele = alt_allele_normalized[base_index]
      else:
        alt_allele = '*'
      if pos == ref_allele_stop or len(ref_allele_normalized) == 1:
        alt_allele = alt_allele_normalized[base_index:]
      allele_map[pos][genotype_index] = alt_allele
    genotype_index += 1
  return allele_map


def get_alleles_from_map(
    allele_map: collections.defaultdict[int, list[str]],
) -> tuple[list[str], list[int]]:
  """Returns the alleles from the allele map."""
  alleles = ['', '', '']
  # First item is position that is unused.
  for _, allele_map_entry in allele_map.items():
    if allele_map_entry[0] != '*':
      alleles[0] += allele_map_entry[0]
    if allele_map_entry[1] != '*':
      alleles[1] += allele_map_entry[1]
    if allele_map_entry[2] != '*':
      alleles[2] += allele_map_entry[2]

  if len(alleles[1]) < len(alleles[2]):
    return alleles, [1, 2]
  else:
    return [alleles[0], alleles[2], alleles[1]], [2, 1]


def set_variant_genotype(
    variant: pysam.VariantRecord, genotype: list[int]
) -> None:
  """Sets the genotype of a variant."""
  for sample in variant.samples:
    variant.samples[sample]['GT'] = genotype
    variant.samples[sample].phased = True


def resolve_overlapping_variants(input_vcf: str, output_vcf: str) -> None:
  """Reads a VCF file and returns a pysam.VCF object."""
  variant_reader = pysam.VariantFile(input_vcf, mode='r')
  bcf_out = pysam.VariantFile(output_vcf, 'w', header=variant_reader.header)
  current_contig = None
  # Group variants by start and end position
  grouped_variants = []
  current_group = []
  last_pos = None

  for rec in variant_reader.fetch():
    # If we are switching contigs.
    if current_contig is None:
      current_contig = rec.contig
    elif rec.contig != current_contig:
      grouped_variants.append(current_group)
      current_group = []
      current_contig = rec.contig
      continue

    # Construct groups of variants that are overlapping.
    if not current_group:
      current_group.append(rec)
      last_pos = rec.stop
    elif (
        rec.contig == current_contig
        and rec.start < last_pos
        and rec.contig not in ['chrX', 'chrY']
    ):
      current_group.append(rec)
      last_pos = max(last_pos, rec.stop)
    else:
      grouped_variants.append(current_group)
      current_group = [rec]
      last_pos = rec.stop
  # If the last group is not empty.
  if current_group:
    grouped_variants.append(current_group)

  for variant_group in grouped_variants:
    if len(variant_group) > 1:
      allele_map = collections.defaultdict(list[str])
      group_stop_position = variant_group[0].stop
      # Generate reference map
      for variant in variant_group:
        group_stop_position = max(group_stop_position, variant.stop)
        # First create a map of position to reference allele
        for pos in range(variant.start, variant.stop):
          if pos not in allele_map:
            allele_map[pos] = ['', '', '']
            allele_map[pos][0] = variant.alleles[0][pos - variant.start]
            allele_map[pos][1] = variant.alleles[0][pos - variant.start]
            allele_map[pos][2] = variant.alleles[0][pos - variant.start]
          elif allele_map[pos][0] != variant.alleles[0][pos - variant.start]:
            raise ValueError(
                'Allele at position %d is different from the reference allele'
                ' for variant %s.' % pos,
                variant,
            )
      # Simple variant consolidation
      for variant in variant_group:
        allele_map = normalize_alleles(variant, allele_map)
      alleles, genotype = get_alleles_from_map(allele_map)
      variant_group[0].alleles = alleles
      set_variant_genotype(variant_group[0], genotype)
      bcf_out.write(variant_group[0])
    else:
      bcf_out.write(variant_group[0])


def main(argv: Sequence[str]) -> None:
  del argv
  logging.info('Truth VCF: %s', _TRUTH_VCF.value)
  resolve_overlapping_variants(_TRUTH_VCF.value, _OUTPUT_VCF.value)


if __name__ == '__main__':
  logging.use_python_logging()
  app.run(main)
