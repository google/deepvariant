# Copyright 2020 Google LLC.
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
"""Functionality for including allele frequencies in DeepVariant pileups."""

import collections
from typing import DefaultDict, Sequence, Optional

from absl import logging

from third_party.nucleus.io import vcf
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from deepvariant.protos import deepvariant_pb2


def get_allele_frequency(variant, index):
  """Gets allele frequency of the index-th alt_base of a Variant proto.

  Args:
    variant: A Variant proto.
    index: The index where we want to query allele frequency.

  Returns:
    A float. The queried allele frequency.

  Raises:
    ValueError: If the Variant proto does not include 'AF' field or the 'AF'
      field encounters an out-of-bound error when querying at position `index`
  """
  if variant.info.get('AF'):
    if index < len(variant.info['AF'].values):
      return variant.info['AF'].values[index].number_value
    else:
      raise ValueError(
          'Invalid index',
          index,
          'for the info[AF] field',
          variant.info['AF'].values,
      )
  raise ValueError('Variant does not have an AF field')


def get_ref_allele_frequency(variant):
  """Gets REF allele frequency for a Variant proto."""
  sum_alt_frequency = 0
  for alt_idx, _ in enumerate(variant.alternate_bases):
    sum_alt_frequency += get_allele_frequency(variant, alt_idx)
  return 1 - sum_alt_frequency


def get_ref_haplotype_and_offset(dv_variant, cohort_variants, ref_reader):
  """Gets reference haplotype that overlaps with given variants and its offset.

  Reference offset is the starting position of the returned haplotype.
  It is the minimal starting positions of dv_variant and cohort_variants.
  Reference haplotype is the reference sequence from reference_offset to the
  maximal ending position of dv_variant and cohort_variants.

  Args:
    dv_variant: A Variant proto.
    cohort_variants: A list of Variant protos.
    ref_reader: A fasta.IndexedFastaReader.

  Returns:
    A tuple of a string and an integer.
      - String: Reference haplotype overlapping with dv_variant and
                cohort_variants.
      - Integer: Offset of the reference haplotype.

  Raises:
    ValueError: If the queried region is invalid.
  """
  # Get the range for haplotypes to be compared.
  min_start_position = min(
      dv_variant.start, min([cv.start for cv in cohort_variants])
  )
  max_end_position = max(
      dv_variant.end, max([cv.end for cv in cohort_variants])
  )

  haplotype_region = ranges.make_range(
      dv_variant.reference_name, min_start_position, max_end_position
  )

  if ref_reader.is_valid(haplotype_region):
    return ref_reader.query(haplotype_region), min_start_position
  else:
    raise ValueError('Invalid IndexedFastaReader region', haplotype_region)


def update_haplotype(variant, reference_haplotype, reference_offset):
  """Updates haplotypes for a variant.

  A list of variant haplotypes are updated given a variant and a reference
  haplotype (this consists of a sequence and an offset wrt to the reference).
  All ALT alleles are updated as independent updated haplotypes.

  Args:
    variant: A Variant proto.
    reference_haplotype: A string extracted from the reference genome.
    reference_offset: An integer. The offset of the starting position of
      reference_haplotype on reference.

  Raises:
    ValueError: Variant.start is smaller than reference_offset.

  Returns:
    A list of haplotype objects. Haplotype objects are stored as dicts:
      {'haplotype': a haplotype (string),
       'alt': an alt allele (string),
       'variant': a Variant proto}
  """
  if variant.start < reference_offset:
    raise ValueError(
        'The starting position of a variant is smaller than its ',
        'corresponding reference offset',
        variant.start,
        reference_offset,
    )
  offset_start = variant.start - reference_offset
  offset_suffix = (
      variant.start + len(variant.reference_bases) - reference_offset
  )
  list_updated_haplotype = []
  for biallelic_variant in variant.alternate_bases:
    updated_haplotype = (
        reference_haplotype[:offset_start]
        + biallelic_variant
        + reference_haplotype[offset_suffix:]
    )
    dict_haplotype = {
        'haplotype': updated_haplotype,
        'alt': biallelic_variant,
        'variant': variant,
    }
    list_updated_haplotype.append(dict_haplotype)

  return list_updated_haplotype


def match_candidate_and_cohort_haplotypes(
    candidate_haps, cohort_haps_and_freqs
):
  """Match candidate haplotypes with cohort haplotypes and update frequency.

  First, we look for exact haplotype matches between candidate and cohorts.
  If there're any matches, the REF allele frequency associated with the matching
  ALT allele is updated as well.

  Second, if no matches are found, we try to find inexact matches, where only
  REF alleles are matched. The inexact matching step is only used to update REF
  allele frequency. If no exact and inexact matches are found, set REF allele
  frequency to 1.

  Args:
    candidate_haps: A list of haplotype objects from a candidate.
    cohort_haps_and_freqs: A list of haplotype objects from cohorts. Haplotype
      objects are stored as dicts: {'haplotype': a haplotype (string), 'alt': an
      alt allele (string), 'variant': a Variant proto}

  Returns:
    A dict with candidate alt alleles as keys, and associated frequencies
    as values.
  """
  dict_allele_frequency = {}
  for candidate_obj in candidate_haps:
    candidate_haplotype = candidate_obj['haplotype']
    candidate_alt = candidate_obj['alt']
    candidate_variant = candidate_obj['variant']

    for cohort_obj in cohort_haps_and_freqs:
      cohort_haplotype = cohort_obj['haplotype']
      # Exact haplotype match.
      if candidate_haplotype == cohort_haplotype:
        cohort_variant = cohort_obj['variant']
        cohort_frequency = get_allele_frequency(
            cohort_variant,
            list(cohort_variant.alternate_bases).index(cohort_obj['alt']),
        )
        dict_allele_frequency[candidate_alt] = cohort_frequency

        # Update REF frequency if it is not in the dictionary.
        if not dict_allele_frequency.get(candidate_variant.reference_bases):
          dict_allele_frequency[candidate_variant.reference_bases] = (
              get_ref_allele_frequency(cohort_variant)
          )

    # For an unmatched alt allele, set the frequency to 0.
    if not dict_allele_frequency.get(candidate_alt):
      dict_allele_frequency[candidate_alt] = 0

  # Calculate REF allele frequency if no exact match was found.
  # It is possible a novel mutation happens at a site where there are other
  # cohort variants. In this case, we cannot simply set REF frequency to 1.
  if sum(dict_allele_frequency.values()) == 0:
    candidate = candidate_haps[0]['variant']
    # Left align variants.
    s_candidate = variant_utils.simplify_variant_alleles(candidate)
    for cohort_obj in cohort_haps_and_freqs:
      s_cohort_variant = variant_utils.simplify_variant_alleles(
          cohort_obj['variant']
      )
      # Try to find inexact matches to set REF allele frequency.
      # Inexact matches here mean only REF alleles match.
      if (
          s_candidate.start == s_cohort_variant.start
          and s_candidate.reference_bases == s_cohort_variant.reference_bases
      ):
        dict_allele_frequency[s_candidate.reference_bases] = (
            get_ref_allele_frequency(s_cohort_variant)
        )

    # If still no match, set REF allele frequency to 1.
    if not dict_allele_frequency.get(candidate.reference_bases):
      dict_allele_frequency[candidate.reference_bases] = 1

  return dict_allele_frequency


def find_matching_allele_frequency(
    variant, population_vcf_reader, ref_reader, padding_bases=0
):
  """Finds the allele frequencies of all the alt alleles for a candidate.

  Args:
    variant: A Variant proto generated by make_examples. Note that it can be
      multi-allelic.
    population_vcf_reader: A VcfReader object that reads associated VCF file for
      a candidate. We want to extract allele frequency information in the VCF.
    ref_reader: A IndexedFastaReader object that reads the reference FASTA.
    padding_bases: An integer that specifies the number of padding bases added
      when performing a VCF query. By default it is set to 0.

  Returns:
    A dict with alleles as keys, and allele frequencies as values
  """
  query_region = ranges.make_range(
      chrom=variant.reference_name,
      start=variant.start - padding_bases,
      end=variant.end + padding_bases,
  )
  # Convert to list because we'll look through cohort_variants more than once.
  try:
    cohort_variants = list(population_vcf_reader.query(query_region))
  except ValueError as e:
    err_msg = str(e)
    if err_msg.startswith("NOT_FOUND: Unknown reference_name '"):
      # Split the message by single quotes (') to get the contig name
      parts = err_msg.split("'")
      contig_name = parts[1]
      logging.warning('population_vcf does not have config %s', contig_name)
      cohort_variants = []
    else:
      raise e

  # Init allele frequency dict using alt bases in the candidate.
  dict_allele_frequency = {}
  for alt_base in variant.alternate_bases:
    dict_allele_frequency[alt_base] = 0

  try:
    reference_haplotype, reference_offset = get_ref_haplotype_and_offset(
        variant, cohort_variants, ref_reader
    )
  except ValueError:
    # If the range associated with variant and cohort_variants is invalid,
    # assume this candidate does not have any ALT alleles.
    dict_allele_frequency = {}
    dict_allele_frequency[variant.reference_bases] = 1
    for alt in variant.alternate_bases:
      dict_allele_frequency[alt] = 0
    return dict_allele_frequency

  candidate_haps = update_haplotype(
      variant, reference_haplotype, reference_offset
  )
  cohort_haps = []
  for cohort_variant in cohort_variants:
    cohort_haps.extend(
        update_haplotype(cohort_variant, reference_haplotype, reference_offset)
    )

  for c in candidate_haps:
    logging.debug('candidate %s, %s', c['haplotype'], c['alt'])
  for c in cohort_haps:
    logging.debug('cohort    %s, %s', c['haplotype'], c['alt'])

  dict_allele_frequency = match_candidate_and_cohort_haplotypes(
      candidate_haps, cohort_haps
  )

  if dict_allele_frequency:
    logging.vlog(
        3,
        'dict_allele_frequency: %s:%d-%d, %s > %s',
        variant.reference_name,
        variant.start,
        variant.end,
        variant.reference_bases,
        dict_allele_frequency,
    )

  return dict_allele_frequency


def make_population_vcf_readers(
    population_vcf_filenames: Sequence[str],
) -> DefaultDict[str, Optional[vcf.VcfReader]]:
  """Creates VcfReaders for the given VCF file paths, organized by reference.

  VcfReaders can be made either from a single VCF that covers all the relevant
  reference sequences or strictly one VCF per reference sequence. By returning
  a defaultdict, any code using the output of this function does not have to
  consider whether there are multiple VCFs or not, it can simply query by
  chromosome and get a reader.

  Args:
    population_vcf_filenames: Paths to files (VCF or VCF.gz) with population
      genotypes.

  Raises:
    ValueError: If there is more than one VCF file containing variants
      from the same chromosome.

  Returns:
    A defaultdict that maps from a reference name to an associated VcfReader.
    If there was only one VCF provided, all references will map to that one
    reader. If more than one VCF was provided, the references will have a
    reader each, while any that were not included will map to None.
  """
  # If only one VCF file is provided.
  if len(population_vcf_filenames) == 1:
    # The DefaultDict allows later code to query for any chromosome and still
    # get the same reader. This is great for compatibility with multi-VCF below.
    return collections.defaultdict(
        lambda: vcf.VcfReader(population_vcf_filenames[0])
    )

  # If more than one VCF files are provided.
  population_vcf_readers = DefaultDict(lambda: None)

  for vcf_filename in population_vcf_filenames:
    population_vcf_reader = vcf.VcfReader(vcf_filename, header=None)

    # Get contig name from the first variant in a file.
    for var in population_vcf_reader:
      reference_name = var.reference_name
      break
    # There should not be more than one VCFs including variants in
    # reference_name.
    if population_vcf_readers.get(reference_name):
      raise ValueError(
          'Variants on %s are included in multiple VCFs' % reference_name
      )
    population_vcf_readers[reference_name] = population_vcf_reader

  return population_vcf_readers


def add_allele_frequencies_to_candidates(
    candidates, population_vcf_reader, ref_reader
):
  """Adds allele frequencies for candidate variants.

  Args:
    candidates: Iterable of DeepVariantCall protos that are the candidates we
      want to process.
    population_vcf_reader: A VcfReader object that reads the associated
      population VCF file for candidates. None if the contig is not found.
    ref_reader: A fasta.IndexedFastaReader.

  Yields:
    DeepVariantCall protos. The same set of input candidates, with field
      allele_frequency filled.
  """
  for candidate in candidates:
    if population_vcf_reader:
      dict_allele_frequency = find_matching_allele_frequency(
          variant=candidate.variant,
          population_vcf_reader=population_vcf_reader,
          ref_reader=ref_reader,
      )
    else:
      # Set ALT frequencies to 0 if population_vcf_reader is None.
      dict_allele_frequency = {}
      dict_allele_frequency[candidate.variant.reference_bases] = 1
      for alt in candidate.variant.alternate_bases:
        dict_allele_frequency[alt] = 0

    yield deepvariant_pb2.DeepVariantCall(
        variant=candidate.variant,
        allele_support=candidate.allele_support,
        allele_frequency=dict_allele_frequency,
    )
