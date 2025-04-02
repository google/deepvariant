# Copyright 2017 Google LLC.
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
"""Core functionality for step one of DeepVariant: Making examples."""

import collections
import itertools
import json
import math
import os
import random
import re
import time
from typing import Any, DefaultDict, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

from absl import flags
from absl import logging
from etils import epath
import numpy as np


from deepvariant import allele_frequency
from deepvariant import calling_regions_utils
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import dv_vcf_constants
from deepvariant import resources
from deepvariant import sample as sample_lib
from deepvariant import variant_caller as vc_base
from deepvariant import vcf_candidate_importer
from deepvariant import very_sensitive_caller
from deepvariant.labeler import customized_classes_labeler
from deepvariant.labeler import haplotype_labeler
from deepvariant.labeler import positional_labeler
from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.python import direct_phasing
from deepvariant.python import make_examples_native as make_examples_native_module
from deepvariant.python import methylation_aware_phasing
from deepvariant.python import pileup_image_native
from deepvariant.realigner import realigner as realigner_module
from deepvariant.small_model import inference as small_model_inference
from deepvariant.small_model import make_small_model_examples
from deepvariant.vendor import timer
from google.protobuf import text_format
from third_party.nucleus.io import fasta
from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import sam
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import struct_utils
from third_party.nucleus.util import utils
from third_party.nucleus.util import variant_utils
# pylint: disable=g-direct-tensorflow-import
from tensorflow.python.lib.io import tf_record
# pylint: enable=g-direct-tensorflow-import

# For --runtime_by_region, these columns will be written out in this order.
RUNTIME_BY_REGION_COLUMNS = (
    'region',
    'get reads',
    'find candidates',
    'make pileup images',
    'write outputs',
    'num reads',
    'num candidates',
    'num examples',
    'small model generate examples',
    'small model call examples',
    'small model write variants',
    'small model total',
)

# For --read_phases_output, these columns will be written out in this order.
READ_PHASES_OUTPUT_COLUMNS = ('fragment_name', 'phase', 'region_order')

# The name used for a sample if one is not specified or present in the reads.
_UNKNOWN_SAMPLE = 'UNKNOWN'

# candidate position -1 designates the end of region.
END_OF_REGION = -1

# candidate position -2 designate the end of partition. This is used to merge
# candidate positions from all shards.
END_OF_PARTITION = -2

# Maximum length of partition in bases. It is limited by available memory.
# TODO: For better flexibility it may be benefitial to expose it as a
# flag.
MAX_PARTITION_LEN = 1000000

# Non DNA regions larger than this value are excluded from processing.
MIN_NON_DNA_REGION = 300000

# Number of regions to use for computing mean coverage.
NUM_REGIONS_FOR_MEAN_COVERAGE = 1500

# Number of loci to use for computing mean coverage.
NUM_LOCI_FOR_MEAN_COVERAGE = 10

# ---------------------------------------------------------------------------
# Selecting variants of specific types (e.g., SNPs)
# ---------------------------------------------------------------------------


def _select_biallelic_snps(v):
  return variant_utils.is_snp(v) and variant_utils.is_biallelic(v)


def _select_biallelic_indels(v):
  return variant_utils.is_indel(v) and variant_utils.is_biallelic(v)


def _select_biallelic_insertions(v):
  return variant_utils.has_insertion(v) and variant_utils.is_biallelic(v)


def _select_biallelic_deletions(v):
  return variant_utils.has_deletion(v) and variant_utils.is_biallelic(v)


VARIANT_TYPE_SELECTORS = {
    'snps': _select_biallelic_snps,
    'indels': _select_biallelic_indels,
    'insertions': _select_biallelic_insertions,
    'deletions': _select_biallelic_deletions,
    'multi-allelics': variant_utils.is_multiallelic,
    'all': lambda v: True,
}

# ---------------------------------------------------------------------------
# Option handling
# ---------------------------------------------------------------------------


def assign_sample_name(sample_name_flag: str, reads_filenames: str) -> str:
  """Returns sample name derived from either sample_name flag or input BAM.

  Function derives sample_name from the flag. If flag is not set then
  sample_name is derived from input BAM.

  Args:
    sample_name_flag: string. sample_name flag value.
    reads_filenames: A list of filenames of an alignments file, e.g. BAM. The
      first of these will be used. May be empty.

  Returns:
    string. Derived sample name.
  """
  if sample_name_flag:
    sample_name = sample_name_flag
  elif reads_filenames:
    with sam.SamReader(reads_filenames.split(',')[0]) as sam_reader:
      sample_name = extract_sample_name_from_sam_reader(sam_reader)
  else:
    sample_name = _UNKNOWN_SAMPLE
  return sample_name


def make_vc_options(
    sample_name: str, flags_obj: flags.FlagValues
) -> deepvariant_pb2.VariantCallerOptions:
  haploid_contigs_str = flags_obj.haploid_contigs or ''
  return deepvariant_pb2.VariantCallerOptions(
      min_count_snps=flags_obj.vsc_min_count_snps,
      min_count_indels=flags_obj.vsc_min_count_indels,
      min_fraction_snps=flags_obj.vsc_min_fraction_snps,
      min_fraction_indels=flags_obj.vsc_min_fraction_indels,
      min_fraction_multiplier=flags_obj.vsc_min_fraction_multiplier,
      max_fraction_indels_for_non_target_sample=flags_obj.vsc_max_fraction_indels_for_non_target_sample,
      max_fraction_snps_for_non_target_sample=flags_obj.vsc_max_fraction_snps_for_non_target_sample,
      # Not specified by default: fraction_reference_sites_to_emit,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=1400605801,
      sample_name=sample_name,
      p_error=flags_obj.p_error,
      max_gq=50,
      gq_resolution=flags_obj.gvcf_gq_binsize,
      ploidy=2,
      skip_uncalled_genotypes=flags_obj.mode == 'training',
      phase_reads_region_padding_pct=dv_constants.PHASE_READS_REGION_PADDING_PCT,
      track_ref_reads=flags_obj.track_ref_reads,
      small_model_vaf_context_window_size=flags_obj.small_model_vaf_context_window_size,
      haploid_contigs=haploid_contigs_str.split(','),
      par_regions_bed=flags_obj.par_regions_bed,
      create_complex_alleles=flags_obj.create_complex_alleles,
      enable_methylation_aware_phasing=flags_obj.enable_methylation_aware_phasing,
      exclude_contigs_for_methylation_phasing=flags_obj.exclude_contigs_for_methylation_phasing,
  )


def parse_proto_enum_flag(
    proto_enum_pb2, flag_value, skip_unspecified_option=True
):
  """Parses a command line flag string value into a protobuf Enum value.

  Args:
    proto_enum_pb2: a enum_type_wrapper.EnumTypeWrapper type containing a proto
      enum definition. For example, this would be
      deepvariant_pb2.MakeExamplesOptions.Mode to get the MakeExamplesOptions
      Mode enum. See:
      https://developers.google.com/protocol-buffers/docs/reference/python-generated#enum
        for more information.
    flag_value: str. The name of the proto enum option from the command line we
      want to convert into the enum value.
    skip_unspecified_option: bool. If True, any enum options that include the
      string 'unspecified' (in any case) will be excluded from the list of
      allowed options in the ValueError raised if flag_value isn't valid.

  Returns:
    The enum value for flag_value in proto_enum_pb2

  Raises:
    ValueError: if flag_value isn't a valid enum name in proto_enum_pb2.
  """
  try:
    return proto_enum_pb2.Value(flag_value)
  except ValueError as exception:
    options = proto_enum_pb2.keys()
    if skip_unspecified_option:
      options = [o for o in options if 'unspecified' not in o.lower()]
    raise ValueError(
        'Unknown enum option "{}". Allowed options are {}'.format(
            flag_value, ','.join(sorted(options))
        )
    ) from exception


def resolve_sam_aux_fields(
    flags_obj: flags.FlagValues,
    provided_channels: List[str],
) -> List[str]:
  """Determines which auxiliary fields to parse based on passed options.

  Previously, options were used to set the aux fields to parse, but the list
  of AUX tags we require can be determined based on the set of flags
  or channels that are passed in.

  Args:
    flags_obj: DeepVariant Flag Options.
    provided_channels: The list of channels provided by the user.

  Returns:
    A list of auxiliary fields to parse.
  """
  aux_fields = set()

  # If the user is not phasing reads natively, but wants to sort by existing HP
  # haplotype tag, then add HP to aux fields to parse.
  if not flags_obj.phase_reads:
    if (
        flags_obj.sort_by_haplotypes
        or flags_obj.reverse_haplotypes
        or flags_obj.hp_tag_for_assembly_polishing
        or 'haplotype' in provided_channels
    ):
      logging.info(
          'Parsing HP AUX tag because --phase_reads=False and'
          ' --sort_by_haplotypes, --reverse_haplotypes,'
          ' --hp_tag_for_assembly_polishing, and/or haplotype channel is'
          ' present.'
      )
      aux_fields.add('HP')

  # Original quality scores require the OQ tag.
  if flags_obj.use_original_quality_scores:
    logging.info(
        'Parsing OQ AUX tag because --use_original_quality_scores is set.'
    )
    aux_fields.add('OQ')

  # Add fields required for channels.
  for base_mod_channel in ['base_methylation', 'base_6ma']:
    if base_mod_channel in provided_channels:
      logging.info(
          'Parsing MM, ML, and MN AUX tags because of base modification'
          ' channel.'
      )
      aux_fields.update(['MM', 'ML', 'MN'])

  if (
      flags_obj.enable_methylation_calling
      or flags_obj.enable_methylation_aware_phasing
  ):
    logging.info(
        'Parsing MM, ML, and MN AUX tags because --enable_methylation_calling'
        ' is set.'
    )
    aux_fields.update(['MM', 'ML', 'MN'])

  return list(aux_fields)


def logging_with_options(
    options: deepvariant_pb2.MakeExamplesOptions, message: str
):
  """If options contain multiple shards, log with task/shard prefix."""
  if options.num_shards > 1:
    prefix = 'Task {}/{}: '.format(options.task_id, options.num_shards)
  else:
    prefix = ''
  logging.info('%s%s', prefix, message)


def log_summary_stats(
    options: deepvariant_pb2.MakeExamplesOptions,
    n_stats: dict[str, int],
) -> None:
  """Prints the summary stats in a neatly-formatted way.

  Example output:
  '''
  Summary stats:
        2869 candidate variants found
        2500 candidate variants phased
         253 examples written
        2493 small model examples called
  '''

  Args:
    options: The MakeExamplesOptions proto.
    n_stats: A dictionary of stats to log.
  """
  stat_to_description = {
      'n_candidates': 'candidate variants found',
      'n_examples': 'examples written',
  }
  if options.phase_reads:
    stat_to_description['n_phased_candidates'] = 'candidate variants phased'
  if options.write_small_model_examples:
    stat_to_description['n_small_model_examples'] = (
        'small model examples written'
    )
  if options.call_small_model_examples:
    stat_to_description['n_small_model_calls'] = 'small model examples called'

  stats_to_log = {k: v for k, v in n_stats.items() if k in stat_to_description}
  longest_stat = max((len(str(x)) for x in stats_to_log.values()))
  message_rows = ['Summary stats:']
  for stat_name, description in stat_to_description.items():
    message = f'\t{str(n_stats[stat_name]).rjust(longest_stat)} {description}'
    message_rows.append(message)
  logging_with_options(options, '\n'.join(message_rows))


# ---------------------------------------------------------------------------
# Simple utilities
# ---------------------------------------------------------------------------


def in_training_mode(options):
  return options.mode == deepvariant_pb2.MakeExamplesOptions.TRAINING


def in_calling_mode(options):
  return options.mode == deepvariant_pb2.MakeExamplesOptions.CALLING


def in_candidate_sweep_mode(options):
  return options.mode == deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP


def gvcf_output_enabled(options):
  """Returns True if we should be generating gVCF output."""
  return bool(options.gvcf_filename)


def only_true(
    *elts: List[reference_pb2.ContigInfo],
) -> List[List[reference_pb2.ContigInfo]]:
  """Returns the sublist of elements that evaluate to True."""
  return [elt for elt in elts if elt]


def extract_sample_name_from_sam_reader(sam_reader):
  """Returns the sample name as derived from the BAM file of reads.

  Args:
    sam_reader: Already opened sam_reader to use to extract the sample names
      from. This sam_reader will not be closed after this function returns.

  Returns:
    The sample ID annotated in the read group.
  """
  samples_list = [
      rg.sample_id for rg in sam_reader.header.read_groups if rg.sample_id
  ]
  samples = set(samples_list)
  if not samples:
    logging.warning(
        (
            'No non-empty sample name found in the input reads. '
            'DeepVariant will use %s as the sample name. You can also '
            'provide a sample name with the --sample_name argument.'
        ),
        dv_constants.DEFAULT_SAMPLE_NAME,
    )
    return dv_constants.DEFAULT_SAMPLE_NAME
  elif len(samples) > 1:
    logging.warning(
        (
            'Multiple samples (%s) were found in the input reads. '
            'Please confirm this is intended. For now, DeepVariant '
            'will use the first sample name %s.'
        ),
        ', '.join(sorted(samples)),
        samples_list[0],
    )
    return samples_list[0]
  return next(iter(samples))


def trim_runtime(seconds: float) -> float:
  """Round seconds (float) to the nearest millisecond."""
  return round(seconds, 3)


# ---------------------------------------------------------------------------
# Utilities for working with labeling metrics
#
# ---------------------------------------------------------------------------


def read_make_examples_run_info(path):
  """Reads a MakeExamplesRunInfo proto in text_format from path."""
  with epath.Path(path).open() as f:
    return text_format.Parse(f.read(), deepvariant_pb2.MakeExamplesRunInfo())


def write_make_examples_run_info(run_info_proto, path):
  """Writes a MakeExamplesRunInfo proto in text_format to path."""
  with epath.Path(path).open('w') as writer:
    writer.write(
        '# proto-file: learning/genomics/deepvariant/protos/deepvariant.proto\n'
        '# proto-message: MakeExamplesRunInfo\n'
    )
    writer.write(text_format.MessageToString(run_info_proto, float_format=''))


# ---------------------------------------------------------------------------
# Region processing
# ---------------------------------------------------------------------------


def _ensure_consistent_contigs(
    ref_contigs: List[reference_pb2.ContigInfo],
    sam_contigs: List[reference_pb2.ContigInfo],
    vcf_contigs: Optional[List[reference_pb2.ContigInfo]],
    exclude_contig_names: Optional[Sequence[str]] = None,
    min_coverage_fraction: float = 1.0,
):
  """Returns the common contigs after ensuring 'enough' overlap.

  Args:
    ref_contigs: list of reference_pb2.ContigInfo protos in the reference
      genome.
    sam_contigs: list of reference_pb2.ContigInfo protos in the SAM/BAM file.
    vcf_contigs: list of reference_pb2.ContigInfo protos in the VCF if in
      training mode, or None otherwise.
    exclude_contig_names: list of strings of contig names to exclude from
      overlap consideration.
    min_coverage_fraction: The fraction of the reference contigs that must be
      shared with all inputs.

  Returns:
    The list of contigs common between all input sources.

  Raises:
    ValueError: The contigs are not sufficiently similar across input sources.
  """
  # Remove any excluded contigs from the ref_contigs, as we want to use the
  # selected contigs for our overlap comparison.
  if exclude_contig_names:
    ref_contigs = [c for c in ref_contigs if c.name not in exclude_contig_names]

  # Compute the common contigs among our inputs, and check that the contigs are
  # sufficiently consistent among each other.
  contigs = common_contigs(only_true(ref_contigs, sam_contigs))
  if vcf_contigs:
    # If VCF contigs exist, we just check the name (not the length).
    vcf_contigs_names = set([x.name for x in vcf_contigs])
    contigs = [x for x in contigs if x.name in vcf_contigs_names]
  validate_reference_contig_coverage(
      ref_contigs, contigs, min_coverage_fraction
  )
  return contigs


def common_contigs(
    contigs_list: List[List[reference_pb2.ContigInfo]],
) -> List[reference_pb2.ContigInfo]:
  """Gets a list of contigs found in all contigs in contigs_list.

  A common contig is considered one where the name and length in basepairs are
  the same.

  Args:
    contigs_list: A sequence of lists of ContigInfo protos.

  Returns:
    A list of ContigInfo protos. Note that the individual protos found in this
    returned list are shared with the ContigInfo protos found in contigs_list,
    so should not be modified.
  """

  def common2(
      contigs1: List[reference_pb2.ContigInfo],
      contigs2: List[reference_pb2.ContigInfo],
  ) -> List[reference_pb2.ContigInfo]:
    """Computes the common contigs between contigs1 and contigs2."""
    map2 = ranges.contigs_dict(contigs2)

    def is_common(contig1: reference_pb2.ContigInfo) -> bool:
      contig2 = map2.get(contig1.name, None)
      return contig2 and contig1.n_bases == contig2.n_bases  # pytype: disable=bad-return-type

    return [c for c in contigs1 if is_common(c)]

  # Compute the common contigs by recursively getting common contigs of our
  # cumulative set of contigs (common) and each contig in other_contigs.
  common = contigs_list[0]
  for other_contigs in contigs_list[1:]:
    common = common2(common, other_contigs)

  return common


def validate_reference_contig_coverage(
    ref_contigs: List[reference_pb2.ContigInfo],
    shared_contigs: List[reference_pb2.ContigInfo],
    min_coverage_fraction: float,
):
  """Validates that shared_contigs spans a sufficient amount of ref_contigs.

  Args:
    ref_contigs: List of ContigInfo protos. All of the contigs from our
      reference genome.
    shared_contigs: The subset of ref_contigs that we found in common with
      ref_contigs and all other genomics data sources.
    min_coverage_fraction: The minimum fraction of basepairs of ref_contigs that
      should be found among the shared_contigs.

  Raises:
    ValueError: If the fraction of covered bases is less than
      min_coverage_fraction.
  """

  def format_contig_matches():
    pieces = []
    common_map = ranges.contigs_dict(shared_contigs)
    for ref_contig in ref_contigs:
      status = 'matched' if ref_contig.name in common_map else 'IS MISSING'
      pieces.append(
          '\n"{}" is {} bp and {}'.format(
              ref_contig.name, ref_contig.n_bases, status
          )
      )
    return ', '.join(pieces)

  ref_bp = ranges.contigs_n_bases(ref_contigs)
  common_bp = ranges.contigs_n_bases(shared_contigs)
  coverage = common_bp / (1.0 * ref_bp)
  if not shared_contigs or coverage < min_coverage_fraction:
    raise ValueError(
        'Reference contigs span {} bases but only {} bases '
        '({:.2%}) were found in common among our input files. '
        'Check that the sources were created on a common genome '
        'reference build. Contig matches were: {}. Here is a '
        'useful article about different human genome reference '
        'builds:\n'
        'https://gatkforums.broadinstitute.org/gatk/discussion/'
        '11010/human-genome-reference-builds-grch38-hg38-b37-hg19'
        '\nPlease make sure the --ref input matches the build '
        'used for the input in --reads.'.format(
            ref_bp, common_bp, coverage, format_contig_matches()
        )
    )


def find_ref_n_regions(
    ref_reader: genomics_reader.GenomicsReader, min_region_len: int
) -> List[range_pb2.Range]:
  """Returns List[nucleus.genomics.v1.Range] regions containing Ns.

  Args:
    ref_reader: genomics_reader.GenomicsReader. Nucleus Fasta reader.
    min_region_len: int. Only regions larger than min_region_len are returned.

  Returns:
    A List of nucleus.genomics.v1.Range containing regions of Ns in the
    reference.
  """
  ref_n_regions = []
  # ref_reader returns tuples of contig_name and vector of bases.
  for ref_name, bases in ref_reader.iterate():
    i = min_region_len - 1
    while i < len(bases):
      b = bases[i]
      if b not in vc_base.CANONICAL_DNA_BASES:
        # Seek back to find the start of N-region.
        j = i
        while j > 0 and bases[j] not in vc_base.CANONICAL_DNA_BASES:
          j -= 1
        start = j if j == 0 else j + 1
        # Seek forward to find the end of N-region.
        j = i
        while j < len(bases) and bases[j] not in vc_base.CANONICAL_DNA_BASES:
          j += 1
        end = j
        if end - start >= min_region_len:
          logging.info('Excluding %s:%d-%d', ref_name, start, end)
          ref_n_regions.append(ranges.make_range(ref_name, start, end))
        i = end
      else:
        i += min_region_len - 1
  return ref_n_regions


def partition_by_candidates(
    regions: ranges.RangeSet, candidate_positions: List[int], max_size: int
) -> List[range_pb2.Range]:
  """Splits our intervals so that none contain more than max_size candidates.

  Slices up the intervals in this RangeSet into a equivalent set of intervals
  (i.e., spanning the same number of candidates), each of which contains at most
  max_size candidates.

  This function does not modify this RangeSet.


  Args:
    regions: All calling regions that we need to process in all shards.
    candidate_positions: List[int]. List of candidate positions. Candidate
      positions are sorted within each genome contig. Candidates in each contig
      are terminated with END_OF_REGION.
    max_size: int > 0. The maximum number of candidates per interval.

  Returns:
    nucleus.genomics.v1.Range protos, in sorted order.

  Raises:
    ValueError: if max_size <= 0.
  """
  if max_size <= 0:
    raise ValueError('max_size must be > 0: {}'.format(max_size))

  partitioned = []
  candidate_it = 0  # Current candidate's index
  # Iterate over regions and cut partitions with equal number of candidates.
  # candidate positions are local to for each region.
  # Candidate position = END_OF_REGION designates the last candidate in a
  # region. (Example: 1, 5, 7, -1, 2, 4, 7, -1)
  for interval in regions:
    num_of_candidates = 0
    refname = interval.reference_name
    partition_start = interval.start
    partition_end = interval.start
    # candidate_positions[candidate_it] == END_OF_REGION designates the last
    # candidate in the region. If it is reached then we need to close the
    # partition.
    while (
        candidate_it < len(candidate_positions)
        and candidate_positions[candidate_it] != END_OF_REGION
        and interval.start <= candidate_positions[candidate_it] < interval.end
    ):
      if (
          num_of_candidates == max_size
          or partition_end - partition_start >= MAX_PARTITION_LEN
      ):
        # It may happen that there are no candidates over the a very long span
        # (For example HG0002 chr1:125,000,000 - 140,000,000). In that case we
        # need to break this interval into smaller partitions. Making allele
        # counters for the very long intervals exhausts the memory quickly.
        for pos in range(partition_start, partition_end, MAX_PARTITION_LEN):
          partitioned.append(
              ranges.make_range(
                  refname, pos, min(partition_end, pos + MAX_PARTITION_LEN)
              )
          )
        partition_start = partition_end
        partition_end = partition_start + 1
        num_of_candidates = 0
      else:
        partition_end = candidate_positions[candidate_it] + 1
        num_of_candidates += 1
      candidate_it += 1

    if (
        candidate_it < len(candidate_positions)
        and candidate_positions[candidate_it] == END_OF_REGION
    ):
      for pos in range(partition_start, interval.end, MAX_PARTITION_LEN):
        partitioned.append(
            ranges.make_range(
                refname, pos, min(interval.end, pos + MAX_PARTITION_LEN)
            )
        )
      candidate_it += 1
    else:
      raise ValueError('Terminating item is missing in candidates list')
  return partitioned


def regions_to_process(
    contigs: Sequence[reference_pb2.ContigInfo],
    partition_size: int,
    calling_regions: Optional[ranges.RangeSet] = None,
    task_id: Optional[int] = None,
    num_shards: Optional[int] = None,
    candidates: Optional[List[int]] = None,
    round_robin_sampling: bool = True,
) -> Iterable[range_pb2.Range]:
  """Determines the regions to process and partitions them into pieces.

  This function divides the genomes into regions we should process by
  intersecting the Ranges spanning all of the contigs with those from
  calling_regions, if provided. These intersected regions are then partitioned
  into pieces no bigger than partition_size bp in length.

  By construction we ensure that the regions are in genomic order, first w.r.t.
  the contigs and then within each contig by start and end of each region.

  This function can further subdivide these regions into a subset appropriate
  for a single task (task_id) among N tasks (num_shards) to process. The
  function ensures that:

    set(all_regions) = union(regions(task_0), ..., regions(task_n))

  when called with task_ids 0 ... N for num_shards = N.

  Args:
    contigs: Sequence of ContigInfo protos. Used to determine the initial ranges
      to process (i.e., all bases of these contigs) and the order of returned
      ranges.
    partition_size: The maximum size to make any region when partitioning.
    calling_regions: None or RangeSet. If provided, we will intersect the
      regions to process so that only those that overlap a region in this set
      are included.
    task_id: int >= 0 or None. The task_id of this job, which will be used to
      subdivide the total set of regions to process into just those that should
      be processed by this job. Must be < num_shards.
    num_shards: int >= 0 or None. The number of shards (i.e., the total number
      of tasks) we are running in parallel. Together with task_id determines the
      subset of regions we want to process.
    candidates: numpy array of int32 containing candidate positions. If
      candidate is provided then partition_by_candidates logic is used.
    round_robin_sampling: If true, sample regions are sampled among tasks in a
      round-robin fashion. This can help balance high-density regions among
      tasks, but results in candidates being output out of order.

  Returns:
    An iterable of nucleus.genomics.v1.Range objects.

  Raises:
    ValueError: if task_id and num_shards are bad or inconsistent.
  """
  if (task_id is None) != (num_shards is None):
    raise ValueError(
        'Both task_id and num_shards must be present if either is',
        task_id,
        num_shards,
    )
  if num_shards:
    if num_shards < 0:
      raise ValueError('num_shards={} must be >= 0'.format(num_shards))
    if task_id < 0 or task_id >= num_shards:
      raise ValueError(
          'task_id={} should be >= 0 and < num_shards={}'.format(
              task_id, num_shards
          )
      )

  regions = ranges.RangeSet.from_contigs(contigs)
  if calling_regions:
    regions = regions.intersection(calling_regions)

  # Depending on candidates parameter we choose the partitioning method.
  if candidates is not None:
    partitioned = list(partition_by_candidates(regions, candidates, 200))
  else:
    # Get number of partitions
    partitioned = list(regions.partition(partition_size))

  if num_shards:
    if round_robin_sampling:
      return (r for i, r in enumerate(partitioned) if i % num_shards == task_id)
    else:
      regions_per_shard = math.ceil(len(partitioned) / num_shards)
      return partitioned[
          task_id * regions_per_shard : (task_id + 1) * regions_per_shard
      ]
  else:
    return partitioned


def fetch_vcf_positions(
    vcf_paths: List[str],
    contigs: Sequence[reference_pb2.ContigInfo],
    calling_regions: Optional[ranges.RangeSet],
) -> List[range_pb2.Range]:
  """Fetches variants present in calling_regions.

  Args:
    vcf_paths: List of paths to VCFs from which to fetch positions.
    contigs: Sequence of ContigInfo protos. Used to determine the initial ranges
      to process (i.e., all bases of these contigs) and the order of returned
      ranges.
    calling_regions: A list of acceptable calling regions.

  Returns:
    Variant positions present in calling_regions.
  """
  # Fetch the set of regions being queried.
  regions = ranges.RangeSet.from_contigs(contigs)
  if calling_regions:
    regions = regions.intersection(calling_regions)

  variant_positions = []
  for vcf_path in vcf_paths:
    with vcf.VcfReader(vcf_path) as vcf_reader:
      for region in regions:
        for variant in vcf_reader.query(region):
          variant_positions.append(variant_utils.variant_position(variant))

  return variant_positions


def filter_regions_by_vcf(
    regions: List[range_pb2.Range], variant_positions: List[range_pb2.Range]
) -> List[range_pb2.Range]:
  """Filter a list of regions to only those that contain variants.

  Args:
    regions: a list of Range objects representing regions to filter on.
    variant_positions: a list of Range objects containing the positions of
      variants.

  Returns:
    filtered_regions: a list of Range objects, each of which appeared in the
        input regions and contains at least one of the input variants.
  """

  def dict_by_chromosome(
      list_of_ranges: List[range_pb2.Range],
  ) -> Dict[str, List[range_pb2.Range]]:
    d = collections.defaultdict(list)
    for r in list_of_ranges:
      d[r.reference_name].append(r)
    for c in d:
      d[c] = sorted(d[c], key=lambda x: (x.start, x.end))
    return d

  region_dict = dict_by_chromosome(regions)
  variant_dict = dict_by_chromosome(variant_positions)
  filtered_regions = []
  for c in region_dict:
    ri = 0
    vi = 0
    if c not in variant_dict:
      # Skip chromosomes with no variants.
      continue
    while ri < len(region_dict[c]) and vi < len(variant_dict[c]):
      region = region_dict[c][ri]
      variant = variant_dict[c][vi]
      if variant.start >= region.start and variant.start < region.end:
        # When the variant falls within the region, then keep the region.
        filtered_regions.append(region)
        # Move both indices because we're already keeping this region, and we
        # don't need to see any more variants inside this same region.
        ri += 1
        vi += 1
      elif region.start < variant.start:
        # Move past this region since the next variant comes later.
        ri += 1
      else:
        # Found another variant in the previous region we already included.
        vi += 1

  return filtered_regions


# ---------------------------------------------------------------------------
# Region processor
# ---------------------------------------------------------------------------


def read_confident_regions(
    options: deepvariant_pb2.MakeExamplesOptions,
    calling_regions: Optional[Sequence[range_pb2.Range]] = None,
) -> Optional[ranges.RangeSet]:
  """Reads in bed file of confident regions.

  Args:
    options: MakeExamplesOptions proto.
    calling_regions: calling regions to intersect with confident regions.

  Returns:
    List of ranges from confident region option or none if option is not set.
  """
  if options.confident_regions_filename:
    confident_regions = ranges.RangeSet.from_bed(
        options.confident_regions_filename,
        intersect_ranges=calling_regions,
    )
    return confident_regions
  else:
    return None


def read_denovo_regions(
    denovo_regions_filename: str,
) -> Optional[ranges.RangeSet]:
  """Read the bedfile provided in options and return a rangeset.

  Args:
    denovo_regions_filename: filename to read denovo regions from.

  Returns:
    List of ranges from denovo region option or none if option is not set.
  """
  if denovo_regions_filename:
    return ranges.RangeSet.from_bed(denovo_regions_filename)
  else:
    return None


def filter_candidates(
    candidates: Iterable[deepvariant_pb2.DeepVariantCall],
    select_variant_types: Sequence[str],
) -> Iterable[deepvariant_pb2.DeepVariantCall]:
  """Yields the candidate variants whose type is one of select_variant_types.

  This function iterates through candidates and yield each candidate in order
  if it satisfies any of the type constraints implied by select_variant_types.
  For example, if select_variant_types = ['snps'] this function will yield
  candidates that are bi-allelic SNPs only. Multiple select types are treated
  as OR'd together, so ['snps', 'indels'] yields candidates that are bi-allelic
  SNPs or indels.

  Args:
    candidates: Iterable of Variant protos. The candidates we want to select
      from.
    select_variant_types: List of str. The names of the variant type selectors
      we want to use to keep/remove variants. Each string must be part of
      VARIANT_TYPE_SELECTORS or an error will be raised.

  Raises:
    ValueError: if any str in select_variant_types isn't present in
      VARIANT_TYPE_SELECTORS.

  Yields:
    Candidates in order.
  """
  if not all(s in VARIANT_TYPE_SELECTORS for s in select_variant_types):
    raise ValueError('Unexpected select variant type', select_variant_types)

  for candidate in candidates:
    v = candidate.variant
    for select_type in select_variant_types:
      selector = VARIANT_TYPE_SELECTORS[select_type]
      if selector(v):
        yield candidate
        break


# ---------------------------------------------------------------------------
# A modified version of reservoir_sample for reads.
# ---------------------------------------------------------------------------


def reservoir_sample_reads(
    iterable_of_reads: Iterator[reads_pb2.Read],
    k: int,
    region: range_pb2.Range,
    max_bases_to_cover: int,
    random_generator: Optional[np.random.RandomState] = None,
) -> List[reads_pb2.Read]:
  """Samples k reads (or cover up to `max_bases_to_cover`) uniformly.

  Args:
    iterable_of_reads: The iterable to sample from.
    k: The number of elements to sample.
    region: The region we're sampling from. This can be used to determine how
      many bases are covered in the region.
    max_bases_to_cover: If this maximum number of bases is reached, the
      samplling will stop.
    random_generator: A random number generator or None.

  Returns:
    A list containing the sample reads.

  Raises:
    ValueError: If k is negative. Or, if k and max_bases_to_cover are both 0.
  """
  # If `max_bases_to_cover` is not set, use the simpler
  # reservoir_sample implementation.
  if not max_bases_to_cover:
    return utils.reservoir_sample(iterable_of_reads, k, random_generator)

  if k < 0:
    raise ValueError('k must be nonnegative, but got {}'.format(k))
  elif k == 0:
    # Because this function is now used both for selecting up to `k` or
    # covering `max_bases_to_cover`, if k is 0, we should set it to a large
    # number (meaning not limiting on that).
    k = float('inf')

  if random_generator is None:
    random_generator = np.random

  sampled_reads = []
  # Keep a list of the number of bases each `sampled_reads` have in the region.
  sampled_reads_overlap_len = []
  bases_covered = 0

  for i, read in enumerate(iterable_of_reads):
    if len(sampled_reads) < k and bases_covered < max_bases_to_cover:
      sampled_reads.append(read)
      overlap_len = ranges.overlap_len(region, utils.read_range(read))
      sampled_reads_overlap_len.append(overlap_len)
      bases_covered += overlap_len
    else:
      j = random_generator.randint(0, i + 1)
      if j < len(sampled_reads):
        # Because this replaces the read at sampled_reads[j], subtract first.
        bases_covered -= sampled_reads_overlap_len[j]
        sampled_reads[j] = read
        overlap_len = ranges.overlap_len(region, utils.read_range(read))
        sampled_reads_overlap_len[j] = overlap_len
        bases_covered += overlap_len

  # At the end, report cases where we covered max_bases_to_cover or more.
  if bases_covered >= max_bases_to_cover:
    # Empirically, bases_covered is likely much more than max_bases_to_cover at
    # this point. Let's do another round of trimming.
    total_bases = 0
    for i, overlap_len in enumerate(sampled_reads_overlap_len):
      total_bases += overlap_len
      if total_bases > max_bases_to_cover:
        sampled_reads = sampled_reads[: i + 1]
        bases_covered = total_bases
        break
    logging.info(
        (
            'In %s:%d-%d: reservoir_sample_reads sampled len(reads)=%s '
            'because bases_covered(%s) > max_bases_to_cover(%s).'
        ),
        region.reference_name,
        region.start,
        region.end,
        len(sampled_reads),
        bases_covered,
        max_bases_to_cover,
    )
  return sampled_reads


class DiagnosticLogger:
  """Writes diagnostic information about the assembler."""

  def __init__(
      self, output_root, normalized_reads_filename='normalized_reads.bam'
  ):
    self.normalized_reads_filename = normalized_reads_filename
    self.output_root = output_root

  def _root_join(self, path, makedirs=True):
    fullpath = os.path.join(self.output_root, path)
    subdir = os.path.dirname(fullpath)
    if makedirs and subdir:
      epath.Path(subdir).mkdir(parents=True, exist_ok=True)
    return fullpath

  def _file_for_region(self, region, basename):
    """Returns the path to a file in a region-specific subdirectory."""
    return self._root_join(os.path.join(ranges.to_literal(region), basename))

  def log_realigned_reads(self, region, reads, shared_header=None):
    """Logs, if enabled, the realigned reads for region."""
    path = self._file_for_region(region, self.normalized_reads_filename)
    logging.warning('writing %d normalized reads to %s', len(reads), path)
    with sam.SamWriter(path, header=shared_header) as writer:
      for read in reads:
        writer.write(read)


class OutputsWriter:
  """Manages all of the outputs of make_examples in a single place."""

  def __init__(self, options, suffix=None):
    outputs = [
        'candidates',
        'examples',
        'gvcfs',
        'runtime',
        'read_phases',
        'sitelist',
        'small_model_examples',
        'call_variant_outputs',
    ]
    self._writers = {k: None for k in outputs}
    self.examples_filename = None

    if options.candidates_filename:
      self._add_writer(
          'candidates',
          dv_utils.get_tf_record_writer(
              self._add_suffix(options.candidates_filename, suffix)
          ),
      )

    if options.examples_filename:
      clean_basename = re.sub(
          r'(\@[0-9]+|-\*?\d*-of-\*?\d*|\.gz)',
          '',
          os.path.basename(options.examples_filename.lower()),
      )
      if not re.fullmatch(r'.*(bagz|tfrecords?)$', clean_basename):
        raise ValueError(
            'Unsupported file extension: %s.\n'
            % os.path.basename(options.examples_filename)
        )
      self.examples_filename = self._add_suffix(
          options.examples_filename, suffix
      )

    if options.gvcf_filename:
      self._add_writer(
          'gvcfs',
          dv_utils.get_tf_record_writer(
              self._add_suffix(options.gvcf_filename, suffix)
          ),
      )

    if options.runtime_by_region:
      self._add_writer(
          'runtime', epath.Path(options.runtime_by_region).open('w')
      )
      writer = self._writers['runtime']
      if writer is not None:
        writer.__enter__()
        writer.write('\t'.join(RUNTIME_BY_REGION_COLUMNS) + '\n')

    if options.call_small_model_examples:
      self._add_writer(
          'call_variant_outputs',
          dv_utils.get_tf_record_writer(
              self._add_suffix(self.examples_filename, 'call_variant_outputs')
          ),
      )

    if options.read_phases_output:
      self._add_writer(
          'read_phases', epath.Path(options.read_phases_output).open('w')
      )
      writer = self._writers['read_phases']
      if writer is not None:
        writer.__enter__()
        writer.write('\t'.join(READ_PHASES_OUTPUT_COLUMNS) + '\n')

    if options.output_sitelist:
      sitelist_fname = options.examples_filename + '.sitelist.tsv'
      self._add_writer('sitelist', epath.Path(sitelist_fname).open('w'))

    if options.write_small_model_examples:
      self._add_writer(
          'small_model_examples',
          dv_utils.get_tf_record_writer(
              self._add_suffix(self.examples_filename, 'small_model')
          ),
      )

    self._deterministic_serialization = options.deterministic_serialization

  def _add_suffix(self, file_path, suffix):
    """Adds suffix to file name if a suffix is given."""
    if not suffix:
      return file_path

    file_dir, file_base = os.path.split(file_path)

    file_split = file_base.split('.')
    file_split[0] = f'{file_split[0]}_{suffix}'
    new_file_base = ('.').join(file_split)

    new_file = os.path.join(file_dir, new_file_base)
    return new_file

  def write_examples(self, *examples):
    self._write('examples', *examples)

  def write_gvcfs(self, *gvcfs):
    self._write('gvcfs', *gvcfs)

  def write_candidates(self, *candidates):
    self._write('candidates', *candidates)

  def write_call_variant_outputs(self, *call_variant_outputs):
    self._write('call_variant_outputs', *call_variant_outputs)

  def write_small_model_examples(self, *examples):
    self._write('small_model_examples', *examples)

  def write_site(
      self,
      call: variants_pb2.Variant,
      label=None,
  ):
    """Writes chrom,pos,ref,alt,label to a sitelist file."""
    chrom_pos_ref_alt = [
        call.reference_name,
        call.start,
        call.reference_bases,
        ','.join(call.alternate_bases),
    ]
    if label:
      label_class = label.features.feature['label'].int64_list.value[0]
      chrom_pos_ref_alt.append(label_class)
    else:
      chrom_pos_ref_alt.append(-1)
    site = '\t'.join(list(map(str, chrom_pos_ref_alt))) + '\n'

    self._write_text('sitelist', site)

  def write_runtime(self, stats_dict: Dict[str, Any]):
    columns = [str(stats_dict.get(k, 'NA')) for k in RUNTIME_BY_REGION_COLUMNS]
    writer = self._writers['runtime']
    if writer is None:
      raise ValueError('Runtime writer unexpectedly found to be None.')
    writer.write('\t'.join(columns) + '\n')

  def write_read_phase(self, read, phase, region_n):
    writer = self._writers['read_phases']
    if writer is not None:
      read_key = read.fragment_name + '/' + str(read.read_number)
      writer.write('\t'.join([read_key, str(phase), str(region_n)]) + '\n')

  def _add_writer(self, name: str, writer: tf_record.TFRecordWriter):
    if name not in self._writers:
      raise ValueError(
          'Expected writer {} to have a None binding in writers.'.format(name)
      )
    if self._writers[name] is not None:
      raise ValueError(
          'Expected writer {} to be bound to None in writers but '
          'saw {} instead'.format(name, self._writers[name])
      )
    self._writers[name] = writer

  def __enter__(self):
    """API function to support with syntax."""
    for writer in self._writers.values():
      if writer is not None:
        writer.__enter__()
    return self

  def __exit__(self, exception_type, exception_value, traceback):
    for writer in self._writers.values():
      if writer is not None:
        writer.__exit__(exception_type, exception_value, traceback)

  def _write(self, writer_name: str, *protos):
    writer = self._writers[writer_name]
    if writer:
      for proto in protos:
        writer.write(
            proto.SerializeToString(
                deterministic=self._deterministic_serialization
            )
        )

  def _write_text(self, writer_name: str, line: str):
    writer = self._writers[writer_name]
    if writer:
      writer.write(line)

  def close_all(self):
    for writer in self._writers.values():
      if writer is not None:
        writer.close()


class RegionProcessor:
  """Creates DeepVariant example protos for a single region on the genome.

  This class helps us to run the very sensitive caller, pileup image creator,
  and variant labeler operations on a single region in parallel across many
  regions using the PoolExecutor API. In order to do this we need three
  separate
  key operations:

  (1) Collect all of the info needed to create our resources (e.g., ref
  reader)
      at construction. We cannot actually initialize those resources in the
      constructor, though, since we actually want different resources in each
      worker process/thread. I.e., we need lazy resource initialization.

  (2) Actually initialize these resources *after* the worker has been forked
      in our process pool. This gives us a fresh resource to use in each
      separate process.

  (3) Process the region to find candidate variants and process those into our
      tf.Example protos.
  """

  def __init__(
      self,
      options: deepvariant_pb2.MakeExamplesOptions,
      calling_regions: Optional[Sequence[range_pb2.Range]] = None,
  ):
    """Creates a new RegionProcess.

    Args:
      options: deepvariant.MakeExamplesOptions proto used to specify our
        resources for calling (e.g., reference_filename).
      calling_regions: A list of ranges to call variants in.
    """
    self.options = options
    self.calling_regions = calling_regions
    self.samples = [
        sample_lib.Sample(options=x) for x in self.options.sample_options
    ]
    self.mean_coverage_per_sample = None
    self.initialized = False
    self._realigner = None
    self._labeler = None
    self._population_vcf_readers = None
    self._ref_reader = None
    self._make_examples_native = None
    self.region_number = 0
    if self.options.phase_reads:
      # One instance of DirectPhasing per lifetime of make_examples.
      self.direct_phasing_cpp = self._make_direct_phasing_obj()
    self.writers_dict = {}
    self.contig_dict = ranges.contigs_dict(
        fasta.IndexedFastaReader(self.options.reference_filename).header.contigs
    )
    self.small_model_example_factory = (
        make_small_model_examples.SmallModelExampleFactory(
            self.options.small_model_vaf_context_window_size,
            sample_names=[sample.options.name for sample in self.samples],
            accept_snps=self.options.small_model_snp_gq_threshold > -1,
            accept_indels=self.options.small_model_indel_gq_threshold > -1,
            accept_multiallelics=self.options.small_model_call_multiallelics,
            expand_by_haplotype=self.options.phase_reads,
        )
    )

  @property
  def realigner(self) -> realigner_module.Realigner:
    if self._realigner is None:
      raise ValueError('realigner is not initialized.')
    return self._realigner

  @realigner.setter
  def realigner(self, realigner: realigner_module.Realigner):
    self._realigner = realigner

  @property
  def labeler(
      self,
  ) -> (
      positional_labeler.PositionalVariantLabeler
      | haplotype_labeler.HaplotypeLabeler
  ):
    if self._labeler is None:
      raise ValueError('labeler is not initialized.')
    return self._labeler

  @labeler.setter
  def labeler(
      self,
      labeler: (
          positional_labeler.PositionalVariantLabeler
          | haplotype_labeler.HaplotypeLabeler
      ),
  ):
    self._labeler = labeler

  @property
  def population_vcf_readers(self) -> DefaultDict[str, Optional[vcf.VcfReader]]:
    if self._population_vcf_readers is None:
      raise ValueError('population_vcf_readers is not initialized.')
    return self._population_vcf_readers

  @population_vcf_readers.setter
  def population_vcf_readers(
      self, population_vcf_readers: DefaultDict[str, Optional[vcf.VcfReader]]
  ):
    self._population_vcf_readers = population_vcf_readers

  @property
  def ref_reader(self) -> fasta.IndexedFastaReader:
    if self._ref_reader is None:
      raise ValueError('ref_reader is not initialized.')
    return self._ref_reader

  @ref_reader.setter
  def ref_reader(self, ref_reader: fasta.IndexedFastaReader):
    self._ref_reader = ref_reader

  @property
  def make_examples_native(
      self,
  ) -> make_examples_native_module.ExamplesGenerator:
    if self._make_examples_native is None:
      raise ValueError('make_examples_native is not initialized.')
    return self._make_examples_native

  @make_examples_native.setter
  def make_examples_native(
      self,
      make_examples_native: make_examples_native_module.ExamplesGenerator,
  ):
    self._make_examples_native = make_examples_native

  def _make_direct_phasing_obj(self) -> direct_phasing.DirectPhasing:
    return direct_phasing.DirectPhasing()

  def _make_allele_counter_for_region(
      self, region: range_pb2.Range, candidate_positions: Iterable[int]
  ) -> allelecounter.AlleleCounter:
    return allelecounter.AlleleCounter(
        self.ref_reader.c_reader,
        region,
        candidate_positions,
        self.options.allele_counter_options,
    )

  def _make_allele_counter_for_read_overlap_region(
      self,
      region: range_pb2.Range,
      full_region: range_pb2.Range,
      candidate_positions: Iterable[int],
  ) -> allelecounter.AlleleCounter:
    return allelecounter.AlleleCounter.Default(
        self.ref_reader.c_reader,
        region,
        full_region,
        candidate_positions,
        self.options.allele_counter_options,
    )

  def _encode_tensor(
      self, image_tensor: np.ndarray
  ) -> Tuple[str, Tuple[int, int, int]]:
    return image_tensor.tostring(), image_tensor.shape

  def _make_sam_readers(
      self,
      reads_filenames: Sequence[str],
      downsample_fraction: float,
  ) -> Optional[List[sam.SamReader]]:
    """Creates a list of SamReaders, one from each filename.

    Args:
      reads_filenames: A list of string read filenames (e.g. for BAM/CRAM
        files). The list may contain empty strings or None, which will be
        skipped.
      downsample_fraction: Fraction by which to downsample. This applies to each
        file in reads_filenames separately.

    Returns:
      A list of sam readers with handles to the files. This may be shorter than
      the input reads_filenames if any of the filenames were empty.
    """
    logging_with_options(
        self.options,
        (
            'Starting from v0.9.0, --use_ref_for_cram is default to true. '
            'If you are using CRAM input, note that we will decode CRAM '
            'using the reference you passed in with --ref'
        ),
    )
    readers = []
    for reads_filename in reads_filenames:
      if reads_filename:
        readers.append(
            sam.SamReader(
                reads_filename,
                ref_name=self.options.ref_name_pangenome,
                context=self.options.allele_counter_options.partition_size,
                chrom_prefix=self.options.ref_chrom_prefix,
                shared_memory_name=self.options.gbz_shared_memory_name,
                create_shared_memory=False,
                use_loaded_shared_memory=self.options.use_loaded_gbz_shared_memory,
                ref_path=self.options.reference_filename
                if self.options.use_ref_for_cram
                else None,
                read_requirements=self.options.read_requirements,
                parse_aux_fields=self.options.parse_sam_aux_fields,
                aux_fields_to_keep=self.options.aux_fields_to_keep,
                hts_block_size=self.options.hts_block_size,
                downsample_fraction=downsample_fraction,
                random_seed=self.options.random_seed,
                use_original_base_quality_scores=self.options.use_original_quality_scores,
            )
        )
    return readers

  def _is_methylated_reference_site(
      self, candidate: deepvariant_pb2.DeepVariantCall
  ) -> bool:
    """Checks if the candidate is a methylated reference site.

    Any reference site that were included as candidate has been checked for
    methylation, hence it's sufficient to only check that it is a reference
    site.

    Args:
      candidate: A DeepVariantCall proto.

    Returns:
      True if the candidate is a methylated reference site.
    """
    return (
        len(candidate.variant.alternate_bases) == 1
        and candidate.variant.alternate_bases[0] == '.'
    )

  def _initialize(self):
    """Initialize the resources needed for this work in the current env."""
    if self.initialized:
      raise ValueError('Cannot initialize this object twice')

    self.ref_reader = fasta.IndexedFastaReader(self.options.reference_filename)

    for sample in self.samples:
      sample.sam_readers = self._make_sam_readers(
          reads_filenames=sample.options.reads_filenames,
          downsample_fraction=sample.options.downsample_fraction,
      )
      sample.in_memory_sam_reader = sam.InMemorySamReader([])
      sample.variant_caller = self._make_variant_caller_from_options(
          sample.options.variant_caller_options,
          sample.options.proposed_variants_filename,
      )
      if (
          self.options.call_small_model_examples
          and sample.options.small_model_path
      ):
        sample.small_model_variant_caller = (
            small_model_inference.SmallModelVariantCaller.from_model_path(
                model_path=sample.options.small_model_path,
                snp_gq_threshold=self.options.small_model_snp_gq_threshold,
                indel_gq_threshold=self.options.small_model_indel_gq_threshold,
                batch_size=self.options.small_model_inference_batch_size,
            )
        )

    if 'allele_frequency' in self.options.pic_options.channels:
      population_vcf_readers = allele_frequency.make_population_vcf_readers(
          self.options.population_vcf_filenames
      )
      self.population_vcf_readers = population_vcf_readers

    if self.options.exclude_variants_vcf_filename:
      self.exclude_variants_vcf_reader = vcf.VcfReader(
          self.options.exclude_variants_vcf_filename
      )
    initialize_realigner = (
        self.options.realigner_enabled
        or self.options.pic_options.alt_aligned_pileup != 'none'
        or self.options.allele_counter_options.track_ref_reads
    )

    if initialize_realigner:
      main_sample = self.samples[self.options.main_sample_index]
      input_bam_header = sam.SamReader(
          main_sample.options.reads_filenames[0]
      ).header
      self.realigner = realigner_module.Realigner(
          self.options.realigner_options,
          self.ref_reader,
          shared_header=input_bam_header,
      )

    self.options.allele_counter_options.enable_methylation_calling = (
        self.options.enable_methylation_calling
    )
    self.options.allele_counter_options.methylation_calling_threshold = (
        self.options.methylation_calling_threshold
    )
    self.options.allele_counter_options.enable_methylation_aware_phasing = (
        self.options.enable_methylation_aware_phasing
    )

    self._channels_enum = pileup_image_native.PileupImageEncoderNative(
        self.options.pic_options
    ).all_channels_enum(self.options.pic_options.alt_aligned_pileup)

    if in_training_mode(self.options):
      self.labeler = self._make_labeler_from_options()

    self.initialized = True

  def initialize(self):
    if not self.initialized:
      self._initialize()

  def _make_labeler_from_options(self):
    """Creates the labeler from options."""
    truth_vcf_reader = vcf.VcfReader(
        self.options.truth_variants_filename,
        excluded_format_fields=['GL', 'GQ', 'PL'],
    )
    confident_regions = read_confident_regions(
        self.options, self.calling_regions
    )

    if (
        self.options.variant_caller
        == deepvariant_pb2.MakeExamplesOptions.VCF_CANDIDATE_IMPORTER
    ):
      logging.info(
          'For --variant_caller=vcf_candidate_importer, we '
          'default the labeler_algorithm to positional_labler.'
      )
      return positional_labeler.PositionalVariantLabeler(
          truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions
      )

    if (
        self.options.labeler_algorithm
        == deepvariant_pb2.MakeExamplesOptions.POSITIONAL_LABELER
    ):
      return positional_labeler.PositionalVariantLabeler(
          truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions
      )
    elif (
        self.options.labeler_algorithm
        == deepvariant_pb2.MakeExamplesOptions.HAPLOTYPE_LABELER
    ):
      return haplotype_labeler.HaplotypeLabeler(
          truth_vcf_reader=truth_vcf_reader,
          ref_reader=self.ref_reader,
          confident_regions=confident_regions,
      )
    elif (
        self.options.labeler_algorithm
        == deepvariant_pb2.MakeExamplesOptions.CUSTOMIZED_CLASSES_LABELER
    ):
      if (
          not self.options.customized_classes_labeler_classes_list
          or not self.options.customized_classes_labeler_info_field_name
      ):
        raise ValueError(
            'For -labeler_algorithm=customized_classes_labeler, '
            'you need to set '
            '-customized_classes_labeler_classes_list and '
            '-customized_classes_labeler_info_field_name.'
        )
      return customized_classes_labeler.CustomizedClassesVariantLabeler(
          truth_vcf_reader=truth_vcf_reader,
          confident_regions=confident_regions,
          classes_list=self.options.customized_classes_labeler_classes_list,
          info_field_name=self.options.customized_classes_labeler_info_field_name,
      )
    else:
      raise ValueError(
          'Unexpected labeler_algorithm', self.options.labeler_algorithm
      )

  def _make_variant_caller_from_options(
      self,
      variant_caller_options: deepvariant_pb2.VariantCallerOptions,
      proposed_variants_filename: str,
  ) -> vc_base.VariantCaller:
    """Creates the variant_caller from options."""
    if (
        self.options.variant_caller
        == deepvariant_pb2.MakeExamplesOptions.VCF_CANDIDATE_IMPORTER
    ):
      if in_training_mode(self.options):
        candidates_vcf = self.options.truth_variants_filename
      else:
        candidates_vcf = proposed_variants_filename
      return vcf_candidate_importer.VcfCandidateImporter(
          variant_caller_options, candidates_vcf
      )
    elif (
        self.options.variant_caller
        == deepvariant_pb2.MakeExamplesOptions.VERY_SENSITIVE_CALLER
    ):
      return very_sensitive_caller.VerySensitiveCaller(variant_caller_options)
    else:
      raise ValueError('Unexpected variant_caller', self.options.variant_caller)

  def writes_examples_in_region(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      region: range_pb2.Range,
      sample_order: List[int],
      n_stats: Dict[str, int],
      runtimes: Dict[str, float],
      role: str,
      denovo_regions: Optional[ranges.RangeSet],
  ) -> Optional[List[int]]:
    """Generates and writes out the examples in a region.

    Args:
      candidates: List of candidates to be processed into examples.
      region: The region to generate examples.
      sample_order: Order of the samples to use when generating examples.
      n_stats: A dictionary that is used to accumulate counts for reporting.
      runtimes: A dictionary that recorded runtime information for reporting.
      role: The role that we make examples for.
      denovo_regions: The regions that contain denovo variants.

    Returns:
      example_shape: a list of 3 integers, representing the example shape in the
        region. If the region contains no examples, return None.
    """
    before_make_pileup_images = time.time()
    example_shape = None
    # Create A tf.Example proto, which includes the candidate variant, the
    # pileup image, and, if in training mode, the truth variants and labels
    # needed for training.
    if in_training_mode(self.options):
      # Initialize labels and types to be updated in the for loop below.
      reads_per_sample = []
      pileup_height = 0
      for sample in self.samples:
        reads_per_sample.append(sample.in_memory_sam_reader.iterate())
        pileup_height += sample.options.pileup_height
      # Unzip list of tuples.
      candidates_list = []
      for candidate, label in self.label_candidates(candidates, region):
        candidates_list.append(candidate)
        is_denovo = False
        if denovo_regions and denovo_regions.variant_overlaps(
            candidate.variant
        ):
          is_denovo = True
        # pylint: disable=unidiomatic-typecheck
        if type(label) is variant_labeler.VariantLabel:
          self.make_examples_native.append_label(
              make_examples_native_module.VariantLabel(
                  label.is_confident, label.variant, label.genotype, is_denovo
              )
          )
        elif (
            # pylint: disable=unidiomatic-typecheck
            type(label)
            is customized_classes_labeler.CustomizedClassesVariantLabel
        ):
          self.make_examples_native.append_label(
              make_examples_native_module.CustomizedClassesLabel(
                  label.is_confident,
                  label.variant,
                  label.truth_variant,
                  label.classes_dict,
                  label.info_field_name,
              )
          )
        else:
          raise ValueError('Unknown VariantLabel type: %s' % (type(label),))

      n_stats_one_region, example_shape_one = (
          self.make_examples_native.write_examples_in_region(
              candidates_list,
              reads_per_sample,
              sample_order,
              role,
              [0.0] * len(self.samples),
          )
      )

      for stat, val in n_stats_one_region.items():
        n_stats[stat] += val
        if stat == 'n_examples':
          if 'num examples' not in runtimes:
            runtimes['num examples'] = 0
          runtimes['num examples'] += val

      if example_shape is None and example_shape_one[0] > 0:
        example_shape = example_shape_one

    else:
      reads_per_sample = []
      pileup_height = 0
      for sample in self.samples:
        reads_per_sample.append(sample.in_memory_sam_reader.iterate())
        pileup_height += sample.options.pileup_height

      n_stats_one_region, example_shape_one = (
          self.make_examples_native.write_examples_in_region(
              candidates,
              reads_per_sample,
              sample_order,
              role,
              [0.0] * len(self.samples),
          )
      )

      for stat, val in n_stats_one_region.items():
        n_stats[stat] += val
        if stat == 'n_examples':
          if 'num examples' not in runtimes:
            runtimes['num examples'] = 0
          runtimes['num examples'] += val

      if example_shape is None and example_shape_one[0] > 0:
        example_shape = example_shape_one

    runtimes['make pileup images'] = trim_runtime(
        time.time() - before_make_pileup_images
    )
    return example_shape

  def write_small_model_examples_in_region(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      read_phases: Dict[str, int],
      sample: sample_lib.Sample,
      region: range_pb2.Range,
      writer: OutputsWriter,
      n_stats: Dict[str, int],
      runtimes: Dict[str, float],
  ) -> None:
    """Writes out the small model training examples in a region.

    Args:
      candidates: List of candidates to be processed into examples.
      read_phases: A dictionary of read names to haplotype phases.
      sample: The sample for which to generate small model examples.
      region: The region to generate examples.
      writer: A OutputsWriter used to write out examples.
      n_stats: A dictionary that is used to accumulate counts for reporting.
      runtimes: A dictionary that recorded runtime information for reporting.
    """
    before_make_summaries = time.time()
    if not in_training_mode(self.options):
      raise ValueError(
          'Writing small model examples is only supported in training mode.'
      )
    training_examples = (
        self.small_model_example_factory.encode_training_examples(
            list(self.label_candidates(candidates, region)),
            read_phases,
            sample.options.order,
        )
    )
    writer.write_small_model_examples(*training_examples)

    n_stats['n_small_model_examples'] += len(training_examples)
    runtimes['make small_model_examples'] = trim_runtime(
        time.time() - before_make_summaries
    )

  def call_small_model_examples_in_region(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      read_phases: Dict[str, int],
      sample: sample_lib.Sample,
      writer: OutputsWriter,
      n_stats: Dict[str, int],
      runtimes: Dict[str, float],
  ) -> Sequence[deepvariant_pb2.DeepVariantCall]:
    """Creates and calls small model examples on candidates in a region.

    Args:
      candidates: List of candidates to be processed into examples.
      read_phases: A dictionary of read names to haplotype phases.
      sample: The sample for which to call small model examples.
      writer: A OutputsWriter used to write out examples.
      n_stats: A dictionary that is used to accumulate counts for reporting.
      runtimes: A dictionary that recorded runtime information for reporting.

    Returns:
      A list of all candidates to be passed to the regular model, either
        because they were skipped or did not pass quality filters.
    """
    before_generate_small_model_examples = time.time()

    inference_example_set = (
        self.small_model_example_factory.encode_inference_examples(
            candidates, read_phases, sample.options.order
        )
    )
    runtimes['small model generate examples'] = trim_runtime(
        time.time() - before_generate_small_model_examples
    )
    if not inference_example_set.candidates_with_alt_allele_indices:
      return inference_example_set.skipped_candidates

    # filtered candidates did not pass the GQ threshold.
    before_call_small_model_examples = time.time()
    assert sample.small_model_variant_caller is not None
    call_variant_outputs, candidates_not_called = (
        sample.small_model_variant_caller.call_variants(
            inference_example_set.candidates_with_alt_allele_indices,
            inference_example_set.inference_examples,
        )
    )
    runtimes['small model call examples'] = trim_runtime(
        time.time() - before_call_small_model_examples
    )
    before_write_variants = time.time()
    writer.write_call_variant_outputs(*call_variant_outputs)
    runtimes['small model write variants'] = trim_runtime(
        time.time() - before_write_variants
    )

    n_stats['n_small_model_calls'] += len(call_variant_outputs)
    runtimes['small model total'] = trim_runtime(
        time.time() - before_generate_small_model_examples
    )
    # pass skipped and filtered candidates to the large model
    inference_example_set.skipped_candidates.extend(candidates_not_called)
    return inference_example_set.skipped_candidates

  def find_candidate_positions(self, region: range_pb2.Range) -> Iterator[int]:
    """Finds all candidate positions within a given region."""
    main_sample = self.samples[self.options.main_sample_index]
    for sample in self.samples:
      reads = itertools.chain()
      for _, sam_reader in enumerate(sample.sam_readers):
        reads = itertools.chain(reads, sam_reader.query(region))
      try:
        sample.in_memory_sam_reader.replace_reads(reads)
        sample.reads = sample.in_memory_sam_reader.query(region)
        max_bases_to_cover = 0
        if self.options.max_reads_for_dynamic_bases_per_region > 0:
          max_bases_to_cover = (
              self.options.max_reads_for_dynamic_bases_per_region
              * (region.end - region.start)
          )
        if self.options.max_reads_per_partition > 0 or max_bases_to_cover > 0:
          random_for_region = np.random.RandomState(self.options.random_seed)
          sample.reads = reservoir_sample_reads(
              sample.reads,
              self.options.max_reads_per_partition,
              region,
              max_bases_to_cover,
              random_for_region,
          )

        sample.allele_counter = self._make_allele_counter_for_region(region, [])

        if sample.options.reads_filenames:
          for read in sample.reads:
            sample.allele_counter.add(read, sample.options.name)
      except ValueError as err:
        error_message = str(err)
        if error_message.startswith('DATA_LOSS:'):
          raise ValueError(
              error_message
              + '\nFailed to parse BAM/CRAM file. '
              'This is often caused by:\n'
              '(1) When using a CRAM file, and setting '
              '--use_ref_for_cram to false (which means you want '
              'to use the embedded ref instead of a ref file), '
              'this error could be because of inability to find '
              'the embedded ref file.\n'
              '(2) Your BAM/CRAM file could be corrupted. Please '
              'check its md5.\n'
              'If you cannot find out the reason why this error '
              'is occurring, please report to '
              'https://github.com/google/deepvariant/issues'
          ) from err
        elif error_message.startswith('NOT_FOUND: Unknown reference_name '):
          raise ValueError(
              '{}\nThe region {} does not exist in {}.'.format(
                  error_message,
                  ranges.to_literal(region),
                  sample.options.reads_filenames,
              )
          ) from err
        else:
          # By default, raise the ValueError as is for now.
          raise err

    # end of self.samples loop:

    allele_counters = {s.options.name: s.allele_counter for s in self.samples}
    # TODO: For phasing we calculate candidates for all samples.
    # If it is done here then we can reuse these results for phasing thus
    # saving runtime.
    candidate_positions = main_sample.variant_caller.get_candidate_positions(
        allele_counters=allele_counters, sample_name=main_sample.options.name
    )
    for pos in candidate_positions:
      yield pos
    # Mark the end of partition
    yield END_OF_PARTITION

  def _only_contains_alts_above_threshold(
      self, variant: variants_pb2.Variant, threshold: float
  ) -> bool:
    """Returns True if the variant only contains alts above the threshold."""
    dict_allele_frequency = allele_frequency.find_matching_allele_frequency(
        variant=variant,
        population_vcf_reader=self.exclude_variants_vcf_reader,
        ref_reader=self.ref_reader,
    )
    alt_afs = [dict_allele_frequency[alt] for alt in variant.alternate_bases]
    if not alt_afs:
      return True
    elif min(alt_afs) >= threshold:
      return True
    return False

  def process(
      self, region: range_pb2.Range, region_n: Optional[int] = None
  ) -> Tuple[
      Dict[str, Sequence[deepvariant_pb2.DeepVariantCall]],
      Dict[str, Sequence[variants_pb2.Variant]],
      Dict[str, Union[float, int]],
      Dict[str, Dict[str, int]],
      int,
  ]:
    """Finds candidates and creates corresponding examples in a region.

    Args:
      region: A nucleus.genomics.v1.Range proto. Specifies the region on the
        genome we should process.
      region_n: Order number of the region being processed by this process.

    Returns:
      (candidates_by_sample, gvcfs_by_sample, runtimes)
      1. candidates_by_sample: A dict keyed by sample role, each a list of
      candidates found, which are deepvariant.DeepVariantCall objects.
      2. gvcfs_by_sample: A dict keyed by sample, each a list of
      nucleus.genomics.v1.Variant protos containing gVCF information for all
      reference sites, if gvcf generation is enabled, otherwise this value is
      [].
      3. runtimes: A dict of runtimes in seconds keyed by stage.
      4. read_phases_by_sample: A dict keyed by sample role, each a dict of
      read_name & read_number to read_phase.
      5. phased_reads_count: The number of reads that were phased.
    """
    runtimes = {}

    if not self.initialized:
      self.initialize()

    # Keep track the number of regions processed.
    self.region_number += 1
    before_get_reads = time.time()
    runtimes['num reads'] = 0
    # Collect reads from multiple BAMs. Each BAM contains a sample.
    sample_reads_list = []
    for sample in self.samples:
      if sample.in_memory_sam_reader is not None:
        # Realigner is called outside region_reads_norealign()
        sample_reads = self.region_reads_norealign(
            region=region,
            sam_readers=sample.sam_readers,
            reads_filenames=sample.options.reads_filenames,
        )
        runtimes['num reads'] += len(sample_reads)
        sample_reads_list.append(sample_reads)
      else:
        sample_reads_list.append([])
    if self.options.joint_realignment:
      sample_reads_list = self.realign_reads_joint_multisample(
          sample_reads_list, region
      )
    else:
      sample_reads_list = self.realign_reads_per_sample_multisample(
          sample_reads_list, region
      )
    for sample_index, sample in enumerate(self.samples):
      sample.in_memory_sam_reader.replace_reads(sample_reads_list[sample_index])

    runtimes['get reads'] = trim_runtime(time.time() - before_get_reads)
    before_find_candidates = time.time()

    # Region is expanded by region_padding number of bases. This functionality
    # is only needed when phase_reads flag is on.
    region_padding_percent = self.options.phase_reads_region_padding_pct
    if self.options.phase_reads and region_padding_percent > 0:
      # When candidate partitioning is used region size is variable. Therefore
      # we need to calculate the padding for each region.
      padding_fraction = int(
          (region.end - region.start) * region_padding_percent / 100
      )
      region_expanded = ranges.expand(
          region, padding_fraction, self.contig_dict
      )

      (
          candidates_by_sample,
          gvcfs_by_sample,
          read_phases_by_sample,
          phased_reads_count,
      ) = self.candidates_in_region(
          region=region, region_n=region_n, padded_region=region_expanded
      )

    else:
      (
          candidates_by_sample,
          gvcfs_by_sample,
          read_phases_by_sample,
          phased_reads_count,
      ) = self.candidates_in_region(region=region, region_n=region_n)

    for sample in self.samples:
      role = sample.options.role
      if sample.options.skip_output_generation:
        continue
      if role not in candidates_by_sample:
        continue
      candidates = candidates_by_sample[role]

      if self.options.select_variant_types:
        candidates = list(
            filter_candidates(candidates, self.options.select_variant_types)
        )

      if (
          hasattr(self, 'exclude_variants_vcf_reader')
          and self.exclude_variants_vcf_reader is not None
      ):
        candidates = [
            candidate
            for candidate in candidates
            if not self._only_contains_alts_above_threshold(
                candidate.variant, self.options.exclude_variants_af_threshold
            )
        ]

      runtimes['find candidates'] = trim_runtime(
          time.time() - before_find_candidates
      )
      before_make_pileup_images = time.time()

      # Get allele frequencies for candidates.
      if 'allele_frequency' in self.options.pic_options.channels:
        candidates = list(
            allele_frequency.add_allele_frequencies_to_candidates(
                candidates=candidates,
                population_vcf_reader=self.population_vcf_readers[
                    region.reference_name
                ],
                ref_reader=self.ref_reader,
            )
        )

      # After any filtering and other changes above, set candidates for sample.
      candidates_by_sample[role] = candidates

      runtimes['make pileup images'] = trim_runtime(
          time.time() - before_make_pileup_images
      )
    runtimes['num candidates'] = sum(
        [len(x) for x in candidates_by_sample.values()]
    )
    return (
        candidates_by_sample,
        gvcfs_by_sample,
        runtimes,
        read_phases_by_sample,
        phased_reads_count,
    )

  def region_reads_norealign(
      self,
      region: range_pb2.Range,
      sam_readers: Optional[Sequence[sam.SamReader]],
      reads_filenames: Optional[Sequence[str]],
  ) -> List[reads_pb2.Read]:
    """Gets reads overlapping the region.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to query reads.
      sam_readers: An iterable of sam.SamReader to query from.
      reads_filenames: Filenames matching sam_readers. This is only used for
        throwing more informative error messages.

    Returns:
      [genomics.deepvariant.core.genomics.Read], reads overlapping the region.
    """
    if sam_readers is None:
      return []

    # reads = itertools.chain([reader.query(region) for reader in sam_readers])
    reads = itertools.chain()
    for sam_reader in sam_readers:
      reads = itertools.chain(reads, sam_reader.query(region))

    try:
      max_bases_to_cover = 0
      if self.options.max_reads_for_dynamic_bases_per_region > 0:
        max_bases_to_cover = (
            self.options.max_reads_for_dynamic_bases_per_region
            * (region.end - region.start)
        )
      if self.options.max_reads_per_partition > 0 or max_bases_to_cover > 0:
        random_for_region = np.random.RandomState(self.options.random_seed)
        reads = reservoir_sample_reads(
            reads,
            self.options.max_reads_per_partition,
            region,
            max_bases_to_cover,
            random_for_region,
        )
      return list(reads)
    except ValueError as err:
      error_message = str(err)
      if error_message.startswith('DATA_LOSS:'):
        raise ValueError(
            error_message
            + '\nFailed to parse BAM/CRAM file. '
            'This is often caused by:\n'
            '(1) When using a CRAM file, and setting '
            '--use_ref_for_cram to false (which means you want '
            'to use the embedded ref instead of a ref file), '
            'this error could be because of inability to find '
            'the embedded ref file.\n'
            '(2) Your BAM/CRAM file could be corrupted. Please '
            'check its md5.\n'
            'If you cannot find out the reason why this error '
            'is occurring, please report to '
            'https://github.com/google/deepvariant/issues'
        ) from err
      elif error_message.startswith('NOT_FOUND: Unknown reference_name '):
        raise ValueError(
            '{}\nThe region {} does not exist in {}.'.format(
                error_message, ranges.to_literal(region), reads_filenames
            )
        ) from err
      else:
        # By default, raise the ValueError as is for now.
        raise err

  def realign_reads(
      self, reads: List[reads_pb2.Read], region: range_pb2.Range
  ) -> List[reads_pb2.Read]:
    """Realign reads overlapping the region.

    Args:
      reads: list of reads.
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to realign reads.

    Returns:
      genomics.deepvariant.core.genomics.Read: realigned reads
    """
    if self.options.realigner_enabled:
      max_read_length_to_realign = 500
      if max_read_length_to_realign > 0:
        long_reads = [
            read
            for read in reads
            if len(read.aligned_sequence) > max_read_length_to_realign
        ]

        short_reads = [
            read
            for read in reads
            if len(read.aligned_sequence) <= max_read_length_to_realign
        ]

        _, realigned_short_reads = self.realigner.realign_reads(
            short_reads, region
        )

        # Long reads will be listed before short reads when both are present.
        # Examples with only short or only long reads will be unaffected.
        return long_reads + realigned_short_reads

      _, reads = self.realigner.realign_reads(reads, region)
    return reads

  def realign_reads_per_sample_multisample(
      self,
      sample_reads_list: List[List[reads_pb2.Read]],
      region: range_pb2.Range,
  ) -> List[List[reads_pb2.Read]]:
    """Realign reads overlapping the region.

    Args:
      sample_reads_list: list of reads-list per sample.
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to realign reads.

    Returns:
      [genomics.deepvariant.core.genomics.Read], realigned reads per sample
    """
    return [
        self.realign_reads(reads_per_sample, region)
        for reads_per_sample in sample_reads_list
    ]

  def realign_reads_joint_multisample(
      self,
      sample_reads_list: List[List[reads_pb2.Read]],
      region: range_pb2.Range,
  ) -> List[List[reads_pb2.Read]]:
    """Realign reads overlapping the region.

    Args:
      sample_reads_list: list of reads-list per sample.
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to realign reads.

    Returns:
    [genomics.deepvariant.core.genomics.Read], realigned reads per sample
    """
    # join reads from all samples
    if len(sample_reads_list) > 1:
      reads = []
      for sample_index, sample_reads in enumerate(sample_reads_list):
        for read in sample_reads:
          read.fragment_name += f'.{sample_index}'
        reads.extend(sample_reads)
    else:
      reads = sample_reads_list[0]

    realigned_reads = self.realign_reads(reads, region)

    sample_realigned_reads_list = [[] for _ in sample_reads_list]

    # demultiplex reads
    if len(sample_reads_list) > 1:
      for read in realigned_reads:
        read.fragment_name, sample_index = read.fragment_name.rsplit('.', 1)
        sample_index = int(sample_index)
        sample_realigned_reads_list[sample_index].append(read)
    else:
      sample_realigned_reads_list = [realigned_reads]
    return sample_realigned_reads_list

  def filter_candidates_by_region(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      region: range_pb2.Range,
  ) -> Sequence[deepvariant_pb2.DeepVariantCall]:
    return [
        candidate
        for candidate in candidates
        if candidate.variant.start >= region.start
        and candidate.variant.start < region.end
    ]

  def _root_join(self, path, makedirs=True):
    fullpath = os.path.join(
        self.options.realigner_options.diagnostics.output_root, path
    )
    subdir = os.path.dirname(fullpath)
    if makedirs and subdir:
      epath.Path(subdir).mkdir(parents=True, exist_ok=True)
    return fullpath

  def _file_for_region(self, region, basename):
    """Returns the path to a file in a region-specific subdirectory."""
    # TODO: This logic currently only works for single sample.
    # Once we extend to multi-sample, we can remove this assert.
    assert len(self.samples) == 1
    return self._root_join(os.path.join(ranges.to_literal(region), basename))

  def log_graph_metrics(self, region, graph):
    """Logs, if enabled, graph construction information for region."""
    if graph:
      dest_file = self._file_for_region(region, 'graph.dot')
      with epath.Path(dest_file).open('w') as f:
        f.write(graph.graphviz())

  def add_phasing_to_candidate(self, candidates, read_id_to_phase):
    """Adds phasing information to candidates.

    Args:
      candidates: A list of DeepVariantCall protos where is written to.
      read_id_to_phase: A dict of read key to phase.

    This function populates the ALT_PS and PS_CONTIG info fields in the
    candidate variant protos. The ALT_PS field is a phased genotype.
    The PS_CONTIG field is a string that identifies the continious contig
    within which the phasing is consistent.

    Returns:
      The number of phased variants.
    """
    # Calling into direct_phasing_cpp to get phased alleles.
    # direct_phasing_cpp preserves the state until the next call to
    # phase().
    phased_variants = self.direct_phasing_cpp.get_phased_variants()
    # phased_variants contains a subset of candidates because not all
    # candidates are phased. We iterate over candidates and
    # phased_variants simultaneously. If candidate is phased, we add phase
    # information to it. If candidate is not phased, we infer phase from
    # supporting reads.
    phased_variants_index = 0
    phase_contig = f'{self.options.task_id}-{self.region_number}'
    for candidate in candidates:
      # if phased_variants_index >= len(phased_variants):
      #   break
      # If variant is phased we can use phasing information from
      # phased_variants.
      if (
          phased_variants_index < len(phased_variants)
          and candidate.variant.start
          == phased_variants[phased_variants_index].position
      ):
        alt_alleles = list(candidate.variant.alternate_bases).copy()
        alt_alleles.insert(0, 'REF')
        phased_genotype = [0] * (len(alt_alleles))
        # phased_variants contains alleles ordered by phase.
        alt_1_order = [
            i
            for i, e in enumerate(alt_alleles)
            if e == phased_variants[phased_variants_index].phase_1_bases
        ]
        alt_2_order = [
            i
            for i, e in enumerate(alt_alleles)
            if e == phased_variants[phased_variants_index].phase_2_bases
        ]
        if alt_1_order and alt_2_order:
          phased_genotype[alt_1_order[0]] = 1
          phased_genotype[alt_2_order[0]] = 2
          variant_utils.set_info(candidate.variant, 'ALT_PS', phased_genotype)
          variant_utils.set_info(candidate.variant, 'PS_CONTIG', phase_contig)
        phased_variants_index += 1
      # If variant is not phased, we infer phase from read phases and supporting
      # reads.
      else:
        # Count the number of reads of each phase supporting the allele.
        # Assign phase to the allele that is supported by more reads. The
        # difference between the number of reads supporting each phase must be
        # at least 2.
        phased_genotype = [0] * (len(candidate.variant.alternate_bases) + 1)
        index = 0
        for allele_index in range(len(candidate.variant.alternate_bases) + 1):
          phases = [0, 0, 0]  # number of reads supporting each phase.
          # allele_index == 0 is REF
          if allele_index == 0:
            # Get REF supporting reads.
            for ref_read_support in candidate.ref_support_ext.read_infos:
              # In multi-sample mode we have reads that support an allele but
              # come from a non-target sample. We ignore these reads.
              if ref_read_support.read_name in read_id_to_phase:
                phases[read_id_to_phase[ref_read_support.read_name]] += 1
          else:
            # Find allele bases, then get supporting reads.
            if (
                allele_index > len(candidate.variant.alternate_bases)
                or allele_index < 1
            ):
              continue
            alt_bases = candidate.variant.alternate_bases[allele_index - 1]
            for read_support in candidate.allele_support_ext[
                alt_bases
            ].read_infos:
              if read_support.read_name in read_id_to_phase:
                phases[read_id_to_phase[read_support.read_name]] += 1

          if phases[1] > phases[2] and phases[1] - phases[2] > 1:
            phased_genotype[index] = 1
          elif phases[2] > phases[1] and phases[2] - phases[1] > 1:
            phased_genotype[index] = 2
          index += 1
        variant_utils.set_info(candidate.variant, 'ALT_PS', phased_genotype)
        variant_utils.set_info(candidate.variant, 'PS_CONTIG', phase_contig)
    return len(phased_variants)

  def candidates_in_region(
      self,
      region: range_pb2.Range,
      region_n: Optional[int] = None,
      padded_region: Optional[range_pb2.Range] = None,
  ) -> Tuple[
      Dict[str, Sequence[deepvariant_pb2.DeepVariantCall]],
      Dict[str, Sequence[variants_pb2.Variant]],
      Dict[str, Dict[str, int]],
      int,
  ]:
    """Finds candidates in the region using the designated variant caller.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to get candidates for.
      region_n: Order number of the region being processed by this process.
      padded_region: A nucleus.genomics.v1.Range object specifying the padded
        region.

    Returns:
      A 3-tuple of (candidates, gvcfs, read_phases).
      The first value, candidates, is a dict keyed by sample role, where each
      item is a list of deepvariant_pb2.DeepVariantCalls objects, in
      coordidate order.
      The second value, gvcfs, is a dict keyed by sample role, where
      each item is a list of nucleus.genomics.v1.Variant protos containing gVCF
      information for all reference sites, if gvcf generation is enabled,
      otherwise the gvcfs value is [].
      The third value, read_phases, is a dict keyed by sample role, where each
      item is a dict keyed by the read name and number, where the value is the
      phase/HP tag of the
      read.
      The fourth value, phased_reads_count, is the number of reads that were
      phased.
    """
    for sample in self.samples:
      sample.reads = sample.in_memory_sam_reader.query(region)

    main_sample = self.samples[self.options.main_sample_index]
    if not main_sample.reads and not gvcf_output_enabled(self.options):
      # If we are generating gVCF output we cannot safely abort early here as
      # we need to return the gVCF records calculated by the caller below.
      return {}, {}, {}, 0

    allele_counters = {}
    candidate_positions = []
    if self.options.allele_counter_options.track_ref_reads:
      # Calculate potential candidate positions from allele counts.
      for sample in self.samples:
        if sample.options.reads_filenames:
          # Calculate potential candidate positions from allele counts
          if padded_region is not None:
            sample.allele_counter = self._make_allele_counter_for_region(
                padded_region, []
            )
          else:
            sample.allele_counter = self._make_allele_counter_for_region(
                region, []
            )

          for read in sample.reads:
            sample.allele_counter.add(read, sample.options.name)
        # Reads iterator needs to be reset since it used in the code below.
        sample.reads = sample.in_memory_sam_reader.query(region)
      allele_counters = {s.options.name: s.allele_counter for s in self.samples}

    for sample in self.samples:
      if self.options.allele_counter_options.track_ref_reads:
        candidate_positions = sample.variant_caller.get_candidate_positions(
            allele_counters=allele_counters, sample_name=sample.options.name
        )
      if sample.options.reads_filenames:
        if (
            self.options.allele_counter_options.normalize_reads
            and not sample.options.skip_normalization
        ):
          reads_start = region.start
          reads_end = region.end
          for read in sample.reads:
            read_last_pos = min(
                self.ref_reader.contig(region.reference_name).n_bases - 1,
                utils.read_end(read),
            )
            if read.alignment.position.position < reads_start:
              reads_start = read.alignment.position.position
            if read_last_pos > reads_end:
              reads_end = read_last_pos
          full_range = range_pb2.Range(
              reference_name=region.reference_name,
              start=reads_start,
              end=reads_end,
          )
          sample.reads = sample.in_memory_sam_reader.query(region)

          if padded_region is not None:
            sample.allele_counter = (
                self._make_allele_counter_for_read_overlap_region(
                    padded_region, full_range, candidate_positions
                )
            )
          else:
            sample.allele_counter = (
                self._make_allele_counter_for_read_overlap_region(
                    region, full_range, candidate_positions
                )
            )
        else:
          if padded_region is not None:
            sample.allele_counter = self._make_allele_counter_for_region(
                padded_region, candidate_positions
            )
          else:
            sample.allele_counter = self._make_allele_counter_for_region(
                region, candidate_positions
            )

        for read in sample.reads:
          if (
              self.options.allele_counter_options.normalize_reads
              and not sample.options.skip_normalization
          ):
            cigar_list, read_shift = sample.allele_counter.normalize_and_add(
                read, sample.options.name
            )
            if cigar_list:
              if read_shift != 0:
                read.alignment.position.position += read_shift
              del read.alignment.cigar[:]
              for el in cigar_list:
                read.alignment.cigar.add(
                    operation=el.operation, operation_length=el.operation_length
                )
          else:
            sample.allele_counter.add(read, sample.options.name)

        allele_counters[sample.options.name] = sample.allele_counter

    candidates = {}
    gvcfs = {}
    read_phases_by_sample = {}
    phased_reads_count = 0
    left_padding = 0
    right_padding = 0
    if padded_region is not None:
      left_padding = region.start - padded_region.start
      right_padding = padded_region.end - region.end
    for sample in self.samples:
      role = sample.options.role
      writer = None
      if role in self.writers_dict:
        writer = self.writers_dict[role]
      if not sample.options.reads_filenames:
        continue
      read_phases_by_sample[role] = {}
      candidates[role], gvcfs[role] = sample.variant_caller.calls_and_gvcfs(
          allele_counters=allele_counters,
          target_sample=sample.options.name,
          include_gvcfs=gvcf_output_enabled(self.options),
          include_med_dp=self.options.include_med_dp,
          left_padding=left_padding,
          right_padding=right_padding,
      )

      # If methylation-aware phasing is enabled, filter for methylated reference
      # sites and SNP candidates.
      # SNP candidates will be phased using direct phasing.
      # Methylated reference sites will be phased using methylation-aware
      # phasing after direct phasing.
      snp_candidates = []
      methylated_ref_sites = []
      if self.options.enable_methylation_aware_phasing:
        for candidate in candidates[role]:
          if self._is_methylated_reference_site(candidate):
            methylated_ref_sites.append(candidate)
          else:
            snp_candidates.append(candidate)

        # Only use SNP candidates for direct phasing
        candidates[role] = snp_candidates

      if self.options.phase_reads and not sample.options.skip_phasing:
        if padded_region is not None:
          reads_to_phase = list(
              sample.in_memory_sam_reader.query(padded_region)
          )
        else:
          reads_to_phase = list(sample.in_memory_sam_reader.query(region))
        for read in reads_to_phase:
          # Remove existing values
          del read.info['HP'].values[:]
        # Skip phasing if number of candidates is over the phase_max_candidates.
        if (
            self.options.phase_max_candidates
            and len(candidates[role]) > self.options.phase_max_candidates
        ):
          logging_with_options(
              self.options,
              'Skip phasing: len(candidates[%s]) is %s.'
              % (role, len(candidates[role])),
          )
        else:
          read_phases = self.direct_phasing_cpp.phase(
              candidates[role], reads_to_phase
          )

          # If methylation-aware phasing is enabled, run it on unphased reads.
          if (
              self.options.enable_methylation_aware_phasing
              and methylated_ref_sites
          ):
            # Perform methylation-aware phasing on unphased reads from direct
            # phasing.
            read_phases = methylation_aware_phasing.phase(
                reads_to_phase, read_phases, methylated_ref_sites
            )

          # Assign phase tag to reads after direct phasing and/or
          # methylation-aware phasing.
          read_id_to_phase = {}
          for read_phase, read in zip(read_phases, reads_to_phase):
            # Remove existing values
            del read.info['HP'].values[:]
            read_key = read.fragment_name + '/' + str(read.read_number)
            read_id_to_phase[read_key] = read_phase

            if self.options.pic_options.reverse_haplotypes:
              if read_phase in [1, 2]:
                read_phase = 1 + (read_phase % 2)
            read.info['HP'].values.add(int_value=read_phase)
            read_phases_by_sample[role][read_key] = read_phase
            if writer and self.options.read_phases_output:
              writer.write_read_phase(read, read_phase, region_n)
          if self.options.output_phase_info:
            phased_reads_count = self.add_phasing_to_candidate(
                candidates[role], read_id_to_phase
            )
          # This logic below will write out the DOT files under the directory
          # specified by the flag --realigner_diagnostics, if phase_reads is
          # set to True.
          # TODO: Extend the logic to work for multi-sample cases.
      if (
          self.options.phase_reads
          and self.options.realigner_options.diagnostics.output_root
          and len(self.samples) == 1  # TODO
      ):
        self.log_graph_metrics(region, self.direct_phasing_cpp)

      if padded_region is not None:
        candidates[role] = self.filter_candidates_by_region(
            candidates[role], region
        )

    return candidates, gvcfs, read_phases_by_sample, phased_reads_count

  def get_channels(self) -> List[int]:
    # All the example would have the same list of channels.
    return self._channels_enum

  def label_candidates(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      region: range_pb2.Range,
  ) -> Iterator[
      Tuple[deepvariant_pb2.DeepVariantCall, variant_labeler.VariantLabel]
  ]:
    """Gets label information for each candidate.

    Args:
      candidates: list[DeepVariantCalls]: The list of candidate variant calls we
        want to label.
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to get candidates for.

    Yields:
      Tuples of (candidate, label_variants.Label objects) for each candidate in
      candidates that could be assigned a label. Candidates that couldn't be
      labeled will not be returned.
    """
    # Set BAM filename (used for training stats).
    for candidate in candidates:
      struct_utils.set_string_field(
          candidate.variant.info, 'BAM_FNAME', self.options.bam_fname
      )

    # Get our list of labels for each candidate variant.
    labels = self.labeler.label_variants(
        [candidate.variant for candidate in candidates], region
    )

    # Remove any candidates we couldn't label, yielding candidate, label pairs.
    for candidate, label in zip(candidates, labels):
      if label.is_confident:
        # Downsample candidates if enabled.
        if (
            self.options.mode
            == deepvariant_pb2.MakeExamplesOptions.Mode.TRAINING
            and self.options.downsample_classes
        ):
          if (
              random.random()
              > self.options.downsample_classes[label.get_class()]
          ):
            continue
        yield candidate, label


def move_to_the_next_non_exhausted_shard(
    shard_index: int, i_th_index: List[int], position_arrays: List[Any]
) -> int:
  """Returns the index of the next non-exhausted shard.

  Args:
    shard_index: int. Index of the current shard being processed.
    i_th_index: List[int]. Current position within i-th shard.
    position_arrays: List[Any]. List of arrays containing candidate positions
      for each shard.

  Returns:
    int. Index of the next shard to be processed.
  """
  i = 0
  while i < len(position_arrays):
    shard_index += 1
    if shard_index >= len(position_arrays):
      shard_index = 0
    if i_th_index[shard_index] < len(position_arrays[shard_index]):
      break
    i += 1
  return shard_index


def merge_ranges_from_files_sequential(position_arrays: List[Any]) -> List[int]:
  """Merges input arrays containing sorted candidate positions.

  positions_array contains all candidate positions for each shart. make_examples
  generates candidate positions in a round robin pattern. So, in order to merge
  all candidate positions from all shards we take candidate positions from the
  first shard, first partition, then second shard first partition and so on.
  Partitions within a shard are separated by END_OF_PARTITION special number.
  <Shard_1 candidates, partition_1>, <Shard_2 candidates, partition_1>, ...
  <Shard_N candidates, partition_1>,
  <Shard_1 candidates, partition_2>, <Shard_2 candidates, partition_2> ...
  <Shard_N candidates, partition_2>,
  ...
  <Shard_N candidates, partition_M>, <Shard_2 candidates, partition_M> ...
  <Shard_N candidates, partition_M>,

  position_arrays is a list of arrays of int32 values. Each list's item contains
  candidate positions for one shard. Candidates positions within each shard are
  not continuous. See regions_to_process() for details.

  Args:
    position_arrays: list of numpy arrays of int32 containing candidate
      positions for each shard.

  Returns:
    List[int] Sorted candidate positions with END_OF_REGION separators.
  """
  candidate_positions_sorted = []
  i_th_index = [0] * len(position_arrays)
  shard_index = 0
  num_arrays_left = len(position_arrays)
  # Iterate until all shards are consumed
  # Add sorted items until -1 is reached. -1 is added as well.
  while num_arrays_left > 0:
    items_added = 0
    # Iterate over positions in one shard.
    while i_th_index[shard_index] < len(position_arrays[shard_index]):
      # Once END_OF_PARTITION is reached we need to move to the next shard.
      if (
          position_arrays[shard_index][i_th_index[shard_index]]
          == END_OF_PARTITION
      ):
        i_th_index[shard_index] += 1
        # If END_OF_REGION is encountered we need to move to the next shard.
        if (
            i_th_index[shard_index] < len(position_arrays[shard_index])
            and position_arrays[shard_index][i_th_index[shard_index]]
            == END_OF_REGION
        ):
          candidate_positions_sorted.append(END_OF_REGION)
          i_th_index[shard_index] += 1
        # Move to the next shard.
        break
      else:
        # Assert that items are sorted
        if candidate_positions_sorted:
          assert (
              position_arrays[shard_index][i_th_index[shard_index]]
              > candidate_positions_sorted[-1]
          )
        candidate_positions_sorted.append(
            position_arrays[shard_index][i_th_index[shard_index]]
        )
        i_th_index[shard_index] += 1
        items_added += 1

    # If all items of the shard are consumed then remove the shard from
    # processing.
    if i_th_index[shard_index] == len(position_arrays[shard_index]):
      num_arrays_left -= 1
    # Move to the next shard
    shard_index = move_to_the_next_non_exhausted_shard(
        shard_index, i_th_index, position_arrays
    )

  logging.info(
      'Total number of candidates: %d', len(candidate_positions_sorted)
  )
  return candidate_positions_sorted


def load_candidate_positions(candidate_path: Any) -> List[int]:
  """Load candidate positions from input file(s)."""
  paths = sharded_file_utils.maybe_generate_sharded_filenames(candidate_path)
  positions = []
  for file_path in paths:
    try:
      with epath.Path(file_path).open('rb') as my_file:
        positions.append(np.frombuffer(my_file.read(), dtype=np.int32))
    except IOError:
      continue
  return merge_ranges_from_files_sequential(positions)


def processing_regions_from_options(
    options: deepvariant_pb2.MakeExamplesOptions,
) -> Tuple[List[range_pb2.Range], Optional[ranges.RangeSet]]:
  """Computes the calling regions from our options.

  This function does all of the work needed to read our input files and region
  specifications to determine the list of regions we should generate examples
  over. It also computes the confident regions needed to label variants.

  Args:
    options: deepvariant.MakeExamplesOptions proto containing information about
      our input data sources.

  Raises:
    ValueError: if the regions to call is empty.

  Returns:
    A tuple of three values, (regions, sample_regions, calling_regions).
    regions: a list of nucleus.genomics.v1.Range protos of the regions we should
    process.
    sample_regions: a list of nucleus.genomics.v1.Range protos of the regions we
    should sample on to calculate global information over the genome such as
    mean coverage.
    calling_regions: a RangeSet containing the calling regions calculated from
    intersection of input regions, confident regions and regions to exclude.
  """
  # Load candidate_positions if the flag is set. Partitioning logic will depend
  # on whether candidate_positions is set.
  candidate_positions = None
  mode_candidate_sweep = deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP
  main_sample_options = options.sample_options[options.main_sample_index]
  if (
      options.mode != mode_candidate_sweep
      and main_sample_options.candidate_positions
  ):
    candidate_positions = load_candidate_positions(
        main_sample_options.candidate_positions
    )

  ref_contigs = fasta.IndexedFastaReader(
      options.reference_filename
  ).header.contigs

  ref_n_regions = None
  if options.discard_non_dna_regions and not options.calling_regions:
    ref_n_regions = find_ref_n_regions(
        fasta.IndexedFastaReader(options.reference_filename), MIN_NON_DNA_REGION
    )

  # Add in confident regions and vcf_contigs if in training mode.
  vcf_contigs = None
  if in_training_mode(options):
    vcf_contigs = vcf.VcfReader(options.truth_variants_filename).header.contigs
    if all([x.n_bases == 0 for x in vcf_contigs]):
      logging.info(
          (
              '%s header does not contain contig lengths. Will skip contig '
              'consistency checking for this file.'
          ),
          options.truth_variants_filename,
      )
      vcf_contigs = None

  main_sample = options.sample_options[options.main_sample_index]
  all_sam_contigs = [
      sam.SamReader(reads_file).header.contigs
      for reads_file in main_sample.reads_filenames
  ]
  sam_contigs = common_contigs(only_true(*all_sam_contigs))

  contigs = _ensure_consistent_contigs(
      ref_contigs,
      sam_contigs,
      vcf_contigs,
      options.exclude_contigs,
      options.min_shared_contigs_basepairs,
  )
  logging_with_options(
      options, 'Common contigs are %s' % [c.name for c in contigs]
  )
  calling_regions = calling_regions_utils.build_calling_regions(
      ref_contigs,
      options.calling_regions,
      options.exclude_calling_regions,
      ref_n_regions,
  )
  if not calling_regions:
    raise ValueError(
        'The regions to call is empty. Check your --regions and '
        '--exclude_regions flags to make sure they are not '
        'resulting in set of empty region to process. This also '
        'happens if you use "chr20" for a BAM where contig names '
        'don\'t have "chr"s (or vice versa).'
    )
  regions = regions_to_process(
      contigs=contigs,
      partition_size=options.allele_counter_options.partition_size,
      calling_regions=calling_regions,
      task_id=options.task_id,
      num_shards=options.num_shards,
      candidates=candidate_positions,
      round_robin_sampling='tfrecord' in options.examples_filename.lower(),
  )

  region_list = list(regions)

  # When using VcfCandidateImporter, it is safe to skip regions without
  # candidates as long as gVCF output is not needed. There is a tradeoff
  # though because it takes time to read the VCF, which is only worth it if
  # there are enough regions.
  if main_sample.proposed_variants_filename and not gvcf_output_enabled(
      options
  ):
    logging_with_options(
        options,
        (
            'Reading VCF to skip processing some regions without '
            'variants in the --proposed_variants VCF.'
        ),
    )
    before = time.time()
    variant_positions = fetch_vcf_positions(
        [
            sample_option.proposed_variants_filename
            for sample_option in options.sample_options
        ],
        contigs,
        calling_regions,
    )
    filtered_regions = filter_regions_by_vcf(region_list, variant_positions)
    time_elapsed = time.time() - before
    logging_with_options(
        options,
        'Filtering regions took {} seconds and reduced the number of '
        'regions to process from {} to {} regions containing variants '
        'from the supplied VCF of proposed variants.'.format(
            trim_runtime(time_elapsed), len(region_list), len(filtered_regions)
        ),
    )
    return filtered_regions, None
  return region_list, calling_regions


def make_examples_runner(options: deepvariant_pb2.MakeExamplesOptions):
  """Runs examples creation stage of deepvariant."""
  resource_monitor = resources.ResourceMonitor().start()
  before_initializing_inputs = time.time()

  logging_with_options(options, 'Preparing inputs')
  (
      regions,
      calling_regions,
  ) = processing_regions_from_options(options)
  main_sample = options.sample_options[options.main_sample_index]
  mode_candidate_sweep = deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP
  candidates_writer = None
  if options.mode == mode_candidate_sweep and main_sample.candidate_positions:
    _, candidate_positions_filename = sharded_file_utils.resolve_filespecs(
        options.task_id, main_sample.candidate_positions
    )
    candidates_writer = epath.Path(candidate_positions_filename).open('wb')

  # Create a processor to create candidates and examples for each region.
  # Replace path in calling regions with the actual calling regions.
  calling_regions = list(calling_regions) if calling_regions else None
  region_processor = RegionProcessor(options, calling_regions)
  region_processor.initialize()

  if options.candidates_filename:
    logging_with_options(
        options, 'Writing candidates to %s' % options.candidates_filename
    )
  if options.gvcf_filename:
    logging_with_options(
        options, 'Writing gvcf records to %s' % options.gvcf_filename
    )

  last_reported = 0

  writers_dict = {}
  samples_that_need_writers = [
      sample
      for sample in region_processor.samples
      if not sample.options.skip_output_generation
  ]
  if not samples_that_need_writers:
    raise ValueError(
        'At least one sample should have skip_output_generation=False.'
    )

  if in_training_mode(options) or len(samples_that_need_writers) == 1:
    writers_dict[options.sample_role_to_train] = OutputsWriter(
        options, suffix=None
    )
  else:
    for sample in samples_that_need_writers:
      if sample.sam_readers is not None:
        writers_dict[sample.options.role] = OutputsWriter(
            options, suffix=sample.options.role
        )
  region_processor.writers_dict = writers_dict
  example_filenames = {
      role: writers_dict[role].examples_filename for role in writers_dict
  }
  region_processor.make_examples_native = (
      make_examples_native_module.ExamplesGenerator(options, example_filenames)
  )

  logging_with_options(
      options,
      'Writing examples to %s'
      % ', '.join(
          [writer.examples_filename for writer in writers_dict.values()]
      ),
  )

  logging_with_options(
      options,
      'Overhead for preparing inputs: %d seconds'
      % (time.time() - before_initializing_inputs),
  )

  running_timer = timer.TimerStart()
  # Ideally this would use dv_constants.NUM_CLASSES, which requires generalizing
  # deepvariant_pb2.MakeExamplesStats to use an array for the class counts.
  n_stats = {
      'n_class_0': 0,
      'n_class_1': 0,
      'n_class_2': 0,
      'n_denovo': 0,
      'n_non_denovo': 0,
      'n_snps': 0,
      'n_indels': 0,
      'n_regions': 0,
      'n_candidates': 0,
      'n_examples': 0,
      'n_small_model_examples': 0,
      'n_small_model_calls': 0,
      'n_phased_candidates': 0,
  }
  example_shape = None
  region_n = 0
  # Get all denovo regions
  denovo_regions = read_denovo_regions(options.denovo_regions_filename)
  for region in regions:
    region_n += 1
    if options.output_debug_info:
      logging_with_options(
          options,
          'Processing %s:%d-%d'
          % (region.reference_name, region.start, region.end),
      )

    if options.mode == mode_candidate_sweep and candidates_writer:
      candidates_in_region = list(
          region_processor.find_candidate_positions(region)
      )
      candidates_writer.write(
          np.array(candidates_in_region, dtype=np.int32).tobytes()
      )
      # Here we mark the end of the calling region
      for cr in calling_regions:
        if cr.reference_name == region.reference_name and cr.end == region.end:
          candidates_writer.write(
              np.array([END_OF_REGION], dtype=np.int32).tobytes()
          )
      continue

    (
        candidates_by_sample,
        gvcfs_by_sample,
        runtimes,
        read_phases_by_sample,
        phased_reads_count,
    ) = region_processor.process(region, region_n)
    for sample in samples_that_need_writers:
      role = sample.options.role
      if role not in candidates_by_sample:
        continue
      if in_training_mode(options) and options.sample_role_to_train != role:
        continue
      if sample.options.skip_output_generation:
        continue
      writer = writers_dict[role]

      if options.write_small_model_examples:
        region_processor.write_small_model_examples_in_region(
            candidates_by_sample[role],
            read_phases_by_sample[role],
            sample,
            region,
            writer,
            n_stats,
            runtimes,
        )
      if options.skip_pileup_image_generation:
        continue

      candidates_for_pileup_images = candidates_by_sample[role]
      if options.call_small_model_examples:
        candidates_not_called_by_small_model = (
            region_processor.call_small_model_examples_in_region(
                candidates_by_sample[role],
                read_phases_by_sample[role],
                sample,
                writer,
                n_stats,
                runtimes,
            )
        )
        candidates_for_pileup_images = candidates_not_called_by_small_model

      region_example_shape = region_processor.writes_examples_in_region(
          candidates_for_pileup_images,
          region,
          sample.options.order,
          n_stats,
          runtimes,
          role,
          denovo_regions,
      )
      if example_shape is None and region_example_shape is not None:
        example_shape = region_example_shape
      gvcfs = gvcfs_by_sample[role]

      n_stats['n_candidates'] += len(candidates_by_sample[role])
      n_stats['n_regions'] += 1
      n_stats['n_phased_candidates'] += phased_reads_count

      before_write_outputs = time.time()
      writer.write_candidates(*candidates_by_sample[role])

      # If we have any gvcf records, write them out. This also serves to
      # protect us from trying to write to the gvcfs output of writer when gvcf
      # generation is turned off. In that case, gvcfs will always be empty and
      # we'll never execute the write.
      if gvcfs:
        writer.write_gvcfs(*gvcfs)

      if options.runtime_by_region:
        runtimes['write outputs'] = runtimes.get('write outputs', 0) + (
            trim_runtime(time.time() - before_write_outputs)
        )
        runtimes['region'] = ranges.to_literal(region)

      # Output timing for every N candidates.
      if (
          int(n_stats['n_candidates'] / options.logging_every_n_candidates)
          > last_reported
          or n_stats['n_regions'] == 1
      ):
        last_reported = int(
            n_stats['n_candidates'] / options.logging_every_n_candidates
        )
        logging_with_options(
            options,
            '%s candidates (%s examples) [%0.2fs elapsed]'
            % (
                n_stats['n_candidates'],
                n_stats['n_examples'],
                running_timer.Stop(),
            ),
        )
        running_timer = timer.TimerStart()
    if options.runtime_by_region:
      # Runtimes are for all samples, so write this only once.
      writers_dict[options.sample_role_to_train].write_runtime(
          stats_dict=runtimes
      )

  for writer in writers_dict.values():
    writer.close_all()
  if options.mode == mode_candidate_sweep and candidates_writer:
    candidates_writer.close()

  # Construct and then write out our MakeExamplesRunInfo proto.
  if options.run_info_filename:
    make_examples_stats = deepvariant_pb2.MakeExamplesStats(
        num_examples=n_stats['n_examples'],
        num_snps=n_stats['n_snps'],
        num_indels=n_stats['n_indels'],
        num_class_0=n_stats['n_class_0'],
        num_class_1=n_stats['n_class_1'],
        num_class_2=n_stats['n_class_2'],
        num_denovo=n_stats['n_denovo'],
        num_nondenovo=n_stats['n_non_denovo'],
    )
    run_info = deepvariant_pb2.MakeExamplesRunInfo(
        options=options,
        resource_metrics=resource_monitor.metrics(),
        stats=make_examples_stats,
    )
    if in_training_mode(options):
      if (
          region_processor.labeler is not None
          and region_processor.labeler.metrics is not None
      ):
        run_info.labeling_metrics.CopyFrom(region_processor.labeler.metrics)
      else:
        logging.warning(
            (
                'Labeling metrics requested but the selected labeling '
                'algorithm %s does not collect metrics; skipping.'
            ),
            options.labeler_algorithm,
        )
    logging_with_options(
        options, 'Writing MakeExamplesRunInfo to %s' % options.run_info_filename
    )
    write_make_examples_run_info(run_info, path=options.run_info_filename)

  # Write to .example_info file. Here we use the examples_filename as prefix.
  # If the examples_filename is sharded, we only write to the first shard.
  # Currently, even in multi-sample scenario, the suffix is not used here
  # because currently all the multiple-sample output will have the same shape
  # and list of channels.
  example_info_filename = dv_utils.get_example_info_json_filename(
      options.examples_filename, options.task_id
  )
  if example_info_filename is not None:
    logging_with_options(
        options, 'Writing example info to %s' % example_info_filename
    )
    example_channels = region_processor.get_channels()
    # example_shape was filled in during the loop above.
    logging.info('example_shape = %s', str(example_shape))
    logging.info('example_channels = %s', str(example_channels))
    with epath.Path(example_info_filename).open('w') as fout:
      json.dump(
          {
              'version': dv_vcf_constants.DEEP_VARIANT_VERSION,
              'shape': example_shape,
              'channels': example_channels,
          },
          fout,
      )

  region_processor.make_examples_native.signal_shard_finished()
  log_summary_stats(options, n_stats)
