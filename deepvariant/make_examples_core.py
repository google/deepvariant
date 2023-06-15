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
import dataclasses
import itertools
import json
import os
import time
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple


from absl import flags
from absl import logging
from etils import epath
import numpy as np

from deepvariant import allele_frequency
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import dv_vcf_constants
from deepvariant import pileup_image
from deepvariant import resources
from deepvariant import variant_caller as vc_base
from deepvariant import vcf_candidate_importer
from deepvariant import very_sensitive_caller
from deepvariant.labeler import customized_classes_labeler
from deepvariant.labeler import haplotype_labeler
from deepvariant.labeler import positional_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.python import direct_phasing
from deepvariant.realigner import realigner
from deepvariant.vendor import timer
from google.protobuf import text_format
from third_party.nucleus.io import fasta
from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io import sam
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import struct_utils
from third_party.nucleus.util import utils
from third_party.nucleus.util import variant_utils

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


def assign_sample_name(sample_name_flag, reads_filenames):
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


def make_vc_options(sample_name, flags_obj: flags.FlagValues):
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


def resolve_sam_aux_fields(flags_obj: flags.FlagValues):
  """Decide value of parse_sam_aux_fields based on other flags."""
  flags_requiring_sam_aux_fields = [
      'sort_by_haplotypes',
      'use_original_quality_scores',
  ]
  flags_using_sam_aux_fields_optionally = ['add_hp_channel']

  parse_sam_aux_fields = flags_obj.parse_sam_aux_fields
  if parse_sam_aux_fields is None:
    # User didn't set the 'parse_sam_aux_fields' flag, so default to False
    # unless a flag is on that would use it.
    parse_sam_aux_fields = False
    for flag_name in (
        flags_requiring_sam_aux_fields + flags_using_sam_aux_fields_optionally
    ):
      if flags_obj[flag_name].value:
        logging.info(
            (
                'Because --%s=true, --parse_sam_aux_fields is set to '
                'true to enable reading auxiliary fields from reads.'
            ),
            flag_name,
        )
      parse_sam_aux_fields = True

  if not parse_sam_aux_fields:
    for flag_name in flags_requiring_sam_aux_fields:
      if flags_obj[flag_name].value:
        raise ValueError(
            f'If --{flag_name} is '
            'set then --parse_sam_aux_fields must be set too.'
        )

    for flag_name in flags_using_sam_aux_fields_optionally:
      if flags_obj[flag_name].value:
        logging.info(
            (
                'Note that --%s is set but --parse_sam_aux_fields is not set.'
                ' This is fine unless you are expecting to use aux fields from'
                ' the alignments file, such as haplotype tags from phasing. If'
                ' you do need to use aux fields, enable --parse_sam_aux_fields.'
            ),
            flag_name,
        )
  return parse_sam_aux_fields


def parse_regions_flag(regions_flag_value):
  if isinstance(regions_flag_value, str):
    regions_flag_value = regions_flag_value.split()
  return regions_flag_value


def logging_with_options(options, message):
  """If options contain multiple shards, log with task/shard prefix."""
  if options.num_shards > 1:
    prefix = 'Task {}/{}: '.format(options.task_id, options.num_shards)
  else:
    prefix = ''
  logging.info('%s%s', prefix, message)


# ---------------------------------------------------------------------------
# Simple utilities
# ---------------------------------------------------------------------------


def in_training_mode(options):
  return options.mode == deepvariant_pb2.MakeExamplesOptions.TRAINING


def in_candidate_sweep_mode(options):
  return options.mode == deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP


def gvcf_output_enabled(options):
  """Returns True if we should be generating gVCF output."""
  return bool(options.gvcf_filename)


def only_true(*elts):
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
    ref_contigs,
    sam_contigs,
    vcf_contigs,
    exclude_contig_names=None,
    min_coverage_fraction=1.0,
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


def common_contigs(contigs_list):
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

  def common2(contigs1, contigs2):
    """Computes the common contigs between contigs1 and contigs2."""
    map2 = ranges.contigs_dict(contigs2)

    def is_common(contig1):
      contig2 = map2.get(contig1.name, None)
      return contig2 and contig1.n_bases == contig2.n_bases

    return [c for c in contigs1 if is_common(c)]

  # Compute the common contigs by recursively getting common contigs of our
  # cumulative set of contigs (common) and each contig in other_contigs.
  common = contigs_list[0]
  for other_contigs in contigs_list[1:]:
    common = common2(common, other_contigs)

  return common


def validate_reference_contig_coverage(
    ref_contigs, shared_contigs, min_coverage_fraction
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


def build_calling_regions(
    contigs, regions_to_include, regions_to_exclude, ref_n_regions
):
  """Builds a RangeSet containing the regions we should call variants in.

  This function intersects the Ranges spanning all of the contigs with those
  from regions_to_include, if not empty, and removes all of the regions in
  regions_to_exclude.

  Args:
    contigs: Sequence of ContigInfo protos. Used to determine the initial ranges
      to process (i.e., all bases of these contigs).
    regions_to_include: RangeSet or iterable that can be converted to a
      RangeSet.
    regions_to_exclude: RangeSet or iterable that can be converted to a
      RangeSet.
    ref_n_regions: List of Range containing non DNA bases to exclude.

  Returns:
    A RangeSet.
  """
  # Initially we are going to call everything in the reference.
  regions = ranges.RangeSet.from_contigs(contigs)

  # If we provided a regions to include, intersect it with all of the regions,
  # producing a common set of regions between the reference and the provided
  # calling regions.
  contig_dict = ranges.contigs_dict(contigs)
  if regions_to_include:
    regions = regions.intersection(
        ranges.RangeSet.from_regions(regions_to_include, contig_dict)
    )

  if ref_n_regions:
    regions.exclude_regions(ranges.RangeSet(ref_n_regions))

  # If we provided regions to exclude, intersect those with the existing calling
  # regions to further refine our set of contigs to process.
  if regions_to_exclude:
    # exclude_regions mutates regions.
    regions.exclude_regions(
        ranges.RangeSet.from_regions(regions_to_exclude, contig_dict)
    )

  return regions


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
    contigs,
    partition_size,
    calling_regions=None,
    task_id=None,
    num_shards=None,
    candidates=None,
):
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

  partitioned = []
  # Depending on candidates parameter we choose the partitioning method.
  if candidates is not None:
    partitioned = partition_by_candidates(regions, candidates, 200)
  else:
    partitioned = regions.partition(partition_size)

  if num_shards:
    return (r for i, r in enumerate(partitioned) if i % num_shards == task_id)
  else:
    return partitioned


def fetch_vcf_positions(vcf_paths, contigs, calling_regions):
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


def filter_regions_by_vcf(regions, variant_positions):
  """Filter a list of regions to only those that contain variants.

  Args:
    regions: a list of Range objects representing regions to filter on.
    variant_positions: a list of Range objects containing the positions of
      variants.

  Returns:
    filtered_regions: a list of Range objects, each of which appeared in the
        input regions and contains at least one of the input variants.
  """

  def dict_by_chromosome(list_of_ranges):
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
# Working with samples
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class Sample(object):
  """Organizes sample-level properties.

  options: A SampleOptions proto containing instructions for how to treat the
      sample, most of which will be set from flags.
  sam_readers: SamReader objects with handles on the `reads_filenames` from the
      options.
  in_memory_sam_reader: InMemorySamReader for this sample, which stores the
      alignments for this sample that have been read into memory from the
      sam_readers.
  reads: A list of reads queried from the sam readers.
  allele_counter: An allele counter object for the sample.
  variant_caller: A variant caller for the sample, should be instantiated using
      the options.variant_caller_options.
  """

  options: deepvariant_pb2.SampleOptions
  sam_readers: Optional[Sequence[sam.SamReader]] = None
  in_memory_sam_reader: Optional[sam.InMemorySamReader] = None
  reads: Optional[List[reads_pb2.Read]] = None
  allele_counter: Optional[allelecounter.AlleleCounter] = None
  variant_caller: Optional[vc_base.VariantCaller] = None

  def __repr__(self):
    return '<Sample {}>'.format(str(self.__dict__))


# ---------------------------------------------------------------------------
# Region processor
# ---------------------------------------------------------------------------


def read_confident_regions(options):
  if options.confident_regions_filename:
    return ranges.RangeSet.from_bed(options.confident_regions_filename)
  else:
    return None


def filter_candidates(candidates, select_variant_types):
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


class DiagnosticLogger(object):
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


class OutputsWriter(object):
  """Manages all of the outputs of make_examples in a single place."""

  def __init__(self, options, suffix=None):
    self._writers = {
        k: None
        for k in ['candidates', 'examples', 'gvcfs', 'runtime', 'read_phases']
    }
    self.examples_filename = None

    if options.candidates_filename:
      self._add_writer(
          'candidates',
          tfrecord.Writer(
              self._add_suffix(options.candidates_filename, suffix)
          ),
      )

    if options.examples_filename:
      self.examples_filename = self._add_suffix(
          options.examples_filename, suffix
      )
      self._add_writer('examples', tfrecord.Writer(self.examples_filename))

    if options.gvcf_filename:
      self._add_writer(
          'gvcfs',
          tfrecord.Writer(self._add_suffix(options.gvcf_filename, suffix)),
      )

    if options.runtime_by_region:
      self._add_writer(
          'runtime', epath.Path(options.runtime_by_region).open('w')
      )
      writer = self._writers['runtime']
      if writer is not None:
        writer.__enter__()
        writer.write('\t'.join(RUNTIME_BY_REGION_COLUMNS) + '\n')

    if options.read_phases_output:
      self._add_writer(
          'read_phases', epath.Path(options.read_phases_output).open('w')
      )
      writer = self._writers['read_phases']
      if writer is not None:
        writer.__enter__()
        writer.write('\t'.join(READ_PHASES_OUTPUT_COLUMNS) + '\n')

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

  def write_runtime(self, stats_dict):
    columns = [str(stats_dict.get(k, 'NA')) for k in RUNTIME_BY_REGION_COLUMNS]
    writer = self._writers['runtime']
    writer.write('\t'.join(columns) + '\n')

  def write_read_phase(self, read, phase, region_n):
    writer = self._writers['read_phases']
    if writer is not None:
      read_key = read.fragment_name + '/' + str(read.read_number)
      writer.write('\t'.join([read_key, str(phase), str(region_n)]) + '\n')

  def _add_writer(self, name, writer):
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

  def _write(self, writer_name, *protos):
    writer = self._writers[writer_name]
    if writer:
      for proto in protos:
        writer.write(proto)

  def close_all(self):
    for writer in self._writers.values():
      if writer is not None:
        writer.close()


class RegionProcessor(object):
  """Creates DeepVariant example protos for a single region on the genome.

  This class helps us to run the very sensitive caller, pileup image creator,
  and variant labeler operations on a single region in parallel across many
  regions using the PoolExecutor API. In order to do this we need separate three
  key operations:

  (1) Collect all of the info needed to create our resources (e.g., ref reader)
      at construction. We cannot actually initialize those resources in the
      constructor, though, since we actually want different resources in each
      worker process/thread. I.e., we need lazy resource initialization.

  (2) Actually initialize these resources *after* the worker has been forked
      in our process pool. This gives us a fresh resource to use in each
      separate process.

  (3) Process the region to find candidate variants and process those into our
      tf.Example protos.
  """

  def __init__(self, options):
    """Creates a new RegionProcess.

    Args:
      options: deepvariant.MakeExamplesOptions proto used to specify our
        resources for calling (e.g., reference_filename).
    """
    self.options = options
    self.samples = [Sample(options=x) for x in self.options.sample_options]
    self.initialized = False
    self.ref_reader = None
    self.realigner = None
    self.pic = None
    self.labeler = None
    self.population_vcf_readers = None
    if self.options.phase_reads:
      # One instance of DirectPhasing per lifetime of make_examples.
      self.direct_phasing_cpp = self._make_direct_phasing_obj()
    self.writers_dict = {}

  def _make_direct_phasing_obj(self):
    return direct_phasing.DirectPhasing()

  def _make_allele_counter_for_region(self, region, candidate_positions):
    return allelecounter.AlleleCounter(
        self.ref_reader.c_reader,
        region,
        candidate_positions,
        self.options.allele_counter_options,
    )

  def _make_allele_counter_for_read_overlap_region(
      self, region, full_region, candidate_positions
  ):
    return allelecounter.AlleleCounter.Default(
        self.ref_reader.c_reader,
        region,
        full_region,
        candidate_positions,
        self.options.allele_counter_options,
    )

  def _encode_tensor(self, image_tensor):
    return image_tensor.tostring(), image_tensor.shape

  def _make_sam_readers(
      self, reads_filenames: Sequence[str], downsample_fraction: float
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

    if self.options.use_allele_frequency:
      population_vcf_readers = allele_frequency.make_population_vcf_readers(
          self.options.population_vcf_filenames
      )
      self.population_vcf_readers = population_vcf_readers

    initialize_raligner = (
        self.options.realigner_enabled
        or self.options.pic_options.alt_aligned_pileup != 'none'
        or self.options.allele_counter_options.track_ref_reads
    )

    if initialize_raligner:
      main_sample = self.samples[self.options.main_sample_index]
      input_bam_header = sam.SamReader(
          main_sample.options.reads_filenames[0]
      ).header
      self.realigner = realigner.Realigner(
          self.options.realigner_options,
          self.ref_reader,
          shared_header=input_bam_header,
      )

    self.pic = pileup_image.PileupImageCreator(
        ref_reader=self.ref_reader,
        options=self.options.pic_options,
        samples=self.samples,
    )

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
    confident_regions = read_confident_regions(self.options)

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
      self, variant_caller_options, proposed_variants_filename
  ):
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
      candidates: List[deepvariant_pb2.DeepVariantCall],
      region: range_pb2.Range,
      sample_order: List[int],
      writer: OutputsWriter,
      n_stats: Dict[str, int],
      runtimes: Dict[str, float],
  ) -> Optional[List[int]]:
    """Generates and writes out the examples in a region.

    Args:
      candidates: List of candidates to be processed into examples.
      region: The region to generate examples.
      sample_order: Order of the samples to use when generating examples.
      writer: A OutputsWriter used to write out examples.
      n_stats: A dictionary that is used to accumulate counts for reporting.
      runtimes: A dictionary that recorded runtime information for reporting.

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
      labels = {i: 0 for i in range(0, dv_constants.NUM_CLASSES)}
      types = {
          dv_utils.EncodedVariantType.SNP: 0,
          dv_utils.EncodedVariantType.INDEL: 0,
          dv_utils.EncodedVariantType.UNKNOWN: 0,
      }
      for candidate, label in self.label_candidates(candidates, region):
        for example in self.create_pileup_examples(
            candidate, sample_order=sample_order
        ):
          self.add_label_to_example(example, label)
          _write_example_and_update_stats(
              example, writer, runtimes, labels, types
          )
          n_stats['n_examples'] += 1
          if example_shape is None:
            example_shape = dv_utils.example_image_shape(example)
      if self.options.run_info_filename:
        n_stats['n_class_0'] += labels[0]
        n_stats['n_class_1'] += labels[1]
        n_stats['n_class_2'] += labels[2]
        n_stats['n_snps'] += types[dv_utils.EncodedVariantType.SNP]
        n_stats['n_indels'] += types[dv_utils.EncodedVariantType.INDEL]
    else:
      for candidate in candidates:
        for example in self.create_pileup_examples(
            candidate, sample_order=sample_order
        ):
          _write_example_and_update_stats(example, writer, runtimes)
          n_stats['n_examples'] += 1
          if example_shape is None:
            example_shape = dv_utils.example_image_shape(example)
    runtimes['make pileup images'] = trim_runtime(
        time.time() - before_make_pileup_images
    )
    return example_shape

  def find_candidate_positions(self, region: range_pb2.Range) -> Iterator[int]:
    """Finds all candidate positions within a given region."""
    main_sample = self.samples[self.options.main_sample_index]
    for sample in self.samples:
      # TODO: Refactor this loop. It is used in other places.
      reads = itertools.chain()
      for _, sam_reader in enumerate(sample.sam_readers):
        reads = itertools.chain(reads, sam_reader.query(region))
      try:
        sample.in_memory_sam_reader.replace_reads(reads)
        sample.reads = sample.in_memory_sam_reader.query(region)
        if self.options.max_reads_per_partition > 0:
          random_for_region = np.random.RandomState(self.options.random_seed)
          sample.reads = utils.reservoir_sample(
              sample.reads,
              self.options.max_reads_per_partition,
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

  def process(self, region, region_n=None):
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
    """
    runtimes = {}

    if not self.initialized:
      self.initialize()

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
      contig_dict = ranges.contigs_dict(
          fasta.IndexedFastaReader(
              self.options.reference_filename
          ).header.contigs
      )
      # When candidate partitioning is used region size is variable. Therefore
      # we need to calculate the padding for each region.
      padding_fraction = int(
          (region.end - region.start) * region_padding_percent / 100
      )
      region_expanded = ranges.expand(region, padding_fraction, contig_dict)

      candidates_by_sample, gvcfs_by_sample = self.candidates_in_region(
          region=region, region_n=region_n, padded_region=region_expanded
      )
    else:
      candidates_by_sample, gvcfs_by_sample = self.candidates_in_region(
          region=region, region_n=region_n
      )

    for sample in self.samples:
      role = sample.options.role
      if role not in candidates_by_sample:
        continue
      candidates = candidates_by_sample[role]

      if self.options.select_variant_types:
        candidates = list(
            filter_candidates(candidates, self.options.select_variant_types)
        )
      runtimes['find candidates'] = trim_runtime(
          time.time() - before_find_candidates
      )
      before_make_pileup_images = time.time()

      # Get allele frequencies for candidates.
      if self.options.use_allele_frequency:
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
    return candidates_by_sample, gvcfs_by_sample, runtimes

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
      if self.options.max_reads_per_partition > 0:
        random_for_region = np.random.RandomState(self.options.random_seed)
        reads = utils.reservoir_sample(
            reads, self.options.max_reads_per_partition, random_for_region
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

  def candidates_in_region(
      self,
      region: range_pb2.Range,
      region_n: Optional[int] = None,
      padded_region: Optional[range_pb2.Range] = None,
  ) -> Tuple[
      Dict[str, Sequence[deepvariant_pb2.DeepVariantCall]],
      Dict[str, Sequence[variants_pb2.Variant]],
  ]:
    """Finds candidates in the region using the designated variant caller.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to get candidates for.
      region_n: Order number of the region being processed by this process.
      padded_region: A nucleus.genomics.v1.Range object specifying the padded
        region.

    Returns:
      A 2-tuple of (candidates, gvcfs).
      The first value, candidates, is a dict keyed by sample role, where each
      item is a list of deepvariant_pb2.DeepVariantCalls objects, in
      coordidate order.
      The second value, gvcfs, is a dict keyed by sample role, where
      each item is a list of nucleus.genomics.v1.Variant protos containing gVCF
      information for all reference sites, if gvcf generation is enabled,
      otherwise the gvcfs value is [].
    """
    for sample in self.samples:
      sample.reads = sample.in_memory_sam_reader.query(region)

    main_sample = self.samples[self.options.main_sample_index]
    if not main_sample.reads and not gvcf_output_enabled(self.options):
      # If we are generating gVCF output we cannot safely abort early here as
      # we need to return the gVCF records calculated by the caller below.
      return {}, {}

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
        if self.options.allele_counter_options.normalize_reads:
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
          if self.options.allele_counter_options.normalize_reads:
            cigar, read_shift = sample.allele_counter.normalize_and_add(
                read, sample.options.name
            )
            if cigar:
              if read_shift != 0:
                read.alignment.position.position += read_shift
              del read.alignment.cigar[:]
              for el in cigar:
                read.alignment.cigar.add(
                    operation=el.operation, operation_length=el.operation_length
                )
          else:
            sample.allele_counter.add(read, sample.options.name)

        allele_counters[sample.options.name] = sample.allele_counter

    candidates = {}
    gvcfs = {}
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
      if (
          in_training_mode(self.options)
          and self.options.sample_role_to_train != role
      ):
        continue
      if not sample.options.reads_filenames:
        continue
      candidates[role], gvcfs[role] = sample.variant_caller.calls_and_gvcfs(
          allele_counters=allele_counters,
          target_sample=sample.options.name,
          include_gvcfs=gvcf_output_enabled(self.options),
          include_med_dp=self.options.include_med_dp,
          left_padding=left_padding,
          right_padding=right_padding,
      )

      if self.options.phase_reads:
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
          # Assign phase tag to reads.
          for read_phase, read in zip(read_phases, reads_to_phase):
            # Remove existing values
            del read.info['HP'].values[:]
            if self.options.pic_options.reverse_haplotypes:
              if read_phase in [1, 2]:
                read_phase = 1 + (read_phase % 2)
            read.info['HP'].values.add(int_value=read_phase)
            if writer and self.options.read_phases_output:
              writer.write_read_phase(read, read_phase, region_n)
        reads_to_phase = None

      if padded_region is not None:
        candidates[role] = self.filter_candidates_by_region(
            candidates[role], region
        )

    return candidates, gvcfs

  def align_to_all_haplotypes(self, variant, reads):
    """For each alternate allele, realign reads to it and get "ref" sequences.

    For alt-aligned pileups, this realigns the reads to each of the alternate
    haplotypes. It also outputs the sequence for each alternate allele, which
    is also needed to build the pileup image.

    Args:
      variant: a nucleus.genomics.v1.Variant containing the alt alleles to align
        against.
      reads: a list of reads (nucleus.genomics.v1.Read) to be realigned around
        the variant.

    Returns:
      dict of alignments keyed by haplotype, dict of window sequences keyed by
          haplotype.
    """

    window_width = self.pic.width
    window_half_width = self.pic.half_width

    alt_alleles = list(variant.alternate_bases)
    contig = variant.reference_name
    ref_start = variant.start
    ref_bases = variant.reference_bases
    ref_end = ref_start + len(ref_bases)

    # Sanity check that the reference_bases in the variant match the reference.
    ref_query_at_variant = self.realigner.ref_reader.query(
        ranges.make_range(contig, ref_start, ref_end)
    )
    if ref_bases != ref_query_at_variant:
      raise ValueError(
          'Error: reference_bases property in variant ({})'
          'does not match the bases in the reference ({}) at that '
          'position.'.format(ref_bases, ref_query_at_variant)
      )

    # Margin must be equal to or more than half the window width.
    # Some extra prefix/suffix can be added to anchor alignments, but currently
    # we don't add extra.
    margin = window_half_width
    valid_end = min(
        self.realigner.ref_reader.contig(contig).n_bases, ref_end + margin
    )
    alignment_region = ranges.make_range(
        contig, max(ref_start - margin, 0), valid_end
    )
    trimmed_reads = [realigner.trim_read(r, alignment_region) for r in reads]
    # Filter reads to a minimum read length of 15 bp after trimming.
    reads = [r for r in trimmed_reads if len(r.aligned_sequence) >= 15]
    prefix = self.realigner.ref_reader.query(
        ranges.make_range(contig, max(ref_start - margin, 0), ref_start)
    )
    suffix = self.realigner.ref_reader.query(
        ranges.make_range(contig, ref_end, valid_end)
    )

    alignments_by_haplotype = {}
    sequences_by_haplotype = {}
    for hap in alt_alleles:
      # Align to each of the alt_alleles:
      alignments_by_haplotype[hap] = self.realigner.align_to_haplotype(
          this_haplotype=hap,
          haplotypes=[hap],
          prefix=prefix,
          suffix=suffix,
          reads=reads,
          contig=contig,
          ref_start=ref_start - len(prefix),
      )
      # Sequence of the alt haplotype in the window:
      end_of_prefix = prefix[-window_half_width:]
      beginning_of_suffix = suffix[: max(window_half_width + 1 - len(hap), 0)]
      sequences_by_haplotype[hap] = end_of_prefix + hap + beginning_of_suffix
      # Long haplotypes can extend past the window, so enforce the width here.
      sequences_by_haplotype[hap] = sequences_by_haplotype[hap][0:window_width]
    return {
        'alt_alignments': alignments_by_haplotype,
        'alt_sequences': sequences_by_haplotype,
    }

  def create_pileup_examples(self, dv_call, sample_order=None):
    """Creates a tf.Example for DeepVariantCall.

    This function calls PileupImageCreator.create_pileup_images on dv_call to
    get raw image tensors for each alt_allele option (see docs for details).
    These tensors are encoded as pngs, and all of the key information is encoded
    as a tf.Example via a call to dv_utils.make_example.

    Args:
      dv_call: A DeepVariantCall.
      sample_order: A list of indices representing the order in which samples
        should be represented in the pileup image. Example: [1,0,2] to swap the
        first and second samples. This is None by default which puts the samples
        in order.

    Returns:
      A list of tf.Example protos.
    """
    reads_for_samples = [
        self.pic.get_reads(
            dv_call.variant, sam_reader=sample.in_memory_sam_reader
        )
        for sample in self.samples
    ]

    logging.vlog(
        3,
        'create_pileup_examples for variant: {}:{}_{}'.format(
            dv_call.variant.reference_name,
            dv_call.variant.start,
            dv_call.variant.reference_bases,
        ),
    )

    # Decide whether each candidate needs ALT-alignment.
    alt_align_this_variant = False
    if self.options.pic_options.alt_aligned_pileup != 'none':
      if self.options.pic_options.types_to_alt_align == 'indels':
        alt_align_this_variant = variant_utils.is_indel(dv_call.variant)
      else:  # types_to_alt_align can only be 'all' or 'indels'.
        alt_align_this_variant = True

    haplotype_alignments_for_samples = None
    haplotype_sequences = None
    if alt_align_this_variant:
      # Align the reads against each alternate allele, saving the sequences of
      # those alleles along with the alignments for pileup images.
      alt_info_for_samples = [
          self.align_to_all_haplotypes(dv_call.variant, reads)
          for reads in reads_for_samples
      ]
      # Each sample has different reads and thus different alt-alignments.
      haplotype_alignments_for_samples = [
          sample['alt_alignments'] for sample in alt_info_for_samples
      ]
      # All samples share the same alt sequences, so select the first one.
      haplotype_sequences = alt_info_for_samples[0]['alt_sequences']


    pileup_images = self.pic.create_pileup_images(
        dv_call=dv_call,
        reads_for_samples=reads_for_samples,
        sample_order=sample_order,
        haplotype_alignments_for_samples=haplotype_alignments_for_samples,
        haplotype_sequences=haplotype_sequences,
    )

    if pileup_images is None:
      # We cannot build a PileupImage for dv_call, issue a warning.
      logging.warning(
          'Could not create PileupImage for candidate at %s:%s',
          dv_call.variant.reference_name,
          dv_call.variant.start,
      )
      return []

    examples = []
    for alt_alleles, image_tensor in pileup_images:
      encoded_tensor, shape = self._encode_tensor(image_tensor)
      examples.append(
          dv_utils.make_example(
              dv_call.variant,
              alt_alleles,
              encoded_tensor,
              shape=shape,
              sequencing_type=self.options.pic_options.sequencing_type,
          )
      )
    return examples

  def get_channels(self) -> List[int]:
    # All the example would have the same list of channels based on `self.pic`.
    return self.pic.get_channels()

  def label_candidates(self, candidates, region):
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
        yield candidate, label

  def add_label_to_example(self, example, label):
    """Adds label information about the assigned label to our example.

    Args:
      example: A tf.Example proto. We will write truth_variant and label into
        this proto.
      label: A variant_labeler.Label object containing the labeling information
        to add to our example.

    Returns:
      The example proto with label fields added.

    Raises:
      ValueError: if label isn't confident.
    """
    if not label.is_confident:
      raise ValueError(
          'Cannot add a non-confident label to an example', example, label
      )
    alt_alleles_indices = dv_utils.example_alt_alleles_indices(example)

    dv_utils.example_set_variant(example, label.variant)

    # Set the label of the example to the # alts given our alt_alleles_indices.
    dv_utils.example_set_label(
        example, label.label_for_alt_alleles(alt_alleles_indices)
    )
    return example


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


def processing_regions_from_options(options):
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
    Two values. The first is a list of nucleus.genomics.v1.Range protos of the
    regions we should process. The second is a RangeSet containing the calling
    regions calculated from intersection of input regions, condident regions
    and regions to exclude.
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
  calling_regions = build_calling_regions(
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


def _write_example_and_update_stats(
    example, writer, runtimes, labels=None, types=None
):
  """Writes out the example using writer; updates labels and types as needed."""
  writer.write_examples(example)
  if runtimes:
    if 'num examples' not in runtimes:
      runtimes['num examples'] = 0
    runtimes['num examples'] += 1
  if labels is not None and types is not None:
    example_label = dv_utils.example_label(example)
    example_type = dv_utils.encoded_variant_type(
        dv_utils.example_variant(example)
    )
    labels[example_label] += 1
    types[example_type] += 1


def make_examples_runner(options):
  """Runs examples creation stage of deepvariant."""
  resource_monitor = resources.ResourceMonitor().start()
  before_initializing_inputs = time.time()

  logging_with_options(options, 'Preparing inputs')
  regions, calling_regions = processing_regions_from_options(options)

  main_sample = options.sample_options[options.main_sample_index]
  mode_candidate_sweep = deepvariant_pb2.MakeExamplesOptions.CANDIDATE_SWEEP
  if options.mode == mode_candidate_sweep and main_sample.candidate_positions:
    _, candidate_positions_filename = sharded_file_utils.resolve_filespecs(
        options.task_id, main_sample.candidate_positions
    )
    candidates_writer = epath.Path(candidate_positions_filename).open('wb')

  # Create a processor to create candidates and examples for each region.
  region_processor = RegionProcessor(options)
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
  if in_training_mode(options) or len(options.sample_options) == 1:
    writers_dict[options.sample_role_to_train] = OutputsWriter(
        options, suffix=None
    )
  else:
    for sample in region_processor.samples:
      if sample.sam_readers is not None:
        writers_dict[sample.options.role] = OutputsWriter(
            options, suffix=sample.options.role
        )
  region_processor.writers_dict = writers_dict

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
      'n_snps': 0,
      'n_indels': 0,
      'n_regions': 0,
      'n_candidates': 0,
      'n_examples': 0,
  }
  example_shape = None
  region_n = 0
  for region in regions:
    region_n += 1

    if options.mode == mode_candidate_sweep:
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

    (candidates_by_sample, gvcfs_by_sample, runtimes) = (
        region_processor.process(region, region_n)
    )
    for sample in region_processor.samples:
      role = sample.options.role
      if role not in candidates_by_sample:
        continue
      writer = writers_dict[role]
      region_example_shape = region_processor.writes_examples_in_region(
          candidates_by_sample[role],
          region,
          sample.options.order,
          writer,
          n_stats,
          runtimes,
      )
      if example_shape is None and region_example_shape is not None:
        example_shape = region_example_shape
      gvcfs = gvcfs_by_sample[role]

      n_stats['n_candidates'] += len(candidates_by_sample[role])
      n_stats['n_regions'] += 1

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

  logging_with_options(
      options, 'Found %s candidate variants' % n_stats['n_candidates']
  )
  logging_with_options(options, 'Created %s examples' % n_stats['n_examples'])
