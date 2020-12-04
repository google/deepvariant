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
"""Step one of DeepVariant: creates tf.Example protos for training/calling."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import os

from absl import app
from absl import flags
from absl import logging
import numpy as np
import tensorflow.compat.v1 as tf

from deeptrio import dt_constants
from deeptrio import very_sensitive_caller
from deeptrio.protos import deeptrio_pb2
from deepvariant import exclude_contigs
from deepvariant import logging_level
from deepvariant import make_examples_utils
from deepvariant import pileup_image
from deepvariant import resources
from deepvariant import tf_utils
from deepvariant.labeler import customized_classes_labeler
from deepvariant.labeler import haplotype_labeler
from deepvariant.labeler import positional_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.realigner import realigner
from deepvariant.vendor import timer
from google.protobuf import text_format
from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils
from third_party.nucleus.util import variant_utils

FLAGS = flags.FLAGS

# Sentinel command line flag value indicating no downsampling should occur.
NO_DOWNSAMPLING = 0.0

# Sentinel command line flag value indicating no random ref sites should be
# emitted.
NO_RANDOM_REF = 0.0

# The name used for a sample if one is not specified or present in the reads.
_UNKNOWN_SAMPLE = 'UNKNOWN'

# The extension we add to our examples path to write our MakeExamplesRunInfo
# protobuf.
_RUN_INFO_FILE_EXTENSION = '.run_info.pbtxt'

# Use a default hts_block_size value of 128 MB (see internal for details) to
# improve SAM/BAM reading throughput, particularly on remote filesystems. Do not
# modify this default parameter without a systematic evaluation of the impact
# across a variety of distributed filesystems!
_DEFAULT_HTS_BLOCK_SIZE = 128 * (1024 * 1024)

flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads.')
flags.DEFINE_string(
    'sample_name_to_call', None,
    'Optional - if not set, default to the value in '
    '--sample_name. The default is set to be backward '
    'compatible. If set, it has to match one of --sample_name, '
    '--sample_name_parent1, or --sample_name_parent2. '
    'This is the sample that we call variants on.')
flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref.')
flags.DEFINE_string(
    'reads_parent1', None,
    'Required. Aligned, sorted, indexed BAM file containing parent 1 reads of '
    'the person we want to call. Should be aligned to a reference genome '
    'compatible with --ref.')
flags.DEFINE_string(
    'reads_parent2', None,
    'Aligned, sorted, indexed BAM file containing parent 2 reads of '
    'the person we want to call. Should be aligned to a reference genome '
    'compatible with --ref.')
flags.DEFINE_bool(
    'use_ref_for_cram', False,
    'If true, use the --ref argument as the reference file for the CRAM '
    'file passed to --reads.  In this case, it is required that the reference '
    'file be located on a local POSIX filesystem.')
flags.DEFINE_string(
    'examples', None,
    'Required. Path to write tf.Example protos in TFRecord format. This is the '
    'root path, as the actual files will be written to this path + suffix, '
    'where suffix corresponds to sample.')
flags.DEFINE_string(
    'candidates', '',
    'Candidate DeepVariantCalls in tfrecord format. For DEBUGGING.')
flags.DEFINE_string('mode', None,
                    'Mode to run. Must be one of calling or training')
flags.DEFINE_string(
    'regions', '',
    'Optional. Space-separated list of regions we want to process. Elements '
    'can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.')
flags.DEFINE_string(
    'exclude_regions', '',
    'Optional. Space-separated list of regions we want to exclude from '
    'processing. Elements can be region literals (e.g., chr20:10-20) or paths '
    'to BED/BEDPE files. Region exclusion happens after processing the '
    '--regions argument, so --region 20 --exclude_regions 20:100 does '
    'everything on chromosome 20 excluding base 100')
flags.DEFINE_string(
    'variant_caller', 'very_sensitive_caller',
    'The caller to use to make examples. Must be one of the VariantCaller enum '
    'values in the DeepTrioOptions proto.')
flags.DEFINE_string(
    'gvcf', '',
    'Optional. Path where we should write gVCF records in TFRecord of Variant '
    'proto format.')
flags.DEFINE_integer(
    'gvcf_gq_binsize', 5,
    'Bin size in which to quantize gVCF genotype qualities. Larger bin size '
    'reduces the number of gVCF records at a loss of quality granularity.')
flags.DEFINE_string(
    'confident_regions', '',
    'Regions that we are confident are hom-ref or a variant in BED format. In '
    'BED or other equivalent format, sorted or unsorted. Contig names must '
    'match those of the reference genome.')
flags.DEFINE_string(
    'truth_variants', '',
    'Tabix-indexed VCF file containing the truth variant calls for this labels '
    'which we use to label our examples.')
flags.DEFINE_integer('task', 0, 'Task ID of this task')
flags.DEFINE_integer(
    'partition_size', 1000,
    'The maximum number of basepairs we will allow in a region before splitting'
    'it into multiple smaller subregions.')
flags.DEFINE_integer(
    'max_reads_per_partition', 1500,
    'The maximum number of reads per partition that we consider before '
    'following processing such as sampling and realigner.')
flags.DEFINE_string(
    'multi_allelic_mode', '',
    'How to handle multi-allelic candidate variants. For DEBUGGING')
flags.DEFINE_bool(
    'realign_reads', True,
    'If True, locally realign reads before calling variants. '
    'Reads longer than 500 bp are never realigned.')
flags.DEFINE_bool(
    'write_run_info', False,
    'If True, write out a MakeExamplesRunInfo proto besides our examples in '
    'text_format.')
# alt_aligned_pileups=rows is unavailable in DeepTrio, since it makes the final
# pileup tensor too tall for InceptionV3.
flags.DEFINE_enum(
    'alt_aligned_pileup', 'none', ['none', 'base_channels', 'diff_channels'],
    'Include alignments of reads against each candidate alternate allele in '
    'the pileup image. "none" turns this feature off. '
    'The default is "none".'
    'Options: "none", "base_channels","diff_channels"')
flags.DEFINE_enum(
    'types_to_alt_align', 'indels', ['indels', 'all'],
    'When --alt_aligned_pileup is not none, this flag determines whether to '
    'align to the alt alleles only for indels or for all variant types '
    'including SNPs. Ignored if --alt_aligned_pileup is "none". This flag is '
    'experimental and is not compatible with the pre-trained release models.')

flags.DEFINE_float(
    'downsample_fraction_child', NO_DOWNSAMPLING,
    'If not ' + str(NO_DOWNSAMPLING) + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input child BAM. This argument makes it easy to create examples '
    'as though the input BAM had less coverage.')
flags.DEFINE_float(
    'downsample_fraction_parents', NO_DOWNSAMPLING,
    'If not ' + str(NO_DOWNSAMPLING) + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input parent BAMs. This argument makes it easy to create examples'
    ' as though the input BAMs had less coverage.')
flags.DEFINE_string(
    'sample_name', '', 'Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads.')
flags.DEFINE_string(
    'sample_name_parent1', '',
    'Parent1 Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads_parent1.')
flags.DEFINE_string(
    'sample_name_parent2', '',
    'Parent2 Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads_parent2.')
flags.DEFINE_string('hts_logging_level',
                    hts_verbose.htsLogLevel.HTS_LOG_WARNING.name,
                    'Sets the htslib logging threshold.')
flags.DEFINE_integer(
    'hts_block_size', _DEFAULT_HTS_BLOCK_SIZE,
    'Sets the htslib block size. Zero or negative uses default htslib setting; '
    'larger values (e.g. 1M) may be beneficial for using remote files. '
    'Currently only applies to SAM/BAM reading.')
flags.DEFINE_integer(
    'min_base_quality', 10,
    'Minimum base quality. This field indicates that we are enforcing a '
    'minimum base quality score for alternate alleles. Alternate alleles will '
    'only be considered if all bases in the allele have a quality greater than '
    'min_base_quality.')
flags.DEFINE_integer(
    'min_mapping_quality', 5,
    'By default, reads with any mapping quality are kept. Setting this field '
    'to a positive integer i will only keep reads that have a MAPQ >= i. Note '
    'this only applies to aligned reads.')
flags.DEFINE_integer(
    'vsc_min_count_snps', 2,
    'SNP alleles occurring at least this many times in our '
    'AlleleCount will be advanced as candidates.')
flags.DEFINE_integer(
    'vsc_min_count_indels', 2,
    'Indel alleles occurring at least this many times in '
    'our AlleleCount will be advanced as candidates.')
flags.DEFINE_float(
    'vsc_min_fraction_snps', 0.12,
    'SNP alleles occurring at least this fraction of all '
    'counts in our AlleleCount will be advanced as '
    'candidates.')
flags.DEFINE_float(
    'vsc_min_fraction_indels', 0.06,
    'Indel alleles occurring at least this fraction of all '
    'counts in our AlleleCount will be advanced as '
    'candidates.')
flags.DEFINE_float(
    'vsc_allele_fraction_trio_coefficient', 0.67,
    'Coefficient that is applied to vsc_min_fraction_snps and '
    'vsc_min_fraction_indels for candidate generation for trio calling.')
flags.DEFINE_float(
    'training_random_emit_ref_sites', NO_RANDOM_REF,
    'If > 0, emit extra random reference examples with this probability.')
flags.DEFINE_integer(
    'pileup_image_height_parent', 0,
    'Height for the parent pileup image. If 0, uses the default height')
flags.DEFINE_integer(
    'pileup_image_height_child', 0,
    'Height for the child pileup image. If 0, uses the default height')
flags.DEFINE_integer(
    'pileup_image_width', 0,
    'Width for the pileup image. If 0, uses the default width')
flags.DEFINE_string(
    'labeler_algorithm', 'haplotype_labeler',
    'Algorithm to use to label examples in training mode. Must be one of the '
    'LabelerAlgorithm enum values in the DeepTrioOptions proto.')
flags.DEFINE_string(
    'customized_classes_labeler_classes_list', '',
    'A comma-separated list of strings that defines customized class labels '
    'for variants. This is only set when labeler_algorithm is '
    'customized_classes_labeler.')
flags.DEFINE_string(
    'customized_classes_labeler_info_field_name', '',
    'The name from the INFO field of VCF where we should get the customized '
    'class labels from. This is only set when labeler_algorithm is '
    'customized_classes_labeler.')
flags.DEFINE_integer(
    'logging_every_n_candidates', 100,
    'Print out the log every n candidates. The smaller the number, the more '
    'frequent the logging information emits.')
flags.DEFINE_bool('keep_duplicates', False, 'If True, keep duplicate reads.')
flags.DEFINE_bool('keep_supplementary_alignments', False,
                  'If True, keep reads marked as supplementary alignments.')
flags.DEFINE_bool('keep_secondary_alignments', False,
                  'If True, keep reads marked as secondary alignments.')
flags.DEFINE_bool(
    'parse_sam_aux_fields', False,
    'If True, auxiliary fields of the SAM/BAM/CRAM records are parsed.')
flags.DEFINE_bool('use_original_quality_scores', False,
                  'If True, base quality scores are read from OQ tag.')
flags.DEFINE_string(
    'select_variant_types', None,
    'If provided, should be a whitespace-separated string of variant types to '
    'keep when generating examples. Permitted values are "snps", "indels", '
    '"multi-allelics", and "all", which select bi-allelic snps, bi-allelic '
    'indels, multi-allelic variants of any type, and all variants, '
    'respectively. Multiple selectors can be specified, so that '
    '--select_variant_types="snps indels" would keep all bi-allelic SNPs and '
    'indels')
flags.DEFINE_string(
    'sequencing_type', None,
    'A string representing input bam file sequencing_type. Permitted values are '
    '"WGS" and "WES", which represent whole genome sequencing and whole exome '
    'sequencing, respectively. This flag is experimental and is not currently '
    'being used.')
flags.DEFINE_bool(
    'sort_by_haplotypes', False,
    'If True, reads are sorted by haplotypes (using HP tag), '
    'parse_sam_aux_fields has to be set for this to work.')

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


_VARIANT_TYPE_SELECTORS = {
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


def parse_proto_enum_flag(proto_enum_pb2,
                          flag_value,
                          skip_unspecified_option=True):
  """Parses a command line flag string value into a protobuf Enum value.

  Args:
    proto_enum_pb2: a enum_type_wrapper.EnumTypeWrapper type containing a proto
      enum definition. For example, this would be
      deeptrio_pb2.DeepTrioOptions.Mode to get the DeepTrioOptions Mode
      enum. See:
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
  except ValueError:
    options = proto_enum_pb2.keys()
    if skip_unspecified_option:
      options = [o for o in options if 'unspecified' not in o.lower()]
    raise ValueError('Unknown enum option "{}". Allowed options are {}'.format(
        flag_value, ','.join(sorted(options))))


def parse_regions_flag(regions_flag_value):
  if isinstance(regions_flag_value, str):
    regions_flag_value = regions_flag_value.split()
  return regions_flag_value


def assign_sample_name(sample_name_flag, reads):
  """Returns sample name derived from either sample_name flag or input BAM.

  Function derives sample_name from the flag. If flag is not set then
  sample_name is derived from input BAM.

  Args:
    sample_name_flag: string. sample_name flag value.
    reads: Iterator of nucleus.genomics.v1.Read reads.

  Returns:
    string. Derived sample name.
  """
  if sample_name_flag:
    sample_name = sample_name_flag
  elif reads and reads is not None:
    with sam.SamReader(reads) as sam_reader:
      sample_name = extract_sample_name_from_sam_reader(sam_reader)
  else:
    sample_name = _UNKNOWN_SAMPLE
  return sample_name


def initialize_variant_caller(the_sample_name, flags_obj):
  return deepvariant_pb2.VariantCallerOptions(
      min_count_snps=flags_obj.vsc_min_count_snps,
      min_count_indels=flags_obj.vsc_min_count_indels,
      min_fraction_snps=flags_obj.vsc_min_fraction_snps,
      min_fraction_indels=flags_obj.vsc_min_fraction_indels,
      vsc_allele_fraction_trio_coefficient=flags_obj
      .vsc_allele_fraction_trio_coefficient,
      # Not specified by default: fraction_reference_sites_to_emit,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=1400605801,
      sample_name=the_sample_name,
      p_error=0.001,
      max_gq=50,
      gq_resolution=flags_obj.gvcf_gq_binsize,
      ploidy=2)


def default_options(add_flags=True, flags_obj=None):
  """Creates a DeepTrioOptions proto populated with reasonable defaults.

  Args:
    add_flags: bool. defaults to True. If True, we will push the value of
      certain FLAGS into our options. If False, those option fields are left
      uninitialized.
    flags_obj: object.  If not None, use as the source of flags, else use global
      FLAGS.

  Returns:
    deeptrio_pb2.DeepTrioOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  if not flags_obj:
    flags_obj = FLAGS

  read_reqs = reads_pb2.ReadRequirements(
      keep_duplicates=flags_obj.keep_duplicates,
      keep_supplementary_alignments=flags_obj.keep_supplementary_alignments,
      keep_secondary_alignments=flags_obj.keep_secondary_alignments,
      min_base_quality=flags_obj.min_base_quality,
      min_mapping_quality=flags_obj.min_mapping_quality,
      min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT)

  logging.info('ReadRequirements are: %s', read_reqs)

  pic_options = pileup_image.default_options(read_requirements=read_reqs)

  # Set DeepTrio-specific pileup constants.
  pic_options.height = dt_constants.PILEUP_DEFAULT_HEIGHT
  pic_options.width = dt_constants.PILEUP_DEFAULT_WIDTH

  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=flags_obj.partition_size, read_requirements=read_reqs)

  sample_name = assign_sample_name(flags_obj.sample_name, flags_obj.reads)
  sample_name_parent1 = assign_sample_name(flags_obj.sample_name_parent1,
                                           flags_obj.reads_parent1)
  sample_name_parent2 = assign_sample_name(flags_obj.sample_name_parent2,
                                           flags_obj.reads_parent2)

  variant_caller_options_child = initialize_variant_caller(
      sample_name, flags_obj)
  variant_caller_options_parent1 = initialize_variant_caller(
      sample_name_parent1, flags_obj)
  variant_caller_options_parent2 = initialize_variant_caller(
      sample_name_parent2, flags_obj)

  options = deeptrio_pb2.DeepTrioOptions(
      exclude_contigs=exclude_contigs.EXCLUDED_HUMAN_CONTIGS,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=609314161,
      # # Not specified by default: calling_regions = 3;
      read_requirements=read_reqs,
      allele_counter_options=allele_counter_options,
      variant_caller_options_child=variant_caller_options_child,
      variant_caller_options_parent1=variant_caller_options_parent1,
      variant_caller_options_parent2=variant_caller_options_parent2,
      pic_options=pic_options,
      n_cores=1,
      task_id=0,
      num_shards=0,
      min_shared_contigs_basepairs=0.9,
  )

  if add_flags:
    options.mode = parse_proto_enum_flag(deeptrio_pb2.DeepTrioOptions.Mode,
                                         flags_obj.mode.upper())

    options.labeler_algorithm = parse_proto_enum_flag(
        deeptrio_pb2.DeepTrioOptions.LabelerAlgorithm,
        flags_obj.labeler_algorithm.upper())

    options.variant_caller = parse_proto_enum_flag(
        deeptrio_pb2.DeepTrioOptions.VariantCaller,
        flags_obj.variant_caller.upper())

    if flags_obj.ref:
      options.reference_filename = flags_obj.ref
    if flags_obj.reads:
      options.reads_filename = flags_obj.reads
    if flags_obj.reads_parent1:
      options.reads_parent1_filename = flags_obj.reads_parent1
    if flags_obj.reads_parent2:
      options.reads_parent2_filename = flags_obj.reads_parent2
    if flags_obj.confident_regions:
      options.confident_regions_filename = flags_obj.confident_regions
    if flags_obj.truth_variants:
      options.truth_variants_filename = flags_obj.truth_variants
    if flags_obj.sequencing_type:
      options.pic_options.sequencing_type = parse_proto_enum_flag(
          deepvariant_pb2.PileupImageOptions.SequencingType,
          flags_obj.sequencing_type)
    if flags_obj.downsample_fraction_child != NO_DOWNSAMPLING:
      options.downsample_fraction_child = flags_obj.downsample_fraction_child
    if flags_obj.downsample_fraction_parents != NO_DOWNSAMPLING:
      options.downsample_fraction_parents = flags_obj.downsample_fraction_parents

    if flags_obj.multi_allelic_mode:
      multi_allelic_enum = {
          'include_het_alt_images':
              deepvariant_pb2.PileupImageOptions.ADD_HET_ALT_IMAGES,
          'exclude_het_alt_images':
              deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES,
      }[flags_obj.multi_allelic_mode]
      options.pic_options.multi_allelic_mode = multi_allelic_enum

    # Pileup heights are assigned to samples later when they are created.
    if flags_obj.pileup_image_height_parent:
      options.height_parent = flags_obj.pileup_image_height_parent
    else:
      options.height_parent = dt_constants.PILEUP_DEFAULT_HEIGHT_PARENT

    if flags_obj.pileup_image_height_child:
      options.height_child = flags_obj.pileup_image_height_child
    else:
      options.height_child = dt_constants.PILEUP_DEFAULT_HEIGHT_CHILD

    if (options.height_parent > 100 or options.height_parent < 10 or
        options.height_child > 100 or options.height_child < 10):
      errors.log_and_raise('Pileup image heights must be between 10 and 100.')

    if flags_obj.pileup_image_width:
      options.pic_options.width = flags_obj.pileup_image_width

    options.pic_options.alt_aligned_pileup = flags_obj.alt_aligned_pileup
    options.pic_options.types_to_alt_align = flags_obj.types_to_alt_align

    if flags_obj.select_variant_types:
      options.select_variant_types[:] = flags_obj.select_variant_types.split()
      for svt in options.select_variant_types:
        if svt not in _VARIANT_TYPE_SELECTORS:
          errors.log_and_raise(
              'Select variant type {} not recognized. Allowed values are {}'
              .format(svt, ', '.join(_VARIANT_TYPE_SELECTORS)),
              errors.CommandLineError)

    num_shards, examples, candidates, gvcf = (
        sharded_file_utils.resolve_filespecs(flags_obj.task,
                                             flags_obj.examples or '',
                                             flags_obj.candidates or '',
                                             flags_obj.gvcf or ''))
    options.examples_filename = examples
    options.candidates_filename = candidates
    options.gvcf_filename = gvcf
    options.task_id = flags_obj.task
    options.num_shards = num_shards
    if flags_obj.use_original_quality_scores and not flags_obj.parse_sam_aux_fields:
      errors.log_and_raise(
          'If use_original_quality_scores is set then parse_sam_aux_fields '
          'must be set too.', errors.CommandLineError)
    options.use_original_quality_scores = flags_obj.use_original_quality_scores

    if flags_obj.sort_by_haplotypes and not flags_obj.parse_sam_aux_fields:
      errors.log_and_raise(
          '--sort_by_haplotypes requires --parse_sam_aux_fields to be set ',
          errors.CommandLineError)
    options.pic_options.sort_by_haplotypes = flags_obj.sort_by_haplotypes

    if flags_obj.write_run_info:
      options.run_info_filename = examples + _RUN_INFO_FILE_EXTENSION

    options.calling_regions.extend(parse_regions_flag(flags_obj.regions))
    options.exclude_calling_regions.extend(
        parse_regions_flag(flags_obj.exclude_regions))

    options.realigner_enabled = flags_obj.realign_reads
    if options.realigner_enabled:
      options.realigner_options.CopyFrom(realigner.realigner_config(flags_obj))

    options.max_reads_per_partition = flags_obj.max_reads_per_partition

    if (options.mode == deeptrio_pb2.DeepTrioOptions.TRAINING and
        flags_obj.training_random_emit_ref_sites != NO_RANDOM_REF):
      options.variant_caller_options_child.fraction_reference_sites_to_emit = (
          flags_obj.training_random_emit_ref_sites)

  return options


# ---------------------------------------------------------------------------
# Simple utilities
# ---------------------------------------------------------------------------


def in_training_mode(options):
  return options.mode == deeptrio_pb2.DeepTrioOptions.TRAINING


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

  Raises:
    ValueError: There is not exactly one unique sample name in the SAM/BAM.
  """
  samples = {
      rg.sample_id for rg in sam_reader.header.read_groups if rg.sample_id
  }
  if not samples:
    raise ValueError(
        'No non-empty sample name found in the input reads. Please provide the '
        'name of the sample with the --sample_name argument.')
  elif len(samples) > 1:
    raise ValueError(
        'Multiple samples ({}) were found in the input reads. DeepVariant can '
        'only call variants from a BAM file containing a single sample.'.format(
            ', '.join(sorted(samples))))
  return next(iter(samples))


# ---------------------------------------------------------------------------
# Utilities for working with labeling metrics
#
# ---------------------------------------------------------------------------


def read_make_examples_run_info(path):
  """Reads a MakeExamplesRunInfo proto in text_format from path."""
  with tf.io.gfile.GFile(path) as f:
    return text_format.Parse(f.read(), deeptrio_pb2.MakeExamplesRunInfo())


def write_make_examples_run_info(run_info_proto, path):
  """Writes a MakeExamplesRunInfo proto in text_format to path."""
  with tf.io.gfile.GFile(path, mode='w') as writer:
    writer.write(text_format.MessageToString(run_info_proto, float_format=''))


# ---------------------------------------------------------------------------
# Region processing
# ---------------------------------------------------------------------------


def _ensure_consistent_contigs(ref_contigs,
                               sam_contigs,
                               vcf_contigs,
                               exclude_contig_names=None,
                               min_coverage_fraction=1.0):
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
  contigs = common_contigs(only_true(ref_contigs, sam_contigs, vcf_contigs))
  validate_reference_contig_coverage(ref_contigs, contigs,
                                     min_coverage_fraction)
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
  # master set of contigs (contigs) and each contig in other_contigs.
  common = contigs_list[0]
  for other_contigs in contigs_list[1:]:
    common = common2(common, other_contigs)

  return common


def validate_reference_contig_coverage(ref_contigs, shared_contigs,
                                       min_coverage_fraction):
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
      pieces.append('"{}" is {} bp and {}'.format(ref_contig.name,
                                                  ref_contig.n_bases, status))
    return ', '.join(pieces)

  ref_bp = ranges.contigs_n_bases(ref_contigs)
  common_bp = ranges.contigs_n_bases(shared_contigs)
  coverage = common_bp / (1. * ref_bp)
  if not shared_contigs or coverage < min_coverage_fraction:
    raise ValueError('Reference contigs span {} bases but only {} bases '
                     '({:.2%}) were found in common among our input files. '
                     'Check that the sources were created on a common genome '
                     'reference build. Contig matches were: {}'.format(
                         ref_bp, common_bp, coverage, format_contig_matches()))


def build_calling_regions(contigs, regions_to_include, regions_to_exclude):
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
        ranges.RangeSet.from_regions(regions_to_include, contig_dict))

  # If we provided regions to exclude, intersect those with the existing calling
  # regions to further refine our set of contigs to process.
  if regions_to_exclude:
    # exclude_regions mutates regions.
    regions.exclude_regions(
        ranges.RangeSet.from_regions(regions_to_exclude, contig_dict))

  return regions


def regions_to_process(contigs,
                       partition_size,
                       calling_regions=None,
                       task_id=None,
                       num_shards=None):
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

  Returns:
    An iterable of nucleus.genomics.v1.Range objects.

  Raises:
    ValueError: if task_id and num_shards are bad or inconsistent.
  """
  if (task_id is None) != (num_shards is None):
    raise ValueError('Both task_id and num_shards must be present if either is',
                     task_id, num_shards)
  if num_shards:
    if num_shards < 0:
      raise ValueError('num_shards={} must be >= 0'.format(num_shards))
    if task_id < 0 or task_id >= num_shards:
      raise ValueError('task_id={} should be >= 0 and < num_shards={}'.format(
          task_id, num_shards))

  regions = ranges.RangeSet.from_contigs(contigs)
  if calling_regions:
    regions = regions.intersection(calling_regions)
  partitioned = regions.partition(partition_size)

  if num_shards:
    return (r for i, r in enumerate(partitioned) if i % num_shards == task_id)
  else:
    return partitioned


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
      _VARIANT_TYPE_SELECTORS or an error will be raised.

  Raises:
    ValueError: if any str in select_variant_types isn't present in
      _VARIANT_TYPE_SELECTORS.

  Yields:
    Candidates in order.
  """
  if not all(s in _VARIANT_TYPE_SELECTORS for s in select_variant_types):
    raise ValueError('Unexpected select variant type', select_variant_types)

  for candidate in candidates:
    v = candidate.variant
    for select_type in select_variant_types:
      selector = _VARIANT_TYPE_SELECTORS[select_type]
      if selector(v):
        yield candidate
        break


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
      options: deepvariant.DeepTrioOptions proto used to specify our resources
        for calling (e.g., reference_filename).
    """
    self.options = options
    self.initialized = False
    self.ref_reader = None
    self.sam_reader = None
    self.in_memory_sam_reader = None
    self.realigner = None
    self.pics = None
    self.labeler = None
    self.variant_caller = None
    self.samples = {}
    self.sample_to_train = None

  def _make_allele_counter_for_region(self, region):
    return allelecounter.AlleleCounter(self.ref_reader.c_reader, region,
                                       self.options.allele_counter_options)

  def _encode_tensor(self, image_tensor):
    return image_tensor.tostring(), image_tensor.shape, 'raw'

  def _make_sam_reader(self,
                       reads_filename,
                       downsample_fraction_param=NO_DOWNSAMPLING):
    return sam.SamReader(
        reads_filename,
        ref_path=FLAGS.ref if FLAGS.use_ref_for_cram else None,
        read_requirements=self.options.read_requirements,
        parse_aux_fields=FLAGS.parse_sam_aux_fields,
        hts_block_size=FLAGS.hts_block_size,
        downsample_fraction=downsample_fraction_param,
        random_seed=self.options.random_seed,
        use_original_base_quality_scores=self.options
        .use_original_quality_scores)

  def _initialize(self):
    """Initialize the resources needed for this work in the current env."""
    if self.initialized:
      raise ValueError('Cannot initialize this object twice')

    self.ref_reader = fasta.IndexedFastaReader(self.options.reference_filename)

    self.sam_reader = self._make_sam_reader(
        self.options.reads_filename, self.options.downsample_fraction_child)
    if self.options.reads_parent1_filename:
      self.sam_reader_parent1 = self._make_sam_reader(
          self.options.reads_parent1_filename,
          self.options.downsample_fraction_parents)
    else:
      self.sam_reader_parent1 = None
    if self.options.reads_parent2_filename:
      self.sam_reader_parent2 = self._make_sam_reader(
          self.options.reads_parent2_filename,
          self.options.downsample_fraction_parents)
    else:
      self.sam_reader_parent2 = None

    # Keep each sample organized with its relevant info.
    child = make_examples_utils.Sample(
        name='child',
        sam_reader=self.sam_reader,
        in_memory_sam_reader=sam.InMemorySamReader([]),
        pileup_height=self.options.height_child)
    parent1 = make_examples_utils.Sample(
        name='parent1',
        sam_reader=self.sam_reader_parent1,
        in_memory_sam_reader=sam.InMemorySamReader([]),
        pileup_height=self.options.height_parent)
    parent2 = make_examples_utils.Sample(
        name='parent2',
        sam_reader=self.sam_reader_parent2,
        in_memory_sam_reader=sam.InMemorySamReader([]),
        pileup_height=self.options.height_parent)

    # Ordering here determines the order of pileups from top to bottom.
    self.samples = {}
    self.samples['child'] = [parent1, child, parent2]
    self.samples['parent1'] = [parent1, child, parent2]
    self.samples['parent2'] = [parent2, child, parent1]

    if self.options.realigner_enabled or self.options.pic_options.alt_aligned_pileup != 'none':
      input_bam_header = sam.SamReader(self.options.reads_filename).header
      self.realigner = realigner.Realigner(
          self.options.realigner_options,
          self.ref_reader,
          shared_header=input_bam_header)

    self.options.pic_options.sequencing_type = deepvariant_pb2.PileupImageOptions.TRIO

    self.pics = {}
    self.pics['child'] = pileup_image.PileupImageCreator(
        ref_reader=self.ref_reader,
        options=self.options.pic_options,
        samples=self.samples['child'])
    self.pics['parent1'] = pileup_image.PileupImageCreator(
        ref_reader=self.ref_reader,
        options=self.options.pic_options,
        samples=self.samples['parent1'])
    self.pics['parent2'] = pileup_image.PileupImageCreator(
        ref_reader=self.ref_reader,
        options=self.options.pic_options,
        samples=self.samples['parent2'])

    if in_training_mode(self.options):
      self.labeler = self._make_labeler_from_options()

    if FLAGS.sample_name_to_call == FLAGS.sample_name:
      self.sample_to_train = 'child'
    elif FLAGS.sample_name_to_call == FLAGS.sample_name_parent1:
      self.sample_to_train = 'parent1'
    else:
      raise ValueError(
          'sample_name_to_call must match either sample_name or sample_name_parent1 '
      )

    self.variant_caller_child = self._make_variant_caller_from_options(
        self.options.variant_caller_options_child)
    self.variant_caller_parent1 = self._make_variant_caller_from_options(
        self.options.variant_caller_options_parent1)
    self.variant_caller_parent2 = self._make_variant_caller_from_options(
        self.options.variant_caller_options_parent2)
    self.initialized = True

  def _make_labeler_from_options(self):
    """Creates the labeler from options."""
    truth_vcf_reader = vcf.VcfReader(
        self.options.truth_variants_filename,
        excluded_format_fields=['GL', 'GQ', 'PL'])
    confident_regions = read_confident_regions(self.options)

    if (self.options.labeler_algorithm ==
        deeptrio_pb2.DeepTrioOptions.POSITIONAL_LABELER):
      return positional_labeler.PositionalVariantLabeler(
          truth_vcf_reader=truth_vcf_reader,
          confident_regions=confident_regions)
    elif (self.options.labeler_algorithm ==
          deeptrio_pb2.DeepTrioOptions.HAPLOTYPE_LABELER):
      return haplotype_labeler.HaplotypeLabeler(
          truth_vcf_reader=truth_vcf_reader,
          ref_reader=self.ref_reader,
          confident_regions=confident_regions)
    elif (self.options.labeler_algorithm ==
          deeptrio_pb2.DeepTrioOptions.CUSTOMIZED_CLASSES_LABELER):
      if (not FLAGS.customized_classes_labeler_classes_list or
          not FLAGS.customized_classes_labeler_info_field_name):
        raise ValueError('For -labeler_algorithm=customized_classes_labeler, '
                         'you need to set '
                         '-customized_classes_labeler_classes_list and '
                         '-customized_classes_labeler_info_field_name.')
      return customized_classes_labeler.CustomizedClassesVariantLabeler(
          truth_vcf_reader=truth_vcf_reader,
          confident_regions=confident_regions,
          classes_list=FLAGS.customized_classes_labeler_classes_list,
          info_field_name=FLAGS.customized_classes_labeler_info_field_name)
    else:
      raise ValueError('Unexpected labeler_algorithm',
                       self.options.labeler_algorithm)

  def _make_variant_caller_from_options(self, variant_caller_options):
    """Creates the variant_caller from options."""
    if (self.options.variant_caller ==
        deeptrio_pb2.DeepTrioOptions.VERY_SENSITIVE_CALLER):
      return very_sensitive_caller.VerySensitiveCaller(variant_caller_options)
    else:
      raise ValueError('Unexpected variant_caller', self.options.variant_caller)

  def initialize(self):
    if not self.initialized:
      self._initialize()

  def process(self, region):
    """Finds candidates and creates corresponding examples in a region.

    Args:
      region: A nucleus.genomics.v1.Range proto. Specifies the region on the
        genome we should process.

    Returns:
      Three values. First is a dict of the found candidates, which are
      deepvariant.DeepVariantCall objects. The second value is a dict of filled
      in tf.Example protos. For example, these will include the candidate
      variant, the pileup image, and, if in training mode, the truth variants
      and labels needed for training. The third value is a dict of
      nucleus.genomics.v1.Variant protos containing gVCF information for all
      reference sites, if gvcf generation is enabled, otherwise returns []. Each
      dict is keyed by sample.
    """
    region_timer = timer.TimerStart()

    # Print some basic information about what we are doing.
    if not self.initialized:
      self._initialize()

    # Get reads in the region for each sample into its in-memory sam reader,
    # optionally realigning reads.
    for sample in self.samples['child']:
      if sample.in_memory_sam_reader is not None:
        sample.in_memory_sam_reader.replace_reads(
            self.region_reads(region, sample.sam_reader))

    # Candidates are created using both parents and child
    candidates_dict, gvcfs_dict = self.candidates_in_region(region)
    examples_dict = {}
    for sample, candidates in candidates_dict.items():
      examples_dict[sample] = []

      if self.options.select_variant_types:
        candidates = list(
            filter_candidates(candidates, self.options.select_variant_types))

      if in_training_mode(self.options):
        for candidate, label in self.label_candidates(candidates, region):
          for example in self.create_pileup_examples(candidate, sample):
            self.add_label_to_example(example, label)
            examples_dict[sample].append(example)
      else:
        for candidate in candidates:
          for example in self.create_pileup_examples(candidate, sample):
            examples_dict[sample].append(example)

      candidates_dict[sample] = candidates
      logging.vlog(
          2, 'Found %s candidates in %s [%d bp, sample %s] '
          '[%0.2fs elapsed]', len(examples_dict[sample]),
          ranges.to_literal(region), ranges.length(region), sample,
          region_timer.Stop())
    return candidates_dict, examples_dict, gvcfs_dict

  def region_reads(self, region, sam_reader):
    """Update in_memory_sam_reader with read alignments overlapping the region.

    If self.options.realigner_enabled is set, uses realigned reads, otherwise
    original reads are returned.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to realign reads.
      sam_reader: sam.SamReader class for querying reads from SAM/BAM.

    Returns:
      [genomics.deepvariant.core.genomics.Read], reads overlapping the region.
    """
    if sam_reader is None:
      return []
    reads = sam_reader.query(region)
    if self.options.max_reads_per_partition > 0:
      random_for_region = np.random.RandomState(self.options.random_seed)
      reads = utils.reservoir_sample(reads,
                                     self.options.max_reads_per_partition,
                                     random_for_region)
    reads = list(reads)
    if self.options.realigner_enabled:
      max_read_length_to_realign = 500
      if max_read_length_to_realign > 0:
        long_reads = [
            read for read in reads
            if len(read.aligned_sequence) > max_read_length_to_realign
        ]

        short_reads = [
            read for read in reads
            if len(read.aligned_sequence) <= max_read_length_to_realign
        ]

        _, realigned_short_reads = self.realigner.realign_reads(
            short_reads, region)

        # Long reads will be listed before short reads when both are present.
        # Examples with only short or only long reads will be unaffected.
        return long_reads + realigned_short_reads

      _, reads = self.realigner.realign_reads(reads, region)
    return reads

  def candidates_in_region(self, region):
    """Finds candidate DeepVariantCall protos in region.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we want
        to get candidates for.

    Returns:
      A 2-tuple. The first value is a dict of deeptrio_pb2.DeepVariantCalls
      objects, in coordidate order. The second value is a dict of
      nucleus.genomics.v1.Variant protos containing gVCF information for all
      reference sites, if gvcf generation is enabled, otherwise returns [].
    """

    # Child is in the middle in self.samples.
    reads_parent1 = self.samples['child'][0].in_memory_sam_reader.query(region)
    reads = self.samples['child'][1].in_memory_sam_reader.query(region)
    reads_parent2 = self.samples['child'][2].in_memory_sam_reader.query(region)

    if not reads and not gvcf_output_enabled(self.options):
      # If we are generating gVCF output we cannot safely abort early here as
      # we need to return the gVCF records calculated by the caller below.
      return {}, {}

    allele_counters = {}

    # Instantiate allele counter for child sample
    allele_counters[self.options.variant_caller_options_child
                    .sample_name] = self._make_allele_counter_for_region(region)

    # Add target sample reads to allele counter
    for read in reads:
      allele_counters[
          self.options.variant_caller_options_child.sample_name].add(
              read, self.options.variant_caller_options_child.sample_name)

    # If parent1 is specified create allele counter and populate it with reads
    if FLAGS.reads_parent1:
      allele_counters[
          FLAGS.sample_name_parent1] = self._make_allele_counter_for_region(
              region)
      for read in reads_parent1:
        allele_counters[FLAGS.sample_name_parent1].add(
            read, FLAGS.sample_name_parent1)

    # If parent2 is specified create allele counter and populate it with reads
    if FLAGS.reads_parent2:
      allele_counters[
          FLAGS.sample_name_parent2] = self._make_allele_counter_for_region(
              region)
      for read in reads_parent2:
        allele_counters[FLAGS.sample_name_parent2].add(
            read, FLAGS.sample_name_parent2)

    candidates = {}
    gvcfs = {}

    if in_training_mode(self.options):
      candidates[self.sample_to_train], gvcfs[self.sample_to_train] = (
          self.variant_caller_child.calls_and_gvcfs(
              allele_counters, gvcf_output_enabled(self.options),
              FLAGS.sample_name_to_call))
      return candidates, gvcfs

    candidates['child'], gvcfs[
        'child'] = self.variant_caller_child.calls_and_gvcfs(
            allele_counters, gvcf_output_enabled(self.options),
            FLAGS.sample_name)
    if FLAGS.reads_parent1:
      candidates['parent1'], gvcfs['parent1'] = (
          self.variant_caller_parent1.calls_and_gvcfs(
              allele_counters, gvcf_output_enabled(self.options),
              FLAGS.sample_name_parent1))
    if FLAGS.reads_parent2:
      candidates['parent2'], gvcfs['parent2'] = (
          self.variant_caller_parent2.calls_and_gvcfs(
              allele_counters, gvcf_output_enabled(self.options),
              FLAGS.sample_name_parent2))
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
    # redacted
    window_width = self.pics['child'].width
    window_half_width = self.pics['child'].half_width

    alt_alleles = list(variant.alternate_bases)
    contig = variant.reference_name
    ref_start = variant.start
    ref_bases = variant.reference_bases
    ref_end = ref_start + len(ref_bases)

    # Sanity check that the reference_bases in the variant match the reference.
    ref_query_at_variant = self.realigner.ref_reader.query(
        ranges.make_range(contig, ref_start, ref_end))
    if ref_bases != ref_query_at_variant:
      raise ValueError('Error: reference_bases property in variant ({})'
                       'does not match the bases in the reference ({}) at that '
                       'position.'.format(ref_bases, ref_query_at_variant))

    # Margin must be equal to or more than half the window width.
    # Some extra prefix/suffix can be added to anchor alignments, but currently
    # we don't add extra.
    margin = window_half_width
    valid_end = min(
        self.realigner.ref_reader.contig(contig).n_bases, ref_end + margin)
    alignment_region = ranges.make_range(contig, max(ref_start - margin, 0),
                                         valid_end)
    trimmed_reads = [realigner.trim_read(r, alignment_region) for r in reads]
    # Filter reads to a minimum read length of 15 bp after trimming.
    reads = [r for r in trimmed_reads if len(r.aligned_sequence) >= 15]
    prefix = self.realigner.ref_reader.query(
        ranges.make_range(contig, max(ref_start - margin, 0), ref_start))
    suffix = self.realigner.ref_reader.query(
        ranges.make_range(contig, ref_end, valid_end))

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
          ref_start=ref_start - len(prefix))
      # Sequence of the alt haplotype in the window:
      end_of_prefix = prefix[-window_half_width:]
      beginning_of_suffix = suffix[:max(window_half_width + 1 - len(hap), 0)]
      sequences_by_haplotype[hap] = end_of_prefix + hap + beginning_of_suffix
      # Long haplotypes can extend past the window, so enforce the width here.
      sequences_by_haplotype[hap] = sequences_by_haplotype[hap][0:window_width]
    return {
        'alt_alignments': alignments_by_haplotype,
        'alt_sequences': sequences_by_haplotype
    }

  def create_pileup_examples(self, dv_call, samples_id):
    """Creates a tf.Example for DeepVariantCall.

    This function calls PileupImageCreator.create_pileup_images on dv_call to
    get raw image tensors for each alt_allele option (see docs for details).
    These tensors are encoded as pngs, and all of the key information is encoded
    as a tf.Example via a call to tf_utils.make_example.

    Args:
      dv_call: A DeepVariantCall.
      samples_id: string Sample id of the sample to call.

    Returns:
      A list of tf.Example protos.
    """
    reads_for_samples = [
        self.pics[samples_id].get_reads(
            dv_call.variant, sam_reader=sample.in_memory_sam_reader)
        for sample in self.samples[samples_id]
    ]
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

    pileup_images = self.pics[samples_id].create_pileup_images(
        dv_call=dv_call,
        reads_for_samples=reads_for_samples,
        haplotype_alignments_for_samples=haplotype_alignments_for_samples,
        haplotype_sequences=haplotype_sequences)

    if pileup_images is None:
      # We cannot build a PileupImage for dv_call, issue a warning.
      logging.warning('Could not create PileupImage for candidate at %s:%s',
                      dv_call.variant.reference_name, dv_call.variant.start)
      return []

    examples = []
    for alt_alleles, image_tensor in pileup_images:
      encoded_tensor, shape, tensor_format = self._encode_tensor(image_tensor)
      examples.append(
          tf_utils.make_example(
              dv_call.variant,
              alt_alleles,
              encoded_tensor,
              shape=shape,
              image_format=tensor_format,
              sequencing_type=deepvariant_pb2.PileupImageOptions.TRIO))
    return examples

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
    # Get our list of labels for each candidate variant.
    labels = self.labeler.label_variants(
        [candidate.variant for candidate in candidates], region)

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
      raise ValueError('Cannot add a non-confident label to an example',
                       example, label)
    alt_alleles_indices = tf_utils.example_alt_alleles_indices(example)

    tf_utils.example_set_variant(example, label.variant)

    # Set the label of the example to the # alts given our alt_alleles_indices.
    tf_utils.example_set_label(example,
                               label.label_for_alt_alleles(alt_alleles_indices))
    return example


def processing_regions_from_options(options):
  """Computes the calling regions from our options.

  This function does all of the work needed to read our input files and region
  specifications to determine the list of regions we should generate examples
  over. It also computes the confident regions needed to label variants.

  Args:
    options: deepvariant.DeepTrioOptions proto containing information about our
      input data sources.

  Raises:
    ValueError: if the regions to call is empty.

  Returns:
    Two values. The first is a list of nucleus.genomics.v1.Range protos of the
    regions we should process. The second is a RangeSet containing the confident
    regions for labeling, or None if we are running in training mode.
  """
  ref_contigs = fasta.IndexedFastaReader(
      options.reference_filename).header.contigs
  sam_contigs = sam.SamReader(options.reads_filename).header.contigs

  # Add in confident regions and vcf_contigs if in training mode.
  vcf_contigs = None
  if in_training_mode(options):
    vcf_contigs = vcf.VcfReader(options.truth_variants_filename).header.contigs

  contigs = _ensure_consistent_contigs(ref_contigs, sam_contigs, vcf_contigs,
                                       options.exclude_contigs,
                                       options.min_shared_contigs_basepairs)
  logging.info('Common contigs are %s', [c.name for c in contigs])
  calling_regions = build_calling_regions(ref_contigs, options.calling_regions,
                                          options.exclude_calling_regions)
  if not calling_regions:
    raise ValueError('The regions to call is empty. Check your --regions and '
                     '--exclude_regions flags to make sure they are not '
                     'resulting in set of empty region to process. This also '
                     'happens if you use "chr20" for a BAM where contig names '
                     'don\'t have "chr"s (or vice versa).')
  regions = regions_to_process(
      contigs=contigs,
      partition_size=options.allele_counter_options.partition_size,
      calling_regions=calling_regions,
      task_id=options.task_id,
      num_shards=options.num_shards)

  return regions


# redacted
class OutputsWriter(object):
  """Manages all of the outputs of make_examples in a single place."""

  def __init__(self, options, suffix=None):
    self._writers = {k: None for k in ['candidates', 'examples', 'gvcfs']}

    if options.candidates_filename:
      self._add_writer(
          'candidates',
          tfrecord.Writer(
              self._add_suffix(options.candidates_filename, suffix)))

    if options.examples_filename:
      self._add_writer(
          'examples',
          tfrecord.Writer(self._add_suffix(options.examples_filename, suffix)))

    if options.gvcf_filename:
      self._add_writer(
          'gvcfs',
          tfrecord.Writer(self._add_suffix(options.gvcf_filename, suffix)))

  def close_all(self):
    for writer in self._writers.values():
      if writer is not None:
        writer.close()

  def _add_suffix(self, file_path, suffix):
    """Adds suffix to file name."""
    if not suffix:
      return file_path

    file_dir, file_base = os.path.split(file_path)

    file_split = file_base.split('.')
    file_split[0] = '%s_%s' % (file_split[0], suffix)
    new_file_base = ('.').join(file_split)

    new_file = os.path.join(file_dir, new_file_base)
    return new_file

  def write_examples(self, *examples):
    self._write('examples', *examples)

  def write_gvcfs(self, *gvcfs):
    self._write('gvcfs', *gvcfs)

  def write_candidates(self, *candidates):
    self._write('candidates', *candidates)

  def _add_writer(self, name, writer):
    if name not in self._writers:
      raise ValueError(
          'Expected writer {} to have a None binding in writers.'.format(name))
    if self._writers[name] is not None:
      raise ValueError('Expected writer {} to be bound to None in writers but '
                       'saw {} instead'.format(name, self._writers[name]))
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


def get_example_counts(examples):
  """Returns a breakdown of examples by categories (label and type)."""
  labels = {0: 0, 1: 0, 2: 0}
  types = {
      tf_utils.EncodedVariantType.SNP: 0,
      tf_utils.EncodedVariantType.INDEL: 0,
      tf_utils.EncodedVariantType.UNKNOWN: 0
  }

  for example in examples:
    example_label = tf_utils.example_label(example)
    example_type = tf_utils.encoded_variant_type(
        tf_utils.example_variant(example))
    labels[example_label] += 1
    types[example_type] += 1
  return labels, types


def make_examples_runner(options):
  """Runs examples creation stage of deepvariant."""
  resource_monitor = resources.ResourceMonitor().start()
  logging.info('Preparing inputs')
  regions = processing_regions_from_options(options)

  # Create a processor to create candidates and examples for each region.
  region_processor = RegionProcessor(options)
  region_processor.initialize()

  logging.info('Writing examples to %s', options.examples_filename)
  if options.candidates_filename:
    logging.info('Writing candidates to %s', options.candidates_filename)
  if options.gvcf_filename:
    logging.info('Writing gvcf records to %s', options.gvcf_filename)

  n_regions, n_candidates, n_examples = 0, 0, 0
  n_class_0, n_class_1, n_class_2 = 0, 0, 0
  n_snps, n_indels = 0, 0
  last_reported = 0
  running_timer = timer.TimerStart()

  writers_dict = {}
  if not in_training_mode(options):
    for sample in region_processor.samples['child']:
      if sample.sam_reader is not None:
        # Only use suffix in calling mode
        suffix = None if in_training_mode(options) else sample.name
        writers_dict[sample.name] = OutputsWriter(options, suffix=suffix)
  else:
    writers_dict['child'] = OutputsWriter(options, suffix=None)

  for region in regions:
    candidates_dict, examples_dict, gvcfs_dict = region_processor.process(
        region)
    for sample in candidates_dict:
      candidates = candidates_dict[sample]
      examples = examples_dict[sample]
      gvcfs = gvcfs_dict[sample]

      writer = writers_dict[sample]

      n_candidates += len(candidates)
      n_examples += len(examples)
      n_regions += 1

      if in_training_mode(options) and options.run_info_filename:
        labels, types = get_example_counts(examples)
        n_class_0 += labels[0]
        n_class_1 += labels[1]
        n_class_2 += labels[2]
        n_snps += types[tf_utils.EncodedVariantType.SNP]
        n_indels += types[tf_utils.EncodedVariantType.INDEL]

      writer.write_candidates(*candidates)

      # If we have any gvcf records, write them out. This if also serves to
      # protect us from trying to write to the gvcfs output of writer when gvcf
      # generation is turned off. In that case, gvcfs will always be empty and
      # we'll never execute the write.
      if gvcfs:
        writer.write_gvcfs(*gvcfs)
      writer.write_examples(*examples)

      if (int(n_candidates / FLAGS.logging_every_n_candidates) > last_reported
          or n_regions == 1):
        last_reported = int(n_candidates / FLAGS.logging_every_n_candidates)
        logging.info('Task %s: %s candidates (%s examples) [%0.2fs elapsed]',
                     options.task_id, n_candidates, n_examples,
                     running_timer.Stop())
        running_timer = timer.TimerStart()

  for writer in writers_dict.values():
    writer.close_all()

  # Construct and then write out our MakeExamplesRunInfo proto.
  if options.run_info_filename:
    make_examples_stats = deepvariant_pb2.MakeExamplesStats(
        num_examples=n_examples,
        num_snps=n_snps,
        num_indels=n_indels,
        num_class_0=n_class_0,
        num_class_1=n_class_1,
        num_class_2=n_class_2)
    run_info = deeptrio_pb2.MakeExamplesRunInfo(
        options=options,
        resource_metrics=resource_monitor.metrics(),
        stats=make_examples_stats)
    if in_training_mode(options):
      if region_processor.labeler.metrics is not None:
        run_info.labeling_metrics.CopyFrom(region_processor.labeler.metrics)
      else:
        logging.warning(
            'Labeling metrics requested but the selected labeling '
            'algorithm %s does not collect metrics; skipping.',
            options.labeler_algorithm)
    logging.info('Writing MakeExamplesRunInfo to %s', options.run_info_filename)
    write_make_examples_run_info(run_info, path=options.run_info_filename)

  logging.info('Found %s candidate variants', n_candidates)
  logging.info('Created %s examples', n_examples)


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: make_examples does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)), errors.CommandLineError)
    del argv  # Unused.

    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()
    hts_verbose.set(hts_verbose.htsLogLevel[FLAGS.hts_logging_level])

    # Set up options; may do I/O.
    options = default_options(add_flags=True, flags_obj=FLAGS)

    # Check arguments that apply to any mode.
    if not options.reference_filename:
      errors.log_and_raise('ref argument is required.', errors.CommandLineError)
    if not options.reads_filename:
      errors.log_and_raise('reads argument is required.',
                           errors.CommandLineError)
    if not options.examples_filename:
      errors.log_and_raise('examples argument is required.',
                           errors.CommandLineError)
    if options.n_cores != 1:
      errors.log_and_raise(
          'Currently only supports n_cores == 1 but got {}.'.format(
              options.n_cores), errors.CommandLineError)

    if in_training_mode(options):
      if not options.truth_variants_filename:
        errors.log_and_raise(
            'truth_variants is required when in training mode.',
            errors.CommandLineError)
      if not options.confident_regions_filename:
        errors.log_and_raise(
            'confident_regions is required when in training mode.',
            errors.CommandLineError)
      if options.gvcf_filename:
        errors.log_and_raise('gvcf is not allowed in training mode.',
                             errors.CommandLineError)
    else:
      # Check for argument issues specific to calling mode.
      if (options.variant_caller_options_child.sample_name == _UNKNOWN_SAMPLE or
          (FLAGS.sample_name_parent1 and
           options.variant_caller_options_parent1.sample_name == _UNKNOWN_SAMPLE
          ) or
          (FLAGS.sample_name_parent2 and
           options.variant_caller_options_parent2.sample_name == _UNKNOWN_SAMPLE
          )):
        errors.log_and_raise(
            'sample_name must be specified for all samples in calling mode.',
            errors.CommandLineError)
      if options.variant_caller_options_child.gq_resolution < 1:
        errors.log_and_raise('gq_resolution must be a non-negative integer.',
                             errors.CommandLineError)

    # Sanity check the sample_names.
    if (options.variant_caller_options_child.sample_name ==
        FLAGS.sample_name_parent1 or
        options.variant_caller_options_child.sample_name ==
        FLAGS.sample_name_parent2):
      errors.log_and_raise(
          'The sample_name of the child is the same as one of '
          'the parents.', errors.CommandLineError)
    # If --sample_name_to_call is not set, use the child's sample_name.
    # This is for backward compatibility.
    if FLAGS.sample_name_to_call is None:
      FLAGS.sample_name_to_call = options.variant_caller_options_child.sample_name
    if FLAGS.sample_name_to_call not in [
        options.variant_caller_options_child.sample_name,
        FLAGS.sample_name_parent1, FLAGS.sample_name_parent2
    ]:
      errors.log_and_raise('--sample_name_to_call has to be one of the trio.',
                           errors.CommandLineError)

    if FLAGS.vsc_allele_fraction_trio_coefficient <= 0 or FLAGS.vsc_allele_fraction_trio_coefficient > 1.0:
      errors.log_and_raise(
          '--vsc_allele_fraction_trio_coefficient must within (0-1] internval.',
          errors.CommandLineError)
    # Run!
    make_examples_runner(options)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'mode',
      'reads',
      'ref',
  ])
  app.run(main)
