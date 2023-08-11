# Copyright 2021 Google LLC.
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
"""Shared flags and option handling for DeepVariant and DeepTrio."""

import re

from absl import flags
from absl import logging

from deepvariant import dv_constants
from deepvariant import exclude_contigs
from deepvariant import make_examples_core
from deepvariant import pileup_image
from deepvariant.protos import deepvariant_pb2
from deepvariant.realigner import realigner
from tensorflow.python.platform import gfile
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import errors

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
    'ref',
    None,
    (
        'Required. Genome reference to use. Must have an associated FAI index'
        ' as well. Supports text or gzipped references. Should match the'
        ' reference used to align the BAM file provided to --reads.'
    ),
)
flags.DEFINE_bool(
    'use_ref_for_cram',
    True,
    (
        'If true, use the --ref argument as the reference file for the CRAM'
        ' file passed to --reads.  In this case, it is required that the'
        ' reference file be located on a local POSIX filesystem. To disable,'
        ' specify --nouse_ref_for_cram.'
    ),
)
flags.DEFINE_string(
    'examples',
    None,
    'Required. Path to write tf.Example protos in TFRecord format.',
)
flags.DEFINE_string(
    'candidates',
    '',
    'Candidate DeepVariantCalls in tfrecord format. For DEBUGGING.',
)
flags.DEFINE_string(
    'mode',
    None,
    'Mode to run. Must be one of calling, training or candidate_sweep',
)
flags.DEFINE_string(
    'regions',
    '',
    (
        'Optional. Space-separated list of regions we want to process. Elements'
        ' can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE'
        ' files.'
    ),
)
flags.DEFINE_string(
    'exclude_regions',
    '',
    (
        'Optional. Space-separated list of regions we want to exclude from'
        ' processing. Elements can be region literals (e.g., chr20:10-20) or'
        ' paths to BED/BEDPE files. Region exclusion happens after processing'
        ' the --regions argument, so --region 20 --exclude_regions 20:100 does'
        ' everything on chromosome 20 excluding base 100'
    ),
)
flags.DEFINE_string(
    'variant_caller',
    'very_sensitive_caller',
    (
        'The caller to use to make examples. Must be one of the VariantCaller'
        ' enum values in the MakeExamplesOptions proto.'
    ),
)
flags.DEFINE_string(
    'gvcf',
    '',
    (
        'Optional. Path where we should write gVCF records in TFRecord of'
        ' Variant proto format.'
    ),
)
flags.DEFINE_integer(
    'gvcf_gq_binsize',
    5,
    (
        'Bin size in which to quantize gVCF genotype qualities. Larger bin size'
        ' reduces the number of gVCF records at a loss of quality granularity.'
        ' Must be a positive integer.'
    ),
)
flags.DEFINE_float(
    'p_error', 0.001, 'Basecalling error for reference confidence model.'
)
flags.DEFINE_bool(
    'include_med_dp',
    False,
    'If true, include MED_DP in the output gVCF records.',
)
flags.DEFINE_string(
    'confident_regions',
    '',
    (
        'Regions that we are confident are hom-ref or a variant in BED format.'
        ' In BED or other equivalent format, sorted or unsorted. Contig names'
        ' must match those of the reference genome.'
    ),
)
flags.DEFINE_string(
    'truth_variants',
    '',
    (
        'Tabix-indexed VCF file containing the truth variant calls for this'
        ' labels which we use to label our examples.'
    ),
)
flags.DEFINE_integer('task', 0, 'Task ID of this task')
flags.DEFINE_integer(
    'partition_size',
    1000,
    (
        'The maximum number of basepairs we will allow in a region before'
        ' splittingit into multiple smaller subregions.'
    ),
)
flags.DEFINE_integer(
    'max_reads_per_partition',
    1500,
    (
        'The maximum number of reads per partition that we consider before '
        'following processing such as sampling and realigner.'
    ),
)
_MAX_READS_FOR_DYNAMIC_BASES_PER_REGION = flags.DEFINE_integer(
    'max_reads_for_dynamic_bases_per_region',
    0,
    (
        'If > 0, set the max number of bases to '
        '(max_reads_for_dynamic_bases_per_region * region length).'
        'This is particularly important for very long reads.'
    ),
)
flags.DEFINE_string(
    'multi_allelic_mode',
    '',
    'How to handle multi-allelic candidate variants. For DEBUGGING',
)
flags.DEFINE_bool(
    'realign_reads',
    True,
    (
        'If True, locally realign reads before calling variants. '
        'Reads longer than 500 bp are never realigned.'
    ),
)
flags.DEFINE_bool(
    'write_run_info',
    False,
    (
        'If True, write out a MakeExamplesRunInfo proto besides our examples in'
        ' text_format.'
    ),
)
flags.DEFINE_enum(
    'alt_aligned_pileup',
    'none',
    ['none', 'base_channels', 'diff_channels', 'rows'],
    (
        'Include alignments of reads against each candidate alternate allele in'
        ' the pileup image. "none" turns this feature off. The default is'
        ' "none".Options: "none", "base_channels","diff_channels", "rows"'
    ),
)
flags.DEFINE_enum(
    'types_to_alt_align',
    'indels',
    ['indels', 'all'],
    (
        'When --alt_aligned_pileup is not none, this flag determines whether to'
        ' align to the alt alleles only for indels or for all variant types'
        ' including SNPs. Ignored if --alt_aligned_pileup is "none". This flag'
        ' is experimental and is not compatible with the pre-trained release'
        ' models.'
    ),
)
flags.DEFINE_string(
    'hts_logging_level',
    hts_verbose.htsLogLevel.HTS_LOG_WARNING.name,
    'Sets the htslib logging threshold.',
)
flags.DEFINE_integer(
    'hts_block_size',
    _DEFAULT_HTS_BLOCK_SIZE,
    (
        'Sets the htslib block size. Zero or negative uses default htslib'
        ' setting; larger values (e.g. 1M) may be beneficial for using remote'
        ' files. Currently only applies to SAM/BAM reading.'
    ),
)
flags.DEFINE_integer(
    'min_base_quality',
    10,
    (
        'Minimum base quality. This field indicates that we are enforcing a'
        ' minimum base quality score for alternate alleles. Alternate alleles'
        ' will only be considered if all bases in the allele have a quality'
        ' greater than min_base_quality.'
    ),
)
flags.DEFINE_integer(
    'min_mapping_quality',
    5,
    (
        'Setting this field to a positive integer i will only keep reads that'
        'have a MAPQ >= i. Note this only applies to aligned reads.'
    ),
)
flags.DEFINE_integer(
    'vsc_min_count_snps',
    2,
    (
        'SNP alleles occurring at least this many times in our '
        'AlleleCount will be advanced as candidates.'
    ),
)
flags.DEFINE_integer(
    'vsc_min_count_indels',
    2,
    (
        'Indel alleles occurring at least this many times in '
        'our AlleleCount will be advanced as candidates.'
    ),
)
flags.DEFINE_float(
    'vsc_min_fraction_snps',
    0.12,
    (
        'SNP alleles occurring at least this fraction of all '
        'counts in our AlleleCount will be advanced as '
        'candidates.'
    ),
)
flags.DEFINE_float(
    'vsc_min_fraction_indels',
    0.06,
    (
        'Indel alleles occurring at least this fraction of all counts in our '
        'AlleleCount will be advanced as candidates.'
    ),
)
_VSC_MIN_FRACTION_MULTIPLIER = flags.DEFINE_float(
    'vsc_min_fraction_multiplier',
    1.0,
    (
        'In candidate generation, this multiplier is applied to the minimum'
        ' allele fraction thresholds (vsc_min_fraction_snps and'
        ' vsc_min_fraction_indels) to adapt thresholds for multi-sample'
        ' calling. This has to in the (0, 1] interval. It can also be set to'
        ' float("inf") programmaticaly to only use candidates from the target'
        ' sample in multi-sample calling.'
    ),
)
_VSC_MAX_FRACTION_INDELS_FOR_NON_TARGET_SAMPLE = flags.DEFINE_float(
    'vsc_max_fraction_indels_for_non_target_sample',
    0.0,
    (
        'In candidate generation, if any non-target sample has more Indels '
        'than this threshold, the candidate will be excluded.'
        'Default is 0.0 which means no max is set.'
    ),
)
_VSC_MAX_FRACTION_SNPS_FOR_NON_TARGET_SAMPLE = flags.DEFINE_float(
    'vsc_max_fraction_snps_for_non_target_sample',
    0.0,
    (
        'In candidate generation, if any non-target sample has more SNPs than '
        'this threshold, the candidate will be excluded.'
        'Default is 0.0 which means no max is set.'
    ),
)
flags.DEFINE_float(
    'training_random_emit_ref_sites',
    NO_RANDOM_REF,
    'If > 0, emit extra random reference examples with this probability.',
)
flags.DEFINE_integer(
    'pileup_image_width',
    0,
    'Width for the pileup image. If 0, uses the default width',
)
flags.DEFINE_string(
    'labeler_algorithm',
    'haplotype_labeler',
    (
        'Algorithm to use to label examples in training mode. Must be one of'
        ' the LabelerAlgorithm enum values in the MakeExamplesOptions proto.'
    ),
)
flags.DEFINE_string(
    'customized_classes_labeler_classes_list',
    '',
    (
        'A comma-separated list of strings that defines customized class labels'
        ' for variants. This is only set when labeler_algorithm is'
        ' customized_classes_labeler.'
    ),
)
flags.DEFINE_string(
    'customized_classes_labeler_info_field_name',
    '',
    (
        'The name from the INFO field of VCF where we should get the customized'
        ' class labels from. This is only set when labeler_algorithm is'
        ' customized_classes_labeler.'
    ),
)
flags.DEFINE_integer(
    'logging_every_n_candidates',
    2000,
    (
        'Print out the log every n candidates. The smaller the number, the more'
        ' frequent the logging information emits.'
    ),
)
flags.DEFINE_bool('keep_duplicates', False, 'If True, keep duplicate reads.')
flags.DEFINE_bool(
    'keep_supplementary_alignments',
    False,
    'If True, keep reads marked as supplementary alignments.',
)
flags.DEFINE_bool(
    'keep_secondary_alignments',
    False,
    'If True, keep reads marked as secondary alignments.',
)
flags.DEFINE_bool(
    'parse_sam_aux_fields',
    None,
    (
        'If True, auxiliary fields of the SAM/BAM/CRAM records are parsed. By'
        ' default this flag is None. This flag will be automatically turned on'
        ' if other flags need it (e.g., sort_by_haplotypes). If it is'
        ' explicitly set by the user (either True or False), the user-specified'
        ' value will be used.'
    ),
)
flags.DEFINE_string(
    'aux_fields_to_keep',
    'HP,OQ',
    (
        'Comma-delimited list of auxiliary BAM fields to keep. '
        'This flag is used only when --parse_sam_aux_fields is '
        'set to true. If set to an empty string, all auxiliary '
        'fields will be kept.'
    ),
)
flags.DEFINE_bool(
    'use_original_quality_scores',
    False,
    'If True, base quality scores are read from OQ tag.',
)
flags.DEFINE_string(
    'select_variant_types',
    None,
    (
        'If provided, should be a whitespace-separated string of variant types'
        ' to keep when generating examples. Permitted values are "snps",'
        ' "indels", "multi-allelics", and "all", which select bi-allelic snps,'
        ' bi-allelic indels, multi-allelic variants of any type, and all'
        ' variants, respectively. Multiple selectors can be specified, so that'
        ' --select_variant_types="snps indels" would keep all bi-allelic SNPs'
        ' and indels'
    ),
)
flags.DEFINE_string(
    'sequencing_type',
    None,
    (
        'A string representing input bam file sequencing_type. Permitted values'
        ' are "WGS" and "WES", which represent whole genome sequencing and'
        ' whole exome sequencing, respectively. This flag is experimental and'
        ' is not currently being used.'
    ),
)
flags.DEFINE_bool(
    'sort_by_haplotypes',
    False,
    (
        'If True, reads are sorted by haplotypes (using HP tag), '
        'parse_sam_aux_fields has to be set for this to work.'
    ),
)
flags.DEFINE_bool(
    'reverse_haplotypes',
    False,
    (
        'If True, reads are sorted by haplotypes (using HP tag) in reverse'
        ' order, parse_sam_aux_fields has to be set for this to work.'
    ),
)
flags.DEFINE_integer(
    'hp_tag_for_assembly_polishing',
    0,
    (
        'If set to > 0, reads with this HP tag will be sorted on top. '
        'sort_by_haplotypes has to be set to True for this to work.'
    ),
)
flags.DEFINE_bool(
    'add_hp_channel',
    False,
    'If true, add another channel to represent HP tags per read.',
)
flags.DEFINE_string(
    'channels',
    None,
    'Comma or space-delimited list of optional channels to add. '
    'Available channels: {}'.format(','.join(dv_constants.OPT_CHANNELS)),
)
flags.DEFINE_bool(
    'add_supporting_other_alt_color',
    False,
    (
        'If True, reads supporting an alt not represented in the '
        'pileup image are colored differently for multiallelics.'
    ),
)
flags.DEFINE_string(
    'population_vcfs',
    None,
    (
        'Optional. Tabix-indexed VCF file (or list of VCFs broken by'
        ' chromosome), separated by comma or space, containing population'
        ' allele frequencies. Each of the item can be a file path, or a'
        ' wildcard pattern.'
    ),
)
flags.DEFINE_bool(
    'use_allele_frequency',
    False,
    (
        'If True, add another channel for pileup images to represent allele '
        'frequency information gathered from population callsets.'
    ),
)
flags.DEFINE_string(
    'runtime_by_region',
    None,
    (
        '[optional] Output filename for a TSV file of runtimes and other stats'
        ' by region. If examples are sharded, this should be sharded into the'
        ' same number of shards as the examples.'
    ),
)
flags.DEFINE_bool(
    'track_ref_reads',
    False,
    (
        'If True, allele counter keeps track of ref supporting reads.'
        'By default allele counter keeps a simple count of number of reads '
        'supporting ref.'
    ),
)
flags.DEFINE_bool(
    'normalize_reads',
    False,
    'If True, allele counter left align INDELs for each read.',
)
flags.DEFINE_bool(
    'keep_legacy_allele_counter_behavior',
    False,
    (
        'If True, the behavior in this commit is reverted: '
        'https://github.com/google/deepvariant/commit/'
        'fbde0674639a28cb9e8004c7a01bbe25240c7d46. '
        'We do not recommend setting this flag to True.'
    ),
)
flags.DEFINE_bool(
    'phase_reads',
    False,
    'Calculate phases and add HP tag to all reads on a fly.',
)
flags.DEFINE_integer(
    'phase_max_candidates',
    5000,
    (
        'Limits the number of candidates for phasing. If number of candidates '
        'exceeds the maximum then phasing is not performed for the window. '
        'This flag is used only when phase_reads is true.'
    ),
)
_ENABLE_JOINT_REALIGNMENT = flags.DEFINE_bool(
    'enable_joint_realignment',
    False,
    (
        'If True, realign reads from all samples together. By default this is '
        'False, which means reads from each sample are realigned per-sample.'
    ),
)
_OUTPUT_LOCAL_READ_PHASING = flags.DEFINE_string(
    'output_local_read_phasing',
    None,
    (
        '[optional] For debugging only. Output filename for a TSV file '
        'containing read phases. If examples are sharded, this should be '
        'sharded into the same number of shards as the examples.'
    ),
)
_DISCARD_NON_DNA_REGIONS = flags.DEFINE_bool(
    'discard_non_dna_regions',
    False,
    (
        'Default is False. If set regions of Ns larger than 300,000bp are'
        'discarded.'
    ),
)

_OUTPUT_SITELIST = flags.DEFINE_bool(
    'output_sitelist',
    False,
    'If True, output a list of sites present in examples output.',
)

_DENOVO_REGIONS = flags.DEFINE_string(
    'denovo_regions',
    '',
    'Regions where variants are de novo. Used to label variants as de novo.',
)


def shared_flags_to_options(
    add_flags,
    flags_obj,
    samples_in_order,
    sample_role_to_train,
    main_sample_index,
) -> deepvariant_pb2.MakeExamplesOptions:
  """Creates options from flags that are shared, along with given samples."""
  read_reqs = reads_pb2.ReadRequirements(
      keep_duplicates=flags_obj.keep_duplicates,
      keep_supplementary_alignments=flags_obj.keep_supplementary_alignments,
      keep_secondary_alignments=flags_obj.keep_secondary_alignments,
      min_base_quality=flags_obj.min_base_quality,
      min_mapping_quality=flags_obj.min_mapping_quality,
      min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT,
  )

  logging.vlog(3, 'ReadRequirements are: %s', read_reqs)

  pic_options = pileup_image.default_options(read_requirements=read_reqs)

  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=flags_obj.partition_size,
      read_requirements=read_reqs,
      track_ref_reads=flags_obj.track_ref_reads,
      normalize_reads=flags_obj.normalize_reads,
      keep_legacy_behavior=flags_obj.keep_legacy_allele_counter_behavior,
  )

  options = deepvariant_pb2.MakeExamplesOptions(
      exclude_contigs=exclude_contigs.EXCLUDED_HUMAN_CONTIGS,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=609314161,
      # # Not specified by default: calling_regions = 3;
      read_requirements=read_reqs,
      allele_counter_options=allele_counter_options,
      pic_options=pic_options,
      n_cores=1,
      task_id=0,
      num_shards=0,
      min_shared_contigs_basepairs=0.9,
      sample_options=samples_in_order,
      main_sample_index=main_sample_index,
      sample_role_to_train=sample_role_to_train,
      output_sitelist=_OUTPUT_SITELIST.value,
  )

  if add_flags:
    options.mode = make_examples_core.parse_proto_enum_flag(
        deepvariant_pb2.MakeExamplesOptions.Mode, flags_obj.mode.upper()
    )

    options.labeler_algorithm = make_examples_core.parse_proto_enum_flag(
        deepvariant_pb2.MakeExamplesOptions.LabelerAlgorithm,
        flags_obj.labeler_algorithm.upper(),
    )

    options.variant_caller = make_examples_core.parse_proto_enum_flag(
        deepvariant_pb2.MakeExamplesOptions.VariantCaller,
        flags_obj.variant_caller.upper(),
    )

    if flags_obj.ref:
      options.reference_filename = flags_obj.ref
    if flags_obj.confident_regions:
      options.confident_regions_filename = flags_obj.confident_regions
    if flags_obj.denovo_regions:
      options.denovo_regions_filename = _DENOVO_REGIONS.value
    if flags_obj.truth_variants:
      options.truth_variants_filename = flags_obj.truth_variants
    if flags_obj.sequencing_type:
      options.pic_options.sequencing_type = (
          make_examples_core.parse_proto_enum_flag(
              deepvariant_pb2.PileupImageOptions.SequencingType,
              flags_obj.sequencing_type,
          )
      )

    if flags_obj.channels:
      channel_set = re.split('[, ]+', flags_obj.channels)
      for channel in channel_set:
        if channel and channel not in dv_constants.OPT_CHANNELS:
          err_msg = (
              'Channel "{}" is not one of the available opt channels: {}'
              .format(channel, ', '.join(dv_constants.OPT_CHANNELS))
          )
          errors.log_and_raise(err_msg, errors.CommandLineError)
      options.pic_options.channels[:] = channel_set
      options.pic_options.num_channels += len(channel_set)

    if flags_obj.multi_allelic_mode:
      multi_allelic_enum = {
          'include_het_alt_images': (
              deepvariant_pb2.PileupImageOptions.ADD_HET_ALT_IMAGES
          ),
          'exclude_het_alt_images': (
              deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES
          ),
      }[flags_obj.multi_allelic_mode]
      options.pic_options.multi_allelic_mode = multi_allelic_enum

    if flags_obj.pileup_image_width:
      options.pic_options.width = flags_obj.pileup_image_width

    # DirectPhasing related flags.
    if flags_obj.phase_reads:
      options.phase_reads = flags_obj.phase_reads
    phase_region_padding = dv_constants.PHASE_READS_REGION_PADDING_PCT
    if phase_region_padding:
      options.phase_reads_region_padding_pct = phase_region_padding
    if flags_obj.phase_max_candidates:
      options.phase_max_candidates = flags_obj.phase_max_candidates

    options.pic_options.alt_aligned_pileup = flags_obj.alt_aligned_pileup
    options.pic_options.types_to_alt_align = flags_obj.types_to_alt_align


    if flags_obj.add_supporting_other_alt_color:
      options.pic_options.other_allele_supporting_read_alpha = 0.3

    if flags_obj.select_variant_types:
      options.select_variant_types[:] = flags_obj.select_variant_types.split()
      for svt in options.select_variant_types:
        if svt not in make_examples_core.VARIANT_TYPE_SELECTORS:
          errors.log_and_raise(
              'Select variant type {} not recognized. Allowed values are {}'
              .format(
                  svt, ', '.join(make_examples_core.VARIANT_TYPE_SELECTORS)
              ),
              errors.CommandLineError,
          )

    (
        num_shards,
        examples,
        candidates,
        gvcf,
        runtime_by_region,
        read_phases_output,
    ) = sharded_file_utils.resolve_filespecs(
        flags_obj.task,
        flags_obj.examples or '',
        flags_obj.candidates or '',
        flags_obj.gvcf or '',
        flags_obj.runtime_by_region or '',
        flags_obj.output_local_read_phasing or '',
    )
    options.examples_filename = examples
    options.candidates_filename = candidates
    options.gvcf_filename = gvcf
    options.include_med_dp = flags_obj.include_med_dp
    options.task_id = flags_obj.task
    options.num_shards = num_shards
    options.runtime_by_region = runtime_by_region
    options.read_phases_output = read_phases_output

    options.parse_sam_aux_fields = make_examples_core.resolve_sam_aux_fields(
        flags_obj=flags_obj
    )
    if flags_obj.aux_fields_to_keep:
      options.aux_fields_to_keep[:] = flags_obj.aux_fields_to_keep.split(',')
    else:
      options.aux_fields_to_keep = None
    options.use_original_quality_scores = flags_obj.use_original_quality_scores

    if flags_obj.add_hp_channel:
      options.pic_options.num_channels += 1
      options.pic_options.add_hp_channel = True

    if flags_obj.hp_tag_for_assembly_polishing < 0:
      errors.log_and_raise(
          '--hp_tag_for_assembly_polishing has to be set to a positive int.',
          errors.CommandLineError,
      )
    if (
        flags_obj.hp_tag_for_assembly_polishing > 0
        and not flags_obj.sort_by_haplotypes
    ):
      errors.log_and_raise(
          (
              '--hp_tag_for_assembly_polishing requires --sort_by_haplotypes to'
              ' be set '
          ),
          errors.CommandLineError,
      )

    options.pic_options.sort_by_haplotypes = flags_obj.sort_by_haplotypes
    options.pic_options.reverse_haplotypes = flags_obj.reverse_haplotypes
    options.pic_options.hp_tag_for_assembly_polishing = (
        flags_obj.hp_tag_for_assembly_polishing
    )

    if flags_obj.write_run_info:
      options.run_info_filename = examples + _RUN_INFO_FILE_EXTENSION

    options.calling_regions.extend(
        make_examples_core.parse_regions_flag(flags_obj.regions)
    )
    options.exclude_calling_regions.extend(
        make_examples_core.parse_regions_flag(flags_obj.exclude_regions)
    )

    options.realigner_enabled = flags_obj.realign_reads
    options.joint_realignment = _ENABLE_JOINT_REALIGNMENT.value
    options.realigner_options.CopyFrom(realigner.realigner_config(flags_obj))

    if (
        options.mode == deepvariant_pb2.MakeExamplesOptions.TRAINING
        and flags_obj.training_random_emit_ref_sites != NO_RANDOM_REF
    ):
      options.sample_options[
          main_sample_index
      ].variant_caller_options.fraction_reference_sites_to_emit = (
          flags_obj.training_random_emit_ref_sites
      )

    if flags_obj.use_allele_frequency and not flags_obj.population_vcfs:
      errors.log_and_raise(
          (
              'If use_allele_frequency is set then population_vcfs '
              'must be provided.'
          ),
          errors.CommandLineError,
      )
    if flags_obj.use_allele_frequency:
      options.use_allele_frequency = flags_obj.use_allele_frequency
      options.pic_options.num_channels += 1
      options.pic_options.use_allele_frequency = True
    if flags_obj.population_vcfs:
      for path in re.split(',| ', flags_obj.population_vcfs):
        options.population_vcf_filenames.extend(gfile.Glob(path))
    options.max_reads_per_partition = flags_obj.max_reads_per_partition
    options.max_reads_for_dynamic_bases_per_region = (
        _MAX_READS_FOR_DYNAMIC_BASES_PER_REGION.value
    )
    options.use_ref_for_cram = flags_obj.use_ref_for_cram
    options.hts_block_size = flags_obj.hts_block_size
    options.logging_every_n_candidates = flags_obj.logging_every_n_candidates
    options.customized_classes_labeler_classes_list = (
        flags_obj.customized_classes_labeler_classes_list
    )
    options.customized_classes_labeler_info_field_name = (
        flags_obj.customized_classes_labeler_info_field_name
    )

  options.discard_non_dna_regions = _DISCARD_NON_DNA_REGIONS.value

  return options


def check_options_are_valid(
    options: deepvariant_pb2.MakeExamplesOptions, main_sample_index: int
):
  """Checks that all the options chosen make sense together."""

  # Check arguments that apply to any mode.
  if not options.reference_filename:
    errors.log_and_raise('ref argument is required.', errors.CommandLineError)

  if not options.examples_filename:
    errors.log_and_raise(
        'examples argument is required.', errors.CommandLineError
    )
  if options.n_cores != 1:
    errors.log_and_raise(
        'Currently only supports n_cores == 1 but got {}.'.format(
            options.n_cores
        ),
        errors.CommandLineError,
    )

  main_sample = options.sample_options[main_sample_index]
  if not main_sample.reads_filenames:
    errors.log_and_raise('reads argument is required.', errors.CommandLineError)

  if make_examples_core.in_candidate_sweep_mode(options):
    # In candidate_sweep mode there is nothing to check here.
    pass
  elif make_examples_core.in_training_mode(options):
    if not options.truth_variants_filename:
      errors.log_and_raise(
          'truth_variants is required when in training mode.',
          errors.CommandLineError,
      )
    if not options.confident_regions_filename:
      if (
          options.variant_caller
          == deepvariant_pb2.MakeExamplesOptions.VCF_CANDIDATE_IMPORTER
      ):
        logging.info(
            'Note: --confident_regions is optional with '
            'vcf_candidate_importer. '
            'You did not specify --confident_regions, which means '
            'examples will be generated for the whole region.'
        )
      else:
        errors.log_and_raise(
            'confident_regions is required when in training mode.',
            errors.CommandLineError,
        )
    if options.gvcf_filename:
      errors.log_and_raise(
          'gvcf is not allowed in training mode.', errors.CommandLineError
      )
    if (
        options.variant_caller
        == deepvariant_pb2.MakeExamplesOptions.VCF_CANDIDATE_IMPORTER
        and main_sample.proposed_variants_filename
    ):
      errors.log_and_raise(
          (
              '--proposed_variants* should not be used with '
              'vcf_candidate_importer in training mode. '
              'Use --truth_variants to pass in the candidates '
              'with correct labels for training.'
          ),
          errors.CommandLineError,
      )
  else:
    # Check for argument issues specific to calling mode.
    for sample in options.sample_options:
      # If there are reads, there must be a sample name too.
      if sample.reads_filenames:
        if sample.variant_caller_options.sample_name == _UNKNOWN_SAMPLE:
          errors.log_and_raise(
              'sample_name must be specified for all samples in calling mode.',
              errors.CommandLineError,
          )
    if main_sample.variant_caller_options.gq_resolution < 1:
      errors.log_and_raise(
          'gq_resolution must be a positive integer.', errors.CommandLineError
      )

    if options.truth_variants_filename:
      errors.log_and_raise(
          'Do not specify --truth_variants in calling mode.',
          errors.CommandLineError,
      )

    if (
        options.variant_caller
        == deepvariant_pb2.MakeExamplesOptions.VCF_CANDIDATE_IMPORTER
    ):
      if any(
          o.proposed_variants_filename is None for o in options.sample_options
      ):
        errors.log_and_raise(
            (
                '--proposed_variants* is required with vcf_candidate_importer'
                ' in calling mode.'
            ),
            errors.CommandLineError,
        )

  multiplier = _VSC_MIN_FRACTION_MULTIPLIER.value
  if (multiplier <= 0 or multiplier > 1.0) and multiplier != float('inf'):
    errors.log_and_raise(
        '--vsc_min_fraction_multiplier must be within (0-1] interval, or set '
        'inf to only use candidates from the target sample. '
        'Currently set to: {}'.format(multiplier),
        errors.CommandLineError,
    )

  total_pileup_height = sum(
      [sample.pileup_height for sample in options.sample_options]
  )
  # Height constraint for Slim InceptionV3 implementation.
  if total_pileup_height < 75 or total_pileup_height > 362:
    errors.log_and_raise('Total pileup image heights must be between 75-362.')
