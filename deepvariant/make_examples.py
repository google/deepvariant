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
"""Step one of DeepVariant: creates tf.Example protos for training/calling."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import re


from absl import app
from absl import flags
from absl import logging

from deepvariant import dv_constants
from deepvariant import exclude_contigs
from deepvariant import logging_level
from deepvariant import make_examples_core
from deepvariant import make_examples_utils
from deepvariant import pileup_image
from deepvariant.protos import deepvariant_pb2
from deepvariant.realigner import realigner
from third_party.nucleus.io import sam
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils

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
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref. '
    'Can provide multiple BAMs (comma-separated).')
flags.DEFINE_bool(
    'use_ref_for_cram', True,
    'If true, use the --ref argument as the reference file for the CRAM '
    'file passed to --reads.  In this case, it is required that the reference '
    'file be located on a local POSIX filesystem. To disable, specify '
    '--nouse_ref_for_cram.')
flags.DEFINE_string(
    'examples', None,
    'Required. Path to write tf.Example protos in TFRecord format.')
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
    'values in the DeepVariantOptions proto.')
flags.DEFINE_string(
    'gvcf', '',
    'Optional. Path where we should write gVCF records in TFRecord of Variant '
    'proto format.')
flags.DEFINE_integer(
    'gvcf_gq_binsize', 5,
    'Bin size in which to quantize gVCF genotype qualities. Larger bin size '
    'reduces the number of gVCF records at a loss of quality granularity.')
flags.DEFINE_bool('include_med_dp', False,
                  'If true, include MED_DP in the output gVCF records.')
flags.DEFINE_string(
    'confident_regions', '',
    'Regions that we are confident are hom-ref or a variant in BED format. In '
    'BED or other equivalent format, sorted or unsorted. Contig names must '
    'match those of the reference genome.')
flags.DEFINE_string(
    'truth_variants', '',
    'Tabix-indexed VCF file containing the truth variant calls for this labels '
    'which we use to label our examples.')
flags.DEFINE_string(
    'proposed_variants', '',
    '(Only used when --variant_caller=vcf_candidate_importer.) '
    'Tabix-indexed VCF file containing the proposed positions and alts for '
    '`vcf_candidate_importer`. The GTs will be ignored.')
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
flags.DEFINE_enum(
    'alt_aligned_pileup', 'none',
    ['none', 'base_channels', 'diff_channels', 'rows'],
    'Include alignments of reads against each candidate alternate allele in '
    'the pileup image. "none" turns this feature off. '
    'The default is "none".'
    'Options: "none", "base_channels","diff_channels", "rows"')
flags.DEFINE_enum(
    'types_to_alt_align', 'indels', ['indels', 'all'],
    'When --alt_aligned_pileup is not none, this flag determines whether to '
    'align to the alt alleles only for indels or for all variant types '
    'including SNPs. Ignored if --alt_aligned_pileup is "none". This flag is '
    'experimental and is not compatible with the pre-trained release models.')
flags.DEFINE_float(
    'downsample_fraction', NO_DOWNSAMPLING,
    'If not ' + str(NO_DOWNSAMPLING) + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input BAM. This argument makes it easy to create examples as '
    'though the input BAM had less coverage.')
flags.DEFINE_string(
    'sample_name', '', 'Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads.')
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
    'Setting this field to a positive integer i will only keep reads that'
    'have a MAPQ >= i. Note this only applies to aligned reads.')
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
    'Indel alleles occurring at least this fraction of all counts in our '
    'AlleleCount will be advanced as candidates.')
flags.DEFINE_float(
    'vsc_min_fraction_multiplier', 1.0,
    'In candidate generation, this multiplier is applied to the minimum allele '
    'fraction thresholds (vsc_min_fraction_snps and vsc_min_fraction_indels) '
    'to adapt thresholds for multi-sample calling.')
flags.DEFINE_float(
    'training_random_emit_ref_sites', NO_RANDOM_REF,
    'If > 0, emit extra random reference examples with this probability.')
flags.DEFINE_integer(
    'pileup_image_height', 0,
    'Height for the pileup image. If 0, uses the default height')
flags.DEFINE_integer(
    'pileup_image_width', 0,
    'Width for the pileup image. If 0, uses the default width')
flags.DEFINE_string(
    'labeler_algorithm', 'haplotype_labeler',
    'Algorithm to use to label examples in training mode. Must be one of the '
    'LabelerAlgorithm enum values in the DeepVariantOptions proto.')
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
    'parse_sam_aux_fields', None,
    'If True, auxiliary fields of the SAM/BAM/CRAM records are parsed. '
    'By default this flag is None. This flag will be automatically turned on '
    'if other flags need it (e.g., sort_by_haplotypes). '
    'If it is explicitly set by the user (either True or False), the '
    'user-specified value will be used.')
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
flags.DEFINE_integer(
    'hp_tag_for_assembly_polishing', 0,
    'If set to > 0, reads with this HP tag will be sorted on top. '
    'sort_by_haplotypes has to be set to True for this to work.')
flags.DEFINE_bool(
    'add_hp_channel', False,
    'If true, add another channel to represent HP tags per read.')
flags.DEFINE_string(
    'channels', None, 'Comma-delimited list of optional channels to add. '
    'Available channels: {}'.format(','.join(dv_constants.OPT_CHANNELS)))
flags.DEFINE_bool(
    'add_supporting_other_alt_color', False,
    'If True, reads supporting an alt not represented in the '
    'pileup image are colored differently for multiallelics.')
flags.DEFINE_string(
    'population_vcfs', None,
    'Optional. Tabix-indexed VCF file (or list of VCFs broken by chromosome),'
    ' separated by comma or space, '
    'containing population allele frequencies.')
flags.DEFINE_bool(
    'use_allele_frequency', False,
    'If True, add another channel for pileup images to represent allele '
    'frequency information gathered from population callsets.')
flags.DEFINE_string(
    'runtime_by_region', None,
    '[optional] Output filename for a TSV file of runtimes and '
    'other stats by region. If examples are sharded, this should be sharded '
    'into the same number of shards as the examples.')


def default_options(add_flags=True, flags_obj=None):
  """Creates a DeepVariantOptions proto populated with reasonable defaults.

  Args:
    add_flags: bool. defaults to True. If True, we will push the value of
      certain FLAGS into our options. If False, those option fields are left
      uninitialized.
    flags_obj: object.  If not None, use as the source of flags, else use global
      FLAGS.

  Returns:
    deepvariant_pb2.DeepVariantOptions protobuf.

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

  logging.vlog(3, 'ReadRequirements are: %s', read_reqs)

  pic_options = pileup_image.default_options(read_requirements=read_reqs)

  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=flags_obj.partition_size, read_requirements=read_reqs)

  if flags_obj.sample_name:
    sample_name = flags_obj.sample_name
  elif flags_obj.reads:
    # If there are multiple BAM files, use the sample name from the first one.
    with sam.SamReader(flags_obj.reads.split(',')[0]) as sam_reader:
      sample_name = make_examples_core.extract_sample_name_from_sam_reader(
          sam_reader)
  else:
    sample_name = _UNKNOWN_SAMPLE

  variant_caller_options = deepvariant_pb2.VariantCallerOptions(
      min_count_snps=flags_obj.vsc_min_count_snps,
      min_count_indels=flags_obj.vsc_min_count_indels,
      min_fraction_snps=flags_obj.vsc_min_fraction_snps,
      min_fraction_indels=flags_obj.vsc_min_fraction_indels,
      min_fraction_multiplier=flags_obj.vsc_min_fraction_multiplier,
      # Not specified by default: fraction_reference_sites_to_emit,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=1400605801,
      sample_name=sample_name,
      p_error=0.001,
      max_gq=50,
      gq_resolution=flags_obj.gvcf_gq_binsize,
      ploidy=2)

  options = deepvariant_pb2.DeepVariantOptions(
      exclude_contigs=exclude_contigs.EXCLUDED_HUMAN_CONTIGS,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=609314161,
      # # Not specified by default: calling_regions = 3;
      read_requirements=read_reqs,
      allele_counter_options=allele_counter_options,
      variant_caller_options=variant_caller_options,
      pic_options=pic_options,
      n_cores=1,
      task_id=0,
      num_shards=0,
      min_shared_contigs_basepairs=0.9)

  if add_flags:
    options.mode = make_examples_core.parse_proto_enum_flag(
        deepvariant_pb2.DeepVariantOptions.Mode, flags_obj.mode.upper())

    options.labeler_algorithm = make_examples_core.parse_proto_enum_flag(
        deepvariant_pb2.DeepVariantOptions.LabelerAlgorithm,
        flags_obj.labeler_algorithm.upper())

    options.variant_caller = make_examples_core.parse_proto_enum_flag(
        deepvariant_pb2.DeepVariantOptions.VariantCaller,
        flags_obj.variant_caller.upper())

    if flags_obj.ref:
      options.reference_filename = flags_obj.ref
    if flags_obj.reads:
      options.reads_filenames.extend(flags_obj.reads.split(','))
    if flags_obj.confident_regions:
      options.confident_regions_filename = flags_obj.confident_regions
    if flags_obj.truth_variants:
      options.truth_variants_filename = flags_obj.truth_variants
    if flags_obj.proposed_variants:
      options.proposed_variants_filename = flags_obj.proposed_variants
    if flags_obj.sequencing_type:
      options.pic_options.sequencing_type = make_examples_core.parse_proto_enum_flag(
          deepvariant_pb2.PileupImageOptions.SequencingType,
          flags_obj.sequencing_type)
    if flags_obj.downsample_fraction != NO_DOWNSAMPLING:
      options.downsample_fraction = flags_obj.downsample_fraction

    if flags_obj.channels:
      channel_set = flags_obj.channels.split(',')
      for channel in channel_set:
        if channel and channel not in dv_constants.OPT_CHANNELS:
          err_msg = 'Channel "{}" is not one of the available opt channels: {}'.format(
              channel, ', '.join(dv_constants.OPT_CHANNELS))
          errors.log_and_raise(err_msg, errors.CommandLineError)
      options.pic_options.channels[:] = channel_set
      options.pic_options.num_channels += len(channel_set)

    if flags_obj.multi_allelic_mode:
      multi_allelic_enum = {
          'include_het_alt_images':
              deepvariant_pb2.PileupImageOptions.ADD_HET_ALT_IMAGES,
          'exclude_het_alt_images':
              deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES,
      }[flags_obj.multi_allelic_mode]
      options.pic_options.multi_allelic_mode = multi_allelic_enum

    if flags_obj.pileup_image_height:
      options.pic_options.height = flags_obj.pileup_image_height
    if flags_obj.pileup_image_width:
      options.pic_options.width = flags_obj.pileup_image_width

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
              .format(svt,
                      ', '.join(make_examples_core.VARIANT_TYPE_SELECTORS)),
              errors.CommandLineError)

    num_shards, examples, candidates, gvcf, runtime_by_region = (
        sharded_file_utils.resolve_filespecs(flags_obj.task,
                                             flags_obj.examples or '',
                                             flags_obj.candidates or '',
                                             flags_obj.gvcf or '',
                                             flags_obj.runtime_by_region or ''))
    options.examples_filename = examples
    options.candidates_filename = candidates
    options.gvcf_filename = gvcf
    options.include_med_dp = flags_obj.include_med_dp
    options.task_id = flags_obj.task
    options.num_shards = num_shards
    options.runtime_by_region = runtime_by_region

    # First, if parse_sam_aux_fields is still None (which means the users didn't
    # set it. We added some logic to decide its value.
    if flags_obj.parse_sam_aux_fields is None:
      flags_that_needs_sam_aux_fields = [
          'add_hp_channel', 'sort_by_haplotypes', 'use_original_quality_scores'
      ]
      flags_obj.parse_sam_aux_fields = make_examples_core.set_parse_sam_aux_fields(
          flags_obj, flags_that_needs_sam_aux_fields)

    for flags_strictly_needs_sam_aux_fields in [
        'sort_by_haplotypes', 'use_original_quality_scores'
    ]:
      if (flags_obj[flags_strictly_needs_sam_aux_fields].value and
          not flags_obj.parse_sam_aux_fields):
        errors.log_and_raise(
            'If --{} is set then parse_sam_aux_fields '
            'must be set too.'.format(flags_strictly_needs_sam_aux_fields),
            errors.CommandLineError)
    options.use_original_quality_scores = flags_obj.use_original_quality_scores

    for flag_optionally_needs_sam_aux_fields in ['add_hp_channel']:
      if (flags_obj[flag_optionally_needs_sam_aux_fields].value and
          not flags_obj.parse_sam_aux_fields):
        logging.warning(
            'WARGNING! --%s is set but --parse_sam_aux_fields is not '
            'set. This will cause aux fields to not be read in. The relevant '
            'values might be zero. For example, for --add_hp_channel, '
            'resulting in an empty HP channel. If this is not what you '
            'intended, please stop and enable --parse_sam_aux_fields.',
            flag_optionally_needs_sam_aux_fields)
    if flags_obj.add_hp_channel:
      options.pic_options.num_channels += 1
      options.pic_options.add_hp_channel = True

    if flags_obj.hp_tag_for_assembly_polishing < 0:
      errors.log_and_raise(
          '--hp_tag_for_assembly_polishing has to be set to a positive int.',
          errors.CommandLineError)
    if (flags_obj.hp_tag_for_assembly_polishing > 0 and
        not flags_obj.sort_by_haplotypes):
      errors.log_and_raise(
          '--hp_tag_for_assembly_polishing requires --sort_by_haplotypes to be '
          'set ', errors.CommandLineError)
    if flags_obj.sort_by_haplotypes and not flags_obj.parse_sam_aux_fields:
      errors.log_and_raise(
          '--sort_by_haplotypes requires --parse_sam_aux_fields to be set ',
          errors.CommandLineError)
    options.pic_options.sort_by_haplotypes = flags_obj.sort_by_haplotypes
    options.pic_options.hp_tag_for_assembly_polishing = flags_obj.hp_tag_for_assembly_polishing

    if flags_obj.write_run_info:
      options.run_info_filename = examples + _RUN_INFO_FILE_EXTENSION

    options.calling_regions.extend(
        make_examples_core.parse_regions_flag(flags_obj.regions))
    options.exclude_calling_regions.extend(
        make_examples_core.parse_regions_flag(flags_obj.exclude_regions))

    options.realigner_enabled = flags_obj.realign_reads
    options.realigner_options.CopyFrom(realigner.realigner_config(flags_obj))

    if (options.mode == deepvariant_pb2.DeepVariantOptions.TRAINING and
        flags_obj.training_random_emit_ref_sites != NO_RANDOM_REF):
      options.variant_caller_options.fraction_reference_sites_to_emit = (
          flags_obj.training_random_emit_ref_sites)

    if (flags_obj.use_allele_frequency and not flags_obj.population_vcfs):
      errors.log_and_raise(
          'If use_allele_frequency is set then population_vcfs '
          'must be provided.', errors.CommandLineError)
    if flags_obj.use_allele_frequency:
      options.use_allele_frequency = flags_obj.use_allele_frequency
      options.pic_options.num_channels += 1
      options.pic_options.use_allele_frequency = True
    if flags_obj.population_vcfs:
      options.population_vcf_filenames.extend(
          re.split(',| ', flags_obj.population_vcfs))

  options.max_reads_per_partition = flags_obj.max_reads_per_partition
  options.use_ref_for_cram = flags_obj.use_ref_for_cram
  options.hts_block_size = flags_obj.hts_block_size
  options.logging_every_n_candidates = flags_obj.logging_every_n_candidates
  options.customized_classes_labeler_classes_list = flags_obj.customized_classes_labeler_classes_list
  options.customized_classes_labeler_info_field_name = flags_obj.customized_classes_labeler_info_field_name

  if flags_obj.parse_sam_aux_fields is not None:
    options.parse_sam_aux_fields = flags_obj.parse_sam_aux_fields

  options.main_sample_index = 0
  return options


def check_options_are_valid(options):
  """Check that DeepVariant options all fit together."""

  # Check arguments that apply to any mode.
  if not options.reference_filename:
    errors.log_and_raise('ref argument is required.', errors.CommandLineError)
  if not options.reads_filenames:
    errors.log_and_raise('reads argument is required.', errors.CommandLineError)
  if not options.examples_filename:
    errors.log_and_raise('examples argument is required.',
                         errors.CommandLineError)
  if options.n_cores != 1:
    errors.log_and_raise(
        'Currently only supports n_cores == 1 but got {}.'.format(
            options.n_cores), errors.CommandLineError)

  # Check for argument issues specific to different modes.
  if make_examples_core.in_training_mode(options):
    if not options.truth_variants_filename:
      errors.log_and_raise('truth_variants is required when in training mode.',
                           errors.CommandLineError)
    if not options.confident_regions_filename:
      if options.variant_caller == \
          deepvariant_pb2.DeepVariantOptions.VCF_CANDIDATE_IMPORTER:
        logging.info('Note: --confident_regions is optional with '
                     'vcf_candidate_importer. '
                     'You did not specify --confident_regions, which means '
                     'examples will be generated for the whole region.')
      else:
        errors.log_and_raise(
            'confident_regions is required when in training mode.',
            errors.CommandLineError)
    if options.gvcf_filename:
      errors.log_and_raise('gvcf is not allowed in training mode.',
                           errors.CommandLineError)
    if (options.variant_caller == \
        deepvariant_pb2.DeepVariantOptions.VCF_CANDIDATE_IMPORTER and
        options.proposed_variants_filename):
      errors.log_and_raise(
          '--proposed_variants should not be used with '
          'vcf_candidate_importer in training mode. '
          'Use --truth_variants to pass in the candidates '
          'with correct labels for training.', errors.CommandLineError)
  else:
    # Check for argument issues specific to calling mode.
    if options.truth_variants_filename:
      errors.log_and_raise('Do not specify --truth_variants in calling mode.',
                           errors.CommandLineError)
    if options.variant_caller_options.sample_name == _UNKNOWN_SAMPLE:
      errors.log_and_raise('sample_name must be specified in calling mode.',
                           errors.CommandLineError)
    if options.variant_caller_options.gq_resolution < 1:
      errors.log_and_raise('gq_resolution must be a non-negative integer.',
                           errors.CommandLineError)
    if options.variant_caller == \
        deepvariant_pb2.DeepVariantOptions.VCF_CANDIDATE_IMPORTER:
      if not options.proposed_variants_filename:
        errors.log_and_raise(
            '--proposed_variants is required with vcf_candidate_importer in '
            'calling mode.', errors.CommandLineError)


def samples_from_options(options):
  """Create an array of one sample from the options given."""
  return [
      make_examples_utils.Sample(
          name=options.variant_caller_options.sample_name,
          nickname='main_sample',
          reads_filenames=options.reads_filenames,
          variant_caller_options=options.variant_caller_options)
  ]


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
    check_options_are_valid(options)

    # Define samples. This is an array of just the one sample for DeepVariant.
    samples = samples_from_options(options)

    # Run!
    make_examples_core.make_examples_runner(options, samples=samples)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'mode',
      'reads',
      'ref',
  ])
  app.run(main)
