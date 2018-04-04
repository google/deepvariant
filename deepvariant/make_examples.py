# Copyright 2017 Google Inc.
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



from tensorflow import flags
import numpy as np
import tensorflow as tf

from absl import logging

from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.io import vcf
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import io_utils
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils
from deepvariant import exclude_contigs
from deepvariant import logging_level
from deepvariant import pileup_image
from deepvariant import tf_utils
from deepvariant import variant_caller
from deepvariant.labeler import haplotype_labeler
from deepvariant.labeler import positional_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.realigner import realigner
from deepvariant.vendor import timer

FLAGS = flags.FLAGS

# Sentinel command line flag value indicating no downsampling should occur.
NO_DOWNSAMPLING = 0.0

# Sentinel command line flag value indicating no random ref sites should be
# emitted.
NO_RANDOM_REF = 0.0

# The name used for a sample if one is not specified or present in the reads.
_UNKNOWN_SAMPLE = 'UNKNOWN'

# Use a default hts_block_size value of 128 MB (see b/69330994 for details) to
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
    'to call. Should be aligned to a reference genome compatible with --ref.')
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
flags.DEFINE_bool('realign_reads', True,
                  'If True, locally realign reads before calling variants.')
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
flags.DEFINE_integer('vsc_min_count_snps', 2,
                     'SNP alleles occurring at least this many times in our '
                     'AlleleCount will be advanced as candidates.')
flags.DEFINE_integer('vsc_min_count_indels', 2,
                     'Indel alleles occurring at least this many times in '
                     'our AlleleCount will be advanced as candidates.')
flags.DEFINE_float('vsc_min_fraction_snps', 0.12,
                   'SNP alleles occurring at least this fraction of all '
                   'counts in our AlleleCount will be advanced as '
                   'candidates.')
flags.DEFINE_float('vsc_min_fraction_indels', 0.12,
                   'Indel alleles occurring at least this fraction of all '
                   'counts in our AlleleCount will be advanced as '
                   'candidates.')
flags.DEFINE_float(
    'training_random_emit_ref_sites', NO_RANDOM_REF,
    'If > 0, emit extra random reference examples with this probability.')
flags.DEFINE_integer(
    'pileup_image_height', 0,
    'Height for the pileup image. If 0, uses the default height')
flags.DEFINE_integer('pileup_image_width', 0,
                     'Width for the pileup image. If 0, uses the default width')
flags.DEFINE_string(
    'labeler_algorithm', 'haplotype_labeler',
    'Algorithm to use to label examples in training mode. Must be one of the '
    'LabelerAlgorithm enum values in the DeepVariantOptions proto.')


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
      deepvariant_pb2.DeepVariantOptions.Mode to get the DeepVariantOptions Mode
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


def default_options(add_flags=True, flags_obj=None):
  """Creates a DeepVariantOptions proto populated with reasonable defaults.

  Args:
    add_flags: bool. defaults to True. If True, we will push the value of
      certain FLAGS into our options. If False, those option fields are left
      uninitialized.
    flags_obj: object.  If not None, use as the source of flags,
      else use global FLAGS.

  Returns:
    deepvariant_pb2.DeepVariantOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  if not flags_obj:
    flags_obj = FLAGS

  read_reqs = reads_pb2.ReadRequirements(
      min_base_quality=10,
      min_mapping_quality=10,
      min_base_quality_mode=reads_pb2.ReadRequirements.ENFORCED_BY_CLIENT)

  pic_options = pileup_image.default_options(read_requirements=read_reqs)

  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=flags_obj.partition_size, read_requirements=read_reqs)

  if flags_obj.sample_name:
    sample_name = flags_obj.sample_name
  elif flags_obj.reads:
    with sam.SamReader(flags_obj.reads) as sam_reader:
      sample_name = extract_sample_name_from_sam_reader(sam_reader)
  else:
    sample_name = _UNKNOWN_SAMPLE

  variant_caller_options = deepvariant_pb2.VariantCallerOptions(
      min_count_snps=flags_obj.vsc_min_count_snps,
      min_count_indels=flags_obj.vsc_min_count_indels,
      min_fraction_snps=flags_obj.vsc_min_fraction_snps,
      min_fraction_indels=flags_obj.vsc_min_fraction_indels,
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
      min_shared_contigs_basepairs=0.9,
  )

  if add_flags:
    options.mode = parse_proto_enum_flag(
        deepvariant_pb2.DeepVariantOptions.Mode, flags_obj.mode.upper())

    options.labeler_algorithm = parse_proto_enum_flag(
        deepvariant_pb2.DeepVariantOptions.LabelerAlgorithm,
        flags_obj.labeler_algorithm.upper())

    if flags_obj.ref:
      options.reference_filename = flags_obj.ref
    if flags_obj.reads:
      options.reads_filename = flags_obj.reads
    if flags_obj.confident_regions:
      options.confident_regions_filename = flags_obj.confident_regions
    if flags_obj.truth_variants:
      options.truth_variants_filename = flags_obj.truth_variants

    if flags_obj.downsample_fraction != NO_DOWNSAMPLING:
      options.downsample_fraction = flags_obj.downsample_fraction

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

    num_shards, examples, candidates, gvcf = io_utils.resolve_filespecs(
        flags_obj.task, flags_obj.examples or '', flags_obj.candidates or '',
        flags_obj.gvcf or '')
    options.examples_filename = examples
    options.candidates_filename = candidates
    options.gvcf_filename = gvcf

    options.calling_regions.extend(parse_regions_flag(flags_obj.regions))
    options.exclude_calling_regions.extend(
        parse_regions_flag(flags_obj.exclude_regions))

    options.task_id = flags_obj.task
    options.num_shards = 0 if num_shards is None else num_shards

    options.realigner_enabled = flags_obj.realign_reads
    if options.realigner_enabled:
      options.realigner_options.CopyFrom(realigner.realigner_config(flags_obj))

    options.max_reads_per_partition = flags_obj.max_reads_per_partition

    if (options.mode == deepvariant_pb2.DeepVariantOptions.TRAINING and
        flags_obj.training_random_emit_ref_sites != NO_RANDOM_REF):
      options.variant_caller_options.fraction_reference_sites_to_emit = (
          flags_obj.training_random_emit_ref_sites)

  return options


# ---------------------------------------------------------------------------
# Simple utilities
# ---------------------------------------------------------------------------


def in_training_mode(options):
  return options.mode == deepvariant_pb2.DeepVariantOptions.TRAINING


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
      rg.sample_id
      for rg in sam_reader.header.read_groups
      if rg.sample_id
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


def _set_variant_genotype(variant, genotype):
  if not variant.calls:
    variant.calls.add(genotype=genotype)
  else:
    variant.calls[0].genotype[:] = genotype


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
  # redacted
  partitioned = regions.partition(partition_size)
  partitioned = ranges.sorted_ranges(partitioned, contigs)

  if num_shards:
    return (r for i, r in enumerate(partitioned) if i % num_shards == task_id)
  else:
    return partitioned


# ---------------------------------------------------------------------------
# Region processor
# ---------------------------------------------------------------------------


def read_confident_regions(options):
  return ranges.RangeSet.from_bed(options.confident_regions_filename)


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
      options: deepvariant.DeepVariantOptions proto used to specify our
        resources for calling (e.g., reference_filename).
    """
    self.options = options
    self.initialized = False
    self.ref_reader = None
    self.sam_reader = None
    self.in_memory_sam_reader = None
    self.realigner = None
    self.pic = None
    self.labeler = None
    self.variant_caller = None

  def _make_allele_counter_for_region(self, region):
    return allelecounter.AlleleCounter(self.ref_reader.get_c_reader(), region,
                                       self.options.allele_counter_options)

  def _encode_tensor(self, image_tensor):
    return image_tensor.tostring(), image_tensor.shape, 'raw'

  def _make_sam_reader(self):
    return sam.SamReader(
        self.options.reads_filename,
        read_requirements=self.options.read_requirements,
        hts_block_size=FLAGS.hts_block_size,
        downsample_fraction=self.options.downsample_fraction,
        random_seed=self.options.random_seed)

  def _initialize(self):
    """Initialize the resources needed for this work in the current env."""
    if self.initialized:
      raise ValueError('Cannot initialize this object twice')

    self.ref_reader = fasta.RefFastaReader(self.options.reference_filename)
    self.sam_reader = self._make_sam_reader()
    self.in_memory_sam_reader = sam.InMemorySamReader([])

    if self.options.realigner_enabled:
      self.realigner = realigner.Realigner(self.options.realigner_options,
                                           self.ref_reader)
    self.pic = pileup_image.PileupImageCreator(
        ref_reader=self.ref_reader,
        sam_reader=self.in_memory_sam_reader,
        options=self.options.pic_options)

    if in_training_mode(self.options):
      self.labeler = self._make_labeler_from_options()

    self.variant_caller = variant_caller.VariantCaller(
        self.options.variant_caller_options)
    self.random = np.random.RandomState(self.options.random_seed)
    self.initialized = True

  def _make_labeler_from_options(self):
    truth_vcf_reader = vcf.VcfReader(
        self.options.truth_variants_filename,
        excluded_format_fields=['GL', 'GQ', 'PL'])
    confident_regions = read_confident_regions(self.options)

    if (self.options.labeler_algorithm ==
        deepvariant_pb2.DeepVariantOptions.POSITIONAL_LABELER):
      return positional_labeler.PositionalVariantLabeler(
          truth_vcf_reader=truth_vcf_reader,
          confident_regions=confident_regions)
    elif (self.options.labeler_algorithm ==
          deepvariant_pb2.DeepVariantOptions.HAPLOTYPE_LABELER):
      return haplotype_labeler.HaplotypeLabeler(
          truth_vcf_reader=truth_vcf_reader,
          ref_reader=self.ref_reader,
          confident_regions=confident_regions)
    else:
      raise ValueError('Unexpected labeler_algorithm',
                       self.options.labeler_algorithm)

  def process(self, region):
    """Finds candidates and creates corresponding examples in a region.

    Args:
      region: A nucleus.genomics.v1.Range proto. Specifies the region on the
        genome we should process.

    Returns:
      Three values. First is a list of the found candidates, which are
      deepvariant.DeepVariantCall objects. The second value is a list of filled
      in tf.Example protos. For example, these will include the candidate
      variant, the pileup image, and, if in training mode, the truth variants
      and labels needed for training. The third value is a list of
      nucleus.genomics.v1.Variant protos containing gVCF information for all
      reference sites, if gvcf generation is enabled, otherwise returns [].
    """
    region_timer = timer.TimerStart()

    # Print some basic information about what we are doing.
    if not self.initialized:
      self._initialize()

    self.in_memory_sam_reader.replace_reads(self.region_reads(region))
    candidates, gvcfs = self.candidates_in_region(region)

    if in_training_mode(self.options):
      examples = [
          self.add_label_to_example(example, label)
          for candidate, label in self.label_candidates(candidates)
          for example in self.create_pileup_examples(candidate)
      ]
    else:
      examples = [
          example for candidate in candidates
          for example in self.create_pileup_examples(candidate)
      ]

    logging.info('Found %s candidates in %s [%0.2fs elapsed]', len(examples),
                 ranges.to_literal(region), region_timer.Stop())
    return candidates, examples, gvcfs

  def region_reads(self, region):
    """Update in_memory_sam_reader with read alignments overlapping the region.

    If self.realigner is set, uses realigned reads, otherwise original reads
    are returned.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we
        want to realign reads.

    Returns:
      [genomics.deepvariant.core.genomics.Read], reads overlapping the region.
    """
    reads = self.sam_reader.query(region)
    if self.options.max_reads_per_partition > 0:
      reads = utils.reservoir_sample(
          reads, self.options.max_reads_per_partition, self.random)
    reads = list(reads)
    if self.realigner:
      _, reads = self.realigner.realign_reads(reads, region)
    return reads

  def candidates_in_region(self, region):
    """Finds candidate DeepVariantCall protos in region.

    Args:
      region: A nucleus.genomics.v1.Range object specifying the region we
      want to get candidates for.

    Returns:
      A 2-tuple. The first value is a list of deepvariant_pb2.DeepVariantCalls
      objects, in coordidate order. The second value is a list of
      nucleus.genomics.v1.Variant protos containing gVCF information for all
      reference sites, if gvcf generation is enabled, otherwise returns [].
    """
    reads = self.in_memory_sam_reader.query(region)
    if not reads and not gvcf_output_enabled(self.options):
      # If we are generating gVCF output we cannot safely abort early here as
      # we need to return the gVCF records calculated by the caller below.
      return [], []

    allele_counter = self._make_allele_counter_for_region(region)
    for read in reads:
      allele_counter.add(read)

    candidates, gvcfs = self.variant_caller.calls_from_allele_counter(
        allele_counter, gvcf_output_enabled(self.options))
    return candidates, gvcfs

  def create_pileup_examples(self, dv_call):
    """Creates a tf.Example for DeepVariantCall.

    This function calls PileupImageCreator.create_pileup_images on dv_call to
    get raw image tensors for each alt_allele option (see docs for details).
    These tensors are encoded as pngs, and all of the key information is encoded
    as a tf.Example via a call to tf_utils.make_example.

    Args:
      dv_call: A DeepVariantCall.

    Returns:
      A list of tf.Example protos.
    """
    pileup_images = self.pic.create_pileup_images(dv_call)
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
              image_format=tensor_format))
    return examples

  def label_candidates(self, candidates):
    """Gets label information for each candidate.

    Args:
      candidates: list[DeepVariantCalls]: The list of candidate variant calls we
        want to label.

    Yields:
      Tuples of (candidate, label_variants.Label objects) for each candidate in
      candidates that could be assigned a label. Candidates that couldn't be
      labeled will not be returned.
    """
    # Get our list of labels for each candidate variant.
    labels = self.labeler.label_variants(
        [candidate.variant for candidate in candidates])

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

    # Set the genotype of the candidate variant to the labeled value.
    candidate = label.variant
    _set_variant_genotype(candidate, label.genotype)
    tf_utils.example_set_variant(example, candidate)

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
    options: deepvariant.DeepVariantOptions proto containing information about
      our input data sources.

  Returns:
    Two values. The first is a list of nucleus.genomics.v1.Range protos of the
    regions we should process. The second is a RangeSet containing the confident
    regions for labeling, or None if we are running in training mode.
  """
  ref_contigs = fasta.RefFastaReader(options.reference_filename).header.contigs
  sam_contigs = sam.SamReader(options.reads_filename).header.contigs

  # Add in confident regions and vcf_contigs if in training mode.
  vcf_contigs = None
  if in_training_mode(options):
    vcf_contigs = vcf.VcfReader(options.truth_variants_filename).header.contigs

  contigs = _ensure_consistent_contigs(ref_contigs, sam_contigs, vcf_contigs,
                                       options.exclude_contigs,
                                       options.min_shared_contigs_basepairs)
  logging.info('Common contigs are %s', [c.name for c in contigs])

  regions = regions_to_process(
      contigs=contigs,
      partition_size=options.allele_counter_options.partition_size,
      calling_regions=build_calling_regions(ref_contigs,
                                            options.calling_regions,
                                            options.exclude_calling_regions),
      task_id=options.task_id,
      num_shards=options.num_shards)

  return regions


# redacted
class OutputsWriter(object):
  """Manages all of the outputs of make_examples in a single place."""

  def __init__(self, options):
    self._writers = {k: None for k in ['candidates', 'examples', 'gvcfs']}

    if options.candidates_filename:
      self._add_writer('candidates',
                       io_utils.RawProtoWriterAdaptor(
                           io_utils.make_tfrecord_writer(
                               options.candidates_filename)))

    if options.examples_filename:
      self._add_writer('examples',
                       io_utils.RawProtoWriterAdaptor(
                           io_utils.make_tfrecord_writer(
                               options.examples_filename)))

    if options.gvcf_filename:
      self._add_writer('gvcfs',
                       io_utils.RawProtoWriterAdaptor(
                           io_utils.make_tfrecord_writer(
                               options.gvcf_filename)))

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
    for writer in self._writers.itervalues():
      if writer is not None:
        writer.__enter__()
    return self

  def __exit__(self, exception_type, exception_value, traceback):
    for writer in self._writers.itervalues():
      if writer is not None:
        writer.__exit__(exception_type, exception_value, traceback)

  def _write(self, writer_name, *protos):
    writer = self._writers[writer_name]
    if writer:
      for proto in protos:
        writer.write(proto)


def make_examples_runner(options):
  """Runs examples creation stage of deepvariant."""
  logging.info('Preparing inputs')
  regions = processing_regions_from_options(options)

  # Create a processor to create candidates and examples for each region.
  region_processor = RegionProcessor(options)

  logging.info('Writing examples to %s', options.examples_filename)
  if options.candidates_filename:
    logging.info('Writing candidates to %s', options.candidates_filename)
  if options.gvcf_filename:
    logging.info('Writing gvcf records to %s', options.gvcf_filename)

  n_regions, n_candidates = 0, 0
  with OutputsWriter(options) as writer:
    for region in regions:
      candidates, examples, gvcfs = region_processor.process(region)
      n_candidates += len(candidates)
      n_regions += 1

      writer.write_candidates(*candidates)

      # If we have any gvcf records, write them out. This if also serves to
      # protect us from trying to write to the gvcfs output of writer when gvcf
      # generation is turned off. In that case, gvcfs will always be empty and
      # we'll never execute the write.
      if gvcfs:
        writer.write_gvcfs(*gvcfs)
      writer.write_examples(*examples)

  logging.info('Found %s candidate variants', n_candidates)


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

    # Check for argument issues specific to train mode.
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
      if options.variant_caller_options.sample_name == _UNKNOWN_SAMPLE:
        errors.log_and_raise('sample_name must be specified in calling mode.',
                             errors.CommandLineError)
      if options.variant_caller_options.gq_resolution < 1:
        errors.log_and_raise('gq_resolution must be a non-negative integer.',
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
  tf.app.run()
