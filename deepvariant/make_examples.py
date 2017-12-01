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



import numpy as np
import tensorflow as tf

from absl import logging

from deepvariant import logging_level
from deepvariant import pileup_image
from deepvariant import tf_utils
from deepvariant import variant_caller
from deepvariant import variant_labeler
from deepvariant.core import errors
from deepvariant.core import genomics_io
from deepvariant.core import htslib_gcp_oauth
from deepvariant.core import io_utils
from deepvariant.core import proto_utils
from deepvariant.core import ranges
from deepvariant.core import utils
from deepvariant.core import variantutils
from deepvariant.core.protos import core_pb2
from deepvariant.core.python import hts_verbose
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.realigner import realigner
from deepvariant.vendor import timer

FLAGS = tf.flags.FLAGS

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

tf.flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference to use. Must have an associated FAI index as '
    'well. Supports text or gzipped references. Should match the reference '
    'used to align the BAM file provided to --reads.')
tf.flags.DEFINE_string(
    'reads', None,
    'Required. Aligned, sorted, indexed BAM file containing the reads we want '
    'to call. Should be aligned to a reference genome compatible with --ref.')
tf.flags.DEFINE_string(
    'examples', None,
    'Required. Path to write tf.Example protos in TFRecord format.')
tf.flags.DEFINE_string(
    'candidates', '',
    'Candidate DeepVariantCalls in tfrecord format. For DEBUGGING.')
tf.flags.DEFINE_string('mode', None,
                       'Mode to run. Must be one of calling or training')
tf.flags.DEFINE_string(
    'regions', '',
    'Optional. Space-separated list of regions we want to process. Elements '
    'can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.')
tf.flags.DEFINE_string(
    'gvcf', '',
    'Optional. Path where we should write gVCF records in TFRecord of Variant '
    'proto format.')
tf.flags.DEFINE_string(
    'confident_regions', '',
    'Regions that we are confident are hom-ref or a variant in BED format. In '
    'BED or other equivalent format, sorted or unsorted. Contig names must '
    'match those of the reference genome.')
tf.flags.DEFINE_string(
    'truth_variants', '',
    'Tabix-indexed VCF file containing the truth variant calls for this labels '
    'which we use to label our examples.')
tf.flags.DEFINE_integer('task', 0, 'Task ID of this task')
tf.flags.DEFINE_integer(
    'partition_size', 1000,
    'The maximum number of basepairs we will allow in a region before splitting'
    'it into multiple smaller subregions.')
tf.flags.DEFINE_integer(
    'max_reads_per_partition', 1500,
    'The maximum number of reads per partition that we consider before '
    'following processing such as sampling and realigner.')
tf.flags.DEFINE_string(
    'multi_allelic_mode', '',
    'How to handle multi-allelic candidate variants. For DEBUGGING')
tf.flags.DEFINE_bool('realign_reads', True,
                     'If True, locally realign reads before calling variants.')
tf.flags.DEFINE_float(
    'downsample_fraction', NO_DOWNSAMPLING,
    'If not ' + str(NO_DOWNSAMPLING) + ' must be a value between 0.0 and 1.0. '
    'Reads will be kept (randomly) with a probability of downsample_fraction '
    'from the input BAM. This argument makes it easy to create examples as '
    'though the input BAM had less coverage.')
tf.flags.DEFINE_string(
    'sample_name', '', 'Sample name to use for our sample_name in the output '
    'Variant/DeepVariantCall protos. If not specified, will be inferred from '
    'the header information from --reads.')
tf.flags.DEFINE_string('hts_logging_level',
                       hts_verbose.htsLogLevel.HTS_LOG_WARNING.name,
                       'Sets the htslib logging threshold.')
tf.flags.DEFINE_integer(
    'hts_block_size', _DEFAULT_HTS_BLOCK_SIZE,
    'Sets the htslib block size. Zero or negative uses default htslib setting; '
    'larger values (e.g. 1M) may be beneficial for using remote files. '
    'Currently only applies to SAM/BAM reading.')
tf.flags.DEFINE_integer('vsc_min_count_snps', 2,
                        'SNP alleles occurring at least this many times in our '
                        'AlleleCount will be advanced as candidates.')
tf.flags.DEFINE_integer('vsc_min_count_indels', 2,
                        'Indel alleles occurring at least this many times in '
                        'our AlleleCount will be advanced as candidates.')
tf.flags.DEFINE_float('vsc_min_fraction_snps', 0.12,
                      'SNP alleles occurring at least this fraction of all '
                      'counts in our AlleleCount will be advanced as '
                      'candidates.')
tf.flags.DEFINE_float('vsc_min_fraction_indels', 0.12,
                      'Indel alleles occurring at least this fraction of all '
                      'counts in our AlleleCount will be advanced as '
                      'candidates.')
tf.flags.DEFINE_float(
    'training_random_emit_ref_sites', NO_RANDOM_REF,
    'If > 0, emit extra random reference examples with this probability.')
tf.flags.DEFINE_integer(
    'pileup_image_height', 0,
    'Height for the pileup image. If 0, uses the default height')
tf.flags.DEFINE_integer(
    'pileup_image_width', 0,
    'Width for the pileup image. If 0, uses the default width')

# ---------------------------------------------------------------------------
# Option handling
# ---------------------------------------------------------------------------


def default_options(add_flags=True, flags=None):
  """Creates a DeepVariantOptions proto populated with reasonable defaults.

  Args:
    add_flags: bool. defaults to True. If True, we will push the value of
      certain FLAGS into our options. If False, those option fields are left
      uninitialized.
    flags: object.  If not None, use as the source of flags,
      else use global FLAGS.

  Returns:
    deepvariant_pb2.DeepVariantOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  if not flags:
    flags = FLAGS

  read_reqs = core_pb2.ReadRequirements(
      min_base_quality=10,
      min_mapping_quality=10,
      min_base_quality_mode=core_pb2.ReadRequirements.ENFORCED_BY_CLIENT)

  pic_options = pileup_image.default_options(read_requirements=read_reqs)

  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      partition_size=flags.partition_size, read_requirements=read_reqs)

  if flags.sample_name:
    sample_name = flags.sample_name
  elif flags.reads:
    sample_name = extract_sample_name_from_reads(flags.reads)
  else:
    sample_name = _UNKNOWN_SAMPLE

  variant_caller_options = deepvariant_pb2.VariantCallerOptions(
      min_count_snps=flags.vsc_min_count_snps,
      min_count_indels=flags.vsc_min_count_indels,
      min_fraction_snps=flags.vsc_min_fraction_snps,
      min_fraction_indels=flags.vsc_min_fraction_indels,
      # Not specified by default: fraction_reference_sites_to_emit,
      # Fixed random seed produced with 'od -vAn -N4 -tu4 < /dev/urandom'.
      random_seed=1400605801,
      sample_name=sample_name,
      p_error=0.001,
      max_gq=50,
      gq_resolution=1,
      ploidy=2)

  options = deepvariant_pb2.DeepVariantOptions(
      exclude_contigs=[
          # The two canonical names for the contig representing the human
          # mitochondrial sequence.
          'chrM',
          'MT',
          # From hs37d5.
          # (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707)  # pylint:disable=line-too-long
          'GL000207.1',
          'GL000226.1',
          'GL000229.1',
          'GL000231.1',
          'GL000210.1',
          'GL000239.1',
          'GL000235.1',
          'GL000201.1',
          'GL000247.1',
          'GL000245.1',
          'GL000197.1',
          'GL000203.1',
          'GL000246.1',
          'GL000249.1',
          'GL000196.1',
          'GL000248.1',
          'GL000244.1',
          'GL000238.1',
          'GL000202.1',
          'GL000234.1',
          'GL000232.1',
          'GL000206.1',
          'GL000240.1',
          'GL000236.1',
          'GL000241.1',
          'GL000243.1',
          'GL000242.1',
          'GL000230.1',
          'GL000237.1',
          'GL000233.1',
          'GL000204.1',
          'GL000198.1',
          'GL000208.1',
          'GL000191.1',
          'GL000227.1',
          'GL000228.1',
          'GL000214.1',
          'GL000221.1',
          'GL000209.1',
          'GL000218.1',
          'GL000220.1',
          'GL000213.1',
          'GL000211.1',
          'GL000199.1',
          'GL000217.1',
          'GL000216.1',
          'GL000215.1',
          'GL000205.1',
          'GL000219.1',
          'GL000224.1',
          'GL000223.1',
          'GL000195.1',
          'GL000212.1',
          'GL000222.1',
          'GL000200.1',
          'GL000193.1',
          'GL000194.1',
          'GL000225.1',
          'GL000192.1',
          'NC_007605',
          'hs37d5',
      ],
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
    if flags.mode == 'training':
      options.mode = deepvariant_pb2.DeepVariantOptions.TRAINING
    elif flags.mode == 'calling':
      options.mode = deepvariant_pb2.DeepVariantOptions.CALLING
    else:
      raise ValueError('Unexpected mode', flags.mode)

    if flags.ref:
      options.reference_filename = flags.ref
    if flags.reads:
      options.reads_filename = flags.reads
    if flags.confident_regions:
      options.confident_regions_filename = flags.confident_regions
    if flags.truth_variants:
      options.truth_variants_filename = flags.truth_variants

    if flags.downsample_fraction != NO_DOWNSAMPLING:
      options.downsample_fraction = flags.downsample_fraction

    if flags.multi_allelic_mode:
      multi_allelic_enum = {
          'include_het_alt_images':
              deepvariant_pb2.PileupImageOptions.ADD_HET_ALT_IMAGES,
          'exclude_het_alt_images':
              deepvariant_pb2.PileupImageOptions.NO_HET_ALT_IMAGES,
      }[flags.multi_allelic_mode]
      options.pic_options.multi_allelic_mode = multi_allelic_enum

    if flags.pileup_image_height:
      options.pic_options.height = flags.pileup_image_height
    if flags.pileup_image_width:
      options.pic_options.width = flags.pileup_image_width

    num_shards, examples, candidates, gvcf = io_utils.resolve_filespecs(
        flags.task, flags.examples or '', flags.candidates or '', flags.gvcf or
        '')
    options.examples_filename = examples
    options.candidates_filename = candidates
    options.gvcf_filename = gvcf

    # redacted
    regions_flag = flags.regions
    if isinstance(regions_flag, str):
      regions_flag = regions_flag.split()
    options.calling_regions.extend(regions_flag)

    options.task_id = flags.task
    options.num_shards = 0 if num_shards is None else num_shards

    if flags.realign_reads:
      options.realigner_enabled = True
      options.realigner_options.CopyFrom(realigner.realigner_config(flags))

    options.max_reads_per_partition = flags.max_reads_per_partition

    if (options.mode == deepvariant_pb2.DeepVariantOptions.TRAINING and
        flags.training_random_emit_ref_sites != NO_RANDOM_REF):
      options.variant_caller_options.fraction_reference_sites_to_emit = (
          flags.training_random_emit_ref_sites)

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


def extract_sample_name_from_reads(reads_path):
  """Returns the sample name as derived from the BAM file of reads.

  Args:
    reads_path: Path to the SAM/BAM file containing a single sample.

  Returns:
    The sample ID annotated in the read group.

  Raises:
    ValueError: There is not exactly one unique sample name in the SAM/BAM.
  """
  with genomics_io.make_sam_reader(reads_path) as sam_reader:
    samples = sam_reader.samples
  if len(samples) != 1:
    raise ValueError('Expected a single sample, found {}'.format(samples))
  sample = next(iter(samples))
  if not sample:
    raise ValueError('Sample name is empty.')
  return sample


# ---------------------------------------------------------------------------
# Region processing
# ---------------------------------------------------------------------------


def common_contigs(contigs_list, exclude_contig_names=None):
  """Gets a list of contigs found in all contigs in contigs_list.

  A common contig is considered one where the name and length in basepairs are
  the same.

  Args:
    contigs_list: A sequence of lists of ContigInfo protos.
    exclude_contig_names: A set/list/etc of str or None. If not None, any contig
      whose name occurs in this sequence of names will be excluded from the list
      of common contigs.

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

  # Remove any excluded contigs from the ref_contigs, as we want to use the
  # selected contigs for our overlap comparison.
  ref_contigs = contigs_list[0]
  if exclude_contig_names:
    ref_contigs = [c for c in ref_contigs if c.name not in exclude_contig_names]

  # Compute the common contigs by recursively getting common contigs of our
  # master set of contigs (contigs) and each contig in other_contigs.
  common = ref_contigs
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
      ranges if randomize_regions==False.
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
    An iterable of learning.genomics.v1.Range objects.

  Raises:
    ValueError: if random is None but randomize_regions is True.
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
# Variant labeler
# ---------------------------------------------------------------------------


class _Counter(object):

  def __init__(self, name, selectp):
    self.name = name
    self.selectp = selectp
    self.n_selected = 0


class VariantCounters(object):
  """Provides stats about the number of variants satisfying pfuncs."""

  def __init__(self, names_and_selectors):
    self.counters = []
    self.n_total = 0
    for name, selector in names_and_selectors:
      self.counters.append(_Counter(name, selector))

  def update(self, variant):
    self.n_total += 1
    for counter in self.counters:
      if counter.selectp(variant):
        counter.n_selected += 1

  def log(self):
    logging.info('----- VariantCounts -----')
    for counter in self.counters:
      percent = (100.0 * counter.n_selected) / (max(self.n_total, 1.0))
      logging.info('%s: %s/%s (%.2f%%)', counter.name, counter.n_selected,
                   self.n_total, percent)


def make_counters():
  """Creates all of the VariantCounters we want to track."""

  def _gt_selector(*gt_types):
    return lambda v: variantutils.genotype_type(v) in gt_types

  return VariantCounters([
      ('All', lambda v: True),
      ('SNPs', variantutils.is_snp),
      ('Indels', variantutils.is_indel),
      ('BiAllelic', variantutils.is_biallelic),
      ('MultiAllelic', variantutils.is_multiallelic),
      ('HomRef', _gt_selector(variantutils.GenotypeType.hom_ref)),
      ('Het', _gt_selector(variantutils.GenotypeType.het)),
      ('HomAlt', _gt_selector(variantutils.GenotypeType.hom_var)),
      ('NonRef',
       _gt_selector(variantutils.GenotypeType.het,
                    variantutils.GenotypeType.hom_var)),
  ])


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
    return allelecounter.AlleleCounter(self.ref_reader, region,
                                       self.options.allele_counter_options)

  def _encode_tensor(self, image_tensor):
    return image_tensor.tostring(), image_tensor.shape, 'raw'

  def _make_sam_reader(self):
    return genomics_io.make_sam_reader(
        self.options.reads_filename,
        self.options.read_requirements,
        hts_block_size=FLAGS.hts_block_size,
        downsample_fraction=self.options.downsample_fraction,
        random_seed=self.options.random_seed)

  def _initialize(self):
    """Initialize the resources needed for this work in the current env."""
    if self.initialized:
      raise ValueError('Cannot initialize this object twice')

    self.ref_reader = genomics_io.make_ref_reader(
        self.options.reference_filename)
    self.sam_reader = self._make_sam_reader()
    self.in_memory_sam_reader = utils.InMemorySamReader([])

    if self.options.realigner_enabled:
      self.realigner = realigner.Realigner(self.options.realigner_options,
                                           self.ref_reader)
    self.pic = pileup_image.PileupImageCreator(
        ref_reader=self.ref_reader,
        sam_reader=self.in_memory_sam_reader,
        options=self.options.pic_options)

    if in_training_mode(self.options):
      self.labeler = variant_labeler.VariantLabeler(
          genomics_io.make_vcf_reader(self.options.truth_variants_filename),
          read_confident_regions(self.options))

    self.variant_caller = variant_caller.VariantCaller(
        self.options.variant_caller_options)
    self.random = np.random.RandomState(self.options.random_seed)
    self.initialized = True

  def process(self, region):
    """Finds candidates and creates corresponding examples in a region.

    Args:
      region: A learning.genomics.v1.Range proto. Specifies the region on the
        genome we should process.

    Returns:
      Three values. First is a list of the found candidates, which are
      deepvariant.DeepVariantCall objects. The second value is a list of filled
      in tf.Example protos. For example, these will include the candidate
      variant, the pileup image, and, if in training mode, the truth variants
      and labels needed for training. The third value is a list of
      learning.genomics.v1.Variant protos containing gVCF information for all
      reference sites, if gvcf generation is enabled, otherwise returns [].
    """
    region_timer = timer.TimerStart()

    # Print some basic information about what we are doing.
    if not self.initialized:
      self._initialize()

    self.in_memory_sam_reader.replace_reads(self.region_reads(region))
    candidates, gvcfs = self.candidates_in_region(region)
    examples = []
    for candidate in candidates:
      for example in self.create_pileup_examples(candidate):
        if in_training_mode(self.options):
          if self.label_variant(example, candidate.variant):
            examples.append(example)
        else:
          examples.append(example)
    logging.info('Found %s candidates in %s [%0.2fs elapsed]', len(examples),
                 ranges.to_literal(region), region_timer.Stop())
    return candidates, examples, gvcfs

  def region_reads(self, region):
    """Update in_memory_sam_reader with read alignments overlapping the region.

    If self.realigner is set, uses realigned reads, otherwise original reads
    are returned.

    Args:
      region: A learning.genomics.v1.Range object specifying the region we
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
      region: A learning.genomics.v1.Range object specifying the region we
      want to get candidates for.

    Returns:
      A 2-tuple. The first value is a list of deepvariant_pb2.DeepVariantCalls
      objects, in coordidate order. The second value is a list of
      learning.genomics.v1.Variant protos containing gVCF information for all
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

  def label_variant(self, example, variant):
    """Adds the truth variant and label for variant to example.

    This function uses VariantLabeler to find a match for variant and writes
    in the correspond truth variant and derived label to our example proto.

    Args:
      example: A tf.Example proto. We will write truth_variant and label into
        this proto.
      variant: A learning.genomics.v1.Variant proto.
        This is the variant we'll use
        to call our VariantLabeler.match to get our truth variant.

    Returns:
      True if the variant was in the confident region (meaning that it could be
        given a label) and False otherwise.
    """
    is_confident, truth_variant = self.labeler.match(variant)
    if not is_confident:
      return False
    alt_alleles = tf_utils.example_alt_alleles(example, variant=variant)
    if variantutils.is_ref(variant):
      label = 0
    else:
      label = self.labeler.match_to_alt_count(variant, truth_variant,
                                              alt_alleles)
    tf_utils.example_set_label(example, label)
    tf_utils.example_set_truth_variant(example, truth_variant)
    return True


def processing_regions_from_options(options):
  """Computes the calling regions from our options.

  This function does all of the work needed to read our input files and region
  specifications to determine the list of regions we should generate examples
  over. It also computes the confident regions need to label variants.

  Args:
    options: deepvariant.DeepVariantOptions proto containing information about
      our input data sources.

  Returns:
    Two values. The first is a list of learning.genomics.v1.Range protos of the
    regions we should process. The second is a RangeSet containing the confident
    regions for labeling, or None if we are running in training mode.
  """
  ref_contigs = genomics_io.make_ref_reader(options.reference_filename).contigs
  sam_contigs = genomics_io.make_sam_reader(options.reads_filename).contigs

  # Add in confident regions and vcf_contigs if in training mode.
  vcf_contigs = None
  if in_training_mode(options):
    vcf_contigs = genomics_io.make_vcf_reader(
        options.truth_variants_filename).contigs

  # Compute the common contigs among our inputs, and check that the contigs are
  # sufficiently consistent among each other.
  contigs = common_contigs(
      only_true(ref_contigs, sam_contigs, vcf_contigs),
      exclude_contig_names=options.exclude_contigs)
  validate_reference_contig_coverage(ref_contigs, contigs,
                                     options.min_shared_contigs_basepairs)
  logging.info('Common contigs are %s', [c.name for c in contigs])

  regions = regions_to_process(
      contigs,
      partition_size=options.allele_counter_options.partition_size,
      calling_regions=ranges.RangeSet.from_regions(
          options.calling_regions, ranges.contigs_dict(ref_contigs)),
      task_id=options.task_id,
      num_shards=options.num_shards)

  return regions


def make_examples_runner(options):
  """Runs examples creation stage of deepvariant."""
  # Counting variants.
  counters = make_counters()

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
  with io_utils.OutputsWriter(options) as writer:
    for region in regions:
      candidates, examples, gvcfs = region_processor.process(region)
      n_candidates += len(candidates)
      n_regions += 1

      writer.write('candidates', *candidates)

      # If we have any gvcf records, write them out. This if also serves to
      # protect us from trying to write to the gvcfs output of writer when gvcf
      # generation is turned off. In that case, gvcfs will always be empty and
      # we'll never execute the write.
      if gvcfs:
        writer.write('gvcfs', *gvcfs)

      for example in examples:
        if in_training_mode(options):
          truth_variant = tf_utils.example_truth_variant(example)
          counters.update(truth_variant)
        writer.write('examples', example)

  logging.info('Found %s candidate variants', n_candidates)
  if in_training_mode(options):
    # This printout is misleading if we are in calling mode.
    counters.log()


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

    # Give htslib authentication access to GCS.
    htslib_gcp_oauth.init()

    # Set up options; may do I/O.
    options = default_options(add_flags=True, flags=FLAGS)

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

    # Run!
    make_examples_runner(options)


if __name__ == '__main__':
  tf.flags.mark_flags_as_required([
      'examples',
      'mode',
      'reads',
      'ref',
  ])
  tf.app.run()
