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
"""Correct read alignment by realigning the read to its most likely haplotype.

This is achieved by constructing de-Bruijn graphs in candidate regions with
potential variations, and determining the mostly likely X haplotypes (where X is
the ploidy).
"""

import csv
import os
import os.path



from absl import flags
from etils import epath

from deepvariant.protos import realigner_pb2
from deepvariant.realigner import window_selector
from deepvariant.realigner.python import debruijn_graph
from deepvariant.realigner.python import fast_pass_aligner
from deepvariant.vendor import timer
from google.protobuf import text_format
from third_party.nucleus.io import sam
from third_party.nucleus.protos import cigar_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import cigar as cigar_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils
from third_party.nucleus.util.cigar import READ_ADVANCING_OPS
from third_party.nucleus.util.cigar import REF_ADVANCING_OPS

_UNSET_WS_INT_FLAG = -1

_WS_USE_WINDOW_SELECTOR_MODEL = flags.DEFINE_bool(
    'ws_use_window_selector_model',
    False,
    'Activate the use of window selector models.',
)
_WS_WINDOW_SELECTOR_MODEL = flags.DEFINE_string(
    'ws_window_selector_model',
    None,
    'Path to a text format proto of the window selector model to use.',
)
_WS_MIN_NUM_SUPPORTING_READS = flags.DEFINE_integer(
    'ws_min_num_supporting_reads',
    _UNSET_WS_INT_FLAG,
    (
        'Minimum number of supporting reads to call a reference position for'
        ' local assembly.'
    ),
)
_WS_MAX_NUM_SUPPORTING_READS = flags.DEFINE_integer(
    'ws_max_num_supporting_reads',
    _UNSET_WS_INT_FLAG,
    (
        'Maximum number of supporting reads to call a reference position for'
        ' local assembly.'
    ),
)
_WS_MIN_MAPQ = flags.DEFINE_integer(
    'ws_min_mapq',
    20,
    (
        'Minimum read alignment quality to consider in calling a reference '
        'position for local assembly.'
    ),
)
_WS_MIN_BASE_QUALITY = flags.DEFINE_integer(
    'ws_min_base_quality',
    20,
    (
        'Minimum base quality to consider in calling a reference position for '
        'local assembly.'
    ),
)
_WS_MIN_WINDOWS_DISTANCE = flags.DEFINE_integer(
    'ws_min_windows_distance',
    80,
    'Minimum distance between candidate windows for local assembly.',
)
_WS_MAX_WINDOW_SIZE = flags.DEFINE_integer(
    'ws_max_window_size',
    1000,
    (
        'Maximum window size to consider for local assembly. Large noisy'
        ' regions are skipped for realignment.'
    ),
)
_WS_REGION_EXPANSION_IN_BP = flags.DEFINE_integer(
    'ws_region_expansion_in_bp',
    20,
    (
        'Number of bases to expand the region when calculating windows; larger'
        ' values add overhead but allow larger nearby events to contribute'
        ' evidence for assembling an region even if they are not contained by'
        ' the region.'
    ),
)
_DBG_MIN_K = flags.DEFINE_integer(
    'dbg_min_k', 10, 'Initial k-mer size to build the graph.'
)
_DBG_MAX_K = flags.DEFINE_integer(
    'dbg_max_k',
    101,
    'Maximum k-mer size. Larger k-mer size is used to resolve graph cycles.',
)
_DBG_STEP_K = flags.DEFINE_integer(
    'dbg_step_k', 1, 'Increment size for k to try in resolving graph cycles.'
)
_DBG_MIN_MAPQ = flags.DEFINE_integer(
    'dbg_min_mapq',
    14,
    'Minimum read alignment quality to consider in building the graph.',
)
_DBG_MIN_BASE_QUALITY = flags.DEFINE_integer(
    'dbg_min_base_quality',
    15,
    (
        'Minimum base quality in a k-mer sequence to consider in building the '
        'graph.'
    ),
)
_DBG_MIN_EDGE_WEIGHT = flags.DEFINE_integer(
    'dbg_min_edge_weight',
    2,
    'Minimum number of supporting reads to keep an edge.',
)
_DBG_MAX_NUM_PATHS = flags.DEFINE_integer(
    'dbg_max_num_paths',
    256,
    (
        'Maximum number of paths within a graph to consider for realignment. '
        'Set max_num_paths to 0 to have unlimited number of paths.'
    ),
)
_ALN_MATCH = flags.DEFINE_integer(
    'aln_match', 4, 'Match score (expected to be a non-negative score).'
)
_ALN_MISMATCH = flags.DEFINE_integer(
    'aln_mismatch', 6, 'Mismatch score (expected to be a non-negative score).'
)
_ALN_GAP_OPEN = flags.DEFINE_integer(
    'aln_gap_open',
    8,
    (
        'Gap open score (expected to be a non-negative score). '
        'Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).'
    ),
)
_ALN_GAP_EXTEND = flags.DEFINE_integer(
    'aln_gap_extend',
    2,
    (
        'Gap extend score (expected to be a non-negative score). '
        'Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).'
    ),
)
_ALN_K = flags.DEFINE_integer(
    'aln_k', 23, 'k-mer size used to index target sequence.'
)
_ALN_ERROR_RATE = flags.DEFINE_float(
    'aln_error_rate', 0.01, 'Estimated sequencing error rate.'
)
_REALIGNER_DIAGNOSTICS = flags.DEFINE_string(
    'realigner_diagnostics',
    '',
    (
        'Root directory where the realigner should place diagnostic output'
        ' (such as a dump of the DeBruijn graph, and a log of metrics'
        ' reflecting the graph and  realignment to the haplotypes).  If empty,'
        ' no diagnostics are output.'
    ),
)
_EMIT_REALIGNED_READS = flags.DEFINE_bool(
    'emit_realigned_reads',
    False,
    (
        'If True, we will emit realigned reads if our realigner_diagnostics are'
        ' also enabled.'
    ),
)
_USE_FAST_PASS_ALIGNER = flags.DEFINE_bool(
    'use_fast_pass_aligner',
    True,
    'This flag should always be True. The old implementation is deprecated.',
)
_MAX_NUM_MISMATCHES = flags.DEFINE_integer(
    'max_num_mismatches',
    2,
    'Num of maximum allowed mismatches for quick read to haplotype alignment.',
)
_REALIGNMENT_SIMILARITY_THRESHOLD = flags.DEFINE_float(
    'realignment_similarity_threshold',
    0.16934,
    'Similarity threshold used in realigner in Smith-Watermanalignment.',
)
_SPLIT_SKIP_READS = flags.DEFINE_bool(
    'split_skip_reads',
    False,
    (
        'If True, splits reads with large SKIP cigar operations '
        'into individual reads. Resulting read parts that are less than '
        '15 bp are filtered out.'
    ),
)
_KMER_SIZE = flags.DEFINE_integer(
    'kmer_size', 32, 'K-mer size for fast pass alinger reads index.'
)

# Margin added to the reference sequence for the aligner module.
_REF_ALIGN_MARGIN = 20

_DEFAULT_MIN_SUPPORTING_READS = 2
_DEFAULT_MAX_SUPPORTING_READS = 300
_ALLELE_COUNT_LINEAR_MODEL_DEFAULT = realigner_pb2.WindowSelectorModel(
    model_type=realigner_pb2.WindowSelectorModel.ALLELE_COUNT_LINEAR,
    allele_count_linear_model=realigner_pb2.WindowSelectorModel.AlleleCountLinearModel(
        bias=-0.683379,
        coeff_soft_clip=2.997000,
        coeff_substitution=-0.086644,
        coeff_insertion=2.493585,
        coeff_deletion=1.795914,
        coeff_reference=-0.059787,
        decision_boundary=3,
    ),
)

# Minimum length of read to retain following splitting with --split_skip_reads.
_MIN_SPLIT_LEN = 15

# ---------------------------------------------------------------------------
# Set configuration settings.
# ---------------------------------------------------------------------------


def window_selector_config(flags_obj):
  """Creates a WindowSelectorOptions proto based on input and default settings.

  Args:
    flags_obj: configuration FLAGS.

  Returns:
    realigner_pb2.WindowSelector protobuf.

  Raises:
    ValueError: If either ws_{min,max}_supporting_reads are set and
      ws_use_window_selector_model is True.
      Or if ws_window_selector_model > ws_max_num_supporting_reads.
      Or if ws_use_window_selector_model is False and
      ws_window_selector_model is not None.
  """
  if not flags_obj.ws_use_window_selector_model:
    if flags_obj.ws_window_selector_model is not None:
      raise ValueError(
          'Cannot specify a ws_window_selector_model '
          'if ws_use_window_selector_model is False.'
      )

    min_num_supporting_reads = (
        _DEFAULT_MIN_SUPPORTING_READS
        if flags_obj.ws_min_num_supporting_reads == _UNSET_WS_INT_FLAG
        else flags_obj.ws_min_num_supporting_reads
    )
    max_num_supporting_reads = (
        _DEFAULT_MAX_SUPPORTING_READS
        if flags_obj.ws_max_num_supporting_reads == _UNSET_WS_INT_FLAG
        else flags_obj.ws_max_num_supporting_reads
    )
    window_selector_model = realigner_pb2.WindowSelectorModel(
        model_type=realigner_pb2.WindowSelectorModel.VARIANT_READS,
        variant_reads_model=realigner_pb2.WindowSelectorModel.VariantReadsThresholdModel(
            min_num_supporting_reads=min_num_supporting_reads,
            max_num_supporting_reads=max_num_supporting_reads,
        ),
    )
  else:
    if flags_obj.ws_min_num_supporting_reads != _UNSET_WS_INT_FLAG:
      raise ValueError(
          'Cannot use both ws_min_num_supporting_reads and '
          'ws_use_window_selector_model flags.'
      )
    if flags_obj.ws_max_num_supporting_reads != _UNSET_WS_INT_FLAG:
      raise ValueError(
          'Cannot use both ws_max_num_supporting_reads and '
          'ws_use_window_selector_model flags.'
      )

    if flags_obj.ws_window_selector_model is None:
      window_selector_model = _ALLELE_COUNT_LINEAR_MODEL_DEFAULT
    else:
      with epath.Path(flags_obj.ws_window_selector_model).open() as f:
        window_selector_model = text_format.Parse(
            f.read(), realigner_pb2.WindowSelectorModel()
        )

  if (
      window_selector_model.model_type
      == realigner_pb2.WindowSelectorModel.VARIANT_READS
  ):
    model = window_selector_model.variant_reads_model
    if model.max_num_supporting_reads < model.min_num_supporting_reads:
      raise ValueError(
          'ws_min_supporting_reads should be smaller than '
          'ws_max_supporting_reads.'
      )

  keep_legacy_behavior = False
  if 'keep_legacy_allele_counter_behavior' in flags_obj:
    keep_legacy_behavior = flags_obj.keep_legacy_allele_counter_behavior
  ws_config = realigner_pb2.WindowSelectorOptions(
      min_mapq=flags_obj.ws_min_mapq,
      min_base_quality=flags_obj.ws_min_base_quality,
      min_windows_distance=flags_obj.ws_min_windows_distance,
      max_window_size=flags_obj.ws_max_window_size,
      region_expansion_in_bp=flags_obj.ws_region_expansion_in_bp,
      window_selector_model=window_selector_model,
      keep_legacy_behavior=keep_legacy_behavior,
  )

  return ws_config


def realigner_config(flags_obj):
  """Creates a RealignerOptions proto based on input and default settings.

  Args:
    flags_obj: configuration FLAGS.

  Returns:
    realigner_pb2.RealignerOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  ws_config = window_selector_config(flags_obj)

  dbg_config = realigner_pb2.DeBruijnGraphOptions(
      min_k=flags_obj.dbg_min_k,
      max_k=flags_obj.dbg_max_k,
      step_k=flags_obj.dbg_step_k,
      min_mapq=flags_obj.dbg_min_mapq,
      min_base_quality=flags_obj.dbg_min_base_quality,
      min_edge_weight=flags_obj.dbg_min_edge_weight,
      max_num_paths=flags_obj.dbg_max_num_paths,
  )

  aln_config = realigner_pb2.AlignerOptions(
      match=flags_obj.aln_match,
      mismatch=flags_obj.aln_mismatch,
      gap_open=flags_obj.aln_gap_open,
      gap_extend=flags_obj.aln_gap_extend,
      k=flags_obj.aln_k,
      error_rate=flags_obj.aln_error_rate,
      max_num_of_mismatches=flags_obj.max_num_mismatches,
      realignment_similarity_threshold=flags_obj.realignment_similarity_threshold,
      kmer_size=flags_obj.kmer_size,
      force_alignment=False,
  )

  diagnostics = realigner_pb2.Diagnostics(
      enabled=bool(flags_obj.realigner_diagnostics),
      output_root=flags_obj.realigner_diagnostics,
      emit_realigned_reads=flags_obj.emit_realigned_reads,
  )

  # The normalize_reads flag could came from the `flags_obj` arg, passed in
  # from make_examples_options.py. It is already part of AlleleCounterOptions in
  # MakeExamplesOptions. Here, we need to set it in RealignerOptions as well
  # because an if statement in fast_pass_aligner.cc needs it to decide whether
  # to run a specific logic.
  # This is not ideal. If there's a way to improve this, please do.
  normalize_reads = False
  if 'normalize_reads' in flags_obj:
    normalize_reads = flags_obj.normalize_reads
  return realigner_pb2.RealignerOptions(
      ws_config=ws_config,
      dbg_config=dbg_config,
      aln_config=aln_config,
      diagnostics=diagnostics,
      split_skip_reads=flags_obj.split_skip_reads,
      normalize_reads=normalize_reads,
  )


class DiagnosticLogger(object):
  """Writes diagnostic information about the assembler."""

  def __init__(
      self,
      config,
      graph_filename='graph.dot',
      metrics_filename='realigner_metrics.csv',
      realigned_reads_filename='realigned_reads.bam',
  ):
    self.config = config
    self.graph_filename = graph_filename
    self.metrics_filename = metrics_filename
    self.realigned_reads_filename = realigned_reads_filename

    # Setup diagnostics outputs if requested.
    if self.enabled:
      self._csv_file = open(self._root_join(self.metrics_filename), 'w')
      self._csv_writer = csv.writer(self._csv_file)
      self._write_csv_line('window', 'k', 'n_haplotypes', 'time')
    else:
      self._csv_file = None
      self._csv_writer = None

  def close(self):
    if self.enabled:
      assert self._csv_file is not None
      self._csv_file.close()

  @property
  def enabled(self):
    return self.config and self.config.enabled

  def _root_join(self, path: str, makedirs: bool = True) -> str:
    fullpath = os.path.join(self.config.output_root, path)
    subdir = os.path.dirname(fullpath)
    if makedirs and subdir:
      epath.Path(subdir).mkdir(parents=True, exist_ok=True)
    return fullpath

  def _write_csv_line(self, *args):
    assert self.enabled, 'only callable when diagnostics are on'
    assert self._csv_writer is not None
    self._csv_writer.writerow(args)

  def _file_for_region(self, region, basename):
    """Returns the path to a file in a region-specific subdirectory."""
    assert self.enabled, 'only callable when diagnostics are on'
    return self._root_join(os.path.join(ranges.to_literal(region), basename))

  def log_realigned_reads(self, region, reads, shared_header=None):
    """Logs, if enabled, the realigned reads for region."""
    if (
        self.enabled
        and self.config.emit_realigned_reads
        and shared_header is not None
    ):
      path = self._file_for_region(region, self.realigned_reads_filename)
      with sam.SamWriter(path, header=shared_header) as writer:
        # For realigned reads, sorting by just looking at starting position is
        # enough.
        for read in sorted(
            reads, key=lambda read: read.alignment.position.position
        ):
          writer.write(read)

  def log_graph_metrics(
      self, region, graph, candidate_haplotypes, graph_building_time
  ):
    """Logs, if enabled, graph construction information for region."""
    if self.enabled:
      if graph:
        dest_file = self._file_for_region(region, self.graph_filename)
        with epath.Path(dest_file).open('w') as f:
          f.write(graph.graphviz())
      self._write_csv_line(
          ranges.to_literal(region),
          graph.kmer_size if graph else 'NA',
          len(candidate_haplotypes),
          graph_building_time,
      )


class AssemblyRegion(object):
  """A region to assemble, holding the region Range and the reads.

  It is not safe to directly modify any of the attributes here. Use the accessor
  functions to add a read to the reads.

  Attributes:
    candidate_haplotypes: realigner.CandidateHaplotypes for this region.
    reads: list[reads_pb2.Read]. Reads for this region.
    region: range_pb2.Range. This is the span of the assembled region on the
      genome.
    read_span: range_pb2.Range. This is the span of reads added to this region.
      The read_span in general is expected to be wider than the region itself,
      since we often include all reads that overlap the region at all. It is
      possible that read_span will be smaller than region, which can happen, for
      example, when we only have reads starts in the middle of the region.
      Here's a picture of when this can happen: ref      : acgtACGTACgtgt region
      :     ------ read1    :       GGa
  read_span:       ---
  """

  def __init__(self, candidate_haplotypes):
    self.candidate_haplotypes = candidate_haplotypes
    self.reads = []
    self._read_span = None

  def __str__(self):
    return (
        'AssemblyRegion(region={}, span={}) with {} haplotypes and {} reads'
    ).format(
        ranges.to_literal(self.region),
        ranges.to_literal(self.read_span),
        len(self.haplotypes),
        len(self.reads),
    )

  @property
  def haplotypes(self):
    """Returns the haplotypes list[str] of our candidate_haplotypes."""
    return self.candidate_haplotypes.haplotypes

  @property
  def region(self):
    return self.candidate_haplotypes.span

  @property
  def read_span(self):
    if self._read_span is None and self.reads:
      spans = [utils.read_range(r) for r in self.reads]
      self._read_span = ranges.make_range(
          spans[0].reference_name,
          min(s.start for s in spans),
          max(s.end for s in spans),
      )
    return self._read_span

  def add_read(self, read):
    self.reads.append(read)
    self._read_span = None  # Adding a read invalidates our _read_span cache.


def assign_reads_to_assembled_regions(assembled_regions, reads):
  """Assign each read to the maximally overlapped window.

  Args:
    assembled_regions: list[AssemblyRegion], list of AssemblyRegion to assign
      reads to. Does not assume AssemblyRegion are sorted.
    reads: iterable[learning.genomics.genomics.Read], to be processed. Does not
      assume the reads are sorted.

  Returns:
    [AssemblyRegion], information on assigned reads for each assembled region.
    list[learning.genomics.genomics.Read], the list of unassigned reads.
  """
  regions = [ar.region for ar in assembled_regions]
  unassigned_reads = []
  for read in reads:
    read_range = utils.read_range(read)
    window_i = ranges.find_max_overlapping(read_range, regions)
    if window_i is not None:
      assembled_regions[window_i].add_read(read)
    else:
      unassigned_reads.append(read)
  return unassigned_reads


def copy_read(read, part):
  """Copies a read proto to create a new read part."""
  new_read = reads_pb2.Read()
  new_read.CopyFrom(read)

  # Reset alignment information.
  # Note: If long reads will be used here, convert
  # to approach used for read trimming.
  new_read.alignment.Clear()
  new_read.aligned_quality[:] = []
  new_read.aligned_sequence = ''
  new_read.alignment.position.reference_name = (
      read.alignment.position.reference_name
  )
  new_read.alignment.position.reverse_strand = (
      read.alignment.position.reverse_strand
  )
  new_read.alignment.mapping_quality = read.alignment.mapping_quality
  new_read.fragment_name = f'{new_read.fragment_name}_p{part}'
  return new_read


def split_reads(reads):
  """Splits reads containing SKIP cigar operations into multiple parts."""

  read_split = []
  for read in reads:
    # Check for SKIP operations within the Cigar String
    if not any(
        [c.operation == cigar_pb2.CigarUnit.SKIP for c in read.alignment.cigar]
    ):
      read_split.append(read)
      continue
    part = 0
    new_read = copy_read(read, part)
    # Get alignment position
    read_start = 0
    read_offset = 0
    reference_offset = 0

    for n, cigar in enumerate(read.alignment.cigar):
      on_last_operation = n + 1 == len(read.alignment.cigar)

      if cigar.operation in REF_ADVANCING_OPS:
        # Set position where the first ref op encountered.
        if not new_read.alignment.position.position:
          new_read.alignment.position.position = (
              read.alignment.position.position + reference_offset
          )
        reference_offset += cigar.operation_length

      if cigar.operation in READ_ADVANCING_OPS:
        read_offset += cigar.operation_length

      if cigar.operation != cigar_pb2.CigarUnit.SKIP:
        new_read.alignment.cigar.append(cigar)

      if cigar.operation == cigar_pb2.CigarUnit.SKIP or on_last_operation:
        new_read.aligned_sequence = read.aligned_sequence[
            read_start:read_offset
        ]
        new_read.aligned_quality.extend(
            read.aligned_quality[read_start:read_offset]
        )
        if len(new_read.aligned_sequence) >= _MIN_SPLIT_LEN:
          read_split.append(new_read)
        if not on_last_operation:
          read_start = read_offset
          part += 1
          new_read = copy_read(read, part)

  return read_split


class Realigner(object):
  """Realign reads in regions to assembled haplotypes.

  This class helps us to realign reads in regions by:

  (1) Create smaller windows in which to operate over the region. These windows
  are created by finding evidence of genetic variation surrounded by stretches
  of reference-matching seqence.

  (2) Build a de-Bruijn assembly graph of the window. Edges are pruned if they
  don't meet the required weight threshold. Every remaining haplotype is listed
  by traversing the graph.

  (3) Realign reads using a Smith-Waterman algorithm to the best candidate
  haplotype and then realign that haplotype to the reference sequence to modify
  the read's alignment.
  """

  def __init__(self, config, ref_reader, shared_header=None):
    """Creates a new Realigner.

    Args:
      config: realigner_pb2.RealignerOptions protobuf.
      ref_reader: GenomeReferenceFai, indexed reference genome to query bases.
      shared_header: header info from the input bam file
    """
    self.config = config
    self.ref_reader = ref_reader
    self.diagnostic_logger = DiagnosticLogger(self.config.diagnostics)
    self.shared_header = shared_header

  def call_debruijn_graph(self, windows, reads):
    """Helper function to call debruijn_graph module."""
    windows_haplotypes = []
    # Build and process de-Bruijn graph for each window.
    sam_reader = sam.InMemorySamReader(reads)

    for window in windows:
      if window.end - window.start > self.config.ws_config.max_window_size:
        continue
      if not self.ref_reader.is_valid(window):
        continue
      ref = self.ref_reader.query(window)
      window_reads = list(sam_reader.query(window))

      with timer.Timer() as t:
        graph = debruijn_graph.build(ref, window_reads, self.config.dbg_config)
      graph_building_time = t.GetDuration()

      if not graph:
        candidate_haplotypes = [ref]
      else:
        candidate_haplotypes = graph.candidate_haplotypes()
      if candidate_haplotypes and candidate_haplotypes != [ref]:
        candidate_haplotypes_info = realigner_pb2.CandidateHaplotypes(
            span=window, haplotypes=candidate_haplotypes
        )
        windows_haplotypes.append(candidate_haplotypes_info)

      self.diagnostic_logger.log_graph_metrics(
          window, graph, candidate_haplotypes, graph_building_time
      )

    return windows_haplotypes

  def call_fast_pass_aligner(self, assembled_region):
    """Helper function to call fast pass aligner module."""
    if not assembled_region.reads:
      return []

    contig = assembled_region.region.reference_name
    ref_start = max(
        0,
        min(assembled_region.read_span.start, assembled_region.region.start)
        - _REF_ALIGN_MARGIN,
    )
    ref_end = min(
        self.ref_reader.contig(contig).n_bases,
        max(assembled_region.read_span.end, assembled_region.region.end)
        + _REF_ALIGN_MARGIN,
    )

    ref_prefix = self.ref_reader.query(
        ranges.make_range(contig, ref_start, assembled_region.region.start)
    )
    ref = self.ref_reader.query(assembled_region.region)

    # If we can't create the ref suffix then return the original alignments.
    if ref_end <= assembled_region.region.end:
      return assembled_region.reads
    else:
      ref_suffix = self.ref_reader.query(
          ranges.make_range(contig, assembled_region.region.end, ref_end)
      )

    ref_seq = ref_prefix + ref + ref_suffix

    fast_pass_realigner = fast_pass_aligner.FastPassAligner()
    # Read sizes may vary. We need this for realigner initialization and sanity
    # checks.
    self.config.aln_config.read_size = len(
        assembled_region.reads[0].aligned_sequence
    )
    self.config.aln_config.force_alignment = False
    fast_pass_realigner.set_normalize_reads(self.config.normalize_reads)
    fast_pass_realigner.set_options(self.config.aln_config)
    fast_pass_realigner.set_reference(ref_seq)
    fast_pass_realigner.set_ref_start(contig, ref_start)
    fast_pass_realigner.set_ref_prefix_len(len(ref_prefix))
    fast_pass_realigner.set_ref_suffix_len(len(ref_suffix))
    fast_pass_realigner.set_haplotypes(
        [
            ref_prefix + target + ref_suffix
            for target in assembled_region.haplotypes
        ]
    )
    return fast_pass_realigner.realign_reads(assembled_region.reads)

  def realign_reads(self, reads, region):
    """Run realigner.

    This is the main function that
      - parses the input reads and reference sequence.
      - select candidate windows for local assembly (WindowSelector (ws)
        module).
        - Windows larger than max_window_size are skipped.
      - build pruned De-Bruijn graph for each candidate window (DeBruijnGraph
        (dbg) module).
        - Graphs with more than max_num_paths candidate haplotypes or
          with reference sequence as the only candidate are skipped.
      - Align reads based on candidate haplotypes (Aligner (aln) module).
      - Output all input reads (whether they required realignment or not).

    Args:
      reads: [`third_party.nucleus.protos.Read` protos]. The list of input reads
        to realign.
      region: A `third_party.nucleus.protos.Range` proto. Specifies the region
        on the genome we should process.

    Returns:
      [realigner_pb2.CandidateHaplotypes]. Information on the list of candidate
        haplotypes.
      [`third_party.nucleus.protos.Read` protos]. The realigned
        reads for the region. NOTE THESE READS MAY NO LONGER BE IN THE SAME
        ORDER AS BEFORE.
    """

    # Split reads containing SKIP operations.
    if self.config.split_skip_reads:
      reads = split_reads(reads)

    # Compute the windows where we need to assemble in the region.
    candidate_windows = window_selector.select_windows(
        self.config.ws_config, self.ref_reader, reads, region
    )

    # Assemble each of those regions.
    candidate_haplotypes = self.call_debruijn_graph(candidate_windows, reads)
    # Create our simple container to store candidate / read mappings.
    assembled_regions = [AssemblyRegion(ch) for ch in candidate_haplotypes]

    # Our realigned_reads start off with all of the unassigned reads.
    realigned_reads = assign_reads_to_assembled_regions(
        assembled_regions, reads
    )

    if not _USE_FAST_PASS_ALIGNER.value:
      raise ValueError(
          '--use_fast_pass_aligner is always true. '
          'The older implementation is deprecated and removed.'
      )
    # Walk over each region and align the reads in that region, adding them to
    # our realigned_reads.
    for assembled_region in assembled_regions:
      realigned_reads_copy = self.call_fast_pass_aligner(assembled_region)
      realigned_reads.extend(realigned_reads_copy)

    self.diagnostic_logger.log_realigned_reads(
        region, realigned_reads, self.shared_header
    )

    return candidate_haplotypes, realigned_reads

  def align_to_haplotype(
      self, this_haplotype, haplotypes, prefix, suffix, reads, contig, ref_start
  ):
    """Align reads to a given haplotype, not necessarily the reference.

    Align reads to a graph of haplotypes, reporting the alignments relative
    to a specific haplotype. This allows treating any alternate allele as
    the reference.

    Args:
      this_haplotype: string. Sequence of the haplotype to treat as reference,
        reporting alignments according to its coordinates.
      haplotypes: list of strings. All haplotypes to use in the graph, including
        this_haplotype.
      prefix: string. Sequence to the left of where the haplotypes differ.
      suffix: string. Sequence to the right of where the haplotypes differ.
      reads: reads to align.
      contig: string. Name of the 'reference' to report in read alignments.
      ref_start: integer. Start position of the region to report in read
        alignments. This should mark the beginning of the prefix sequence.

    Returns:
      Reads. Realigned and reported relative to the chosen haplotype.
    """
    if not reads:
      return []
    fast_pass_realigner = fast_pass_aligner.FastPassAligner()
    aln_config = self.config.aln_config
    aln_config.read_size = len(reads[0].aligned_sequence)
    aln_config.force_alignment = True
    fast_pass_realigner.set_options(aln_config)
    fast_pass_realigner.set_reference(prefix + this_haplotype + suffix)
    fast_pass_realigner.set_ref_start(contig, ref_start)

    # Testing found that when the prefix and suffix both go right up to the
    # ref/alt variants, the alignment does not work well, so a margin of 100
    # bases on each side of the variant are used here to pad each
    # haplotype with enough sequence to align against. While some further
    # testing showed this could be reduced, 100 is the only value that has been
    # tested with a full training experiment.
    central_allele_margin = min(len(prefix), len(suffix), 100)
    fast_pass_realigner.set_ref_prefix_len(len(prefix) - central_allele_margin)
    fast_pass_realigner.set_ref_suffix_len(len(suffix) - central_allele_margin)
    extended_haplotypes = [prefix + target + suffix for target in haplotypes]
    fast_pass_realigner.set_haplotypes(extended_haplotypes)
    return fast_pass_realigner.realign_reads(reads)


def trim_cigar(cigar, ref_trim, ref_length):
  """Trim a cigar string to a certain reference length.

  Args:
    cigar: list of `nucleus.protos.CigarUnit`s of the original read alignment.
    ref_trim: integer. Number of reference bases to trim off the beginning of
      the read.
    ref_length: integer. Number of reference bases to cover with the read, the
      middle part that is not trimmed from the start or end of the read.

  Returns:
    new_cigar: list of `nucleus.protos.CigarUnit`s of the final read alignment,
        after the left and/or right have been trimmed off.
    read_trim: The number of bases of the read that are trimmed off.
    new_read_length: The number of bases of the read that remain after trimming.
        This is different from the final number of reference bases; for example,
        an insertion makes the read longer without affecting the reference.
  """
  # First consume the ref until the trim is covered.
  trim_remaining = ref_trim
  # Then consume the ref until the ref_length is covered.
  ref_to_cover_remaining = ref_length
  read_trim = 0
  new_cigar = []
  new_read_length = 0
  for cigar_unit in cigar:
    c_operation_length = cigar_unit.operation_length
    # Each operation moves forward in the ref, the read, or both.
    advances_ref = cigar_unit.operation in cigar_utils.REF_ADVANCING_OPS
    advances_read = cigar_unit.operation in cigar_utils.READ_ADVANCING_OPS
    ref_step = c_operation_length if advances_ref else 0
    # First, use up each operation until the trimmed area is covered.
    if trim_remaining > 0:
      if ref_step <= trim_remaining:
        # Fully apply to the trim.
        trim_remaining -= ref_step
        read_trim += c_operation_length if advances_read else 0
        continue
      else:
        # Partially apply to finish the trim.
        ref_step -= trim_remaining
        read_trim += trim_remaining if advances_read else 0
        # If trim finishes here, the rest of the ref_step can apply to the
        # next stage and count towards covering the given ref window.
        c_operation_length = ref_step
        trim_remaining = 0

    # Once the trim is done, start applying cigar entries to covering the ref
    # window.
    if trim_remaining == 0:
      if ref_step <= ref_to_cover_remaining:
        # Fully apply to the window.
        new_cigar.append(
            cigar_pb2.CigarUnit(
                operation=cigar_unit.operation,
                operation_length=c_operation_length,
            )
        )
        ref_to_cover_remaining -= ref_step
        new_read_length += c_operation_length if advances_read else 0
      else:
        # Partially apply to finish the window.
        c_operation_length = ref_to_cover_remaining
        new_cigar.append(
            cigar_pb2.CigarUnit(
                operation=cigar_unit.operation,
                operation_length=c_operation_length,
            )
        )
        new_read_length += c_operation_length if advances_read else 0
        ref_to_cover_remaining = 0
        break

  return new_cigar, read_trim, new_read_length


def trim_read(read, region):
  """Trim a read down to the part that aligns within a given region.

  The following properties of the read are updated, trimming on both sides as
  necessary to save only the parts of the read that fit fully within the
  region, potentially starting and ending at the region's boundaries:
  - The alignment position (read.alignment.position.position).
  - The read sequence (read.aligned_sequence).
  - Base qualities (read.aligned_quality).
  - The cigar string of the alignment (read.alignment.cigar)

  Args:
    read: A `nucleus.protos.Read` that is aligned to the region.
    region: A `nucleus.protos.Range` region.

  Returns:
    a new `nucleus.protos.Read` trimmed to the region.
  """
  if not read.alignment:
    raise ValueError('Read must already be aligned.')

  read_start = read.alignment.position.position

  trim_left = max(region.start - read_start, 0)

  ref_length = region.end - max(region.start, read_start)
  new_cigar, read_trim, new_read_length = trim_cigar(
      read.alignment.cigar, trim_left, ref_length
  )

  # Copy everything but aligned_sequence and aligned_quality fields of the read
  # to get all recursive properties and prevent mutating the original.
  new_read = reads_pb2.Read()
  new_read.fragment_name = read.fragment_name
  new_read.id = read.id
  new_read.read_group_id = read.read_group_id
  new_read.read_group_set_id = read.read_group_set_id
  new_read.read_number = read.read_number
  new_read.fragment_length = read.fragment_length
  new_read.number_reads = read.number_reads
  for each_info_key in read.info:
    new_read.info[each_info_key].CopyFrom(read.info[each_info_key])
  new_read.alignment.position.position = read.alignment.position.position
  new_read.alignment.position.reference_name = (
      read.alignment.position.reference_name
  )
  new_read.alignment.position.reverse_strand = (
      read.alignment.position.reverse_strand
  )
  new_read.alignment.mapping_quality = read.alignment.mapping_quality
  # Following fields are not needed but we copy them for consistency:
  new_read.next_mate_position.CopyFrom(read.next_mate_position)
  new_read.proper_placement = read.proper_placement
  new_read.duplicate_fragment = read.duplicate_fragment
  new_read.failed_vendor_quality_checks = read.failed_vendor_quality_checks
  new_read.secondary_alignment = read.secondary_alignment
  new_read.supplementary_alignment = read.supplementary_alignment

  if trim_left != 0:
    new_read.alignment.position.position = region.start
  # Set aligned_sequence, a string:
  new_read.aligned_sequence = read.aligned_sequence[
      read_trim : read_trim + new_read_length
  ]
  # Set aligned_quality, a repeated integer:
  new_read.aligned_quality[:] = read.aligned_quality[
      read_trim : read_trim + new_read_length
  ]

  # Direct assignment on a repeated message field is not allowed, so setting
  # the cigar by using 'extend'.
  new_read.alignment.cigar.extend(new_cigar)

  return new_read
