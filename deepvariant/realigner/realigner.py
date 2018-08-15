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
"""Correct read alignment by realigning the read to its most likely haplotype.

This is achieved by constructing de-Bruijn graphs in candidate regions with
potential variations, and determining the mostly likely X haplotypes (where X is
the ploidy).
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import csv
import os
import os.path

from absl import flags
import tensorflow as tf

from third_party.nucleus.io import sam
from third_party.nucleus.util import ranges
from third_party.nucleus.util import utils
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import aligner
from deepvariant.realigner import window_selector
from deepvariant.realigner.python import debruijn_graph
from deepvariant.realigner.python import fast_pass_aligner
from deepvariant.vendor import timer
from google.protobuf import text_format

_UNSET_WS_INT_FLAG = -1

flags.DEFINE_bool('ws_use_window_selector_model', True,
                  'Activate the use of window selector models.')
flags.DEFINE_string(
    'ws_window_selector_model', None,
    'Path to a text format proto of the window selector model to use.')
flags.DEFINE_integer(
    'ws_min_num_supporting_reads', _UNSET_WS_INT_FLAG,
    'Minimum number of supporting reads to call a reference position for local '
    'assembly.')
flags.DEFINE_integer(
    'ws_max_num_supporting_reads', _UNSET_WS_INT_FLAG,
    'Maximum number of supporting reads to call a reference position for local '
    'assembly.')
flags.DEFINE_integer(
    'ws_min_mapq', 20,
    'Minimum read alignment quality to consider in calling a reference '
    'position for local assembly.')
flags.DEFINE_integer(
    'ws_min_base_quality', 20,
    'Minimum base quality to consider in calling a reference position for '
    'local assembly.')
flags.DEFINE_integer(
    'ws_min_windows_distance', 80,
    'Minimum distance between candidate windows for local assembly.')
flags.DEFINE_integer(
    'ws_max_window_size', 1000,
    'Maximum window size to consider for local assembly. Large noisy regions '
    'are skipped for realignment.')
flags.DEFINE_integer(
    'ws_region_expansion_in_bp', 20,
    'Number of bases to expand the region when calculating windows; larger '
    'values add overhead but allow larger nearby events to contribute evidence '
    'for assembling an region even if they are not contained by the region.')
flags.DEFINE_integer('dbg_min_k', 10, 'Initial k-mer size to build the graph.')
flags.DEFINE_integer(
    'dbg_max_k', 101,
    'Maximum k-mer size. Larger k-mer size is used to resolve graph cycles.')
flags.DEFINE_integer('dbg_step_k', 1,
                     'Increment size for k to try in resolving graph cycles.')
flags.DEFINE_integer(
    'dbg_min_mapq', 14,
    'Minimum read alignment quality to consider in building the graph.')
flags.DEFINE_integer(
    'dbg_min_base_quality', 15,
    'Minimum base quality in a k-mer sequence to consider in building the '
    'graph.')
flags.DEFINE_integer('dbg_min_edge_weight', 2,
                     'Minimum number of supporting reads to keep an edge.')
flags.DEFINE_integer(
    'dbg_max_num_paths', 256,
    'Maximum number of paths within a graph to consider for realignment. '
    'Set max_num_paths to 0 to have unlimited number of paths.')
flags.DEFINE_integer('aln_match', 4,
                     'Match score (expected to be a non-negative score).')
flags.DEFINE_integer('aln_mismatch', 6,
                     'Mismatch score (expected to be a non-negative score).')
flags.DEFINE_integer(
    'aln_gap_open', 8, 'Gap open score (expected to be a non-negative score). '
    'Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).')
flags.DEFINE_integer(
    'aln_gap_extend', 2,
    'Gap extend score (expected to be a non-negative score). '
    'Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).')
flags.DEFINE_integer('aln_k', 23, 'k-mer size used to index target sequence.')
flags.DEFINE_float('aln_error_rate', .01, 'Estimated sequencing error rate.')
flags.DEFINE_string(
    'realigner_diagnostics', '',
    'Root directory where the realigner should place diagnostic output (such as'
    ' a dump of the DeBruijn graph, and a log of metrics reflecting the graph '
    'and  realignment to the haplotypes).  If empty, no diagnostics are output.'
)
flags.DEFINE_bool(
    'emit_realigned_reads', False,
    'If True, we will emit realigned reads if our realigner_diagnostics are '
    'also enabled.')
flags.DEFINE_bool(
    'use_fast_pass_aligner', True,
    'If True, fast_pass_aligner (improved performance) implementation is used ')
flags.DEFINE_integer(
    'max_num_mismatches', 2,
    'Num of maximum allowed mismatches for quick read to '
    'haplotype alignment.')
flags.DEFINE_float(
    'realignment_similarity_threshold', 0.16934,
    'Similarity threshold used in realigner in Smith-Waterman'
    'alignment.')
flags.DEFINE_integer('kmer_size', 32,
                     'K-mer size for fast pass alinger reads index.')


# Margin added to the reference sequence for the aligner module.
_REF_ALIGN_MARGIN = 20

_DEFAULT_MIN_SUPPORTING_READS = 2
_DEFAULT_MAX_SUPPORTING_READS = 300
_ALLELE_COUNT_LINEAR_MODEL_DEFAULT = realigner_pb2.WindowSelectorModel(
    model_type=realigner_pb2.WindowSelectorModel.ALLELE_COUNT_LINEAR,
    allele_count_linear_model=realigner_pb2.WindowSelectorModel.
    AlleleCountLinearModel(
        bias=-0.683379,
        coeff_soft_clip=2.997000,
        coeff_substitution=-0.086644,
        coeff_insertion=2.493585,
        coeff_deletion=1.795914,
        coeff_reference=-0.059787,
        decision_boundary=3))

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
      raise ValueError('Cannot specify a ws_window_selector_model '
                       'if ws_use_window_selector_model is False.')

    min_num_supporting_reads = (
        _DEFAULT_MIN_SUPPORTING_READS
        if flags_obj.ws_min_num_supporting_reads == _UNSET_WS_INT_FLAG else
        flags_obj.ws_min_num_supporting_reads)
    max_num_supporting_reads = (
        _DEFAULT_MAX_SUPPORTING_READS
        if flags_obj.ws_max_num_supporting_reads == _UNSET_WS_INT_FLAG else
        flags_obj.ws_max_num_supporting_reads)
    window_selector_model = realigner_pb2.WindowSelectorModel(
        model_type=realigner_pb2.WindowSelectorModel.VARIANT_READS,
        variant_reads_model=realigner_pb2.WindowSelectorModel.
        VariantReadsThresholdModel(
            min_num_supporting_reads=min_num_supporting_reads,
            max_num_supporting_reads=max_num_supporting_reads))
  else:
    if flags_obj.ws_min_num_supporting_reads != _UNSET_WS_INT_FLAG:
      raise ValueError('Cannot use both ws_min_num_supporting_reads and '
                       'ws_use_window_selector_model flags.')
    if flags_obj.ws_max_num_supporting_reads != _UNSET_WS_INT_FLAG:
      raise ValueError('Cannot use both ws_max_num_supporting_reads and '
                       'ws_use_window_selector_model flags.')

    if flags_obj.ws_window_selector_model is None:
      window_selector_model = _ALLELE_COUNT_LINEAR_MODEL_DEFAULT
    else:
      with tf.gfile.GFile(flags_obj.ws_window_selector_model) as f:
        window_selector_model = text_format.Parse(
            f.read(), realigner_pb2.WindowSelectorModel())

  if (window_selector_model.model_type ==
      realigner_pb2.WindowSelectorModel.VARIANT_READS):
    model = window_selector_model.variant_reads_model
    if model.max_num_supporting_reads < model.min_num_supporting_reads:
      raise ValueError('ws_min_supporting_reads should be smaller than '
                       'ws_max_supporting_reads.')

  ws_config = realigner_pb2.WindowSelectorOptions(
      min_mapq=flags_obj.ws_min_mapq,
      min_base_quality=flags_obj.ws_min_base_quality,
      min_windows_distance=flags_obj.ws_min_windows_distance,
      max_window_size=flags_obj.ws_max_window_size,
      region_expansion_in_bp=flags_obj.ws_region_expansion_in_bp,
      window_selector_model=window_selector_model)

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
      max_num_paths=flags_obj.dbg_max_num_paths)

  aln_config = realigner_pb2.AlignerOptions(
      match=flags_obj.aln_match,
      mismatch=flags_obj.aln_mismatch,
      gap_open=flags_obj.aln_gap_open,
      gap_extend=flags_obj.aln_gap_extend,
      k=flags_obj.aln_k,
      error_rate=flags_obj.aln_error_rate,
      max_num_of_mismatches=flags_obj.max_num_mismatches,
      realignment_similarity_threshold=flags_obj.
      realignment_similarity_threshold,
      kmer_size=flags_obj.kmer_size)

  diagnostics = realigner_pb2.Diagnostics(
      enabled=bool(flags_obj.realigner_diagnostics),
      output_root=flags_obj.realigner_diagnostics,
      emit_realigned_reads=flags_obj.emit_realigned_reads)

  return realigner_pb2.RealignerOptions(
      ws_config=ws_config,
      dbg_config=dbg_config,
      aln_config=aln_config,
      diagnostics=diagnostics)


class DiagnosticLogger(object):
  """Writes diagnostic information about the assembler."""

  def __init__(self,
               config,
               graph_filename='graph.dot',
               metrics_filename='realigner_metrics.csv',
               realigned_reads_filename='realigned_reads.tfrecord'):
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
      self._csv_file.close()

  @property
  def enabled(self):
    return self.config and self.config.enabled

  def _root_join(self, path, makedirs=True):
    fullpath = os.path.join(self.config.output_root, path)
    subdir = os.path.dirname(fullpath)
    if makedirs and subdir:
      tf.gfile.MakeDirs(subdir)
    return fullpath

  def _write_csv_line(self, *args):
    assert self.enabled, 'only callable when diagnostics are on'
    self._csv_writer.writerow(args)

  def _file_for_region(self, region, basename):
    """Returns the path to a file in a region-specific subdirectory."""
    assert self.enabled, 'only callable when diagnostics are on'
    return self._root_join(os.path.join(ranges.to_literal(region), basename))

  def log_realigned_reads(self, region, reads):
    """Logs, if enabled, the realigned reads for region."""
    if self.enabled and self.config.emit_realigned_reads:
      path = self._file_for_region(region, self.realigned_reads_filename)
      with sam.SamWriter(path) as writer:
        for read in reads:
          writer.write(read)

  def log_graph_metrics(self, region, graph, candidate_haplotypes,
                        graph_building_time):
    """Logs, if enabled, graph construction information for region."""
    if self.enabled:
      if graph:
        dest_file = self._file_for_region(region, self.graph_filename)
        with tf.gfile.FastGFile(dest_file, 'w') as f:
          f.write(graph.graphviz())
      self._write_csv_line(
          ranges.to_literal(region), graph.kmer_size if graph else 'NA',
          len(candidate_haplotypes), graph_building_time)


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
  possible that read_span will be smaller than region, which can happen,
  for example, when we only have reads starts in the middle of the region.
  Here's a picture of when this can happen:

  ref      : acgtACGTACgtgt
  region   :     ------
  read1    :       GGa
  read_span:       ---
  """

  def __init__(self, candidate_haplotypes):
    self.candidate_haplotypes = candidate_haplotypes
    self.reads = []
    self._read_span = None

  def __str__(self):
    return ('AssemblyRegion(region={}, span={}) with {} haplotypes and {} '
            'reads').format(
                ranges.to_literal(self.region),
                ranges.to_literal(self.read_span), len(self.haplotypes),
                len(self.reads))

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
      self._read_span = ranges.make_range(spans[0].reference_name,
                                          min(s.start for s in spans),
                                          max(s.end for s in spans))
    return self._read_span

  def add_read(self, read):
    self.reads.append(read)
    self._read_span = None  # Adding a read invalidates our _read_span cache.


def assign_reads_to_assembled_regions(assembled_regions, reads):
  """Assign each read to the maximally overlapped window.

  Args:
    assembled_regions: list[AssemblyRegion], list of AssemblyRegion to assign
      reads to. Does not assume AssemblyRegion are sorted.
    reads: iterable[learning.genomics.genomics.Read], to be processed. Does
      not assume the reads are sorted.

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

  def __init__(self, config, ref_reader):
    """Creates a new Realigner.

    Args:
      config: realigner_pb2.RealignerOptions protobuf.
      ref_reader: GenomeReferenceFai, indexed reference genome to query bases.
    """
    self.config = config
    self.ref_reader = ref_reader
    self.diagnostic_logger = DiagnosticLogger(self.config.diagnostics)

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
            span=window, haplotypes=candidate_haplotypes)
        windows_haplotypes.append(candidate_haplotypes_info)

      self.diagnostic_logger.log_graph_metrics(
          window, graph, candidate_haplotypes, graph_building_time)

    return windows_haplotypes

  def call_aligner(self, assembled_region):
    """Helper function to call aligner module."""
    if not assembled_region.reads:
      return []

    contig = assembled_region.region.reference_name
    ref_start = max(
        0,
        min(assembled_region.read_span.start, assembled_region.region.start) -
        _REF_ALIGN_MARGIN)
    ref_end = min(
        self.ref_reader.contig(contig).n_bases,
        max(assembled_region.read_span.end, assembled_region.region.end) +
        _REF_ALIGN_MARGIN)

    ref_prefix = self.ref_reader.query(
        ranges.make_range(contig, ref_start, assembled_region.region.start))
    ref = self.ref_reader.query(assembled_region.region)

    # If we can't create the ref suffix then return the original alignments.
    if ref_end <= assembled_region.region.end:
      return assembled_region.reads
    else:
      ref_suffix = self.ref_reader.query(
          ranges.make_range(contig, assembled_region.region.end, ref_end))

    ref_region = ranges.make_range(contig, ref_start, ref_end)
    ref_seq = ref_prefix + ref + ref_suffix
    reads_aligner = aligner.Aligner(self.config.aln_config, ref_region, ref_seq)
    return reads_aligner.align_reads([
        ref_prefix + target + ref_suffix
        for target in assembled_region.haplotypes
    ], assembled_region.reads)

  def call_fast_pass_aligner(self, assembled_region):
    """Helper function to call fast pass aligner module."""
    if not assembled_region.reads:
      return []

    contig = assembled_region.region.reference_name
    ref_start = max(
        0,
        min(assembled_region.read_span.start, assembled_region.region.start) -
        _REF_ALIGN_MARGIN)
    ref_end = min(
        self.ref_reader.contig(contig).n_bases,
        max(assembled_region.read_span.end, assembled_region.region.end) +
        _REF_ALIGN_MARGIN)

    ref_prefix = self.ref_reader.query(
        ranges.make_range(contig, ref_start, assembled_region.region.start))
    ref = self.ref_reader.query(assembled_region.region)

    # If we can't create the ref suffix then return the original alignments.
    if ref_end <= assembled_region.region.end:
      return assembled_region.reads
    else:
      ref_suffix = self.ref_reader.query(
          ranges.make_range(contig, assembled_region.region.end, ref_end))

    ref_seq = ref_prefix + ref + ref_suffix

    fast_pass_realigner = fast_pass_aligner.FastPassAligner()
    # Read sizes may vary. We need this for realigner initialization and sanity
    # checks.
    self.config.aln_config.read_size = len(
        assembled_region.reads[0].aligned_sequence)
    fast_pass_realigner.set_options(self.config.aln_config)
    fast_pass_realigner.set_reference(ref_seq)
    fast_pass_realigner.set_ref_start(contig, ref_start)
    fast_pass_realigner.set_haplotypes([
        ref_prefix + target + ref_suffix
        for target in assembled_region.haplotypes
    ])
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
      reads: [`third_party.nucleus.protos.Read` protos]. The
        list of input reads to realign.
      region: A `third_party.nucleus.protos.Range` proto.
        Specifies the region on the genome we should process.

    Returns:
      [realigner_pb2.CandidateHaplotypes]. Information on the list of candidate
        haplotypes.
      [`third_party.nucleus.protos.Read` protos]. The realigned
        reads for the region. NOTE THESE READS MAY NO LONGER BE IN THE SAME
        ORDER AS BEFORE.
    """
    # Compute the windows where we need to assemble in the region.
    candidate_windows = window_selector.select_windows(
        self.config.ws_config, self.ref_reader, reads, region)

    # Assemble each of those regions.
    candidate_haplotypes = self.call_debruijn_graph(candidate_windows, reads)
    # Create our simple container to store candidate / read mappings.
    assembled_regions = [AssemblyRegion(ch) for ch in candidate_haplotypes]

    # Our realigned_reads start off with all of the unassigned reads.
    realigned_reads = assign_reads_to_assembled_regions(assembled_regions,
                                                        reads)

    # Walk over each region and align the reads in that region, adding them to
    # our realigned_reads.
    for assembled_region in assembled_regions:
      if flags.FLAGS.use_fast_pass_aligner:
        realigned_reads_copy = self.call_fast_pass_aligner(assembled_region)
      else:
        realigned_reads_copy = self.call_aligner(assembled_region)

      realigned_reads.extend(realigned_reads_copy)

    self.diagnostic_logger.log_realigned_reads(region, realigned_reads)

    return candidate_haplotypes, realigned_reads
