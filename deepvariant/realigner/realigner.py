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


import tensorflow as tf

from deepvariant.core import genomics_io
from deepvariant.core import ranges
from deepvariant.core import utils
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import aligner
from deepvariant.realigner import window_selector
from deepvariant.realigner.python import debruijn_graph
from deepvariant.vendor import timer

tf.flags.DEFINE_integer(
    'ws_min_num_supporting_reads', 3,
    'Minimum number of supporting reads to call a reference position for local '
    'assembly.')
tf.flags.DEFINE_integer(
    'ws_max_num_supporting_reads', 300,
    'Maximum number of supporting reads to call a reference position for local '
    'assembly.')
tf.flags.DEFINE_integer(
    'ws_min_mapq', 20,
    'Minimum read alignment quality to consider in calling a reference '
    'position for local assembly.')
tf.flags.DEFINE_integer(
    'ws_min_base_quality', 20,
    'Minimum base quality to consider in calling a reference position for '
    'local assembly.')
tf.flags.DEFINE_integer(
    'ws_min_windows_distance', 70,
    'Minimum distance between candidate windows for local assembly.')
tf.flags.DEFINE_integer(
    'ws_max_window_size', 1000,
    'Maximum window size to consider for local assembly. Large noisy regions '
    'are skipped for realignment.')
tf.flags.DEFINE_integer('dbg_min_k', 10,
                        'Initial k-mer size to build the graph.')
tf.flags.DEFINE_integer(
    'dbg_max_k', 100,
    'Maximum k-mer size. Larger k-mer size is used to resolve graph cycles.')
tf.flags.DEFINE_integer(
    'dbg_step_k', 1, 'Increment size for k to try in resolving graph cycles.')
tf.flags.DEFINE_integer(
    'dbg_min_mapq', 14,
    'Minimum read alignment quality to consider in building the graph.')
tf.flags.DEFINE_integer(
    'dbg_min_base_quality', 17,
    'Minimum base quality in a k-mer sequence to consider in building the '
    'graph.')
tf.flags.DEFINE_integer('dbg_min_edge_weight', 2,
                        'Minimum number of supporting reads to keep an edge.')
tf.flags.DEFINE_integer(
    'dbg_max_num_paths', 256,
    'Maximum number of paths within a graph to consider for realignment. '
    'Set max_num_paths to 0 to have unlimited number of paths.')
tf.flags.DEFINE_integer('aln_match', 4,
                        'Match score (expected to be a non-negative score).')
tf.flags.DEFINE_integer('aln_mismatch', 6,
                        'Mismatch score (expected to be a non-negative score).')
tf.flags.DEFINE_integer(
    'aln_gap_open', 8, 'Gap open score (expected to be a non-negative score). '
    'Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).')
tf.flags.DEFINE_integer(
    'aln_gap_extend', 1,
    'Gap extend score (expected to be a non-negative score). '
    'Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).')
tf.flags.DEFINE_integer('aln_k', 23,
                        'k-mer size used to index target sequence.')
tf.flags.DEFINE_float('aln_error_rate', .01, 'Estimated sequencing error rate.')
tf.flags.DEFINE_string(
    'realigner_diagnostics', '',
    'Root directory where the realigner should place diagnostic output (such as'
    ' a dump of the DeBruijn graph, and a log of metrics reflecting the graph '
    'and  realignment to the haplotypes).  If empty, no diagnostics are output.'
)
tf.flags.DEFINE_bool(
    'emit_realigned_reads', False,
    'If True, we will emit realigned reads if our realigner_diagnostics are '
    'also enabled.')

# Margin added to the reference sequence for the aligner module.
_REF_ALIGN_MARGIN = 20

# ---------------------------------------------------------------------------
# Set configuration settings.
# ---------------------------------------------------------------------------


def realigner_config(flags):
  """Creates a RealignerOptions proto based on input and default settings.

  Args:
    flags: configuration FLAGS.

  Returns:
    realigner_pb2.RealignerOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  ws_config = realigner_pb2.RealignerOptions.WindowSelectorOptions(
      min_num_supporting_reads=flags.ws_min_num_supporting_reads,
      max_num_supporting_reads=flags.ws_max_num_supporting_reads,
      min_mapq=flags.ws_min_mapq,
      min_base_quality=flags.ws_min_base_quality,
      min_windows_distance=flags.ws_min_windows_distance,
      max_window_size=flags.ws_max_window_size)

  dbg_config = realigner_pb2.RealignerOptions.DeBruijnGraphOptions(
      min_k=flags.dbg_min_k,
      max_k=flags.dbg_max_k,
      step_k=flags.dbg_step_k,
      min_mapq=flags.dbg_min_mapq,
      min_base_quality=flags.dbg_min_base_quality,
      min_edge_weight=flags.dbg_min_edge_weight,
      max_num_paths=flags.dbg_max_num_paths)

  aln_config = realigner_pb2.RealignerOptions.AlignerOptions(
      match=flags.aln_match,
      mismatch=flags.aln_mismatch,
      gap_open=flags.aln_gap_open,
      gap_extend=flags.aln_gap_extend,
      k=flags.aln_k,
      error_rate=flags.aln_error_rate)

  diagnostics = realigner_pb2.RealignerOptions.Diagnostics(
      enabled=bool(flags.realigner_diagnostics),
      output_root=flags.realigner_diagnostics,
      emit_realigned_reads=flags.emit_realigned_reads)

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
      with genomics_io.make_read_writer(path) as writer:
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
    self.window_sel = window_selector.WindowSelector(self.config.ws_config)
    self.diagnostic_logger = DiagnosticLogger(self.config.diagnostics)

  def call_window_selector(self, region, reads):
    """Helper function to call window_selector module."""
    return sorted(
        self.window_sel.process_reads(
            self.ref_reader.bases(region), reads, region.reference_name,
            region.start),
        key=ranges.as_tuple)

  def call_debruijn_graph(self, windows, reads):
    """Helper function to call debruijn_graph module."""
    windows_haplotypes = []
    # Build and process de-Bruijn graph for each window.
    for window in windows:
      if window.end - window.start > self.config.ws_config.max_window_size:
        continue
      if not self.ref_reader.is_valid_interval(window):
        continue
      ref = self.ref_reader.bases(window)
      # redacted
      dbg_reads = [
          read for read in reads
          if ranges.ranges_overlap(window, utils.read_range(read))
      ]

      with timer.Timer() as t:
        graph = debruijn_graph.build(ref, dbg_reads, self.config.dbg_config)
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

    ref_prefix = self.ref_reader.bases(
        ranges.make_range(contig, ref_start, assembled_region.region.start))
    ref = self.ref_reader.bases(assembled_region.region)

    # If we can't create the ref suffix then return the original alignments.
    if ref_end <= assembled_region.region.end:
      return assembled_region.reads
    else:
      ref_suffix = self.ref_reader.bases(
          ranges.make_range(contig, assembled_region.region.end, ref_end))

    ref_region = ranges.make_range(contig, ref_start, ref_end)
    ref_seq = ref_prefix + ref + ref_suffix
    reads_aligner = aligner.Aligner(self.config.aln_config, ref_region, ref_seq)
    return reads_aligner.align_reads([
        ref_prefix + target + ref_suffix
        for target in assembled_region.haplotypes
    ], assembled_region.reads)

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
      reads: [`learning.genomics.deepvariant.core.genomics.Read` protos]. The
        list of input reads to realign.
      region: A `learning.genomics.deepvariant.core.genomics.Range` proto.
        Specifies the region on the genome we should process.

    Returns:
      [realigner_pb2.CandidateHaplotypes]. Information on the list of candidate
        haplotypes.
      [`learning.genomics.deepvariant.core.genomics.Read` protos]. The realigned
        reads for the region. NOTE THESE READS MAY NO LONGER BE IN THE SAME
        ORDER AS BEFORE.
    """
    # Compute the windows where we need to assemble in the region.
    candidate_windows = self.call_window_selector(region, reads)
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
      realigned_reads.extend(self.call_aligner(assembled_region))

    self.diagnostic_logger.log_realigned_reads(region, realigned_reads)

    return candidate_haplotypes, realigned_reads
