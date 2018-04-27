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
"""Determine genomic ranges to perform local assembly."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import ranges
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.realigner.python import window_selector as cpp_window_selector


def _candidates_from_reads(config, ref_reader, reads, region):
  """Returns a dictionary mapping positions to non-reference counts.

  Following cigar operations generate candidate position:
    - ALIGNMENT_MATCH, SEQUENCE_MISMATCH, SEQUENCE_MATCH: at mismatch positions
      in the read when compared to the reference sequence.
    - DELETE: at positions within [cigar_start, cigar_start + cigar_len)
    - INSERT, CLIP_SOFT: at positions within
        [cigar_start - cigar_len, cigar_start + cigar_len)

  Following filters are applied:
   - A read with low-quality alignment score is ignored.
   - A variation position corresponding to a low quality base is ignored. For
     multi-base events (e.g., insertions) variant positions are generated for
     all affected bases provided *all* of the read qualities of the
     corresponding bases are high-quality.

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref_reader: GenomeReference. Indexed reference genome to query bases.
    reads: list[nucleus.protos.Read]. The reads we are processing into
      candidate positions.
    region: nucleus.protos.Range. The region we are processing.

  Returns:
    A dictionary. The keys are reference positions within region. The value of
    each key is the number of candidate observations (see above) seen at that
    position summed up over all reads. The dictionary doesn't contain any
    position keys that would have a value of 0 (i.e., it's sparse).
  """
  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      read_requirements=reads_pb2.ReadRequirements(
          min_mapping_quality=config.min_mapq,
          min_base_quality=config.min_base_quality))

  expanded_region = ranges.expand(
      region,
      config.region_expansion_in_bp,
      contig_map=ranges.contigs_dict(ref_reader.header.contigs))
  allele_counter = allelecounter.AlleleCounter(
      ref_reader.get_c_reader(), expanded_region, allele_counter_options)

  for read in reads:
    allele_counter.add(read)

  counts_vec = cpp_window_selector.candidates_from_allele_counter(
      allele_counter)

  return {
      expanded_region.start + i: count
      for i, count in enumerate(counts_vec)
      if count > 0 and ranges.position_overlaps(
          region.reference_name, expanded_region.start + i, region)
  }


def _candidates_to_windows(config, candidate_pos, ref_name):
  """"Process candidate positions to determine windows for local assembly.

  Following filters are applied:
    - Candidate position with low number of supporting reads is ignored.
    - Candidate position with too many of supporting reads is ignored.

  Windows are within range of
    [min(pos) - config.min_windows_distance,
     max(pos) + config.min_windows_distance)

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    candidate_pos: A dictionary with ref_pos as key and number of supporting
      reads as its value.
    ref_name: Reference name, used in setting the output
      genomics.range.reference_name value.

  Returns:
    A sorted list of nucleus.protos.Range protos for all windows in this region.
  """
  windows = []

  def _add_window(start_pos, end_pos):
    windows.append(
        ranges.make_range(ref_name, start_pos - config.min_windows_distance,
                          end_pos + config.min_windows_distance))

  start_pos, end_pos = None, None
  for pos, count in sorted(candidate_pos.iteritems()):
    if count < config.min_num_supporting_reads:
      continue
    if count > config.max_num_supporting_reads:
      continue
    if start_pos is None:
      start_pos = pos
      end_pos = pos
    elif pos > end_pos + config.min_windows_distance:
      _add_window(start_pos, end_pos)
      start_pos = pos
      end_pos = pos
    else:
      end_pos = pos
  if start_pos is not None:
    _add_window(start_pos, end_pos)

  return sorted(windows, key=ranges.as_tuple)


def select_windows(config, ref_reader, reads, region):
  """"Process reads to determine candidate windows for local assembly.

  Windows are within range of
    [0 - config.min_windows_distance, ref_len + config.min_windows_distance)

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref_reader: GenomeReference. Indexed reference genome to query bases.
    reads: A list of genomics.Read records.
    region: nucleus.protos.Range. The region we are processing.

  Returns:
    A list of nucleus.protos.Range protos sorted by their genomic position.
  """
  # This is a fast path for the case where we have no reads, so we have no
  # windows to assemble.
  if not reads:
    return []

  candidates = _candidates_from_reads(config, ref_reader, reads, region)
  return _candidates_to_windows(config, candidates, region.reference_name)
