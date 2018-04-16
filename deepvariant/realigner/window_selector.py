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

from collections import defaultdict

from third_party.nucleus.protos import cigar_pb2
from third_party.nucleus.util import ranges
from deepvariant.realigner import utils


def _process_align_match(config, cigar, ref, read, ref_pos, read_pos):
  for _ in range(cigar.operation_length):
    # Check for ref_pos ranges, if the read is partially overlapped with
    # ref sequence, it might be out-of-range and result in index error.
    if ref_pos >= len(ref):
      break
    if ref_pos >= 0:
      if ref[ref_pos] != read.aligned_sequence[read_pos]:
        if read.aligned_quality[read_pos] >= config.min_base_quality:
          yield ref_pos
    read_pos += 1
    ref_pos += 1


def _process_seq_mismatch(config, cigar, ref, read, ref_pos, read_pos):
  del ref  # Unused in processing cigar seq_mismatch operation.
  for _ in range(cigar.operation_length):
    if read.aligned_quality[read_pos] >= config.min_base_quality:
      yield ref_pos
    read_pos += 1
    ref_pos += 1


def _process_insert(config, cigar, ref, read, ref_pos, read_pos):
  del ref  # Unused in processing cigar insert operation.
  for i in range(cigar.operation_length):
    if read.aligned_quality[read_pos] >= config.min_base_quality:
      yield ref_pos + i
      yield ref_pos - cigar.operation_length + i
    read_pos += 1


def _process_soft_clip(config, cigar, ref, read, ref_pos, read_pos):
  del ref  # Unused in processing cigar soft_clip operation.
  if ref_pos == read.alignment.position.position:
    offset = -cigar.operation_length
  else:
    offset = 0
  for i in range(cigar.operation_length):
    if read.aligned_quality[read_pos] >= config.min_base_quality:
      yield ref_pos + offset + i
    read_pos += 1


def _process_delete(config, cigar, ref, read, ref_pos, read_pos):
  # Unused in processing cigar delete operation.
  del config, ref, read, read_pos
  for _ in range(cigar.operation_length):
    yield ref_pos
    ref_pos += 1


def _process_read(config, ref, read, ref_offset):
  """Yields reference positions corresponding to read's variations.

  See _candidates_from_reads for more details on this algorithm.

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref: reference sequence
    read: A genomics.Read record.
    ref_offset: Start offset for reference position.

  Yields:
    reference positions within range of [0, seq_len], corresponding
    to read's variations.

  Raises:
    ValueError: for unsupported cigar operation.
  """

  if read.alignment.mapping_quality < config.min_mapq:
    return

  ref_pos = read.alignment.position.position - ref_offset
  read_pos = 0
  # Use set(), as some cigar operations might generate duplicated positions,
  # E.g. for insertions, it extends the candidate positions to
  # [ins_pos - ins_len, ins_pos + ins_len] which might overlap with some
  # nearby mismatches.
  positions = set()
  for cigar in read.alignment.cigar:
    # Break if it reached the end of reference sequence.
    if ref_pos >= len(ref):
      break
    if cigar.operation not in utils.CIGAR_OPS:
      raise ValueError('Unexpected CIGAR operation', cigar, read)

    if cigar.operation == cigar_pb2.CigarUnit.ALIGNMENT_MATCH:
      positions.update(
          _process_align_match(config, cigar, ref, read, ref_pos, read_pos))
      read_pos += cigar.operation_length
      ref_pos += cigar.operation_length
    elif cigar.operation == cigar_pb2.CigarUnit.SEQUENCE_MISMATCH:
      positions.update(
          _process_seq_mismatch(config, cigar, ref, read, ref_pos, read_pos))
      read_pos += cigar.operation_length
      ref_pos += cigar.operation_length
    elif cigar.operation == cigar_pb2.CigarUnit.INSERT:
      positions.update(
          _process_insert(config, cigar, ref, read, ref_pos, read_pos))
      read_pos += cigar.operation_length
    elif cigar.operation == cigar_pb2.CigarUnit.CLIP_SOFT:
      positions.update(
          _process_soft_clip(config, cigar, ref, read, ref_pos, read_pos))
      read_pos += cigar.operation_length
    # redacted
    elif (cigar.operation == cigar_pb2.CigarUnit.DELETE or
          cigar.operation == cigar_pb2.CigarUnit.SKIP):
      positions.update(
          _process_delete(config, cigar, ref, read, ref_pos, read_pos))
      ref_pos += cigar.operation_length
    elif cigar.operation == cigar_pb2.CigarUnit.SEQUENCE_MATCH:
      ref_pos += cigar.operation_length
      read_pos += cigar.operation_length
    elif (cigar.operation == cigar_pb2.CigarUnit.CLIP_HARD or
          cigar.operation == cigar_pb2.CigarUnit.PAD):
      pass

  # Yield positions within the range
  for pos in sorted(positions):
    if pos >= 0 and pos < len(ref):
      yield pos + ref_offset


def _candidates_from_reads(config, ref, reads, region):
  """Returns a dictionary mapping positions to non-reference counts.

  Following cigar operations generate candidate position:
    - ALIGNMENT_MATCH: at mismatch positions
    - SEQUENCE_MISMATCH: at positions within
        [cigar_start, cigar_start + cigar_len)
    - DELETE: at positions within [cigar_start, cigar_start + cigar_len)
    - INSERT: at positions within
        [cigar_start - cigar_len, cigar_start + cigar_len)
    - CLIP_SOFT: at positions within [cigar_start, cigar_start + cigar_len)
    - SKIP: at positions within [cigar_start, cigar_start + cigar_len)

  Following filters are applied:
   - A read with low-quality alignment score is ignored.
   - A variation position corresponding to a low quality base is ignored.

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref: str. The reference bases as a string, where the first character of
      ref starts at region.start and the last position in region.end.
    reads: list[nucleus.protos.Read]. The reads we are processing into
      candidate positions.
    region: nucleus.protos.Range. The region we are processing.

  Returns:
    A dictionary. The keys are reference positions within region. The value of
    each key is the number of candidate observations (see above) seen at that
    position summed up over all reads. The dictionary doesn't contain any
    position keys that would have a value of 0 (i.e., it's sparse).
  """
  # A list of candidate positions mapping to their number of supporting reads.
  candidates = defaultdict(int)

  for read in reads:
    for ref_pos in _process_read(config, ref, read, region.start):
      candidates[ref_pos] += 1

  return candidates


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


# redacted
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
  ref = ref_reader.query(region)
  candidates = _candidates_from_reads(config, ref, reads, region)
  return _candidates_to_windows(config, candidates, region.reference_name)
