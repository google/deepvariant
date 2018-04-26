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

import collections

from absl import logging

from third_party.nucleus.protos import cigar_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util import ranges
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter
from deepvariant.realigner import utils
from deepvariant.realigner.python import window_selector as cpp_window_selector
from deepvariant.vendor import timer

# redacted
# If True, we will run both the C++ and Python versions of the window selector
# and compare the results.
_COMPARE_CPP_PY = False

# redacted
# If True, we will use the C++ version of the window selector code. If False,
# we will use the original Python version.
_USE_CPP_IMPLEMENTATION = True


def _process_aligned_bases(config, cigar, ref, read, ref_pos, read_pos):
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


def _process_insert_soft_clip(config, cigar, ref, read, ref_pos, read_pos):
  del ref  # Unused in processing cigar insert operation.
  if all(read.aligned_quality[read_pos + i] >= config.min_base_quality
         for i in range(cigar.operation_length)):
    for i in range(cigar.operation_length):
      yield ref_pos + i
      yield ref_pos - cigar.operation_length + i


def _process_delete(config, cigar, ref, read, ref_pos, read_pos):
  # Unused in processing cigar delete operation.
  del config, ref, read, read_pos
  for _ in range(cigar.operation_length):
    yield ref_pos
    ref_pos += 1


def _process_read(config, ref, read, region):
  """Yields reference positions corresponding to read's variations.

  See _candidates_from_reads for more details on this algorithm.

  Args:
    config: learning.genomics.deepvariant.realigner.WindowSelectorOptions
      options determining the behavior of this window selector.
    ref: reference sequence
    read: A genomics.Read record.
    region: nucleus.protos.Range. The region we are constructing windows in.

  Yields:
    A tuple of (reference_position, count). The reference positions is within
    the range of [0, seq_len], corresponding to read's variations, and count is
    the number of "non-reference" events that occurred at that position.

  Raises:
    ValueError: for unsupported cigar operation.
  """

  if read.alignment.mapping_quality < config.min_mapq:
    return

  region_size = region.end - region.start
  ref_pos = read.alignment.position.position - region.start
  read_pos = 0
  positions = []
  for cigar in read.alignment.cigar:
    # Break if it reached the end of reference sequence.
    if ref_pos >= region_size:
      break
    if cigar.operation not in utils.CIGAR_OPS:
      raise ValueError('Unexpected CIGAR operation', cigar, read)

    if cigar.operation in {
        cigar_pb2.CigarUnit.ALIGNMENT_MATCH, cigar_pb2.CigarUnit.SEQUENCE_MATCH,
        cigar_pb2.CigarUnit.SEQUENCE_MISMATCH
    }:
      positions.extend(
          _process_aligned_bases(config, cigar, ref, read, ref_pos, read_pos))
      read_pos += cigar.operation_length
      ref_pos += cigar.operation_length
    elif cigar.operation in {
        cigar_pb2.CigarUnit.INSERT, cigar_pb2.CigarUnit.CLIP_SOFT
    }:
      positions.extend(
          _process_insert_soft_clip(config, cigar, ref, read, ref_pos,
                                    read_pos))
      read_pos += cigar.operation_length
    elif cigar.operation == cigar_pb2.CigarUnit.DELETE:
      positions.extend(
          _process_delete(config, cigar, ref, read, ref_pos, read_pos))
      ref_pos += cigar.operation_length
    elif (cigar.operation == cigar_pb2.CigarUnit.CLIP_HARD or
          cigar.operation == cigar_pb2.CigarUnit.PAD or
          cigar.operation == cigar_pb2.CigarUnit.SKIP):
      pass

  # Yield positions within the range.
  for pos, count in collections.Counter(positions).iteritems():
    if pos >= 0 and pos < region_size:
      yield pos + region.start, count


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
  if _COMPARE_CPP_PY or _USE_CPP_IMPLEMENTATION:
    with timer.Timer() as t2:
      cpp_candidates = _candidates_from_reads_cpp(config, ref_reader, reads,
                                                  region)
      candidates = cpp_candidates

  if _COMPARE_CPP_PY or not _USE_CPP_IMPLEMENTATION:
    with timer.Timer() as t1:
      python_candidates = _candidates_from_reads_python(config, ref_reader,
                                                        reads, region)
      if not _USE_CPP_IMPLEMENTATION:
        candidates = python_candidates

  if _COMPARE_CPP_PY:
    logging.info('Region %s', ranges.to_literal(region))
    logging.info('N. reads %s', len(reads))
    logging.info('Python candidates %s', python_candidates)
    logging.info('C++    candidates %s', cpp_candidates)
    logging.info('Speed up: python: %5f, c++ %5f, relative cost %2f%%',
                 t1.GetDuration(), t2.GetDuration(),
                 (100.0 * t2.GetDuration()) / t1.GetDuration())
    if python_candidates != cpp_candidates:
      s = _show_diff(python_candidates, cpp_candidates)
      # raise ValueError('Python and C++ candidates are different: %s' % s)
      logging.warn('Python and C++ candidates are different: %s', s)

  return candidates


# redacted
def _show_diff(py_candidates, cpp_candidates):
  keys = set(py_candidates) | set(cpp_candidates)
  result = {}
  for key in keys:
    pyr = py_candidates.get(key)
    cppr = cpp_candidates.get(key)
    result[key] = (pyr, cppr, '*** BAD' if pyr != cppr else '')
  return '\n'.join('@{}: py={} cpp={} {}'.format(k, p, c, d)
                   for k, (p, c, d) in sorted(result.iteritems()))


# redacted
def _candidates_from_reads_python(config, ref_reader, reads, region):
  """See candidates_from_reads."""
  # A list of candidate positions mapping to their number of supporting reads.
  candidates = collections.defaultdict(int)

  ref = ref_reader.query(region)
  for read in reads:
    for ref_pos, count in _process_read(config, ref, read, region):
      candidates[ref_pos] += count

  return candidates


# redacted
def _candidates_from_reads_cpp(config, ref_reader, reads, region):
  """See candidates_from_reads."""
  allele_counter_options = deepvariant_pb2.AlleleCounterOptions(
      read_requirements=reads_pb2.ReadRequirements(
          min_mapping_quality=config.min_mapq,
          min_base_quality=config.min_base_quality))

  allele_counter = allelecounter.AlleleCounter(ref_reader.get_c_reader(),
                                               region, allele_counter_options)

  for read in reads:
    allele_counter.add(read)

  counts_vec = cpp_window_selector.candidates_from_allele_counter(
      allele_counter)
  return {
      region.start + i: count for i, count in enumerate(counts_vec) if count > 0
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
