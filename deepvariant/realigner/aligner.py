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
"""Perform detailed alignment of a query sequence against multiple targets."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections

import copy
import math
import re


from absl import logging

from deepvariant.core.genomics import cigar_pb2
from deepvariant.realigner import utils
from deepvariant.realigner.python import ssw as libssw

_CIGAR_CLIP_OPS = frozenset(
    [cigar_pb2.CigarUnit.CLIP_SOFT, cigar_pb2.CigarUnit.CLIP_HARD])

_SSW_CIGAR_CHARS = {
    'M': cigar_pb2.CigarUnit.ALIGNMENT_MATCH,
    'I': cigar_pb2.CigarUnit.INSERT,
    'D': cigar_pb2.CigarUnit.DELETE
}


class LibSSWPairwiseAligner(object):
  """Class for performing pairwise sequence alignment."""

  # redacted
  # align multiple queries against a single target. So we should swap query and
  # target usage throughout this codebase.
  def __init__(self, query_seq, match, mismatch, gap_open_penalty,
               gap_extend_penalty):
    """Creates a new LibSSWPairwiseAligner, aligning against query_seq.

    Arguments:
      query_seq: str, the string we'll align our targets against.
      match: int > 0, alignment match score.
      mismatch: int > 0, alignment mismatch score.
      gap_open_penalty: int > 0, alignment gap_open penalty score.
      gap_extend_penalty: int > 0, alignment gap_extension penalty score.
    """
    self.filter = libssw.Filter()
    self.aligner = libssw.Aligner.construct(
        match_score=match,
        mismatch_penalty=mismatch,
        gap_opening_penalty=gap_open_penalty,
        gap_extending_penalty=gap_extend_penalty)
    self.query = query_seq

  def align(self, target):
    self.aligner.set_reference_sequence(target)
    return LibSSWAlignmentFacade(self.aligner.align(self.query, self.filter))


class LibSSWAlignmentFacade(object):
  """Facade class for libssw alignments giving same interface as skbio."""

  def __init__(self, libssw_alignment):
    self.libssw_alignment = libssw_alignment
    self._simplified_cigar = None

  @property
  def query_begin(self):
    return self.libssw_alignment.query_begin

  @property
  def query_end(self):
    return self.libssw_alignment.query_end

  @property
  def target_begin(self):
    return self.libssw_alignment.ref_begin

  @property
  def target_end_optimal(self):
    return self.libssw_alignment.ref_end

  @staticmethod
  def _simplify_cigar_string(cigar_string):
    """Simplifies a CIGAR string containing X/=/S codes.

    Args:
      cigar_string: a CIGAR string possibly containing X/=/S codes

    Returns:
      an "equivalent" CIGAR string where X/= ops have been replaced by M, S ops
      have been elided, and successive concordant ops have then been merged.
      Note that S ops are only legal at the beginning and end of the input CIGAR
      string, but we do not assert that condition here.
    """
    simplified_cigar_ops = []  # list of (int, str) tuples
    for op in re.finditer('([0-9]+)([A-Z=]+)', cigar_string):
      op_len = int(op.group(1))
      op_code = op.group(2)
      # Elide soft clips.
      if op_code == 'S':
        continue
      # Replace X/= with M
      if op_code in ['X', '=']:
        op_code = 'M'
      # If same op code as previous op (if there was one), coalesce.
      if simplified_cigar_ops:
        prev_op_len, prev_op_code = simplified_cigar_ops[-1]
        if op_code == prev_op_code:
          simplified_cigar_ops[-1] = (op_len + prev_op_len, prev_op_code)
          continue
      simplified_cigar_ops.append((op_len, op_code))
    return ''.join('{}{}'.format(*op) for op in simplified_cigar_ops)

  @property
  def cigar(self):
    # Not calculated and cached yet?
    if self._simplified_cigar is None:
      self._simplified_cigar = self._simplify_cigar_string(
          self.libssw_alignment.cigar_string)
    return self._simplified_cigar

  @property
  def optimal_alignment_score(self):
    return self.libssw_alignment.sw_score

  def __str__(self):
    return ('LibSSWAlignmentFacade(query_begin={}, query_end={}, '
            'target_begin={}, target_end_optimal={}, cigar={}, '
            'optimal_alignment_score={}'
            .format(self.query_begin, self.query_end, self.target_begin,
                    self.target_end_optimal, self.cigar,
                    self.optimal_alignment_score))


_PAIRWISE_ALIGNERS = {
    'libssw': LibSSWPairwiseAligner,
}


def make_pairwise_aligner(query_seq,
                          match,
                          mismatch,
                          gap_open_penalty,
                          gap_extend_penalty,
                          implementation='libssw'):
  """Factory method for getting a pairwise aligner by implementation name."""
  cls = _PAIRWISE_ALIGNERS[implementation]
  return cls(query_seq, match, mismatch, gap_open_penalty, gap_extend_penalty)


def _sw_start_offsets(target_kmer_index, read_kmers):
  """Returns target start offsets for Smith-Waterman alignment.

  Args:
    target_kmer_index: {str: [int]}, target kmer index.
    read_kmers: [str], read kmers.

  Returns:
    A sorted list of start offsets.
  """
  return sorted({(s - i)
                 for i, kmer in enumerate(read_kmers)
                 for s in target_kmer_index[kmer]})


class SingleAlnOp(object):
  """Encapsulates alignment operation in basepair resolution.

  Each instance specifies the position offset and cigar operation type.
  """

  def __init__(self, pos, cigar_op):
    """Initializes the alignment information at a single position.

    Args:
      pos: int, operation start position.
      cigar_op: cigar_pb2.CigarUnit.Operation, cigar operation.

    Returns:
      None
    """
    self.pos = pos
    self.cigar_op = cigar_op

  def __repr__(self):
    cigar_op_name = cigar_pb2.CigarUnit.Operation.Name(self.cigar_op)
    return 'SingleAlnOp(%i, %s)' % (self.pos, cigar_op_name)

  def end_pos(self):
    """Determines alignment operation end offset w.r.t. the query sequence.

    For example, deletion uses 0 basepair from the sequence so its end
    offset is the same as the start offset.

    Returns:
      int, alignment operation end offset.

    Raises:
      ValueError: if cigar_op isn't one of the expected ops.
    """
    if (self.cigar_op in utils.CIGAR_ALIGN_OPS or
        self.cigar_op in utils.CIGAR_INSERT_OPS):
      return self.pos + 1
    elif (self.cigar_op in utils.CIGAR_DELETE_OPS or
          self.cigar_op in utils.CIGAR_NO_OPS):
      return self.pos
    raise ValueError('Unexpected cigar_op', self.cigar_op)


class Target(object):
  """Target alignment information.

  Attributes:
    sequence: Target sequence.
    k: int, index k-mer size.
    kmer_index: {str: [int]}, mapping each k-mer to its offsets.
    gaps: [SingleAlnOp], target gaps (i.e. indels) with respect to the
      reference.
    ref_pos_mapping: [int], maps target offsets to reference positions.
  """

  def __init__(self, sequence):
    """Initialize the target information.

    Args:
     sequence: str, target sequence.

    Returns:
      None
    """
    self.sequence = str(sequence).upper()
    self.k = None
    self.kmer_index = None
    self.gaps = None
    self.ref_pos_mapping = None

  def build_target_index(self, k):
    """Build a target index, mapping its k-mers to their offsets.

    Indexing fails if the size of sequence is smaller than k.

    Args:
      k: int, indexing size.

    Returns:
      False if indexing fails, otherwise True.
    """
    self.k = k
    self.kmer_index = collections.defaultdict(list)
    if k > len(self.sequence):
      return False

    kmers = [self.sequence[i:i + k] for i in range(len(self.sequence) - k + 1)]
    for i, kmer in enumerate(kmers):
      self.kmer_index[kmer].append(i)
    return True


class Read(object):
  """Read alignment information."""

  def __init__(self, read):
    """Initialize the read information.

    Encapsulates information on
      - original learning.genomics.genomics.core.Read information.
      - read kmers.
      - best realignment info (Target instance, target offset and alignment).

    Args:
      read: learning.genomics.genomics.core.Read, read's original alignment
        record.

    Returns:
      None
    """
    self.original_readalignment = read
    self.kmers = None
    self.target = None
    self.target_offset = None
    self.alignment = None
    self.readalignment = None

  def aligned_seq(self):
    return str(self.original_readalignment.aligned_sequence)

  def set_read_kmers(self, k):
    seq = self.aligned_seq().upper()
    self.kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]

  def update_best_alignment(self, target, target_offset, alignment):
    if (not self.alignment or alignment.optimal_alignment_score >
        self.alignment.optimal_alignment_score):
      self.target = target
      self.target_offset = target_offset
      self.alignment = alignment


class Aligner(object):
  """Perform detailed alignment of a query sequence against multiple targets."""

  def __init__(self, config, ref_region, ref_seq):
    """Initialize aligner information.

    It encapsulates information on:
      - config: realigner configuration.
      - ssw_kwargs: Smith-Waterman algorithm configuration.
      - ref_region: reference name and genomic range.
      - ref_seq: reference sequence.
      - targets: [Aligner.Target], targets information.
      - reads: [Aligner.read], reads information.

    Args:
      config: Realigner.AlignerOptions record.
      ref_region: range_pb2.Range, specifying the reference name and genomic
        range.
      ref_seq: str, reference sequence.

    Returns:
      None
    """
    self.config = config
    self.ref_region = ref_region
    self.ref_seq = ref_seq.upper()
    self.targets = None

  def _make_pairwise_aligner(self, query_seq, implementation='libssw'):
    return make_pairwise_aligner(
        implementation=implementation,
        query_seq=query_seq,
        match=self.config.match,
        mismatch=self.config.mismatch,
        gap_open_penalty=self.config.gap_open,
        gap_extend_penalty=self.config.gap_extend)

  def set_targets(self, target_sequences):
    """Set targets information.

    Note: Target sequences are stored in uppercase.
    Note: The target sequence matching the reference is stored at the end, since
          update_best_alignment() always keeps the first best alignment. We want
          the reference sequence to have the lowest priority, since it is
          included in the candidate haplotypes, even if it is not supported
          by reads.
    Note: Aligns target sequences to the reference and store gaps information.

    Args:
      target_sequences: [str], a list of target sequences.
    """
    # redacted
    if self.ref_seq is not None and self.ref_seq in target_sequences:
      target_sequences.remove(self.ref_seq)
      target_sequences.append(self.ref_seq)

    self.targets = []
    for target_seq in target_sequences:
      target = Target(target_seq)
      if target.build_target_index(self.config.k):
        # redacted
        # replace local aligner ssw.StripedSmithWaterman().
        # redacted
        # single aligner against ref_seq as the target, and then aligning each
        # target to it. But I'm not doing that right now since
        # _set_target_alignment_info would need to be updated.
        alignment = self._make_pairwise_aligner(target.sequence).align(
            self.ref_seq)
        if self._set_target_alignment_info(target, alignment):
          self.targets.append(target)

  def _set_target_alignment_info(self, target, alignment):
    """Set target alignment gaps and positions with respect to the reference.

    Args:
      target: Target, to be processed.
      alignment: ssw.AlignmentStructure, target alignment result against the
        reference.

    Returns:
      False, if setting target alignment fails.

    Raises:
      ValueError: if cigar_op isn't one of the expected ops.
    """
    # Note: query_end is inclusive.
    if (alignment.query_begin != 0 or
        alignment.query_end != len(target.sequence) - 1):
      logging.warning('aligner: Target alignment should be end-to-end.')
      return False
    target.gaps = self._cigar_str_to_gaps(alignment.cigar, 0)

    target.ref_pos_mapping = [None] * len(target.sequence)
    ref_pos = self.ref_region.start
    target_pos = 0
    for gap in target.gaps:
      pre_gap_len = gap.pos - target_pos
      target.ref_pos_mapping[target_pos:gap.pos] = range(
          ref_pos, ref_pos + pre_gap_len)
      ref_pos += pre_gap_len
      target_pos += pre_gap_len
      if gap.cigar_op in utils.CIGAR_INSERT_OPS:
        target.ref_pos_mapping[target_pos] = ref_pos
        target_pos += 1
      elif gap.cigar_op in utils.CIGAR_DELETE_OPS:
        ref_pos += 1
      elif gap.cigar_op in utils.CIGAR_ALIGN_OPS:
        target.ref_pos_mapping[target_pos] = ref_pos
        ref_pos += 1
        target_pos += 1
      else:
        raise ValueError('Unexpected cigar_op', gap.cigar_op)
    pre_gap_len = len(target.sequence) - target_pos
    target.ref_pos_mapping[target_pos:] = range(ref_pos, ref_pos + pre_gap_len)
    assert target.ref_pos_mapping[0] >= self.ref_region.start
    assert target.ref_pos_mapping[-1] <= self.ref_region.end
    return True

  def _indel_to_single_aln_ops(self, start_pos, cigar_op, size):
    if cigar_op == cigar_pb2.CigarUnit.INSERT:
      return [SingleAlnOp(start_pos + i, cigar_op) for i in range(size)]
    elif cigar_op == cigar_pb2.CigarUnit.DELETE:
      return [SingleAlnOp(start_pos, cigar_op) for i in range(size)]
    else:
      raise ValueError('Unexpected cigar_op', cigar_op)

  def _cigar_str_to_gaps(self, cigar, query_offset):
    """Parses alignment cigar string to a list of gap positions.

    Note: starting and ending clipped bases are not included.

    Args:
      cigar: str, cigar string.
      query_offset: int, query start position offset.

    Returns:
     [SingleAlnOp], a sorted list by SingleAlnOp.pos field.
    """
    gaps = []
    cigar_ls = re.split('(' + '|'.join(_SSW_CIGAR_CHARS.keys()) + ')', cigar)
    cigar_ops = [_SSW_CIGAR_CHARS[c] for c in cigar_ls[1:-1:2]]
    cigar_lens = [int(i) for i in cigar_ls[:-1:2]]
    query_pos = 0
    query_pos += query_offset
    for cigar_op, cigar_len in zip(cigar_ops, cigar_lens):
      if cigar_op == cigar_pb2.CigarUnit.INSERT:
        gaps += self._indel_to_single_aln_ops(query_pos, cigar_op, cigar_len)
        query_pos += cigar_len
      elif cigar_op == cigar_pb2.CigarUnit.DELETE:
        gaps += self._indel_to_single_aln_ops(query_pos, cigar_op, cigar_len)
      elif cigar_op == cigar_pb2.CigarUnit.ALIGNMENT_MATCH:
        query_pos += cigar_len
      else:
        assert cigar_op in _SSW_CIGAR_CHARS.values()
    return gaps

  def _read_gaps(self, read):
    """Returns read alignment gaps (i.e. indels) with respect to the reference.

    It merges alignment gaps from read -> target and target -> reference
    and returns a list of read -> reference gaps.

    Args:
      read: Aligner.Read, to be processed.

    Returns:
      [SingleAlnOp], a sorted list by SingleAlnOps.pos field.

    Raises:
      ValueError: if cigar_op isn't one of the expected ops.
    """

    def _read_pos_offset_adjusment(single_aln_op):
      if (single_aln_op.cigar_op in utils.CIGAR_ALIGN_OPS or
          single_aln_op.cigar_op in utils.CIGAR_INSERT_OPS):
        return 1
      elif single_aln_op.cigar_op in utils.CIGAR_DELETE_OPS:
        return -1
      else:
        raise ValueError('Unexpected cigar op', single_aln_op.cigar_op)

    read_gaps = self._cigar_str_to_gaps(read.alignment.cigar,
                                        read.alignment.query_begin)
    target_gaps = read.target.gaps
    target_offset = read.alignment.query_begin - (
        read.alignment.target_begin + read.target_offset)
    gaps = []
    read_i = 0
    target_i = 0
    while True:
      read_gap = (read_gaps[read_i] if read_i < len(read_gaps) else None)
      target_gap = (
          target_gaps[target_i] if target_i < len(target_gaps) else None)
      if not (read_gap or target_gap):
        break
      elif not target_gap:
        gaps.append(read_gap)
        target_offset += _read_pos_offset_adjusment(read_gap)
        read_i += 1
      elif not read_gap:
        if target_gap.pos < read.alignment.target_begin + read.target_offset:
          target_i += 1
        elif (target_gap.pos <=
              read.alignment.target_end_optimal + read.target_offset):
          gaps.append(
              SingleAlnOp(target_gap.pos + target_offset, target_gap.cigar_op))
          target_i += 1
        else:
          target_i = len(target_gaps)
      elif target_gap.pos < read.alignment.target_begin + read.target_offset:
        target_i += 1
      elif (target_gap.pos >
            read.alignment.target_end_optimal + read.target_offset):
        target_i = len(target_gaps)
      elif read_gap.pos < target_gap.pos + target_offset:
        gaps.append(read_gap)
        target_offset += _read_pos_offset_adjusment(read_gap)
        read_i += 1
      elif target_gap.pos + target_offset < read_gap.pos:
        gaps.append(
            SingleAlnOp(target_gap.pos + target_offset, target_gap.cigar_op))
        target_i += 1
      elif read_gap.cigar_op != target_gap.cigar_op:
        assert target_gap.pos + target_offset == read_gap.pos
        target_offset += _read_pos_offset_adjusment(read_gap)
        read_i += 1
        target_i += 1
      elif target_gap.cigar_op == cigar_pb2.CigarUnit.DELETE:
        assert target_gap.pos + target_offset == read_gap.pos
        gaps.append(
            SingleAlnOp(target_gap.pos + target_offset, target_gap.cigar_op))
        target_i += 1
      elif read_gap.cigar_op == cigar_pb2.CigarUnit.INSERT:
        assert target_gap.pos + target_offset == read_gap.pos
        gaps.append(read_gap)
        target_offset += _read_pos_offset_adjusment(read_gap)
        read_i += 1
    return gaps

  def _add_cigar_unit(self, cigar, cigar_op, cigar_len):
    if cigar_len <= 0:
      return
    if cigar and cigar[-1].operation == cigar_op:
      cigar[-1].operation_length += cigar_len
    else:
      cigar_unit = cigar_pb2.CigarUnit(
          operation=cigar_op, operation_length=cigar_len)
      cigar.extend([cigar_unit])

  def _candidate_clipped_bases(self, cigar):
    """Returns starting cigar units that could be treated as clipped bases.

    Args:
      cigar: str, cigar string.

    Returns:
      (int, int), number of cigar units, and read bases associated with clipped
      bases.

    Raises:
      ValueError: if cigar_op isn't one of the expected ops.
    """
    read_offset = 0
    index = 0
    for cigar_unit in cigar:
      if cigar_unit.operation in utils.CIGAR_DELETE_OPS:
        index += 1
      elif cigar_unit.operation in utils.CIGAR_INSERT_OPS:
        read_offset += cigar_unit.operation_length
        index += 1
      elif cigar_unit.operation in utils.CIGAR_ALIGN_OPS:
        break
      else:
        raise ValueError('Unexpected cigar_op', cigar_unit.operation)
    return index, read_offset

  def _gaps_to_cigar_ops(self, gaps, seq_start_offset, seq_end_offset):
    """Converts a list of gaps to a list of cigar units.

    Args:
      gaps: [SingleAlnOp], a sorted list of gaps.
      seq_start_offset: int, sequence start alignment position.
      seq_end_offset: int, sequence end alignment position.

    Returns:
      [cigar_pb2.CigarUnit], a list of cigar units.
    """
    cigar_ops = []
    if not gaps:
      self._add_cigar_unit(cigar_ops, cigar_pb2.CigarUnit.ALIGNMENT_MATCH,
                           seq_end_offset - seq_start_offset)
    else:
      self._add_cigar_unit(cigar_ops, cigar_pb2.CigarUnit.ALIGNMENT_MATCH,
                           gaps[0].pos - seq_start_offset)
      for i, read_gap in enumerate(gaps):
        self._add_cigar_unit(cigar_ops, read_gap.cigar_op, 1)
        try:
          next_gap_pos = gaps[i + 1].pos
        except IndexError:
          next_gap_pos = seq_end_offset
        self._add_cigar_unit(cigar_ops, cigar_pb2.CigarUnit.ALIGNMENT_MATCH,
                             next_gap_pos - read_gap.end_pos())
    return cigar_ops

  def sanity_check_readalignment(self, readalignment):
    """Sanity check readalignment information.

    Validates:
      - alignment position is within the reference range (self.ref_region).
      - cigar matches the read length.

    Args:
      readalignment: learning.genomics.genomics.core.Read, to be validated.

    Raises:
      ValueError: if sanity check fails.

    """
    if not readalignment.HasField('alignment'):
      return

    if (readalignment.alignment.position.reference_name !=
        self.ref_region.reference_name):
      raise ValueError(
          'readalignment validation: '
          'read reference name is inconsistent with reference information.')
    aln_pos = readalignment.alignment.position.position
    if aln_pos < self.ref_region.start:
      raise ValueError('readalignment validation: '
                       'read start position is out of reference genomic range.')

    query_len = 0
    ref_len = 0
    for cigar_op in readalignment.alignment.cigar:
      if cigar_op.operation in utils.CIGAR_ALIGN_OPS:
        ref_len += cigar_op.operation_length
        query_len += cigar_op.operation_length
      elif cigar_op.operation in utils.CIGAR_INSERT_OPS:
        query_len += cigar_op.operation_length
      elif cigar_op.operation in utils.CIGAR_DELETE_OPS:
        ref_len += cigar_op.operation_length
      elif cigar_op.operation in utils.CIGAR_NO_OPS:
        pass
      else:
        raise ValueError('readalignment validation: Unexpected cigar_op %i',
                         self.cigar_op)

    if query_len != len(readalignment.aligned_sequence):
      raise ValueError('readalignment validation: '
                       'cigar is inconsistent with the read length.')
    if aln_pos + ref_len > self.ref_region.end:
      raise ValueError('readalignment validation: '
                       'read end position is out of reference genomic range.')

  def set_read_cigar(self, read):
    """Set read's cigar information, given alignment result.

    Final realignment cigar string is computed by merging read -> target and
    target -> reference alignment gaps.

    Args:
      read: Aligner.Read, to be processed.

    Updates:
      read's alignment.cigar field.
    """
    read.readalignment.alignment.ClearField('cigar')

    read_gaps = self._read_gaps(read)
    cigar = self._gaps_to_cigar_ops(
        read_gaps,
        read.alignment.query_begin,
        # Note: ssw.TabularMSA.query_end is inclusive.
        read.alignment.query_end + 1)

    cigar_start_i, start_read_clipped_bases = self._candidate_clipped_bases(
        cigar)
    if cigar_start_i < len(cigar):
      tmp_end_i, end_read_clipped_bases = self._candidate_clipped_bases(
          cigar[::-1])
      cigar_end_i = len(cigar) - tmp_end_i
      assert cigar_start_i < cigar_end_i
    else:
      cigar_end_i, end_read_clipped_bases = len(cigar), 0

    # Construct the cigar string, including hard and soft clipped bases.
    if (read.original_readalignment.alignment.cigar[0].operation ==
        cigar_pb2.CigarUnit.CLIP_HARD):
      self._add_cigar_unit(
          read.readalignment.alignment.cigar, cigar_pb2.CigarUnit.CLIP_HARD,
          read.original_readalignment.alignment.cigar[0].operation_length)
    self._add_cigar_unit(read.readalignment.alignment.cigar,
                         cigar_pb2.CigarUnit.CLIP_SOFT,
                         read.alignment.query_begin)
    self._add_cigar_unit(read.readalignment.alignment.cigar,
                         cigar_pb2.CigarUnit.CLIP_SOFT,
                         start_read_clipped_bases)
    read.readalignment.alignment.cigar.extend(cigar[cigar_start_i:cigar_end_i])
    self._add_cigar_unit(read.readalignment.alignment.cigar,
                         cigar_pb2.CigarUnit.CLIP_SOFT, end_read_clipped_bases)
    self._add_cigar_unit(read.readalignment.alignment.cigar,
                         cigar_pb2.CigarUnit.CLIP_SOFT,
                         len(read.aligned_seq()) - read.alignment.query_end - 1)
    if (read.original_readalignment.alignment.cigar[-1].operation ==
        cigar_pb2.CigarUnit.CLIP_HARD):
      self._add_cigar_unit(
          read.readalignment.alignment.cigar, cigar_pb2.CigarUnit.CLIP_HARD,
          read.original_readalignment.alignment.cigar[-1].operation_length)

  def set_readalignment(self, read):
    """Set read's readalignment information, given alignment results.

    Args:
      read: Aligner.Read, to be processed.

    Updates:
      Read.readalignment field.
    """
    if not read.alignment:
      return

    read.readalignment = copy.deepcopy(read.original_readalignment)
    read.readalignment.alignment.position.position = (
        read.target.ref_pos_mapping[read.target_offset
                                    + read.alignment.target_begin])
    self.set_read_cigar(read)
    self.sanity_check_readalignment(read.readalignment)

    # If the whole read sequence is clipped, reset the readalignment to None.
    if all(cigar_op.operation in _CIGAR_CLIP_OPS
           for cigar_op in read.readalignment.alignment.cigar):
      read.readalignment = None

  # redacted
  # functions that serve only to figure out candidate start/end positions in
  # target_seq where we should consider aligning read_seq. This logic could be
  # cleaned up.
  def _ssw_alignment(self, aligner, read_seq, target_seq, target_offset):
    """Perform Smith-Waterman alignment.

    This function aligns read_seq to a subsequence of target_seq starting at
    target_offset. The function continues to grow the length of the subsequence,
    until an end-to-end alignment of the read sequence can be found w.r.t. the
    target sequence.

    Args:
      aligner: PairwiseAligner to use for alignment.
      read_seq: str, read sequence.
      target_seq: str, target sequence.
      target_offset: int, target start offset.

    Returns:
      Tuple of (target_start_offset, ssw.AlignmentStructure) or None if an
      alignment is not found.
    """
    # redacted
    # really doesn't make sense to call this an error rate.
    terminal_seq_len = int(math.ceil(len(read_seq) * self.config.error_rate))
    while True:
      target_start = target_offset - terminal_seq_len
      target_end = target_offset + len(read_seq) + terminal_seq_len
      if target_end > len(target_seq) or target_start < 0:
        return None, None
      alignment = aligner.align(target_seq[target_start:target_end])
      if (alignment.target_end_optimal == target_end - target_start - 1 or
          alignment.target_begin == 0):
        terminal_seq_len *= 2
      else:
        return target_start, alignment

  def realign_read(self, read, implementation='libssw'):
    """Realign read to best scoring target sequence.

    Align read to all targets and find the one with the highest
    alignment score. Final read cigar is generated based on read's alignment
    to the target sequence and the target's variations with respect to
    the reference.

    Args:
      read: Read, to be processed.
      implementation: the aligner implementation to use.  'libssw' currently
        only available option.
    Returns:
      Updated read with realignment information.
    """
    read.set_read_kmers(self.config.k)
    aligner = self._make_pairwise_aligner(
        read.aligned_seq(), implementation=implementation)
    for target in self.targets:
      read_seq = read.aligned_seq()
      target_offsets = _sw_start_offsets(target.kmer_index, read.kmers)
      for target_offset in target_offsets:
        alignment_offset, alignment = self._ssw_alignment(
            aligner, read_seq, target.sequence, target_offset)
        if alignment:
          read.update_best_alignment(target, alignment_offset, alignment)
    self.set_readalignment(read)
    return read

  def align_reads(self, targets, reads):
    """Re-align reads to reference given candidate targets.

    If a read's realignment is not successful, its original alignments is
    returned.

    Note: Reference and target sequences are stored in uppercase.

    Args:
      targets: [str], a list of target sequences.
      reads: [learning.genomics.genomics.core.Read], a list of reads to
        perform realignment.

    Returns:
      [learning.genomics.genomics.core.Read], a list of read realignments.
      If there is no target sequence or the only target sequence is the
      reference, original input reads are returned.
    """
    # Remove any duplicates but make sure it's sorted.
    targets = sorted(set(targets))
    if not targets or targets == [self.ref_seq]:
      return reads
    self.set_targets(targets)
    realigned_reads = []
    for read in reads:
      realigned_read = self.realign_read(Read(read))
      if realigned_read.readalignment:
        realigned_reads.append(realigned_read.readalignment)
      else:
        realigned_reads.append(realigned_read.original_readalignment)

    return realigned_reads
