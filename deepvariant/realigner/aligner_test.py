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
"""Tests for deepvariant .realigner.aligner."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import test_utils
from deepvariant.core import cigar as _cigar
from deepvariant.core import ranges
from deepvariant.core.genomics import cigar_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import aligner


class PairwiseAlignerTest(parameterized.TestCase):

  def test_sw_start_offsets(self):
    """Test Aligner._sw_start_offsets()."""
    k = 3
    read = aligner.Read(
        test_utils.make_read(
            'AaGAt', start=0, cigar=[(5, 'M')], quals=[64] * 5, name='read_1'))
    read.set_read_kmers(k)
    target = aligner.Target('TgATCAGATAAG')
    target.build_target_index(k)
    self.assertEqual([-1, 4, 9],
                     aligner._sw_start_offsets(target.kmer_index, read.kmers))

  def _align_seqs(self,
                  seq1,
                  seq2,
                  match=1,
                  mismatch=2,
                  gap_open=4,
                  gap_extend=1):
    pw_aligner = aligner.make_pairwise_aligner(
        seq1,
        match=match,
        mismatch=mismatch,
        gap_open_penalty=gap_open,
        gap_extend_penalty=gap_extend)
    return pw_aligner.align(seq2)

  def assertAlignmentEqual(self, alignment, expected_dict):
    for k, v in expected_dict.iteritems():
      self.assertEqual(
          getattr(alignment, k), v,
          'Expected field {} of alignment [{}] to be {} but got {}'.format(
              k, alignment, v, getattr(alignment, k)))

  @parameterized.parameters(
      # Try out two exact matches of different lengths.
      dict(
          seq1='ttAtt',
          seq2='ttAtt',
          params=dict(match=1, mismatch=2),
          expected_score=5,
          expected_cigar='5M'),
      dict(
          seq1='tttAttt',
          seq2='tttAttt',
          params=dict(match=1, mismatch=2),
          expected_score=7,
          expected_cigar='7M'),

      # One mismatch in the middle.
      dict(
          seq1='tttAttt',
          seq2='tttCttt',
          params=dict(match=1, mismatch=2),
          expected_score=6 - 2,
          expected_cigar='7M'),
      # Check that match and mismatch scores are respected.
      dict(
          seq1='tttAttt',
          seq2='tttCttt',
          params=dict(match=2, mismatch=3),
          expected_score=6 * 2 - 3,
          expected_cigar='7M'),
      dict(
          seq1='tttAttt',
          seq2='tttCttt',
          params=dict(match=4, mismatch=2),
          expected_score=6 * 4 - 2,
          expected_cigar='7M'),

      # Now for some insertion/deletions.
      dict(
          seq1='ttAtt',
          seq2='tttt',
          params=dict(match=4, mismatch=2, gap_open=4, gap_extend=2),
          expected_score=4 * 4 - 4,
          expected_cigar='2M1I2M'),

      # Same as above sequence, but reversed to the cigar is different.
      dict(
          seq1='tttt',
          seq2='ttAtt',
          params=dict(match=4, mismatch=2, gap_open=4, gap_extend=2),
          expected_score=4 * 4 - 4,
          expected_cigar='2M1D2M'),

      # Gap extension is respected.
      dict(
          seq1='ttAAtt',
          seq2='tttt',
          params=dict(match=4, mismatch=2, gap_open=2, gap_extend=1),
          expected_score=4 * 4 - 2 * 1 - 1 * 1,
          expected_cigar='2M2I2M'),
      dict(
          seq1='ttAAAtt',
          seq2='tttt',
          params=dict(match=4, mismatch=2, gap_open=2, gap_extend=1),
          expected_score=4 * 4 - 2 * 1 - 2 * 1,
          expected_cigar='2M3I2M'),
  )
  def test_pairwise_alignment(self, seq1, seq2, params, expected_cigar,
                              expected_score):
    alignment = self._align_seqs(seq1, seq2, **params)
    self.assertAlignmentEqual(alignment,
                              dict(
                                  query_begin=0,
                                  query_end=len(seq1) - 1,
                                  target_begin=0,
                                  cigar=expected_cigar,
                                  target_end_optimal=len(seq2) - 1,
                                  optimal_alignment_score=expected_score))

  @parameterized.parameters(
      dict(
          query='ttACT',
          target='ttACTtt',
          expected=dict(
              query_begin=0,
              query_end=4,
              target_begin=0,
              target_end_optimal=4,
              cigar='5M')),
      dict(
          query='ACTtt',
          target='ttACTtt',
          expected=dict(
              query_begin=0,
              query_end=4,
              target_begin=2,
              target_end_optimal=6,
              cigar='5M')),
      dict(
          query='ACT',
          target='ttACTtt',
          expected=dict(
              query_begin=0,
              query_end=2,
              target_begin=2,
              target_end_optimal=4,
              cigar='3M')),
      dict(
          query='ACTtt',
          target='ACT',
          expected=dict(
              query_begin=0,
              query_end=2,
              target_begin=0,
              target_end_optimal=2,
              cigar='3M')),
      dict(
          query='ttACT',
          target='ACT',
          expected=dict(
              query_begin=2,
              query_end=4,
              target_begin=0,
              target_end_optimal=2,
              cigar='3M')),
  )
  def test_start_ends(self, query, target, expected):
    alignment = self._align_seqs(query, target)
    self.assertAlignmentEqual(alignment, expected)


class AlignerTest(parameterized.TestCase):

  def make_test_aligner(self, ref_seq=None, region=None):
    config = realigner_pb2.RealignerOptions.AlignerOptions(
        match=1, mismatch=1, gap_open=2, gap_extend=1, k=3, error_rate=.02)
    ref_seq = ref_seq or 'AAAAAAA'
    region = region or ranges.make_range('ref', 10, 10 + len(ref_seq))
    return aligner.Aligner(config, region, ref_seq)

  @parameterized.parameters(
      ('AATA', {
          'AAT': [0],
          'ATA': [1]
      }),
      ('ATcATCA', {
          'ATC': [0, 3],
          'TCA': [1, 4],
          'CAT': [2]
      }),
      ('AAtA', {
          'AAT': [0],
          'ATA': [1]
      }),
      ('AT', None),
  )
  def test_build_target_index(self, seq, expected_kmers_or_none):
    """Test Aligner.set_targets()."""
    target = aligner.Target(seq)
    result = target.build_target_index(k=3)
    self.assertEqual(result, bool(expected_kmers_or_none))
    if expected_kmers_or_none is None:
      self.assertEqual(len(target.kmer_index), 0)
    else:
      self.assertEqual(expected_kmers_or_none, target.kmer_index)

  @parameterized.parameters(
      ('ATCATCA', 'ATcATCA', [], []),
      ('ATCAAATTTCA', 'ATCAATTTCA', [3], [cigar_pb2.CigarUnit.DELETE]),
      ('ATCAAATTTCA', 'ATCTAAATTTCA', [3], [cigar_pb2.CigarUnit.INSERT]),
      # redacted
      # global alignment.
      ('ATCATCA', 'ATCCA', None, None),
      ('ATCAAATTTCA', 'ATAAATTTCTTA', None, None))
  def test_set_targets(self, ref_seq, target_seq, expected_pos,
                       expected_cigar_op):
    """Test Aligner.align_targets()."""
    align_reads = self.make_test_aligner(ref_seq)
    align_reads.set_targets([target_seq])

    if expected_pos is not None:
      self.assertEqual(expected_pos,
                       [var.pos for var in align_reads.targets[0].gaps])
      self.assertEqual(expected_cigar_op,
                       [var.cigar_op for var in align_reads.targets[0].gaps])
    else:
      self.assertEqual([], align_reads.targets)

  @parameterized.parameters((2), (3))
  def test_ssw_alignment(self, target_offset):
    """Test Aligner._ssw_alignment()."""
    align_reads = self.make_test_aligner()
    read_seq = 'AaGAt'
    target_seq = 'TgAAGATCAGA'
    pw_aligner = align_reads._make_pairwise_aligner(read_seq)

    start_offset, alignment = align_reads._ssw_alignment(
        pw_aligner, read_seq, target_seq, target_offset)
    self.assertEqual(2, alignment.target_begin + start_offset)
    self.assertEqual(5, alignment.optimal_alignment_score)
    self.assertEqual('5M', alignment.cigar)

  @parameterized.parameters(
      ('AaGAt', 'AATA', None, None, 'Read has no common k-mer with target.'),
      ('AaGAt', 'TTAAGAtA', 2, '5M', 'Read has a perfect match.'),
      ('AAAAAAATAAA', 'AAGAAAAAAAA', 0, '2M1D5M1I3M',
       'Read has one insertion and one deletion.'),
      ('TTCAAAGTC', 'AGTCAAAGTCC', 2, '8M',
       'Read starts with a mismatch which should be clipped.'))
  def test_realign_read(self, read_seq, target_seq, expected_align_start,
                        expected_cigar, comment):
    """Test Aligner.test_align_read_to_target()."""
    read = aligner.Read(
        test_utils.make_read(
            read_seq,
            chrom='ref',
            start=0,
            cigar=[(len(read_seq), 'M')],
            quals=[64] * len(read_seq),
            name='read'))
    align_reads = self.make_test_aligner(ref_seq=target_seq)
    align_reads.set_targets([target_seq])

    align_reads.realign_read(read)

    if expected_align_start:
      self.assertEqual(align_reads.targets[0], read.target, comment)
      self.assertEqual(expected_align_start,
                       read.target_offset + read.alignment.target_begin,
                       comment)
      self.assertEqual(expected_cigar, read.alignment.cigar, comment)
    else:
      self.assertIsNone(read.target, comment)
      self.assertIsNone(read.target_offset, comment)
      self.assertIsNone(read.alignment, comment)

  @parameterized.parameters(
      ('10M', 0, [], 'Cigar with no indels.'),
      ('10M', 2, [], 'Cigar with no indels with starting soft-clipped bases.'),
      ('5M2D5M', 2, ['SingleAlnOp(7, DELETE)', 'SingleAlnOp(7, DELETE)'],
       'Cigar has one deletion and starting soft-clipped bases.'),
      ('5M2D4M1I3M', 0, [
          'SingleAlnOp(5, DELETE)', 'SingleAlnOp(5, DELETE)',
          'SingleAlnOp(9, INSERT)'
      ], 'Read sequence has one deletion and one insertion.'),
  )
  def test_cigar_str_to_gaps(self, cigar, query_offset, expected_variants_repr,
                             comment):
    """Test Aligner._cigar_str_to_gaps()."""
    align_reads = self.make_test_aligner()

    self.assertEqual(expected_variants_repr, [
        str(var) for var in align_reads._cigar_str_to_gaps(cigar, query_offset)
    ], comment)

  @parameterized.parameters(
      ('chr1', 10, 100, 'chr1', 20, 20, [(20, 'M')], ''),
      ('chr1', 10, 100, 'chr2', 20, 20, [(20, 'M')],
       'readalignment validation: read reference name is inconsistent with '
       'reference information.'),
      ('chr1', 10, 100, 'chr1', 90, 20, [(20, 'M')],
       'readalignment validation: read end position is out of reference '
       'genomic range.'),
      ('chr1', 10, 100, 'chr1', 75, 20, [(10, 'M'), (10, 'D'), (10, 'M')],
       'readalignment validation: read end position is out of reference '
       'genomic range.'),
      ('chr1', 10, 100, 'chr1', 70, 20, [(10, 'M'), (10, 'D'), (10, 'M')], ''),
      ('chr1', 10, 100, 'chr1', 80, 20, [(5, 'M'), (10, 'I'), (5, 'M')], ''),
      ('chr1', 10, 100, 'chr1', 20, 10, [(5, 'M'), (10, 'I'), (5, 'M')],
       'readalignment validation: cigar is inconsistent with the read length.'),
      ('chr1', 10, 100, 'chr1', 20, 20, [(10, 'H'), (5, 'M'), (10, 'I'),
                                         (5, 'S')], ''),
      ('chr1', 10, 100, 'chr1', 20, 10, [(10, 'H'), (5, 'M'), (10, 'D'),
                                         (5, 'S')], ''),
  )
  def test_sanity_check_readalignment(self, ref_name, ref_start, ref_end,
                                      read_chrom, read_start, read_len,
                                      read_cigar, exception_msg):
    """Test Aligner.sanity_check_readalignment()."""
    region = ranges.make_range(ref_name, ref_start, ref_end)
    ref_seq = 'A' * (ref_end - ref_start)
    align_reads = self.make_test_aligner(ref_seq, region)
    read = test_utils.make_read(
        'A' * read_len,
        chrom=read_chrom,
        start=read_start,
        cigar=read_cigar,
        quals=[64] * read_len,
        name='read')
    if exception_msg:
      with self.assertRaisesRegexp(ValueError, exception_msg):
        align_reads.sanity_check_readalignment(read)
    else:
      align_reads.sanity_check_readalignment(read)

  @parameterized.parameters(
      ('TGCATGG', 23, [(7, 'M')], 'Read is a prefect match to the reference.'),
      ('TGCAAGG', 23, [(7, 'M')],
       'Read has one mismatch w.r.t. the reference.'), ('TAAGCAGGG', 23, [
           (1, 'M'), (2, 'I'), (3, 'M'), (1, 'D'), (3, 'M')
       ], 'Read is a perfect match to the 2nd target.'),
      ('CAGGGGG', 25, [(2, 'M'), (1, 'D'),
                       (5, 'M')], 'Read is a perfect match to the 2nd target.'),
      ('TAAGCAGGAGG', 23, [(1, 'M'), (2, 'I'), (3, 'M'), (1, 'D'), (5, 'M')],
       'Read has one mismatch w.r.t. the 2nd target.'), ('AATAAAGCAGGG', 21, [
           (3, 'M'), (3, 'I'), (3, 'M'), (1, 'D'), (3, 'M')
       ], 'Read has one insertion w.r.t. the 2nd target.'), ('AATAGCAGGG', 21, [
           (3, 'M'), (1, 'I'), (3, 'M'), (1, 'D'), (3, 'M')
       ], 'Read has one deletion w.r.t. the 2nd target.'),
      ('AATAAAGCGGGGGA', 21, [(3, 'M'), (3, 'I'), (2, 'M'), (2, 'D'), (6, 'M')],
       'Read has one insertion and one deletion w.r.t. the 2nd target.'),
      ('GCAAGGGGGA', 24, [(10, 'M')],
       'Read insertion overlaps with the deletion in the 2nd target.'),
      ('TTTAAGCAGGGGGC', 23, [(2, 'S'), (1, 'M'), (2, 'I'), (3, 'M'), (1, 'D'),
                              (5, 'M'), (1, 'S')],
       'Read has clipped bases w.r.t. the 2nd target.'), ('AAGCAGGGGGC', 24, [
           (2, 'S'), (3, 'M'), (1, 'D'), (5, 'M'), (1, 'S')
       ], 'Read starts in an insertion within the 2nd target.'),
      ('AAAGCAGGGGGC', 24, [(3, 'S'), (3, 'M'), (1, 'D'), (5, 'M'), (1, 'S')],
       'Read starts in an insertion within the 2nd target, followed by an '
       'insertion in read to target alignment.'), ('GGGGG', 28, [
           (5, 'M')
       ], 'Read starts after a deletion within the 2nd target.'))
  def test_align_reads_simple(self, read_seq, expected_align_pos,
                              expected_cigar, comment):
    """Test Aligner.align_reads(). Simple tests.

    Targets consist of
      - original reference sequence.
      - a sequence with 'AA' insertion at position 14 and
      -                 'T' deletion at position 19.

    Args:
      read_seq: str, read sequence.
      expected_align_pos: int, expected aligned position
      expected_cigar: [(int, str)], expected cigar information.
      comment: str, test comment.
    """
    ref_seq = 'AAAAAAAAAAAAATGCATGGGGGATTTTTTTTTTT'
    region = ranges.make_range('ref', 10, 10 + len(ref_seq))
    align_reads = self.make_test_aligner(ref_seq, region)
    # redacted
    # implemented. For local alignment, it ensures that there are enough exact
    # matches between the reference and target for end-to-end alignment.
    targets = [ref_seq, 'AAAAAAAAAAAAATAAGCAGGGGGATTTTTTTTTTT']
    read = test_utils.make_read(
        read_seq,
        chrom='ref',
        start=0,
        cigar=[(len(read_seq), 'M')],
        quals=[64] * len(read_seq),
        name='read')
    aligned_reads = align_reads.align_reads(targets, [read])
    self.assertEqual(expected_align_pos,
                     aligned_reads[0].alignment.position.position, comment)
    self.assertEqual(
        _cigar.to_cigar_units(expected_cigar),
        list(aligned_reads[0].alignment.cigar), comment)

    read = test_utils.make_read(
        read_seq,
        chrom='ref',
        start=0,
        cigar=[(2, 'H'), (len(read_seq), 'M'), (1, 'H')],
        quals=[64] * len(read_seq),
        name='read')
    aligned_reads = align_reads.align_reads(targets, [read])
    expected_cigar_w_hard_clip = [(2, 'H')] + expected_cigar + [(1, 'H')]
    self.assertEqual(
        _cigar.to_cigar_units(expected_cigar_w_hard_clip),
        list(aligned_reads[0].alignment.cigar), comment)

  def test_align_read_with_whole_clippd_seq(self):
    """Test Aligner.align_reads() when the whole read sequence is clipped."""
    ref_seq = ('TTTGTTTGTTTGTGTTTGTGTTTTTGTTTGTTTGTGTTTGTGTTTGTTTGTGGTTTGTGT'
               'GTTTGTGTTTGTGTTGGTTTG')
    ref_len = len(ref_seq)
    align_reads = self.make_test_aligner(ref_seq)
    target_ins = 'AAAAAGTGGGGGGGAAGTGGGGAAAAA'
    targets = [
        ref_seq,
        ref_seq[:int(ref_len / 2)] + target_ins + ref_seq[int(ref_len / 2):]
    ]
    read_seq = 'CCC' + target_ins + 'CCC'
    read = test_utils.make_read(
        read_seq,
        chrom='ref',
        start=10,
        cigar=[(len(read_seq), 'M')],
        quals=[64] * len(read_seq),
        name='read')
    aligned_reads = align_reads.align_reads(targets, [read])
    self.assertEqual(read, aligned_reads[0],
                     'Read should have its original alignment.')

  def test_no_bad_soft_clipping(self):
    self.skipTest('Enable when b/63143285 global alignment is fixed')
    common = 'CTA'
    read_seq = common + 'GA'
    ref_seq = 'N' + common + 'CA' + 'N'
    alt_seq = 'A' + ref_seq
    targets = [ref_seq, alt_seq]

    region = ranges.make_range('ref', 0, len(ref_seq))
    align_reads = self.make_test_aligner(ref_seq, region)

    read = test_utils.make_read(
        read_seq,
        chrom='ref',
        start=0,
        cigar=[(len(read_seq), 'M')],
        quals=[35] * len(read_seq),
        name='read')
    realigned = align_reads.align_reads(targets, [read])[0]

    # redacted
    # 5M as we'd expect for this read:
    # read_seq: -CTAGA-
    # ref_seq : NCGTCAN
    # But the current algorithm produces a local alignment of the read against
    # the haplotypes, and the G <=> C mismatch causes the local aligner to
    # simply skip those bases instead of incurring the mismatch penalty for it,
    # resulting in a 3M2S read (GA clipped off) instead of the better 5M result.
    self.assertEqual([_cigar.to_cigar_unit(len(read_seq), 'M')],
                     list(realigned.alignment.cigar))


class LibSSWAlignmentFacadeTest(parameterized.TestCase):
  """Tests for special logic in the wrapper class for libssw alignments."""

  @parameterized.parameters(('5M', '5M'), ('2I2D2M', '2I2D2M'),
                            ('2X1=2X', '5M'), ('2S5M2S', '5M'))
  def test_cigar_simplification(self, cigar, expected_simplified_cigar):
    self.assertEqual(
        expected_simplified_cigar,
        aligner.LibSSWAlignmentFacade._simplify_cigar_string(cigar))


if __name__ == '__main__':
  absltest.main()
