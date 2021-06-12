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



from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from deepvariant import testdata
from deepvariant.protos import realigner_pb2
from deepvariant.realigner.python import debruijn_graph


def setUpModule():
  testdata.init()


class DeBruijnGraphWrapTest(parameterized.TestCase):
  """Basic tests for the wrapped DeBruijnGraph class."""

  def dbg_options(self):
    return realigner_pb2.DeBruijnGraphOptions(
        min_k=12,
        max_k=50,
        step_k=2,
        min_mapq=20,
        min_base_quality=20,
        min_edge_weight=2,
        max_num_paths=10)

  def single_k_dbg_options(self, k):
    """Get a DeBruijnGraphOptions allowing us to try a single kmer size."""
    test_options = self.dbg_options()
    test_options.min_k = k
    test_options.max_k = k
    test_options.step_k = 1
    return test_options

  def assertGraphEqual(self, graphviz_string, dbg):
    """Assert that the DeBruijn has the given graphviz representation.

    Args:
      graphviz_string: the graphviz representation, potentially including common
        leading whitespace.
      dbg: the DeBruijn graph object.
    """
    # Remove all whitespace before comparison to avoid failing over trivial
    # indentation / newline differences.
    self.assertEqual(''.join(graphviz_string.split()),
                     ''.join(dbg.graphviz().split()))

  def test_basics(self):
    """Basic example."""
    ref_str = 'GATTACA'
    read_str = 'GATGACA'
    read = test_utils.make_read(
        read_str,
        chrom='chr20',
        start=1,
        cigar=[(len(read_str), 'M')],
        quals=[30] * len(read_str),
        name='read')

    self.assertEqual(self.single_k_dbg_options(3).min_k, 3)
    # Use two reads so read path doesn't get pruned.
    dbg = debruijn_graph.build(ref_str, [read, read],
                               self.single_k_dbg_options(3))

    self.assertItemsEqual([ref_str, read_str], dbg.candidate_haplotypes())

    self.assertGraphEqual(
        """\
          digraph G {
          0[label=GAT];
          1[label=ATT];
          2[label=TTA];
          3[label=TAC];
          4[label=ACA];
          5[label=ATG];
          6[label=TGA];
          7[label=GAC];
          0->1 [label=1 color=red];
          1->2 [label=1 color=red];
          2->3 [label=1 color=red];
          3->4 [label=1 color=red];
          0->5 [label=2];
          5->6 [label=2];
          6->7 [label=2];
          7->4 [label=2];
          }
          """, dbg)

  def test_pruning_1(self):
    """Test that pruning removes a path traced by only one read."""
    ref_str = 'GATTACA'
    read_str = 'GATGACA'
    read = test_utils.make_read(
        read_str,
        chrom='chr20',
        start=1,
        cigar=[(len(read_str), 'M')],
        quals=[30] * len(read_str),
        name='read')
    dbg = debruijn_graph.build(ref_str, [read], self.single_k_dbg_options(3))
    self.assertGraphEqual(
        """\
        digraph G {
        0[label=GAT];
        1[label=ATT];
        2[label=TTA];
        3[label=TAC];
        4[label=ACA];
        0->1 [label=1 color=red];
        1->2 [label=1 color=red];
        2->3 [label=1 color=red];
        3->4 [label=1 color=red];
        }
        """, dbg)

  def test_pruning_2(self):
    """Test that pruning removes edges not between source and sink."""
    ref_str = 'GATTACA'
    read_str = 'CCGATGACACC'
    read = test_utils.make_read(
        read_str,
        chrom='chr20',
        start=1,
        cigar=[(len(read_str), 'M')],
        quals=[30] * len(read_str),
        name='read')
    # Use two reads so read path doesn't get pruned.
    dbg = debruijn_graph.build(ref_str, [read, read],
                               self.single_k_dbg_options(3))

    self.assertGraphEqual(
        """\
        digraph G {
        0[label=GAT];
        1[label=ATT];
        2[label=TTA];
        3[label=TAC];
        4[label=ACA];
        5[label=ATG];
        6[label=TGA];
        7[label=GAC];
        0->1 [label=1 color=red];
        1->2 [label=1 color=red];
        2->3 [label=1 color=red];
        3->4 [label=1 color=red];
        0->5 [label=2];
        5->6 [label=2];
        6->7 [label=2];
        7->4 [label=2];
        }
        """, dbg)

  @parameterized.parameters(
      # No bad positions => all edges get +1 to counts.
      dict(bad_position=None, dropped_edges={}),

      # Ref and read are   GATTACA
      # Bad position: 0 => *
      # Breaks kmers:      GA->AT
      dict(bad_position=0, dropped_edges={'GA->AT'}),

      # Ref and read are   GATTACA
      # Bad position: 1 =>  *
      # Breaks kmers:      GA->AT, AT->TT
      dict(bad_position=1, dropped_edges={'GA->AT', 'AT->TT'}),

      # Ref and read are   GATTACA
      # Bad position: 2 =>   *
      # Breaks kmers:      GA->AT, AT->TT, TT->TA
      dict(bad_position=2, dropped_edges={'GA->AT', 'AT->TT', 'TT->TA'}),

      # Ref and read are   GATTACA
      # Bad position: 3 =>    *
      # Breaks kmers:      AT->TT, TT->TA, TA->AC
      dict(bad_position=3, dropped_edges={'AT->TT', 'TT->TA', 'TA->AC'}),

      # Ref and read are   GATTACA
      # Bad position: 4 =>     *
      # Breaks kmers:      TT->TA, TA->AC, AC->CA
      dict(bad_position=4, dropped_edges={'TT->TA', 'TA->AC', 'AC->CA'}),

      # Ref and read are   GATTACA
      # Bad position: 5 =>      *
      # Breaks kmers:      TA->AC, AC->CA
      dict(bad_position=5, dropped_edges={'TA->AC', 'AC->CA'}),

      # Ref and read are   GATTACA
      # Bad position: 6 =>       *
      # Breaks kmers:      AC->CA
      dict(bad_position=6, dropped_edges={'AC->CA'}),
  )
  def test_adding_edges_with_bad_positions(self, bad_position, dropped_edges):
    """Test that we filter out edges containing low-quality basecalls."""
    ref_str = 'GATTACA'
    read_str = 'GATTACA'

    kmer_indices = {
        'GA': 0,
        'AT': 1,
        'TT': 2,
        'TA': 3,
        'AC': 4,
        'CA': 5,
    }

    def kmer_to_index_edge(kmer_edge):
      k1, k2 = kmer_edge.split('->')
      return '{}->{}'.format(kmer_indices[k1], kmer_indices[k2])

    dropped_edges = {kmer_to_index_edge(edge) for edge in dropped_edges}

    for bad_type in ['qual', 'base']:
      bases = list(read_str)
      quals = [30] * len(bases)
      cigar = [(len(bases), 'M')]
      if bad_position is not None:
        if bad_type == 'qual':
          quals[bad_position] = 1
        elif bad_type == 'base':
          bases[bad_position] = 'N'
        else:
          raise ValueError('Unexpected base type')

      read = test_utils.make_read(
          ''.join(bases), start=0, cigar=cigar, quals=quals)

      # Use two reads so read path doesn't get pruned.
      dbg = debruijn_graph.build(ref_str, [read, read],
                                 self.single_k_dbg_options(2))

      expected_edges = '\n'.join(
          '{} [label={} color=red];'.format(edge, 1 if edge in
                                            dropped_edges else 3)
          for edge in ['0->1', '1->2', '2->3', '3->4', '4->5'])

      self.assertGraphEqual(
          """\
            digraph G {
            0[label=GA];
            1[label=AT];
            2[label=TT];
            3[label=TA];
            4[label=AC];
            5[label=CA];
            %s
            }
            """ % expected_edges, dbg)

  def test_straightforward_region(self):
    ref_reader = fasta.IndexedFastaReader(testdata.CHR20_FASTA)
    bam_reader = sam.SamReader(testdata.CHR20_BAM)
    region = ranges.parse_literal('chr20:10,000,000-10,000,100')
    ref_seq = ref_reader.query(region)

    all_reads = list(bam_reader.query(region))
    dbg30 = debruijn_graph.build(ref_seq, all_reads,
                                 self.single_k_dbg_options(30))
    self.assertIsNotNone(dbg30)
    self.assertEqual([ref_seq], dbg30.candidate_haplotypes())

  def test_complex_region(self):
    # There is a heterozygous 9 bp deletion of tandem TGA repeat.
    # "chr20:10,095,379-10,095,500"
    ref_reader = fasta.IndexedFastaReader(testdata.CHR20_FASTA)
    bam_reader = sam.SamReader(testdata.CHR20_BAM)
    region = ranges.parse_literal('chr20:10,095,379-10,095,500')
    ref_seq = ref_reader.query(region)
    reads = list(bam_reader.query(region))
    dbg = debruijn_graph.build(ref_seq, reads, self.dbg_options())
    self.assertIsNotNone(dbg)
    self.assertEqual(44, dbg.kmer_size)
    self.assertLen(dbg.candidate_haplotypes(), 2)
    self.assertIn(ref_seq, dbg.candidate_haplotypes())

  def test_k_exceeds_read_length(self):
    """This is a regression test for internal."""
    # If k > read length, no edges will go into the graph from this read.
    # This crashed prior to the bugfix.
    ref_str = 'GATTACATG'
    read_str = 'GATGACA'
    read = test_utils.make_read(
        read_str,
        chrom='chr20',
        start=1,
        cigar=[(len(read_str), 'M')],
        quals=[30] * len(read_str),
        name='read')
    dbg = debruijn_graph.build(ref_str, [read, read],
                               self.single_k_dbg_options(8))
    self.assertIsNotNone(dbg)

  def test_k_exceeds_ref_length(self):
    """This is a regression test for internal."""
    # We don't allow a k >= ref length.  This crashed prior to the bugfix.
    ref_str = 'GATTACA'
    dbg = debruijn_graph.build(ref_str, [], self.single_k_dbg_options(7))
    self.assertIsNone(dbg)
    dbg = debruijn_graph.build(ref_str, [], self.single_k_dbg_options(8))
    self.assertIsNone(dbg)

  @parameterized.parameters(
      dict(ref='ACGTACGT', smallest_good_k=5),
      dict(ref='ACGTAAACGT', smallest_good_k=5),
      dict(ref='ACGTAAACGTAAA', smallest_good_k=8),
      dict(ref='AAACGTAAACGT', smallest_good_k=7),
      dict(ref='AAACGTAAACGTAAA', smallest_good_k=10),
      # Actual example where the cycle detector failed because the cycle only
      # occurs with the last kmer in the reference.
      dict(
          ref=(
              'TGGTAAGTTTATAAGGTTATAAGCTGAGAGGTTTTGCTGATCTTGGCTGAGCTCAGCTGGGCAGGTC'
              'TTCCGGTCTTGGCTGGGGTTCACTGACACACAAGCAGCTGACAGTTGGCTGATCTAGGATGGCCTCA'
              'GCTGGG'),
          smallest_good_k=11),
  )
  def test_ref_cycle_detector(self, ref, smallest_good_k):
    min_k = max(smallest_good_k - 5, 1)
    max_k = min(smallest_good_k + 5, len(ref))
    for k in range(min_k, max_k):
      # The build fails, returning None, with a k < smallest_good_k. If
      # k >= smallest_good_k, then we expect a real non-None instance.
      result = debruijn_graph.build(ref, [], self.single_k_dbg_options(k))
      if k < smallest_good_k:
        self.assertIsNone(result, 'Cycle not detected for k={}'.format(k))
      else:
        self.assertIsNotNone(result, 'False cycle detected for k={}'.format(k))


if __name__ == '__main__':
  absltest.main()
