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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import textwrap



from absl.testing import absltest

from deepvariant import test_utils
from deepvariant.core import genomics_io
from deepvariant.core import ranges
from deepvariant.protos import realigner_pb2
from deepvariant.realigner.python import debruijn_graph


def setUpModule():
  test_utils.init()


class DeBruijnGraphWrapTest(absltest.TestCase):
  """Basic tests for the wrapped DeBruijnGraph class."""

  def dbg_options(self):
    return realigner_pb2.RealignerOptions.DeBruijnGraphOptions(
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
    test_options.max_k = k + 1
    test_options.step_k = 1
    return test_options

  def assertGraphEqual(self, graphviz_string, dbg):
    """Assert that the DeBruijn has the given graphviz representation.

    Args:
      graphviz_string: the graphviz representation, potentially including common
        leading whitespace.
      dbg: the DeBruijn graph object.
    """
    self.assertEqual(textwrap.dedent(graphviz_string), dbg.graphviz())

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

    self.assertGraphEqual("""\
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
    self.assertGraphEqual("""\
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

    self.assertGraphEqual("""\
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

  def test_filtering_by_base(self):
    """Test that we filter out edges containing non-canonical bases."""
    ref_str = 'GATTACA'
    read_str = 'GATNTACA'
    read = test_utils.make_read(
        read_str,
        chrom='chr20',
        start=1,
        cigar=[(len(read_str), 'M')],
        quals=[30] * len(read_str),
        name='read')

    # Use two reads so read path doesn't get pruned.
    dbg = debruijn_graph.build(ref_str, [read, read],
                               self.single_k_dbg_options(2))

    self.assertGraphEqual("""\
        digraph G {
        0[label=GA];
        1[label=AT];
        2[label=TT];
        3[label=TA];
        4[label=AC];
        5[label=CA];
        0->1 [label=3 color=red];
        1->2 [label=1 color=red];
        2->3 [label=1 color=red];
        3->4 [label=3 color=red];
        4->5 [label=3 color=red];
        }
        """, dbg)

  def test_filtering_by_qual(self):
    """Test that we filter out edges containing low-quality basecalls."""
    ref_str = 'GATTACA'
    read_str = 'GATGTACA'
    read = test_utils.make_read(
        read_str,
        chrom='chr20',
        start=1,
        cigar=[(len(read_str), 'M')],
        quals=[30, 30, 30, 1, 30, 30, 30, 30],
        name='read')

    # Use two reads so read path doesn't get pruned.
    dbg = debruijn_graph.build(ref_str, [read, read],
                               self.single_k_dbg_options(2))

    self.assertGraphEqual("""\
        digraph G {
        0[label=GA];
        1[label=AT];
        2[label=TT];
        3[label=TA];
        4[label=AC];
        5[label=CA];
        0->1 [label=3 color=red];
        1->2 [label=1 color=red];
        2->3 [label=1 color=red];
        3->4 [label=3 color=red];
        4->5 [label=3 color=red];
        }
        """, dbg)

  def test_straightforward_region(self):
    ref_reader = genomics_io.make_ref_reader(test_utils.CHR20_FASTA)
    bam_reader = genomics_io.make_sam_reader(test_utils.CHR20_BAM)
    region = ranges.parse_literal('chr20:10,000,000-10,000,100')
    ref_seq = ref_reader.bases(region)

    all_reads = list(bam_reader.query(region))
    dbg30 = debruijn_graph.build(ref_seq, all_reads,
                                 self.single_k_dbg_options(30))
    self.assertIsNotNone(dbg30)
    self.assertEqual([ref_seq], dbg30.candidate_haplotypes())

  def test_complex_region(self):
    # There is a heterozygous 9 bp deletion of tandem TGA repeat.
    # "chr20:10,095,379-10,095,500"
    ref_reader = genomics_io.make_ref_reader(test_utils.CHR20_FASTA)
    bam_reader = genomics_io.make_sam_reader(test_utils.CHR20_BAM)
    region = ranges.parse_literal('chr20:10,095,379-10,095,500')
    ref_seq = ref_reader.bases(region)
    reads = list(bam_reader.query(region))
    dbg = debruijn_graph.build(ref_seq, reads, self.dbg_options())
    self.assertIsNotNone(dbg)
    self.assertEqual(44, dbg.kmer_size)
    self.assertEqual(2, len(dbg.candidate_haplotypes()))
    self.assertIn(ref_seq, dbg.candidate_haplotypes())

  def test_k_exceeds_read_length(self):
    """This is a regression test for b/64564513."""
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
    """This is a regression test for b/64564513."""
    # We don't allow a k >= ref length.  This crashed prior to the bugfix.
    ref_str = 'GATTACA'
    dbg = debruijn_graph.build(ref_str, [], self.single_k_dbg_options(7))
    self.assertIsNone(dbg)
    dbg = debruijn_graph.build(ref_str, [], self.single_k_dbg_options(8))
    self.assertIsNone(dbg)


if __name__ == '__main__':
  absltest.main()
