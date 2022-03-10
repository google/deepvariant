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
"""Tests for deepvariant .realigner.realigner."""

import csv
import itertools
import os


from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
import numpy as np
import six
import tensorflow as tf

from deepvariant import testdata
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import realigner
from deepvariant.realigner import utils
from third_party.nucleus.io import fasta
from third_party.nucleus.io import sam
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import cigar as cigar_utils
from third_party.nucleus.util import ranges


FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


def _get_reads(region):
  with sam.SamReader(testdata.CHR20_BAM) as in_sam_reader:
    return list(in_sam_reader.query(region))


def _get_reads_and_header(region):
  with sam.SamReader(testdata.CHR20_BAM) as in_sam_reader:
    return list(in_sam_reader.query(region)), in_sam_reader.header


def _test_assembled_region(region_str, haplotypes=None):
  return realigner.AssemblyRegion(
      realigner_pb2.CandidateHaplotypes(
          span=ranges.parse_literal(region_str), haplotypes=haplotypes or []))


class ReadAssignmentTests(parameterized.TestCase):

  def setUp(self):
    reads = [
        test_utils.make_read('ACG', start=1, cigar='3M', name='read1'),
        test_utils.make_read('ACG', start=6, cigar='3M', name='read2'),
        test_utils.make_read('ACG', start=9, cigar='3M', name='read3'),
        test_utils.make_read('ACG', start=28, cigar='3M', name='read4'),
        test_utils.make_read('A' * 10, start=3, cigar='10M', name='read5'),
    ]
    self.reads = {read.fragment_name: read for read in reads}
    self.regions = {
        'r1': _test_assembled_region('chr1:1-5'),
        'r2': _test_assembled_region('chr1:10-15'),
        'r3': _test_assembled_region('chr1:20-30'),
    }
    self.assembled_regions = [self.regions[r] for r in sorted(self.regions)]

  def get_reads_by_name(self, names):
    return [self.reads[name] for name in names]

  def test_construction(self):
    aregion = _test_assembled_region('chr1:1-5', haplotypes=['A', 'C'])
    self.assertEqual(aregion.region, ranges.parse_literal('chr1:1-5'))
    self.assertEqual(aregion.haplotypes, ['A', 'C'])
    self.assertEqual(aregion.reads, [])

  def test_adding_reads(self):
    aregion = _test_assembled_region('chr1:3-15')

    # We haven't added any reads, so reads is empty and the span is None.
    self.assertEqual(aregion.reads, [])
    self.assertIsNone(aregion.read_span)

    # Add read2, giving us a real read span and a read in our region's reads.
    read_to_add = self.get_reads_by_name(['read2'])[0]
    expected_reads = [read_to_add]
    aregion.add_read(read_to_add)
    self.assertEqual(aregion.reads, expected_reads)
    self.assertEqual(aregion.read_span, ranges.parse_literal('chr1:7-9'))

    # Add read1, increasing the span on the left.
    read_to_add = self.get_reads_by_name(['read1'])[0]
    expected_reads += [read_to_add]
    aregion.add_read(read_to_add)
    self.assertEqual(aregion.reads, expected_reads)
    self.assertEqual(aregion.read_span, ranges.parse_literal('chr1:2-9'))

    # Finally, add in all of the reads.
    reads_to_add = self.get_reads_by_name(['read3', 'read4', 'read5'])
    expected_reads += reads_to_add
    for read in reads_to_add:
      aregion.add_read(read)
    self.assertEqual(aregion.reads, expected_reads)
    self.assertEqual(aregion.read_span, ranges.parse_literal('chr1:2-31'))

  @parameterized.parameters(
      # Single read tests.
      # read1 overlaps r1.
      dict(read_name='read1', expected_region='r1'),
      # read2 falls between r1 and r2, should be unassigned.
      dict(read_name='read2', expected_region=None),
      # read3 starts before r2 but overlaps it.
      dict(read_name='read3', expected_region='r2'),
      # read4 starts in r3 but extends beyond it.
      dict(read_name='read4', expected_region='r3'),
      # read5 overlaps r1 and r2 but is more in r2 then r1.
      dict(read_name='read5', expected_region='r2'),
  )
  def test_assign_reads_to_assembled_regions_single_read(
      self, read_name, expected_region):
    assignment = {expected_region: [read_name]} if expected_region else {}
    self.assertReadsGoToCorrectRegions(
        reads=self.get_reads_by_name([read_name]),
        expected_assignments=assignment)

  @parameterized.parameters(
      # Let's make sure adding all of the reads together results in the correct
      # assignment across all regions.
      dict(
          read_names=names,
          expected_assignments={
              'r1': ['read1'],
              'r2': ['read3', 'read5'],
              'r3': ['read4'],
          }) for names in itertools.permutations(
              ['read1', 'read2', 'read3', 'read4', 'read5']))
  def test_assign_reads_to_assembled_regions_multiple_reads(
      self, read_names, expected_assignments):
    self.assertReadsGoToCorrectRegions(
        self.get_reads_by_name(read_names), expected_assignments)

  def assertReadsGoToCorrectRegions(self, reads, expected_assignments):
    unassigned = realigner.assign_reads_to_assembled_regions(
        self.assembled_regions, reads)

    # Every read should be in the assembled regions or unassigned.
    six.assertCountEqual(
        self,
        [r for ar in self.assembled_regions for r in ar.reads] + unassigned,
        reads)

    # Go through each region and make sure the reads that are supposed to
    # appear in each region do in appear there.
    for region_name, region in self.regions.items():
      expected_reads = self.get_reads_by_name(
          expected_assignments.get(region_name, []))
      six.assertCountEqual(self, region.reads, expected_reads)


class RealignerTest(parameterized.TestCase):

  def setUp(self):
    self.ref_reader = fasta.IndexedFastaReader(testdata.CHR20_FASTA)
    # TODO: Update the tests to reflect the new default (False).
    FLAGS.ws_use_window_selector_model = True
    self.config = realigner.realigner_config(FLAGS)
    self.reads_realigner = realigner.Realigner(self.config, self.ref_reader)

  @parameterized.parameters(
      # Arguments passed by ws_{min,max}_supporting_reads.
      dict(
          model=None, min_supporting=2, max_supporting=300, use_ws_model=False),
      # No flags passed for the window_selection.
      dict(
          model=None, min_supporting=-1, max_supporting=-1, use_ws_model=False),
      # VariantReadsThresholdModel.
      dict(
          model='VARIANT_READS_THRESHOLD',
          min_supporting=-1,
          max_supporting=-1,
          use_ws_model=True),
      # AlleleCountLinearModel.
      dict(
          model='ALLELE_COUNT_LINEAR',
          min_supporting=-1,
          max_supporting=-1,
          use_ws_model=True),
      # Use the default AlleleCountLinearModel.
      dict(model=None, min_supporting=-1, max_supporting=-1, use_ws_model=True))
  @flagsaver.flagsaver
  def test_window_selector_model_flags(self, model, min_supporting,
                                       max_supporting, use_ws_model):
    # This indirection is needed because the symbols in testdata are not set
    # when the @parameterized decorator is called.
    symbol_to_testdata = {
        None: None,
        'VARIANT_READS_THRESHOLD': testdata.WS_VARIANT_READS_THRESHOLD_MODEL,
        'ALLELE_COUNT_LINEAR': testdata.WS_ALLELE_COUNT_LINEAR_MODEL
    }
    FLAGS.ws_max_num_supporting_reads = max_supporting
    FLAGS.ws_min_num_supporting_reads = min_supporting
    FLAGS.ws_window_selector_model = symbol_to_testdata[model]
    FLAGS.ws_use_window_selector_model = use_ws_model
    # We only make sure that reading the model does not crash or raise
    # exceptions.
    _ = realigner.realigner_config(FLAGS)

  @flagsaver.flagsaver
  def test_window_selector_model_flags_failures(self):
    with six.assertRaisesRegex(
        self, ValueError, 'ws_min_supporting_reads should be smaller than ws_'
        'max_supporting_reads.'):
      FLAGS.ws_max_num_supporting_reads = 1
      FLAGS.ws_min_num_supporting_reads = 2
      FLAGS.ws_window_selector_model = None
      FLAGS.ws_use_window_selector_model = False
      _ = realigner.realigner_config(FLAGS)

    with six.assertRaisesRegex(
        self, ValueError, 'Cannot specify a ws_window_selector_model '
        'if ws_use_window_selector_model is False.'):
      FLAGS.ws_max_num_supporting_reads = -1
      FLAGS.ws_min_num_supporting_reads = -1
      FLAGS.ws_window_selector_model = testdata.WS_ALLELE_COUNT_LINEAR_MODEL
      FLAGS.ws_use_window_selector_model = False
      _ = realigner.realigner_config(FLAGS)

    with six.assertRaisesRegex(
        self, ValueError, 'Cannot use both ws_min_num_supporting_reads and '
        'ws_use_window_selector_model flags.'):
      FLAGS.ws_max_num_supporting_reads = -1
      FLAGS.ws_min_num_supporting_reads = 1
      FLAGS.ws_window_selector_model = None
      FLAGS.ws_use_window_selector_model = True
      _ = realigner.realigner_config(FLAGS)

    with six.assertRaisesRegex(
        self, ValueError, 'Cannot use both ws_max_num_supporting_reads and '
        'ws_use_window_selector_model flags.'):
      FLAGS.ws_max_num_supporting_reads = 1
      FLAGS.ws_min_num_supporting_reads = -1
      FLAGS.ws_window_selector_model = None
      FLAGS.ws_use_window_selector_model = True
      _ = realigner.realigner_config(FLAGS)

  @parameterized.parameters(
      dict(
          region_literal='chr20:10,095,379-10,095,500',
          expected_window_literal='chr20:10,095,352-10,095,553',
          expected_haplotypes={
              'TAGTGATCTAGTCCTTTTTGTTGTGCAAAAGGAAGTGCTAAAATCAGAATGAGAACCATGGTCA'
              'CCTGACATAGACACAAGTGATGATGATGATGATGATGATGATGATGATGATGATATCCATGTTC'
              'AAGTACTAATTCTGGGCAAGACACTGTTCTAAGTGCTATGAATATATTACCTCATTTAATCATC'
              'T',
              'TAGTGATCTAGTCCTTTTTGTTGTGCAAAAGGAAGTGCTAAAATCAGAATGAGAACCATGGTCA'
              'CCTGACATAGACACAAGTGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATA'
              'TCCATGTTCAAGTACTAATTCTGGGCAAGACACTGTTCTAAGTGCTATGAATATATTACCTCAT'
              'TTAATCATCT'
          },
          comment='There is a heterozygous 9 bp deletion of tandem TGA repeat.'
      ),
      dict(
          region_literal='chr20:10,046,080-10,046,307',
          expected_window_literal='chr20:10,046,096-10,046,267',
          expected_haplotypes={
              'CCCAAAAAAAGAGTTAGGGATGCTGGAAAGGCAGAAAGAAAAGGGAAGGGAAGAGGAAGGGGAA'
              'AAGGAAAGAAAAAAAAGAAAGAAAGAAAGAGAAAGAAAGAGAAAGAGAAAGAAAGAGGAAAGAG'
              'AGAAAGAGAAAGAGAAGGAAAGAGAAAGAAAGAGAAGGAAAGAG',
              'CCCAAAAAAAGAGTTAGGGATGCTGGAAAGGCAGAAAGAAAAGGGAAGGGAAGAGGAAGGGGAA'
              'AAGGAAAGAAAAAAAAGAAAGAAAGAAAGAGAAAGAGAAAGAAAGAGGAAAGAGAGAAAGAGAA'
              'AGAGAAGGAAAGAGAAAGAAAGAGAAGGAAAGAG'
          },
          comment='There is a heterozygous 10 bp deletion.'),
  )
  def test_realigner_example_region(self, region_literal,
                                    expected_window_literal,
                                    expected_haplotypes, comment):
    region = ranges.parse_literal(region_literal)
    reads = _get_reads(region)
    windows_haplotypes, realigned_reads = self.reads_realigner.realign_reads(
        reads, region)

    self.assertEqual(len(reads), len(realigned_reads))
    self.assertEqual(
        ranges.parse_literal(expected_window_literal),
        windows_haplotypes[0].span)
    self.assertEqual(expected_haplotypes, set(windows_haplotypes[0].haplotypes))

  @parameterized.parameters(
      [
          dict(
              region_literal='chr20:10,046,080-10,046,307',
              variant_literal='chr20:10,046,179-10,046,188')
      ],)
  def test_realigner_example_variant(self, region_literal, variant_literal):
    """All overlapping reads should include 10bp deletion at chr20:10046178."""
    region = ranges.parse_literal(region_literal)
    variant = ranges.parse_literal(variant_literal)

    reads = _get_reads(region)
    _, realigned_reads = self.reads_realigner.realign_reads(reads, region)

    for read in realigned_reads:
      has_variant = False
      self.assertTrue(read.HasField('alignment'))
      self.assertEqual(variant.reference_name,
                       read.alignment.position.reference_name)
      ref_pos = read.alignment.position.position
      for cigar in read.alignment.cigar:
        self.assertIn(cigar.operation, utils.CIGAR_OPS)
        if cigar.operation in utils.CIGAR_ALIGN_OPS:
          ref_pos += cigar.operation_length
        elif cigar.operation in utils.CIGAR_DELETE_OPS:
          if (ref_pos == variant.start and
              cigar.operation_length == variant.end - ref_pos):
            has_variant = True
          ref_pos += cigar.operation_length
      if (read.alignment.position.position <= variant.start and
          ref_pos >= variant.end):
        self.assertTrue(has_variant)

  def test_realigner_doesnt_create_invalid_intervals(self):
    """Tests that read sets don't result in a crash in reference_fai.cc."""
    region = ranges.parse_literal('chr20:63,025,320-63,025,520')

    # pylint: disable=g-complex-comprehension
    reads = [
        test_utils.make_read(
            'ACCGT' * 50,
            start=63025520 - 250,
            cigar='250M',
            quals=list(np.tile(range(30, 35), 50))) for _ in range(20)
    ]
    # pylint: enable=g-complex-comprehension
    self.reads_realigner.realign_reads(reads, region)

    # These reads are aligned off the edge of the contig. Note that the
    # reference bases in this interval are all Ns as well.
    # pylint: disable=g-complex-comprehension
    reads = [
        test_utils.make_read(
            'TTATA' * 50,
            start=63025520 - 200,
            cigar='200M50S',
            quals=list(np.tile(range(30, 35), 50))) for _ in range(20)
    ]
    # pylint: enable=g-complex-comprehension
    self.reads_realigner.realign_reads(reads, region)

  @parameterized.parameters(
      dict(enabled=False, emit_reads=False),
      dict(enabled=True, emit_reads=False),
      dict(enabled=True, emit_reads=True),
  )
  def test_realigner_diagnostics(self, enabled, emit_reads):
    # Make sure that by default we aren't emitting any diagnostic outputs.
    dx_dir = test_utils.test_tmpfile('dx_enabled{}_emitreads_{}'.format(
        enabled, emit_reads))
    region_str = 'chr20:10046178-10046188'
    region = ranges.parse_literal(region_str)
    assembled_region_str = 'chr20:10046096-10046267'
    reads, header = _get_reads_and_header(region)
    self.config = realigner.realigner_config(FLAGS)
    self.config.diagnostics.enabled = enabled
    self.config.diagnostics.output_root = dx_dir
    self.config.diagnostics.emit_realigned_reads = emit_reads
    self.reads_realigner = realigner.Realigner(self.config, self.ref_reader,
                                               header)
    _, _ = self.reads_realigner.realign_reads(reads, region)
    self.reads_realigner.diagnostic_logger.close()  # Force close all resources.

    if not enabled:
      # Make sure our diagnostic output isn't emitted.
      self.assertFalse(tf.io.gfile.exists(dx_dir))
    else:
      # Our root directory exists.
      self.assertTrue(tf.io.gfile.isdir(dx_dir))

      # We expect a realigner_metrics.csv in our rootdir with 1 entry in it.
      metrics_file = os.path.join(
          dx_dir, self.reads_realigner.diagnostic_logger.metrics_filename)
      self.assertTrue(tf.io.gfile.exists(metrics_file))
      with tf.io.gfile.GFile(metrics_file) as fin:
        rows = list(csv.DictReader(fin))
        self.assertLen(rows, 1)
        self.assertEqual(
            set(rows[0].keys()), {'window', 'k', 'n_haplotypes', 'time'})
        self.assertEqual(rows[0]['window'], assembled_region_str)
        self.assertEqual(int(rows[0]['k']), 25)
        self.assertTrue(int(rows[0]['n_haplotypes']), 2)
        # Check that our runtime is reasonable (greater than 0, less than 10 s).
        self.assertTrue(0.0 < float(rows[0]['time']) < 10.0)

      # As does the subdirectory for this region.
      region_subdir = os.path.join(dx_dir, assembled_region_str)
      self.assertTrue(tf.io.gfile.isdir(region_subdir))

      # We always have a graph.dot
      self.assertTrue(
          tf.io.gfile.exists(
              os.path.join(
                  region_subdir,
                  self.reads_realigner.diagnostic_logger.graph_filename)))

      reads_file = os.path.join(
          dx_dir, region_str,
          self.reads_realigner.diagnostic_logger.realigned_reads_filename)

      # if emit_reads=False then file should not exist and vice versa.
      self.assertEqual(emit_reads, tf.io.gfile.exists(reads_file))

  @parameterized.parameters(
      dict(
          read_seq='AAGGAAGTGCTAAAATCAGAATGAGAACCATGG'
          'ATCCATGTTCAAGTACTAATTCTGGGC',
          prefix='AGTGATCTAGTCCTTTTTGTTGTGCAAAAGGAAGTGCTAAAATCAGAATGAGAACCATGG',
          suffix='ATCCATGTTCAAGTACTAATTCTGGGCAAGACACTGTTCTAAGTGCTATGAATATATTACC',
          haplotypes=['CATCATCAT', ''],
          expected_cigars=['33M9D27M', '60M']),
      dict(
          read_seq='TTGCCCGGGCATAAGGTGTTTCGGAGAAGCCTAG'
          'TATATATA'
          'CTCCGGTTTTTAAGTAGGGTCGTAGCAG',
          prefix='AACGGGTCTACAAGTCTCTGCGTGTTGCCCGGGCATAAGGTGTTTCGGAGAAGCCTAG',
          suffix='CTCCGGTTTTTAAGTAGGGTCGTAGCAGCAAAGTAAGAGTGGAACGCGTGGGCGACTA',
          haplotypes=['', 'TATATATA'],
          expected_cigars=['34M8I28M', '70M']),
      dict(
          read_seq='AAAAAAAAAAGGGGGGGGGGATTTTTTTTTTTTTCCCCCCCCCCCCCCC',
          prefix='AAAAAAAAAAGGGGGGGGGG',
          suffix='TTTTTTTTTTTTTCCCCCCCCCCCCCCC',
          haplotypes=['A', ''],
          expected_cigars=['49M', '20M1I28M']),
  )
  def test_align_to_haplotype(self, read_seq, prefix, suffix, haplotypes,
                              expected_cigars):
    test_read = test_utils.make_read(read_seq, start=1)
    reads = [test_read]
    # Align to each haplotype in turn.
    for i in range(len(haplotypes)):
      aligned_reads = self.reads_realigner.align_to_haplotype(
          haplotypes[i], haplotypes, prefix, suffix, reads, 'test', 1)
      self.assertEqual(len(reads), len(aligned_reads))
      self.assertEqual(
          cigar_utils.format_cigar_units(aligned_reads[0].alignment.cigar),
          expected_cigars[i])

  @parameterized.parameters(
      dict(alt_allele='CATTACA', ref_buffer_length=70, read_buffer_length=20),
      dict(alt_allele='CATTACA', ref_buffer_length=20, read_buffer_length=20),
      dict(alt_allele='G', ref_buffer_length=70, read_buffer_length=20),
      # At or below read_buffer_length=15 the reads start to come back
      # unaligned, but this depends on the specific ref and alt alleles, so
      # this does not include exhaustive tests for how low these values can go.
  )
  def test_align_to_haplotype_stress_tests(self, alt_allele, ref_buffer_length,
                                           read_buffer_length):
    """Testing what happens when read and reference sequences are shorter."""
    # Start with long prefix and suffix to enable cutting it down as necessary
    whole_prefix = 'AGTGATCTAGTCCTTTTTGTTGTGCAAAAGGAAGTGCTAAAATCAGAATGAGAACCATGGTCACCTGACATAGAC'
    whole_suffix = 'ATCCATGTTCAAGTACTAATTCTGGGCAAGACACTGTTCTAAGTGCTATGAATATATTACCTCATTTAATCATCT'

    ref_prefix = whole_prefix[-ref_buffer_length:]
    ref_suffix = whole_suffix[:ref_buffer_length]

    # Make two haplotypes.
    ref_allele = ''
    haplotypes = [ref_allele, alt_allele]

    # Simulate one read from the reference and one from the alt haplotype.
    read_prefix = ref_prefix[-read_buffer_length:]
    read_suffix = ref_suffix[:read_buffer_length]

    expected_cigars = [
        # Aligning to ref haplotype: Insertion.
        '{}M{}I{}M'.format(len(read_prefix), len(alt_allele), len(read_suffix)),
        # Aligning to alt haplotype: All matching.
        '{}M'.format(len(read_prefix) + len(alt_allele) + len(read_suffix))
    ]

    reads = [
        test_utils.make_read(read_prefix + alt_allele + read_suffix, start=1)
    ]
    # Align to each haplotype in turn.
    for i in range(len(haplotypes)):
      aligned_reads = self.reads_realigner.align_to_haplotype(
          haplotypes[i], haplotypes, ref_prefix, ref_suffix, reads, 'test', 1)
      self.assertEqual(len(reads), len(aligned_reads))
      self.assertEqual(
          cigar_utils.format_cigar_units(aligned_reads[0].alignment.cigar),
          expected_cigars[i])

  def test_align_to_haplotype_empty_reads(self):
    # Empty reads as input should return empty reads as output.
    aligned_reads = self.reads_realigner.align_to_haplotype(
        this_haplotype='G',
        haplotypes=['G', ''],
        prefix='AAA',
        suffix='AAA',
        reads=[],
        contig='test',
        ref_start=1)
    self.assertEqual(aligned_reads, [])

  @parameterized.parameters(
      dict(
          # No change.
          read_seq='AAGGAAGTGCTAAAATCAGAATGAGAACCA',
          cigar='30M',
          expected_cigars=['30M'],
          expected_sequences=['AAGGAAGTGCTAAAATCAGAATGAGAACCA'],
          expected_positions=[1]),
      dict(
          # Basic split.
          read_seq='AAGGAAGTGCTAAAATCAGAATGAGAACCA',
          cigar='15M5000N15M',
          expected_cigars=['15M', '15M'],
          expected_sequences=['AAGGAAGTGCTAAAA', 'TCAGAATGAGAACCA'],
          expected_positions=[1, 5016]),
      dict(
          # Split with 15bp filter.
          read_seq='AAGGAAGTGCTAAAATCAGAATGAGAACCA',
          cigar='10M10N20M',
          expected_cigars=['20M'],
          expected_sequences=['TAAAATCAGAATGAGAACCA'],
          expected_positions=[21]),
      dict(
          # Many small splits filtered out.
          read_seq='AAGGAAGTGCTAAAATCAGAATGAGAACCA',
          cigar='5M5N5M5N5M5N5M5N5M5N5M',
          expected_cigars=[],
          expected_sequences=[],
          expected_positions=[]),
      dict(
          # Large split.
          read_seq='AAGGAAGTGCTAAAATCAGAATGAGAACCA',
          cigar='2M5000N28M',
          expected_cigars=['28M'],
          expected_sequences=['GGAAGTGCTAAAATCAGAATGAGAACCA'],
          expected_positions=[5003]),
      dict(
          # Insertion.
          read_seq='AAGGAAGTGCTAATTTTTAATCAGAATGAGAACCA',
          cigar='15M5I15M',
          expected_cigars=['15M5I15M'],
          expected_sequences=['AAGGAAGTGCTAATTTTTAATCAGAATGAGAACCA'],
          expected_positions=[1]),
      dict(
          # Insertion + Split.
          read_seq='AAGGAAGTGCTAAAAGGGGGTCAGAATGAGAACCA',
          cigar='15M5I50N15M',
          expected_cigars=['15M5I', '15M'],
          expected_sequences=['AAGGAAGTGCTAAAAGGGGG', 'TCAGAATGAGAACCA'],
          expected_positions=[1, 66]),
      dict(
          # Deletion.
          read_seq='AAGGAAGTGCTAATTTTTAATCAGAATGAGAACCA',
          cigar='15M5D15M',
          expected_cigars=['15M5D15M'],
          expected_sequences=['AAGGAAGTGCTAATTTTTAATCAGAATGAGAACCA'],
          expected_positions=[1]),
      dict(
          # Deletion + Split.
          read_seq='AAGGAAGTGCTAATTTCAGAATGAGAACCA',
          cigar='15M5D50N15M',
          expected_cigars=['15M5D', '15M'],
          expected_sequences=['AAGGAAGTGCTAATT', 'TCAGAATGAGAACCA'],
          expected_positions=[1, 71]),
      dict(
          # Sequence Match/Mismatch + Split.
          read_seq='CCCCGGACACTTCTAGTTTGTCGGAGCGAGTC',
          cigar='15=1X1=20N15=',
          expected_cigars=['15=1X1=', '15='],
          expected_sequences=['CCCCGGACACTTCTAGT', 'TTGTCGGAGCGAGTC'],
          expected_positions=[1, 38]),
      dict(
          # Soft Clip + Split.
          read_seq='TGAGCTAGTAGAATTTAGGGAGAAAGATTAATGCG',
          cigar='15S5M50N15M',
          expected_cigars=['15S5M', '15M'],
          expected_sequences=['TGAGCTAGTAGAATTTAGGG', 'AGAAAGATTAATGCG'],
          expected_positions=[1, 56]),
      dict(
          # Hard Clip + Split.
          read_seq='ATCCCGGCCACGTTAATCCCGGCCACGTTA',
          cigar='15H15M50N15M15H',
          expected_cigars=['15H15M', '15M15H'],
          expected_sequences=['ATCCCGGCCACGTTA', 'ATCCCGGCCACGTTA'],
          expected_positions=[1, 66]),
  )
  def test_split_reads(self, read_seq, cigar, expected_cigars,
                       expected_sequences, expected_positions):
    test_read = test_utils.make_read(read_seq, cigar=cigar, start=1)
    reads = realigner.split_reads([test_read])
    for i in range(len(reads)):
      # Check sequences
      self.assertEqual(reads[i].aligned_sequence, expected_sequences[i])
      # Check cigars
      self.assertEqual(
          cigar_utils.format_cigar_units(reads[i].alignment.cigar),
          expected_cigars[i])
      # Check reference positions
      self.assertEqual(reads[i].alignment.position.position,
                       expected_positions[i])
    self.assertLen(reads, len(expected_sequences))


class RealignerIntegrationTest(absltest.TestCase):

  def test_realigner_end2end(self):
    ref_reader = fasta.IndexedFastaReader(testdata.CHR20_FASTA)
    config = realigner.realigner_config(FLAGS)
    reads_realigner = realigner.Realigner(config, ref_reader)
    region_str = 'chr20:10,000,000-10,009,999'
    windows_count = 0

    regions = ranges.RangeSet.from_regions([region_str])
    for region in regions.partition(1000):
      with sam.SamReader(
          testdata.CHR20_BAM,
          read_requirements=reads_pb2.ReadRequirements()) as sam_reader:
        in_reads = list(sam_reader.query(region))
      windows, out_reads = reads_realigner.realign_reads(in_reads, region)

      # We should always get back all of the reads we sent in. Instead of just
      # checking the lengths are the same, make sure all the read names are the
      # same.
      six.assertCountEqual(self, [r.fragment_name for r in in_reads],
                           [r.fragment_name for r in out_reads])

      # Check each window to make sure it's reasonable.
      for window in windows:
        # We always expect the reference sequence to be one of our haplotypes.
        ref_seq = ref_reader.query(window.span)
        self.assertIn(ref_seq, set(window.haplotypes))
      windows_count += len(windows)

    self.assertGreater(windows_count, 0)


class TrimTest(parameterized.TestCase):

  @parameterized.parameters(
      dict(
          cigar='3M2D5M3I10M',
          ref_trim=6,
          ref_length=9,
          expected_cigar='4M3I5M',
          expected_read_trim=4,
          expected_read_length=12,
          comment='Start and end window in different match operations.'),
      dict(
          cigar='30M',
          ref_trim=5,
          ref_length=10,
          expected_cigar='10M',
          expected_read_trim=5,
          expected_read_length=10,
          comment='Start and end window in the same cigar entry'),
      dict(
          cigar='10D10M',
          ref_trim=5,
          ref_length=10,
          expected_cigar='5D5M',
          expected_read_trim=0,
          expected_read_length=5,
          comment='Start window in a deletion'),
      dict(
          cigar='10I10M',
          ref_trim=5,
          ref_length=5,
          expected_cigar='5M',
          expected_read_trim=15,
          expected_read_length=5,
          comment='Start window in an insertion'),
      dict(
          cigar='10M',
          ref_trim=5,
          ref_length=10,
          expected_cigar='5M',
          expected_read_trim=5,
          expected_read_length=5,
          comment='Read ends before the window'),
      dict(
          cigar='10M',
          ref_trim=20,
          ref_length=10,
          expected_cigar='',
          expected_read_trim=10,
          expected_read_length=0,
          comment='Read ends before the trim'),
      dict(
          cigar='10M20D10M',
          ref_trim=12,
          ref_length=5,
          expected_cigar='5D',
          expected_read_trim=10,
          expected_read_length=0,
          comment='Deletion covers the whole window'),
      dict(
          cigar='10M20I10M',
          ref_trim=10,
          ref_length=20,
          expected_cigar='20I10M',
          expected_read_trim=10,
          expected_read_length=30,
          comment='Trim to edge of an insertion'),
      dict(
          cigar='10M2I10M',
          ref_trim=0,
          ref_length=20,
          expected_cigar='10M2I10M',
          expected_read_trim=0,
          expected_read_length=22,
          comment='Zero trim'),
  )
  def test_trim_cigar(self, cigar, ref_trim, ref_length, expected_cigar,
                      expected_read_trim, expected_read_length, comment):
    read = test_utils.make_read('AAAATAAAATAAAATAAAATA', start=100, cigar=cigar)
    output_cigar, output_read_trim, output_read_length = realigner.trim_cigar(
        read.alignment.cigar, ref_trim, ref_length)
    self.assertEqual(
        cigar_utils.format_cigar_units(output_cigar),
        expected_cigar,
        msg='Wrong cigar for: {}'.format(comment))
    self.assertEqual(
        output_read_trim,
        expected_read_trim,
        msg='Wrong read trim for: {}'.format(comment))
    self.assertEqual(
        output_read_length,
        expected_read_length,
        msg='Wrong read length for: {}'.format(comment))
    self.assertEqual(
        cigar_utils.format_cigar_units(read.alignment.cigar),
        cigar,
        msg='Cigar in original read was mutated.')

  @parameterized.parameters([
      # Window region literals are 1-based, but all other coordinates are
      # 0-based: chr1:11-20 means start at 10 and end at 20 (exclusive).
      dict(
          window='chr1:11-20',
          cigar='9M',
          start=8,
          read_length=9,
          expected_cigar='7M',
          expected_position=10,
          expected_read_length=7,
          comment='Trim first 2 bases'),
      dict(
          window='chr1:11-20',
          cigar='9M',
          start=13,
          read_length=9,
          expected_cigar='7M',
          expected_position=13,
          expected_read_length=7,
          comment='Trim last 2 bases'),
      dict(
          window='chr1:11-20',
          cigar='5M',
          start=12,
          read_length=5,
          expected_cigar='5M',
          expected_position=12,
          expected_read_length=5,
          comment='Read fits entirely inside window'),
      dict(
          window='chr1:11-20',
          cigar='9M',
          start=10,
          read_length=9,
          expected_cigar='9M',
          expected_position=10,
          expected_read_length=9,
          comment='Read starts and ends at window edges'),
  ])
  def test_trim_read(self, window, cigar, start, read_length, expected_cigar,
                     expected_position, expected_read_length, comment):
    read = test_utils.make_read(
        'A' * read_length, start=start, cigar=cigar, quals=[30] * read_length)
    region = ranges.parse_literal(window)
    output = realigner.trim_read(read, region)
    self.assertEqual(
        expected_cigar,
        cigar_utils.format_cigar_units(output.alignment.cigar),
        msg='Wrong cigar for case: {}'.format(comment))
    # Start position of the alignment.
    self.assertEqual(
        output.alignment.position.position,
        expected_position,
        msg='Wrong position for case: {}'.format(comment))
    # Read sequence.
    self.assertLen(
        output.aligned_sequence,
        expected_read_length,
        msg='Wrong length of aligned_sequence for case: {}'.format(comment))
    # Base quality scores.
    self.assertLen(
        output.aligned_quality,
        expected_read_length,
        msg='Wrong  length of aligned_quality for case: {}'.format(comment))



if __name__ == '__main__':
  absltest.main()
