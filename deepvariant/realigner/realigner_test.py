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
"""Tests for deepvariant .realigner.realigner."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import csv
import itertools
import os


from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant import test_utils
from deepvariant.core import genomics_io
from deepvariant.core import io_utils
from deepvariant.core import ranges
from deepvariant.core.genomics import reads_pb2
from deepvariant.core.protos import core_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.realigner import realigner
from deepvariant.realigner import utils

FLAGS = tf.flags.FLAGS


def setUpModule():
  test_utils.init()


def _get_reads(region):
  with genomics_io.make_sam_reader(test_utils.CHR20_BAM) as in_sam_reader:
    return list(in_sam_reader.query(region))


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
          })
      for names in itertools.permutations(
          ['read1', 'read2', 'read3', 'read4', 'read5']))
  def test_assign_reads_to_assembled_regions_multiple_reads(
      self, read_names, expected_assignments):
    self.assertReadsGoToCorrectRegions(
        self.get_reads_by_name(read_names), expected_assignments)

  def assertReadsGoToCorrectRegions(self, reads, expected_assignments):
    unassigned = realigner.assign_reads_to_assembled_regions(
        self.assembled_regions, reads)

    # Every read should be in the assembled regions or unassigned.
    self.assertCountEqual(
        [r for ar in self.assembled_regions for r in ar.reads] + unassigned,
        reads)

    # Go through each region and make sure the reads that are supposed to
    # appear in each region do in appear there.
    for region_name, region in self.regions.iteritems():
      expected_reads = self.get_reads_by_name(
          expected_assignments.get(region_name, []))
      self.assertCountEqual(region.reads, expected_reads)


class RealignerTest(parameterized.TestCase):

  def setUp(self):
    self.ref_reader = genomics_io.make_ref_reader(test_utils.CHR20_FASTA)
    self.config = realigner.realigner_config(FLAGS)
    self.reads_realigner = realigner.Realigner(self.config, self.ref_reader)

  @parameterized.parameters((
      'chr20:10,095,379-10,095,500', 'chr20:10095363-10095543', {
          'TCCTTTTTGTTGTGCAAAAGGAAGTGCTAAAATCAGAATGAGAACCATGGTCACCTGACATAGACACA'
          'AGTGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA'
          'TATCCATGTTCAAGTACTAATTCTGGGCAAGACACTGTTCTAAGTGCTATGAATATATTACCTCAT',
          'TCCTTTTTGTTGTGCAAAAGGAAGTGCTAAAATCAGAATGAGAACCATGGTCACCTGACATAGACACA'
          'AGTGATGATGATGATGATGATGATGATGATGATGATGA'
          'TATCCATGTTCAAGTACTAATTCTGGGCAAGACACTGTTCTAAGTGCTATGAATATATTACCTCAT'
      }, 'There is a heterozygous 9 bp deletion of tandem TGA repeat.'
  ), ('chr20:10,046,080-10,046,307', 'chr20:10046107-10046276', {
      'AGTTAGGGATGCTGGAAAGGCAGAAAGAAAAGGGAAGGGAAGAGGAAGGGGAAAAGGAAAGAAAAAAAAG'
      'AAAGAAAGAAAGAGAAAGAAAGAGAAAGAGAAAGAAAGAGGAAAGAGAGA'
      'AAGAGAAAGAGAAGGAAAGAGAAAGAAAGAGAAGGAAAGAGAGAAAGAGA',
      'AGTTAGGGATGCTGGAAAGGCAGAAAGAAAAGGGAAGGGAAGAGGAAGGGGAAAAGGAAAGAAAAAAAAG'
      'AAAGAAAGAAAGAGAAAGAGAAAGAAAGAGGAAAGAGAGAAAGAGAAAGA'
      'GAAGGAAAGAGAAAGAAAGAGAAGGAAAGAGAGAAAGAGA'
  }, 'There is a heterozygous 10 bp deletion.'))
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
        windows_haplotypes[0].span, comment)
    self.assertEqual(expected_haplotypes, set(windows_haplotypes[0].haplotypes),
                     comment)

  @parameterized.parameters(('chr20:10,046,080-10,046,307',
                             'chr20:10,046,179-10,046,188'))
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
    read = test_utils.make_read(
        'ACCGT' * 50,
        start=63025520 - 250,
        cigar='250M',
        quals=range(30, 35) * 50,
        name='read1')
    reads = [read] * 20
    region = ranges.parse_literal('chr20:63,025,320-63,025,520')
    self.reads_realigner.realign_reads(reads, region)

    # These reads are aligned off the edge of the contig.
    read = test_utils.make_read(
        'TTATA' * 50,
        start=63025520 - 200,
        cigar='200M50S',
        quals=range(30, 35) * 50,
        name='read1')
    reads = [read] * 20
    self.reads_realigner.realign_reads(reads, region)

  @parameterized.parameters(
      dict(enabled=False, emit_reads=False),
      dict(enabled=True, emit_reads=False),
      dict(enabled=True, emit_reads=True),
  )
  def test_realigner_diagnostics(self, enabled, emit_reads):
    # Make sure that by default we aren't emitting any diagnostic outputs.
    dx_dir = test_utils.test_tmpfile('dx')
    region_str = 'chr20:10046179-10046188'
    region = ranges.parse_literal(region_str)
    assembled_region_str = 'chr20:10046109-10046257'
    reads = _get_reads(region)
    self.config = realigner.realigner_config(FLAGS)
    self.config.diagnostics.enabled = enabled
    self.config.diagnostics.output_root = dx_dir
    self.config.diagnostics.emit_realigned_reads = emit_reads
    self.reads_realigner = realigner.Realigner(self.config, self.ref_reader)
    _, realigned_reads = self.reads_realigner.realign_reads(reads, region)
    self.reads_realigner.diagnostic_logger.close()  # Force close all resources.

    if not enabled:
      # Make sure our diagnostic output isn't emitted.
      self.assertFalse(tf.gfile.Exists(dx_dir))
    else:
      # Our root directory exists.
      self.assertTrue(tf.gfile.IsDirectory(dx_dir))

      # We expect a realigner_metrics.csv in our rootdir with 1 entry in it.
      metrics_file = os.path.join(
          dx_dir, self.reads_realigner.diagnostic_logger.metrics_filename)
      self.assertTrue(tf.gfile.Exists(metrics_file))
      with tf.gfile.FastGFile(metrics_file) as fin:
        rows = list(csv.DictReader(fin))
        self.assertEqual(len(rows), 1)
        self.assertEqual(
            set(rows[0].keys()), {'window', 'k', 'n_haplotypes', 'time'})
        self.assertEqual(rows[0]['window'], assembled_region_str)
        self.assertEqual(int(rows[0]['k']), 25)
        self.assertTrue(int(rows[0]['n_haplotypes']), 2)
        # Check that our runtime is reasonable (greater than 0, less than 10 s).
        self.assertTrue(0.0 < float(rows[0]['time']) < 10.0)

      # As does the subdirectory for this region.
      region_subdir = os.path.join(dx_dir, assembled_region_str)
      self.assertTrue(tf.gfile.IsDirectory(region_subdir))

      # We always have a graph.dot
      self.assertTrue(
          tf.gfile.Exists(
              os.path.join(
                  region_subdir,
                  self.reads_realigner.diagnostic_logger.graph_filename)))

      reads_file = os.path.join(
          dx_dir, region_str,
          self.reads_realigner.diagnostic_logger.realigned_reads_filename)
      if emit_reads:
        self.assertTrue(tf.gfile.Exists(reads_file))
        reads_from_dx = io_utils.read_tfrecords(reads_file, reads_pb2.Read)
        self.assertCountEqual(reads_from_dx, realigned_reads)
      else:
        self.assertFalse(tf.gfile.Exists(reads_file))


class RealignerIntegrationTest(absltest.TestCase):

  def test_realigner_end2end(self):
    ref_reader = genomics_io.make_ref_reader(test_utils.CHR20_FASTA)
    config = realigner.realigner_config(FLAGS)
    reads_realigner = realigner.Realigner(config, ref_reader)
    region_str = 'chr20:10,000,000-10,009,999'

    regions = ranges.RangeSet.from_regions([region_str])
    for region in regions.partition(1000):
      with genomics_io.make_sam_reader(
          test_utils.CHR20_BAM, core_pb2.ReadRequirements()) as sam_reader:
        in_reads = list(sam_reader.query(region))
      windows, out_reads = reads_realigner.realign_reads(in_reads, region)

      # We should always get back all of the reads we sent in. Instead of just
      # checking the lengths are the same, make sure all the read names are the
      # same.
      self.assertCountEqual([r.fragment_name for r in in_reads],
                            [r.fragment_name for r in out_reads])

      # Make sure we assembled at least one windows in the region.
      self.assertNotEqual(0, len(windows))

      # Check each window to make sure it's reasonable.
      for window in windows:
        # We always expect the reference sequence to be one of our haplotypes.
        ref_seq = ref_reader.bases(window.span)
        self.assertIn(ref_seq, set(window.haplotypes))


if __name__ == '__main__':
  absltest.main()
