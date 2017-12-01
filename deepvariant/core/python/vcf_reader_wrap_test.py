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
"""Tests for vcf_reader CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest

from deepvariant.core import ranges
from deepvariant.core import test_utils
from deepvariant.core.protos import core_pb2
from deepvariant.core.python import vcf_reader

expected_sites_contigs = [
    core_pb2.ContigInfo(name='chr1', pos_in_fasta=0, n_bases=248956422),
    core_pb2.ContigInfo(name='chr2', pos_in_fasta=1, n_bases=242193529),
    core_pb2.ContigInfo(name='chr3', pos_in_fasta=2, n_bases=198295559),
    core_pb2.ContigInfo(name='chr4', pos_in_fasta=3, n_bases=190214555),
    core_pb2.ContigInfo(name='chr5', pos_in_fasta=4, n_bases=181538259),
    core_pb2.ContigInfo(name='chr6', pos_in_fasta=5, n_bases=170805979),
    core_pb2.ContigInfo(name='chr7', pos_in_fasta=6, n_bases=159345973),
    core_pb2.ContigInfo(name='chr8', pos_in_fasta=7, n_bases=145138636),
    core_pb2.ContigInfo(name='chr9', pos_in_fasta=8, n_bases=138394717),
    core_pb2.ContigInfo(name='chr10', pos_in_fasta=9, n_bases=133797422),
    core_pb2.ContigInfo(name='chr11', pos_in_fasta=10, n_bases=135086622),
    core_pb2.ContigInfo(name='chr12', pos_in_fasta=11, n_bases=133275309),
    core_pb2.ContigInfo(name='chr13', pos_in_fasta=12, n_bases=114364328),
    core_pb2.ContigInfo(name='chr14', pos_in_fasta=13, n_bases=107043718),
    core_pb2.ContigInfo(name='chr15', pos_in_fasta=14, n_bases=101991189),
    core_pb2.ContigInfo(name='chr16', pos_in_fasta=15, n_bases=90338345),
    core_pb2.ContigInfo(name='chr17', pos_in_fasta=16, n_bases=83257441),
    core_pb2.ContigInfo(name='chr18', pos_in_fasta=17, n_bases=80373285),
    core_pb2.ContigInfo(name='chr19', pos_in_fasta=18, n_bases=58617616),
    core_pb2.ContigInfo(name='chr20', pos_in_fasta=19, n_bases=64444167),
    core_pb2.ContigInfo(name='chr21', pos_in_fasta=20, n_bases=46709983),
    core_pb2.ContigInfo(name='chr22', pos_in_fasta=21, n_bases=50818468),
    core_pb2.ContigInfo(name='chrX', pos_in_fasta=22, n_bases=156040895),
    core_pb2.ContigInfo(name='chrY', pos_in_fasta=23, n_bases=57227415),
    core_pb2.ContigInfo(name='chrM', pos_in_fasta=24, n_bases=16569),
]

# pylint: disable=line-too-long
expected_samples_filters = [
    core_pb2.VcfFilterInfo(id='PASS', description='All filters passed'),
    core_pb2.VcfFilterInfo(id='LowQual', description='Low	quality'),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL95.00to96.00',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	0.9364	<=	x	<	1.0415'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL96.00to97.00',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	0.8135	<=	x	<	0.9364'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL97.00to99.00',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	0.323	<=	x	<	0.8135'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL99.00to99.50',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	-0.1071	<=	x	<	0.323'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL99.50to99.90',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	-1.845	<=	x	<	-0.1071'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL99.90to99.95',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	-3.2441	<=	x	<	-1.845'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL99.95to100.00+',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod	<	-57172.0693'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheINDEL99.95to100.00',
        description=
        'Truth	sensitivity	tranche	level	for	INDEL	model	at	VQS	Lod:	-57172.0693	<=	x	<	-3.2441'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheSNP99.50to99.60',
        description=
        'Truth	sensitivity	tranche	level	for	SNP	model	at	VQS	Lod:	-0.751	<=	x	<	-0.6681'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheSNP99.60to99.80',
        description=
        'Truth	sensitivity	tranche	level	for	SNP	model	at	VQS	Lod:	-1.0839	<=	x	<	-0.751'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheSNP99.80to99.90',
        description=
        'Truth	sensitivity	tranche	level	for	SNP	model	at	VQS	Lod:	-1.7082	<=	x	<	-1.0839'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheSNP99.90to99.95',
        description=
        'Truth	sensitivity	tranche	level	for	SNP	model	at	VQS	Lod:	-3.0342	<=	x	<	-1.7082'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheSNP99.95to100.00+',
        description=
        'Truth	sensitivity	tranche	level	for	SNP	model	at	VQS	Lod	<	-40235.9641'
    ),
    core_pb2.VcfFilterInfo(
        id='VQSRTrancheSNP99.95to100.00',
        description=
        'Truth	sensitivity	tranche	level	for	SNP	model	at	VQS	Lod:	-40235.9641	<=	x	<	-3.0342'
    )
]

# pylint: enable=line-too-long


class WrapVcfReaderTests(absltest.TestCase):

  def setUp(self):
    self.unindexed_options = core_pb2.VcfReaderOptions()
    self.indexed_options = core_pb2.VcfReaderOptions(
        index_mode=core_pb2.INDEX_BASED_ON_FILENAME)
    self.sites_vcf = test_utils.genomics_core_testdata('test_sites.vcf')
    self.samples_vcf = test_utils.genomics_core_testdata('test_samples.vcf.gz')
    self.sites_reader = vcf_reader.VcfReader.from_file(self.sites_vcf,
                                                       self.unindexed_options)
    self.samples_reader = vcf_reader.VcfReader.from_file(
        self.samples_vcf, self.indexed_options)

  def test_vcf_iterate(self):
    iterable = self.sites_reader.iterate()
    self.assertEqual(test_utils.iterable_len(iterable), 5)

  def test_vcf_contigs(self):
    self.assertEqual(expected_sites_contigs, self.sites_reader.contigs)

  def test_vcf_filters(self):
    self.assertEqual(expected_samples_filters, self.samples_reader.filters)

  def test_vcf_samples(self):
    self.assertEqual(self.sites_reader.samples, [])
    self.assertEqual(self.samples_reader.samples, ['NA12878_18_99'])

  def test_vcf_query(self):
    range1 = ranges.parse_literal('chr3:100,000-500,000')
    iterable = self.samples_reader.query(range1)
    self.assertEqual(test_utils.iterable_len(iterable), 4)

  def test_from_file_raises_with_missing_source(self):
    with self.assertRaisesRegexp(ValueError,
                                 'Not found: Could not open missing.vcf'):
      vcf_reader.VcfReader.from_file('missing.vcf', self.unindexed_options)

  def test_from_file_raises_with_missing_index(self):
    with self.assertRaisesRegexp(ValueError, 'Not found: No index found for'):
      vcf_reader.VcfReader.from_file(
          test_utils.genomics_core_testdata('test_sites.vcf'),
          self.indexed_options)

  def test_ops_on_closed_reader_raise(self):
    with self.samples_reader:
      pass
    # At this point the reader is closed.
    with self.assertRaisesRegexp(ValueError, 'Cannot Iterate a closed'):
      self.samples_reader.iterate()
    with self.assertRaisesRegexp(ValueError, 'Cannot Query a closed'):
      self.samples_reader.query(
          ranges.parse_literal('chr1:10,000,000-10,000,100'))

  def test_query_on_unindexed_reader_raises(self):
    with vcf_reader.VcfReader.from_file(self.samples_vcf,
                                        self.unindexed_options) as reader:
      with self.assertRaisesRegexp(ValueError, 'Cannot query without an index'):
        reader.query(ranges.parse_literal('chr1:10,000,000-10,000,100'))

  def test_query_raises_with_bad_range(self):
    with self.assertRaisesRegexp(ValueError, 'Unknown reference_name'):
      self.samples_reader.query(ranges.parse_literal('XXX:1-10'))
    with self.assertRaisesRegexp(ValueError, 'Malformed region'):
      self.samples_reader.query(ranges.parse_literal('chr1:0-5'))
    with self.assertRaisesRegexp(ValueError, 'Malformed region'):
      self.samples_reader.query(ranges.parse_literal('chr1:6-5'))
    with self.assertRaisesRegexp(ValueError, 'Malformed region'):
      self.samples_reader.query(ranges.parse_literal('chr1:10-5'))

  def test_context_manager(self):
    with vcf_reader.VcfReader.from_file(self.sites_vcf,
                                        self.unindexed_options) as f:
      self.assertEqual(expected_sites_contigs, f.contigs)

  # Commented out because we in fact don't detect the malformed VCF yet. It is
  # unclear if it's even possible to detect the issue with the API provided by
  # htslib.
  # def test_vcf_iterate_raises_on_malformed_record(self):
  #   malformed = test_utils.genomics_core_testdata('malformed.vcf')
  #   reader = vcf_reader.VcfReader.from_file(malformed, self.unindexed_options)
  #   iterable = iter(reader.iterate())
  #   self.assertIsNotNone(next(iterable))
  #   with self.assertRaises(ValueError):
  #     print(list(iterable))


if __name__ == '__main__':
  absltest.main()
