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
"""Tests for VcfWriter CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant.core import genomics_io
from deepvariant.core import io_utils
from deepvariant.core import test_utils
from deepvariant.core.genomics import variants_pb2
from deepvariant.core.protos import core_pb2
from deepvariant.core.python import vcf_writer

_DOUBLE_CLOSE_ERROR = 'Cannot close an already closed VcfWriter'
_WRITE_TO_CLOSED_ERROR = 'Cannot write to closed VCF stream'
_WRONG_NUMBER_OF_SAMPLES = 'Variant call count must match number of samples'
_DISCORDANT_SAMPLE_NAMES_ERROR = (
    'Out-of-order call set names, or unrecognized call set name, with respect '
    'to samples declared in VCF header.')
_UNKNOWN_CONTIG_ERROR = "Record's reference name is not available in VCF header"
_FILTER_NOT_FOUND_ERROR = 'Filter must be found in header'


class WrapVcfWriterTest(parameterized.TestCase):

  def setUp(self):
    self.out_fname = test_utils.test_tmpfile('output.vcf')
    self.options = core_pb2.VcfWriterOptions(
        contigs=[
            core_pb2.ContigInfo(name='Chr1', n_bases=50, pos_in_fasta=0),
            core_pb2.ContigInfo(name='Chr2', n_bases=25, pos_in_fasta=1),
        ],
        sample_names=['Fido', 'Spot'],
        filters=[])
    self.writer = vcf_writer.VcfWriter.to_file(self.out_fname, self.options)
    self.variant = test_utils.make_variant(
        chrom='Chr1', start=10, alleles=['A', 'C'])
    self.variant.calls.add(genotype=[0, 0], call_set_name='Fido')
    self.variant.calls.add(genotype=[0, 1], call_set_name='Spot')

  def test_writing_canned_variants(self):
    """Tests writing all the variants that are 'canned' in our tfrecord file."""

    # This file is in the TF record format
    tfrecord_file = test_utils.genomics_core_testdata(
        'test_samples.vcf.golden.tfrecord')

    writer_options = core_pb2.VcfWriterOptions(
        contigs=[
            core_pb2.ContigInfo(name='chr1', n_bases=248956422),
            core_pb2.ContigInfo(name='chr2', n_bases=242193529),
            core_pb2.ContigInfo(name='chr3', n_bases=198295559),
            core_pb2.ContigInfo(name='chrX', n_bases=156040895)
        ],
        sample_names=['NA12878_18_99'],
        filters=[
            core_pb2.VcfFilterInfo(id='LowQual'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL95.00to96.00'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL96.00to97.00'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL97.00to99.00'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL99.00to99.50'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL99.50to99.90'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL99.90to99.95'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL99.95to100.00+'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheINDEL99.95to100.00'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheSNP99.50to99.60'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheSNP99.60to99.80'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheSNP99.80to99.90'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheSNP99.90to99.95'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheSNP99.95to100.00+'),
            core_pb2.VcfFilterInfo(id='VQSRTrancheSNP99.95to100.00'),
        ])

    variant_records = list(
        io_utils.read_tfrecords(tfrecord_file, proto=variants_pb2.Variant))
    out_fname = test_utils.test_tmpfile('output.vcf')
    with vcf_writer.VcfWriter.to_file(out_fname, writer_options) as writer:
      for record in variant_records[:5]:
        writer.write(record)

    # Check: are the variants written as expected?
    # pylint: disable=line-too-long
    expected_vcf_content = [
        '##fileformat=VCFv4.2\n',
        '##FILTER=<ID=PASS,Description="All filters passed">\n',
        '##FILTER=<ID=LowQual,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL95.00to96.00,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL96.00to97.00,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL97.00to99.00,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL99.00to99.50,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL99.50to99.90,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL99.90to99.95,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL99.95to100.00+,Description="">\n',
        '##FILTER=<ID=VQSRTrancheINDEL99.95to100.00,Description="">\n',
        '##FILTER=<ID=VQSRTrancheSNP99.50to99.60,Description="">\n',
        '##FILTER=<ID=VQSRTrancheSNP99.60to99.80,Description="">\n',
        '##FILTER=<ID=VQSRTrancheSNP99.80to99.90,Description="">\n',
        '##FILTER=<ID=VQSRTrancheSNP99.90to99.95,Description="">\n',
        '##FILTER=<ID=VQSRTrancheSNP99.95to100.00+,Description="">\n',
        '##FILTER=<ID=VQSRTrancheSNP99.95to100.00,Description="">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth of all '
        'passing filters reads.">\n',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth of all '
        'passing filters reads for each allele.">\n',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele '
        'fractions.">\n',
        '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype '
        'likelihoods, log10 encoded">\n',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype '
        'likelihoods, Phred encoded">\n',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of '
        'the interval">\n', '##contig=<ID=chr1,length=248956422>\n',
        '##contig=<ID=chr2,length=242193529>\n',
        '##contig=<ID=chr3,length=198295559>\n',
        '##contig=<ID=chrX,length=156040895>\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878_18_99\n',
        'chr1\t13613\t.\tT\tA\t39.88\tVQSRTrancheSNP99.90to99.95\t.\tGT:GQ:DP:AD:PL\t0/1:16:4:1,3:68,0,16\n',
        'chr1\t13813\t.\tT\tG\t90.28\tPASS\t.\tGT:GQ:DP:AD:PL\t1/1:9:3:0,3:118,9,0\n',
        'chr1\t13838\trs28428499\tC\tT\t62.74\tPASS\t.\tGT:GQ:DP:AD:PL\t1/1:6:2:0,2:90,6,0\n',
        'chr1\t14397\trs756427959\tCTGT\tC\t37.73\tPASS\t.\tGT:GQ:DP:AD:PL\t0/1:75:5:3,2:75,0,152\n',
        'chr1\t14522\t.\tG\tA\t49.77\tVQSRTrancheSNP99.60to99.80\t.\tGT:GQ:DP:AD:PL\t0/1:78:10:6,4:78,0,118\n'
    ]
    # pylint: enable=line-too-long

    with tf.gfile.GFile(out_fname, 'r') as f:
      self.assertEqual(f.readlines(), expected_vcf_content)

  def test_write_variant_is_ok(self):
    self.assertIsNone(self.writer.write(self.variant))

  def test_write_raises_with_unknown_contig(self):
    with self.assertRaisesRegexp(ValueError, _UNKNOWN_CONTIG_ERROR):
      self.variant.reference_name = 'BadChrom'
      self.writer.write(self.variant)

  def test_write_raises_with_unknown_filter(self):
    with self.assertRaisesRegexp(ValueError, _FILTER_NOT_FOUND_ERROR):
      self.variant.filter[:] = ['BadFilter']
      self.writer.write(self.variant)

  @parameterized.parameters(
      ([], _WRONG_NUMBER_OF_SAMPLES),
      (['Spot'], _WRONG_NUMBER_OF_SAMPLES),
      (['Fido'], _WRONG_NUMBER_OF_SAMPLES),
      (['Unknown', 'Fido'], _DISCORDANT_SAMPLE_NAMES_ERROR),
      (['Spot', 'Unknown'], _DISCORDANT_SAMPLE_NAMES_ERROR),
      (['Spot', 'Fido'], _DISCORDANT_SAMPLE_NAMES_ERROR),  # Out of order.
      (['Fido', 'Spot', 'Extra'], _WRONG_NUMBER_OF_SAMPLES),
  )
  def test_write_raises_with_unknown_sample(self, sample_names, message):
    with self.assertRaisesRegexp(ValueError, message):
      del self.variant.calls[:]
      for sample_name in sample_names:
        self.variant.calls.add(genotype=[0, 0], call_set_name=sample_name)
      self.writer.write(self.variant)

  def test_context_manager(self):
    with self.writer:
      # Writing within the context manager succeeds.
      self.assertIsNone(self.writer.write(self.variant))

    # self.writer should be closed, so writing again will fail.
    with self.assertRaisesRegexp(ValueError, _WRITE_TO_CLOSED_ERROR):
      self.writer.write(self.variant)

  def test_double_context_manager(self):
    with self.writer:
      # Writing within the context manager succeeds.
      self.assertIsNone(self.writer.write(self.variant))

    with self.assertRaisesRegexp(ValueError, _DOUBLE_CLOSE_ERROR):
      # Entering the closed writer should be fine.
      with self.writer:
        pass  # We want to raise an error on exit, so nothing to do in context.


class WrapVcfWriterRoundTripTests(parameterized.TestCase):

  @parameterized.parameters(('test_samples.vcf',), ('test_samples.vcf.gz',),
                            ('test_sites.vcf',))
  def test_round_trip_vcf(self, test_datum_name):
    # Round-trip variants through writing and reading:
    # 1. Read variants v1 from VcfReader;
    # 2. Write v1 to vcf using our VcfWriter;
    # 3. Read back in using VcfReader -- v2;
    # 4. compare v1 and v2.
    in_file = test_utils.genomics_core_testdata(test_datum_name)
    out_file = test_utils.test_tmpfile('output_' + test_datum_name)

    v1_reader = genomics_io.make_vcf_reader(in_file, use_index=False)
    v1_records = list(v1_reader.iterate())
    self.assertTrue(v1_records, 'Reader failed to find records')

    writer_options = core_pb2.VcfWriterOptions(
        contigs=v1_reader.contigs,
        sample_names=v1_reader.samples,
        filters=v1_reader.filters)

    with vcf_writer.VcfWriter.to_file(out_file, writer_options) as writer:
      for record in v1_records:
        writer.write(record)

    v2_reader = genomics_io.make_vcf_reader(out_file, use_index=False)
    v2_records = list(v2_reader.iterate())

    self.assertEqual(v1_records, v2_records,
                     'Round-tripped variants not as expected')


if __name__ == '__main__':
  absltest.main()
