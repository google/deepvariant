# Copyright 2018 Google LLC.
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
"""Tests for third_party.nucleus.io.vcf."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from tensorflow.python.platform import gfile
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges


class VcfHeaderCacheTests(parameterized.TestCase):
  """Test the functionality of the VcfHeaderCache class."""

  def setUp(self):
    self.vcf_reader = vcf.VcfReader(
        test_utils.genomics_core_testdata('test_sites.vcf'))
    self.cache = self.vcf_reader.field_access_cache

  @parameterized.parameters(
      'DP',
      'AF',
      'END',
      'ExcessHet',
      'culprit',
  )
  def test_valid_info_get_funcs(self, field_name):
    fn = self.cache.info_field_get_fn(field_name)
    self.assertTrue(callable(fn))

  @parameterized.parameters(
      'DP',
      'AF',
      'END',
      'ExcessHet',
      'culprit',
      'HaplotypeScore',
      'InbreedingCoeff',
  )
  def test_valid_info_set_funcs(self, field_name):
    fn = self.cache.info_field_set_fn(field_name)
    self.assertTrue(callable(fn))

  def test_invalid_info_funcs(self):
    with self.assertRaises(KeyError):
      self.cache.info_field_get_fn('RGQ')
    with self.assertRaises(KeyError):
      self.cache.info_field_set_fn('PID')

  @parameterized.parameters(
      'AD',
      'DP',
      'PID',
      'RGQ',
  )
  def test_valid_format_get_funcs(self, field_name):
    fn = self.cache.format_field_get_fn(field_name)
    self.assertTrue(callable(fn))

  @parameterized.parameters(
      'AD',
      'DP',
      'PID',
      'RGQ',
  )
  def test_valid_format_set_funcs(self, field_name):
    fn = self.cache.format_field_set_fn(field_name)
    self.assertTrue(callable(fn))

  def test_invalid_format_funcs(self):
    with self.assertRaises(KeyError):
      self.cache.format_field_get_fn('culprit')
    with self.assertRaises(KeyError):
      self.cache.format_field_set_fn('ExcessHet')


class VcfReaderTests(absltest.TestCase):
  """Test the iteration functionality provided by vcf.VcfReader."""

  def setUp(self):
    self.sites_reader = vcf.VcfReader(
        test_utils.genomics_core_testdata('test_sites.vcf'))

    self.samples_reader = vcf.VcfReader(
        test_utils.genomics_core_testdata('test_samples.vcf.gz'))

  def test_vcf_iterate(self):
    self.assertEqual(test_utils.iterable_len(self.sites_reader.iterate()), 5)

  def test_vcf_query(self):
    range1 = ranges.parse_literal('chr3:100,000-500,000')
    self.assertEqual(
        test_utils.iterable_len(self.samples_reader.query(range1)), 4)

  def test_vcf_iter(self):
    n = 0
    for _ in self.sites_reader:
      n += 1
    self.assertEqual(n, 5)

  def test_fail_multiple_concurrent_iterations(self):
    range1 = ranges.parse_literal('chr3:100,000-500,000')
    reads = self.samples_reader.query(range1)
    for read in reads:
      pass

    r2 = self.samples_reader.query(range1)
    with self.assertRaisesRegexp(ValueError, 'No underlying iterable. This '):
      next(r2)

  def test_c_reader(self):
    self.assertNotEqual(self.sites_reader.c_reader, 0)
    self.assertNotEqual(self.samples_reader.c_reader, 0)

    tfrecord_reader = vcf.VcfReader(
        test_utils.genomics_core_testdata('test_samples.vcf.golden.tfrecord'))
    self.assertNotEqual(tfrecord_reader.c_reader, 0)


class VcfReaderInputTests(absltest.TestCase):
  """Tests VcfReader behavior on specific inputs."""

  def test_header_format_mixed_order(self):
    """Tests reading a VCF with unconventional FORMAT field definition.

    Tests reading a VCF in which the properties of the format
    fields are defined in mixed order in the header. For example,

    ##FORMAT=<ID=GT,Type=String,Number=1,Description="GT description">

    (In normal VCFs "Number" should come before "Type".)
    """
    with vcf.VcfReader(
        test_utils.genomics_core_testdata(
            'header_format_mixed_order.vcf')) as vreader:
      formats = vreader.header.formats
      variants = list(vreader)
    self.assertLen(formats, 1)
    self.assertEqual(formats[0].id, 'GT')
    self.assertEqual(formats[0].number, '1')
    self.assertEqual(formats[0].type, 'String')
    self.assertEqual(formats[0].description, 'GT description')
    self.assertLen(variants, 2)
    self.assertEqual(variants[0].calls[0].genotype, [0, 1])
    self.assertEqual(variants[1].calls[0].genotype, [1, 1])


def _format_expected_variant(ref, alts, format_spec, *samples):
  base = ['20', 1, '.', ref, alts, 0, '.', '.', format_spec]
  return base + list(samples)


def _format_test_variant(alleles, call_infos):
  variant = test_utils.make_variant(chrom='20', start=0, alleles=alleles)
  for i, call_info in enumerate(call_infos):
    call = variant.calls.add(call_set_name='sample' + str(i))
    for key, value in call_info.items():
      if not isinstance(value, (list, tuple)):
        value = [value]
      call.info[key].values.extend(
          [struct_pb2.Value(int_value=v) for v in value])
  return variant


class VcfWriterTests(parameterized.TestCase):
  """Tests for VcfWriter."""

  def assertWrittenVCFRecordsEqual(self, path, expected_lines):

    def cleanup_line(line):
      if isinstance(line, (list, tuple)):
        return '\t'.join(str(x) for x in line)
      else:
        return line

    expected_lines = [cleanup_line(line) for line in expected_lines]
    with gfile.Open(path, 'r') as fin:
      self.assertEqual([
          line.strip() for line in fin.readlines() if not line.startswith('#')
      ], expected_lines)

  def write_variant_to_tempfile(self, variant):
    output_path = test_utils.test_tmpfile('test.vcf')
    header = variants_pb2.VcfHeader(
        contigs=[reference_pb2.ContigInfo(name='20')],
        sample_names=[call.call_set_name for call in variant.calls],
        formats=[
            variants_pb2.VcfFormatInfo(
                id='DP', number='1', type='Integer', description='Read depth'),
            variants_pb2.VcfFormatInfo(
                id='AD',
                number='R',
                type='Integer',
                description='Read depth for each allele')
        ])
    writer = vcf.VcfWriter(output_path, header=header)
    with writer:
      writer.write(variant)
    return output_path

  @parameterized.parameters(
      # Check that our DP field is getting written out properly.
      (_format_test_variant(['A', 'T'], [{
          'DP': 1
      }, {
          'DP': 2
      }]), _format_expected_variant('A', 'T', 'DP', '1', '2')),
      # Checks that we get the missing value when DP is missing in some samples.
      (_format_test_variant(['A', 'T'], [{
          'DP': 1
      }, {}]), _format_expected_variant('A', 'T', 'DP', '1', '.')),
      (_format_test_variant(['A', 'T'], [{}, {
          'DP': 2
      }]), _format_expected_variant('A', 'T', 'DP', '.', '2')),
  )
  def test_single_value_format_field(self, variant, expected_vcf_line):
    self.assertWrittenVCFRecordsEqual(
        self.write_variant_to_tempfile(variant), [expected_vcf_line])

  @parameterized.parameters(
      # Check that our AD field is getting written correctly.
      (_format_test_variant(['A', 'T'], [{
          'AD': [0, 1]
      }, {
          'AD': [2, 3]
      }]), _format_expected_variant('A', 'T', 'AD', '0,1', '2,3')),
      (_format_test_variant(['A', 'T'], [{}, {
          'AD': [2, 3]
      }]), _format_expected_variant('A', 'T', 'AD', '.', '2,3')),
      (_format_test_variant(['A', 'T'], [{
          'AD': [0, 1]
      }, {}]), _format_expected_variant('A', 'T', 'AD', '0,1', '.')),
      # Let's try a tri-allelic site where we have 3 AD values / sample.
      (_format_test_variant(['A', 'T', 'C'], [{
          'AD': [0, 1, 2]
      }, {
          'AD': [4, 5, 6]
      }]), _format_expected_variant('A', 'T,C', 'AD', '0,1,2', '4,5,6')),
      # Check that we handle missing values properly.
      (_format_test_variant(['A', 'T', 'C'], [{
          'AD': [0, 1, 2]
      }, {}]), _format_expected_variant('A', 'T,C', 'AD', '0,1,2', '.')),
      (_format_test_variant(['A', 'T', 'C'], [{}, {
          'AD': [4, 5, 6]
      }]), _format_expected_variant('A', 'T,C', 'AD', '.', '4,5,6')),
  )
  def test_multi_value_format_field(self, variant, expected_vcf_line):
    self.assertWrittenVCFRecordsEqual(
        self.write_variant_to_tempfile(variant), [expected_vcf_line])

  @parameterized.parameters(
      # Now let's combine some AD and DP fields.
      (_format_test_variant(['A', 'T', 'C'], [{
          'DP': 3,
          'AD': [0, 1, 2]
      }, {
          'DP': 12,
          'AD': [3, 4, 5]
      }]), _format_expected_variant('A', 'T,C', 'DP:AD', '3:0,1,2', '12:3,4,5')
      ),
      (_format_test_variant(['A', 'T', 'C'], [{
          'DP': 3
      }, {
          'AD': [3, 4, 5]
      }]), _format_expected_variant('A', 'T,C', 'DP:AD', '3:.', '.:3,4,5')),
  )
  def test_multiple_format_fields(self, variant, expected_vcf_line):
    self.assertWrittenVCFRecordsEqual(
        self.write_variant_to_tempfile(variant), [expected_vcf_line])


class VcfWriterHeaderlessTests(absltest.TestCase):
  """Tests for VcfWriter with exclude_header=True."""

  def test_headerless_vcf(self):
    """Writes a headerless vcf and reads it back out."""
    test_vcf = test_utils.genomics_core_testdata('test_sites.vcf')
    output_vcf = test_utils.test_tmpfile('output.vcf')
    expected_variants = []
    with vcf.VcfReader(test_vcf) as reader:
      with vcf.VcfWriter(
          output_vcf, header=reader.header, exclude_header=True) as writer:
        for record in reader:
          expected_variants.append(record)
          writer.write(record)

      with vcf.VcfReader(output_vcf, header=reader.header) as actual_reader:
        self.assertEqual(expected_variants, list(actual_reader))


class VcfRoundtripTests(parameterized.TestCase):
  """Test the ability to round-trip VCF files."""

  def setUp(self):
    self.header = (
        '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters '
        'passed">\n##INFO=<ID=DB,Number=0,Type=Flag,Description="In '
        'dbSNP">\n##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Min '
        'DP">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic'
        ' depths">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read '
        'depth">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype '
        'likelihood,Phred-encoded">\n##contig=<ID=chr1,length=248956422>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n')
    self.record_format_strings = [
        'chr1\t13613\t.\tT\tA\t39.88\tPASS\t{info}\t{fmt}\t0/1{efmts1}\t1/1{efmts2}\n',
        'chr1\t13813\trs1\tT\tG\t90.28\tPASS\t{info}\t{fmt}\t1/1{efmts1}\t0|1{efmts2}\n',
        'chr1\t13838\t.\tC\tT\t62.74\tPASS\t{info}\t{fmt}\t0/1{efmts1}\t0/1{efmts2}\n',
    ]

  @parameterized.parameters(
      dict(
          expected_infos=['DB;MIN_DP=4', 'MIN_DP=15', 'DB;MIN_DP=10'],
          expected_fmt='GT:AD:DP:PL',
          expected_fmt1=[
              ':1,3:4:10,5,0', ':11,13:24:55,0,50', ':5,5:10:20,0,20'
          ],
          expected_fmt2=[
              ':1,19:20:100,90,0', ':7,8:15:15,0,12', ':.:10:0,0,50'
          ],
      ),
      dict(
          expected_infos=['DB', '.', 'DB'],
          expected_fmt='GT:AD:DP:PL',
          expected_fmt1=[
              ':1,3:4:10,5,0', ':11,13:24:55,0,50', ':5,5:10:20,0,20'
          ],
          expected_fmt2=[
              ':1,19:20:100,90,0', ':7,8:15:15,0,12', ':.:10:0,0,50'
          ],
          reader_excluded_info=['MIN_DP'],
      ),
      dict(
          expected_infos=['DB', '.', 'DB'],
          expected_fmt='GT',
          expected_fmt1=['', '', ''],
          expected_fmt2=['', '', ''],
          reader_excluded_info=['MIN_DP'],
          reader_excluded_format=['AD', 'DP', 'PL'],
      ),
      dict(
          expected_infos=['DB', '.', 'DB'],
          expected_fmt='GT',
          expected_fmt1=['', '', ''],
          expected_fmt2=['', '', ''],
          writer_excluded_info=['MIN_DP'],
          writer_excluded_format=['AD', 'DP', 'PL'],
      ),
      dict(
          expected_infos=['DB', '.', 'DB'],
          expected_fmt='GT',
          expected_fmt1=['', '', ''],
          expected_fmt2=['', '', ''],
          reader_excluded_info=['MIN_DP'],
          reader_excluded_format=['AD'],
          writer_excluded_info=['MIN_DP'],
          writer_excluded_format=['DP', 'PL'],
      ),
  )
  def test_roundtrip(self,
                     expected_infos,
                     expected_fmt,
                     expected_fmt1,
                     expected_fmt2,
                     reader_excluded_info=None,
                     reader_excluded_format=None,
                     writer_excluded_info=None,
                     writer_excluded_format=None):
    expected_records = [
        record.format(info=info, fmt=expected_fmt, efmts1=e1,
                      efmts2=e2) for record, info, e1, e2 in zip(
                          self.record_format_strings, expected_infos,
                          expected_fmt1, expected_fmt2)
    ]
    expected = self.header + ''.join(expected_records)
    for info_map_pl in [False, True]:
      with vcf.VcfReader(
          test_utils.genomics_core_testdata('test_py_roundtrip.vcf'),
          excluded_info_fields=reader_excluded_info,
          excluded_format_fields=reader_excluded_format,
          store_gl_and_pl_in_info_map=info_map_pl) as reader:
        records = list(reader.iterate())
        output_path = test_utils.test_tmpfile(
            'test_roundtrip_tmpfile_{}.vcf'.format(info_map_pl))
        with vcf.VcfWriter(
            output_path,
            header=reader.header,
            excluded_info_fields=writer_excluded_info,
            excluded_format_fields=writer_excluded_format,
            retrieve_gl_and_pl_from_info_map=info_map_pl) as writer:
          for record in records:
            writer.write(record)

      with open(output_path) as f:
        actual = f.read()
      self.assertEqual(actual, expected)


class InMemoryVcfReaderTests(parameterized.TestCase):
  """Test the functionality provided by vcf.InMemoryVcfReader."""

  def setUp(self):
    self.variants = [
        test_utils.make_variant(chrom='1', start=10),
        test_utils.make_variant(chrom='1', start=20),
        test_utils.make_variant(chrom='1', start=30),
        test_utils.make_variant(chrom='2', start=25),
        test_utils.make_variant(chrom='2', start=55),
        test_utils.make_variant(chrom='3', start=10),
    ]
    self.header = variants_pb2.VcfHeader(
        contigs=[
            reference_pb2.ContigInfo(name='1', n_bases=100),
            reference_pb2.ContigInfo(name='2', n_bases=100),
            reference_pb2.ContigInfo(name='3', n_bases=100),
            reference_pb2.ContigInfo(name='4', n_bases=100),
        ],
        filters=[],
        sample_names=['NA12878'])
    self.reader = vcf.InMemoryVcfReader(
        self.variants, self.header)

  def test_iterate(self):
    """Tests that iterate returns an iterable containing our variants."""
    self.assertEqual(list(self.reader.iterate()), self.variants)

  def test_header(self):
    """Tests that the reader provides us back the header we gave it."""
    self.assertEqual(self.reader.header, self.header)

  @parameterized.parameters(
      dict(query='1', expected_variant_indices=[0, 1, 2]),
      dict(query='2', expected_variant_indices=[3, 4]),
      dict(query='3', expected_variant_indices=[5]),
      dict(query='4', expected_variant_indices=[]),
      dict(query='1:1-15', expected_variant_indices=[0]),
      dict(query='1:1-25', expected_variant_indices=[0, 1]),
      dict(query='1:1-35', expected_variant_indices=[0, 1, 2]),
      dict(query='1:15-35', expected_variant_indices=[1, 2]),
      dict(query='1:25-35', expected_variant_indices=[2]),
  )
  def test_query(self, query, expected_variant_indices):
    range1 = ranges.parse_literal(query, ranges.contigs_dict(
        self.header.contigs))
    self.assertEqual(
        list(self.reader.query(range1)),
        [self.variants[i] for i in expected_variant_indices])


if __name__ == '__main__':
  absltest.main()
