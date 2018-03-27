# Copyright 2018 Google Inc.
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
"""Tests for deepvariant.util.io.vcf."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.util.io import vcf
from deepvariant.util.genomics import reference_pb2
from deepvariant.util.genomics import struct_pb2
from deepvariant.util.genomics import variants_pb2
from deepvariant.util import ranges
from deepvariant.util import test_utils
from tensorflow.python.platform import gfile


class VcfHeaderCacheTests(parameterized.TestCase):
  """Test the functionality of the VcfHeaderCache class."""

  def setUp(self):
    self.vcf_reader = vcf.VcfReader(
        test_utils.genomics_core_testdata('test_sites.vcf'), use_index=False)
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
        test_utils.genomics_core_testdata('test_sites.vcf'), use_index=False)

    self.samples_reader = vcf.VcfReader(
        test_utils.genomics_core_testdata('test_samples.vcf.gz'),
        use_index=True)

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


def _format_expected_variant(ref, alts, format_spec, *samples):
  base = ['20', 1, '.', ref, alts, 0, '.', '.', format_spec]
  return base + list(samples)


def _format_test_variant(alleles, call_infos):
  variant = test_utils.make_variant(chrom='20', start=0, alleles=alleles)
  for i, call_info in enumerate(call_infos):
    call = variant.calls.add(call_set_name='sample' + str(i))
    for key, value in call_info.iteritems():
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
    with gfile.FastGFile(path, 'r') as fin:
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
