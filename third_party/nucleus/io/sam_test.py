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
"""Tests for third_party.nucleus.util.io."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import itertools

from absl.testing import absltest
from absl.testing import parameterized

import six

from tensorflow.python.platform import gfile
from third_party.nucleus.io import sam
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges


class SamReaderTests(parameterized.TestCase):
  """Test the iteration functionality provided by io.SamReader."""

  def test_sam_iterate(self):
    reader = sam.SamReader(test_utils.genomics_core_testdata('test.sam'))
    with reader:
      self.assertEqual(test_utils.iterable_len(reader.iterate()), 6)

  def test_bam_iterate(self):
    reader = sam.SamReader(test_utils.genomics_core_testdata('test.bam'))
    with reader:
      self.assertEqual(test_utils.iterable_len(reader.iterate()), 106)

  def test_bam_iterate_partially(self):
    """Verify that iteration provides results incrementally, not all at once."""
    reader = sam.SamReader(test_utils.genomics_core_testdata('test.bam'))
    with reader:
      iterable = reader.iterate()
      # We expect 106 records in total.
      for _ in range(10):
        results = list(itertools.islice(iterable, 10))
        self.assertEqual(len(results), 10)
      results = list(itertools.islice(iterable, 10))
      self.assertEqual(len(results), 6)

  def test_sam_query(self):
    reader = sam.SamReader(test_utils.genomics_core_testdata('test.bam'))
    expected = [(ranges.parse_literal('chr20:10,000,000-10,000,100'), 106),
                (ranges.parse_literal('chr20:10,000,000-10,000,000'), 45)]
    with reader:
      for interval, n_expected in expected:
        with reader.query(interval) as iterable:
          self.assertEqual(test_utils.iterable_len(iterable), n_expected)

  def test_sam_query_alternate_index_name(self):
    reader = sam.SamReader(
        test_utils.genomics_core_testdata('test_alternate_index.bam'))
    expected = [(ranges.parse_literal('chr20:10,000,000-10,000,100'), 106),
                (ranges.parse_literal('chr20:10,000,000-10,000,000'), 45)]
    with reader:
      for interval, n_expected in expected:
        with reader.query(interval) as iterable:
          self.assertEqual(test_utils.iterable_len(iterable), n_expected)

  @parameterized.parameters(('\t'.join(x[0] for x in items), {
      k: v for t in items for k, v in t[1].items()
  }) for r in [1, 2] for items in itertools.permutations(
      [
          ('X1:i:0', {
              'X1': 0
          }),
          ('X2:i:127', {
              'X2': 127
          }),
          ('X3:i:255', {
              'X3': 255
          }),
          ('X4:A:a', {
              'X4': 'a'
          }),
          ('X5:f:1.234', {
              'X5': 1.234
          }),
          ('X6:Z:string', {
              'X6': 'string'
          }),
          ('X7:Z:with spaces', {
              'X7': 'with spaces'
          }),
          ('ZJ:Z:', {
              'ZJ': ''
          }),  # Empty string.
          # We skip H hex byte-array tags as they appear deprecated.
          ('X8:H:1AE301', {}),
          ('X9:B:i,1', {
              'X9': [1]
          }),
          ('XA:B:i,1,2', {
              'XA': [1, 2]
          }),
          ('XB:B:i,1,2,3', {
              'XB': [1, 2, 3]
          }),
          ('XC:B:i,1,2,3,4,5,6,7,8,9,10', {
              'XC': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
          }),
          ('XD:B:I,1,2,3', {
              'XD': [1, 2, 3]
          }),
          ('XE:B:c,1,2,3', {
              'XE': [1, 2, 3]
          }),
          ('XF:B:C,1,2,3', {
              'XF': [1, 2, 3]
          }),
          ('XG:B:f,0.12,0.34', {
              'XG': [0.12, 0.34]
          }),
          ('XH:B:s,1,2,3', {
              'XH': [1, 2, 3]
          }),
          ('XI:B:S,1,2,3', {
              'XI': [1, 2, 3]
          }),
      ],
      r=r))
  def test_parsing_aux_tags(self, tag_string, expected_info):
    # Minimal header line to create a valid SAM file.
    reads = self._parse_read_with_aux_tags(tag_string)
    self.assertLen(reads, 1)
    self.assertInfoMapEqual(reads[0].info, expected_info)

  @parameterized.parameters(
      '\t'.join(tags) for r in [1, 2, 3] for tags in itertools.permutations(
          [
              'X2:i:x',  # Integer with character value.
              'X3:f:string',  # A string instead of the expected float.
              'X4:A:',  # Supposed to be single char, but we none here.
              'X5:A:ab',  # Supposed to be single char, but we have two here.
              'X6:B:i',  # Empty byte array.
              'X7:B:i,1.23',  # Integer byte array with a float.
              'X8:B:i,string',  # Integer byte array with a string.
              'X8:B:f,string',  # Float byte array with a string.
              'X8:B:z,1,2,3',  # z is not a valid subtype.
          ],
          r=r))
  def test_survives_malformed_lines(self, tag_string):
    try:
      reads = self._parse_read_with_aux_tags(tag_string)
      # If we didn't detect the error, make sure we actually still parsed the
      # read itself.
      self.assertLen(reads, 1)
      self.assertEqual(reads[0].fragment_name, 'read_name')
      self.assertEqual(reads[0].aligned_sequence, 'CCC')
      self.assertEqual(reads[0].alignment.position.reference_name, 'chr1')
      self.assertEqual(reads[0].alignment.position.position, 0)
    except ValueError as e:
      if 'Failed to parse SAM record' not in str(e):
        self.fail('Parsing failed but unexpected exception was seen: ' + str(e))

  def _parse_read_with_aux_tags(self, tag_string):
    # Minimal header line to create a valid SAM file.
    header_lines = '@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n'
    # A single stock read we'll add our AUX fields to.
    read = 'read_name\t0\tchr1\t1\t0\t3M\t*\t0\t0\tCCC\tAAA\t' + tag_string
    path = test_utils.test_tmpfile('aux_tags.bam')
    with gfile.Open(path, 'w') as fout:
      fout.write(header_lines)
      fout.write(read + '\n')
    with sam.SamReader(path, parse_aux_fields=True) as reader:
      return list(reader.iterate())

  def assertInfoMapEqual(self, info_map, expected_info):
    self.assertCountEqual(
        info_map.keys(), expected_info.keys(),
        'info has {} keys but we expected {}'.format(info_map.keys(),
                                                     expected_info.keys()))
    for key, expected_values in expected_info.items():
      if not isinstance(expected_values, list):
        expected_values = [expected_values]
      for actual_value, expected_value in zip(info_map[key].values,
                                              expected_values):
        if isinstance(expected_value, float):
          self.assertAlmostEqual(actual_value.number_value, expected_value)
        elif isinstance(expected_value, six.integer_types):
          self.assertEqual(actual_value.int_value, expected_value)
        elif isinstance(expected_value, str):
          self.assertEqual(actual_value.string_value, expected_value)
        else:
          self.fail('Unsupported expected_value type {}'.format(expected_value))

  @parameterized.parameters(
      # These expected counts are deterministic because we always set the random
      # seed in each test.
      # There are 106 total reads if we iterate.
      ('iterate', None, 1.0, 106),
      ('iterate', None, 0.5, 59),
      ('iterate', None, 0.25, 31),
      # There are 45 total reads if we don't downsample.
      ('query', 'chr20:10,000,000-10,000,000', 1.0, 45),
      ('query', 'chr20:10,000,000-10,000,000', 0.5, 25),
      ('query', 'chr20:10,000,000-10,000,000', 0.25, 13),
  )
  def test_downsampling(self, method, maybe_range, fraction, expected_n_reads):
    reader = sam.SamReader(
        test_utils.genomics_core_testdata('test.bam'),
        downsample_fraction=fraction,
        random_seed=12345)
    with reader:
      if method == 'iterate':
        reads_iter = reader.iterate()
      elif method == 'query':
        reads_iter = reader.query(ranges.parse_literal(maybe_range))
      else:
        self.fail('Unexpected method ' + str(method))
      self.assertEqual(test_utils.iterable_len(reads_iter), expected_n_reads)


# Note that CRAM version 2.1 files work with Nucleus but they cannot be used in
# our test here because CRAM 2.1 embeds an exact path to the reference file
# which LEAKR flags as leaking internal google paths.
@parameterized.parameters(
    dict(
        filename='test_cram.embed_ref_0_version_3.0.cram',
        has_embedded_ref=False),
    dict(
        filename='test_cram.embed_ref_1_version_3.0.cram',
        has_embedded_ref=True),
)
class CramReaderTests(parameterized.TestCase):
  """Test io.SamReader on CRAM formatted files."""

  def _make_reader(self, filename, has_embedded_ref):
    if has_embedded_ref:
      # If we have an embedded reference, force the reader to use it by not
      # providing an argument for ref_path.
      return sam.SamReader(test_utils.genomics_core_testdata(filename))
    else:
      # Otherwise we need to explicitly override the reference encoded in the UR
      # of the CRAM file to use the path provided to our test.fasta.
      return sam.SamReader(
          test_utils.genomics_core_testdata(filename),
          ref_path=test_utils.genomics_core_testdata('test.fasta'))

  def test_header(self, filename, has_embedded_ref):
    with self._make_reader(filename, has_embedded_ref) as reader:
      self.assertEqual(reader.header.format_version, '1.3')
      self.assertEqual([contig.name for contig in reader.header.contigs],
                       ['chrM', 'chr1', 'chr2'])

  def test_iterate(self, filename, has_embedded_ref):
    with self._make_reader(filename, has_embedded_ref) as reader:
      reads = list(reader.iterate())
      self.assertEqual(len(reads), 3)
      self.assertEqual([read.fragment_name for read in reads],
                       ['cram1', 'cram2', 'cram3'])
      self.assertEqual([read.aligned_sequence for read in reads], [
          'CCCTAACCCTAACCCTAACCCTAACCCTANNNNNN',
          ('TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCAAAACGAATCAAAAAAGAAAAACGAA'
           'AAAAAAA'),
          'CACAGACGCTT'
      ])

  def test_query(self, filename, has_embedded_ref):
    with self._make_reader(filename, has_embedded_ref) as reader:
      for interval, n_expected in [('chr1:1-100', 3), ('chr2:1-121', 0)]:
        with reader.query(ranges.parse_literal(interval)) as iterable:
          self.assertEqual(test_utils.iterable_len(iterable), n_expected)


class ReadWriterTests(parameterized.TestCase):
  """Tests for sam.SamWriter."""

  def setUp(self):
    self.read1 = test_utils.make_read(
        bases='ACCGT',
        chrom='chr1',
        start=10,
        cigar='5M',
        mapq=50,
        quals=range(30, 35),
        name='read1')
    self.read2 = test_utils.make_read(
        bases='AACCTT',
        chrom='chr2',
        start=15,
        cigar='7M',
        mapq=40,
        quals=range(20, 26),
        name='read2')
    self.contigs = [
        reference_pb2.ContigInfo(name='chr1'),
        reference_pb2.ContigInfo(name='chr2'),
    ]
    self.header = reads_pb2.SamHeader()

  def test_make_read_writer_tfrecords(self):
    outfile = test_utils.test_tmpfile('test.tfrecord')
    writer = sam.SamWriter(outfile, header=self.header)

    # Test that the writer is a context manager and that we can write a read to
    # it.
    with writer:
      writer.write(self.read1)
      writer.write(self.read2)

    # Our output should have exactly one read in it.
    self.assertEqual([self.read1, self.read2],
                     list(
                         tfrecord.read_tfrecords(outfile,
                                                 proto=reads_pb2.Read)))

  @parameterized.parameters('test.bam', 'test.sam')
  def test_roundtrip_writer(self, filename):
    output_path = test_utils.test_tmpfile(filename)
    original_reader = sam.SamReader(test_utils.genomics_core_testdata(filename))
    original_records = list(original_reader.iterate())
    with sam.SamWriter(output_path, header=original_reader.header) as writer:
      for record in original_records:
        writer.write(record)
    with sam.SamReader(output_path) as new_reader:
      self.assertEqual(original_records, list(new_reader.iterate()))

  @parameterized.parameters(
      dict(
          filename='test_cram.embed_ref_0_version_3.0.cram',
          has_embedded_ref=False),
      dict(
          filename='test_cram.embed_ref_1_version_3.0.cram',
          has_embedded_ref=True))
  def test_roundtrip_cram_writer(self, filename, has_embedded_ref):
    output_path = test_utils.test_tmpfile(filename)
    writer_ref_path = test_utils.genomics_core_testdata('test.fasta')
    reader_ref_path = ''
    if not has_embedded_ref:
      reader_ref_path = writer_ref_path
    original_reader = sam.SamReader(
        test_utils.genomics_core_testdata(filename), ref_path=reader_ref_path)
    original_records = list(original_reader.iterate())
    with sam.SamWriter(
        output_path,
        header=original_reader.header,
        ref_path=writer_ref_path,
        embed_ref=has_embedded_ref) as writer:
      for record in original_records:
        writer.write(record)
    with sam.SamReader(output_path, ref_path=reader_ref_path) as new_reader:
      self.assertEqual(original_records, list(new_reader.iterate()))


if __name__ == '__main__':
  absltest.main()
