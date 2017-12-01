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
"""Smoke test for authenticated GCS file I/O."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

from absl.testing import absltest

from deepvariant.core import genomics_io
from deepvariant.core import ranges

# At present we  can only test accessing *public* files on GCS as blaze
# sanitizes the environment, removing our access to gcloud credentials.
REMOTE_BAM = 'gs://genomics-public-data/platinum-genomes/bam/NA12877_S1.bam'
REMOTE_FASTA = 'gs://genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa'
QUERY_WINDOW = 'chr20:10000001-10000100'
EXPECTED_REFERENCE_SEQUENCE = 'AACAAGTTCCAGAAGATAGCTAGAGGATGGGAGCACATGAAGAGCAGATCACAACCATCCCTGGGGAGCCCAGCCTGGACCAGCTACAGCCACCCGTCTC'
EXPECTED_READS_IN_WINDOW = 95


# cd to a temporary directory so that htslib can freely write index files.
def setUpModule():
  os.chdir(absltest.get_default_test_tmpdir())


class GenomicsIoGCSTest(absltest.TestCase):

  def setUp(self):
    self.query_window = ranges.parse_literal(QUERY_WINDOW)

  def test_remote_fasta(self):
    reader = genomics_io.make_ref_reader(REMOTE_FASTA)
    reference_segment = reader.bases(self.query_window)
    self.assertEqual(EXPECTED_REFERENCE_SEQUENCE, reference_segment)

  def test_remote_bam(self):
    reader = genomics_io.make_sam_reader(REMOTE_BAM)
    reads = list(reader.query(self.query_window))
    self.assertEqual(EXPECTED_READS_IN_WINDOW, len(reads))


if __name__ == '__main__':
  absltest.main()
