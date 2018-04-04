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
"""Utilities to help with testing DeepVariant code."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os



from third_party.nucleus.testing import test_utils as nucleus_test_utils

DEEPVARIANT_DATADIR = ''


def deepvariant_testdata(filename):
  """Gets the path to filename in genomics/deepvariant/testdata.

  These paths are only known at runtime, after flag parsing
  has occurred.

  Args:
    filename: The name of a testdata file in the core genomics testdata
      directory. For example, if you have a test file in
      "learning/genomics/deepvariant/testdata/foo.txt", filename should be
      "foo.txt" to get a path to it.

  Returns:
    The absolute path to a testdata file.
  """
  return nucleus_test_utils.genomics_testdata(
      os.path.join('deepvariant/testdata', filename), DEEPVARIANT_DATADIR)


CHR20_FASTA = None
CHR20_BAM = None
GOLDEN_TRAINING_EXAMPLES = None
GOLDEN_CALLING_CANDIDATES = None
GOLDEN_CALLING_EXAMPLES = None
CONFIDENT_REGIONS_BED = None
TRUTH_VARIANTS_VCF = None
GOLDEN_POSTPROCESS_INPUT = None
GOLDEN_POSTPROCESS_OUTPUT = None
GOLDEN_POSTPROCESS_GVCF_INPUT = None
GOLDEN_POSTPROCESS_GVCF_OUTPUT = None

N_GOLDEN_TRAINING_EXAMPLES = 49
N_GOLDEN_CALLING_EXAMPLES = 82


def init():
  """Initialize global variables from flag values."""
  global CHR20_FASTA
  global CHR20_BAM
  global GOLDEN_TRAINING_EXAMPLES
  global GOLDEN_CALLING_CANDIDATES
  global GOLDEN_CALLING_EXAMPLES
  global CONFIDENT_REGIONS_BED
  global TRUTH_VARIANTS_VCF
  global GOLDEN_POSTPROCESS_INPUT
  global GOLDEN_POSTPROCESS_OUTPUT
  global GOLDEN_POSTPROCESS_GVCF_INPUT
  global GOLDEN_POSTPROCESS_GVCF_OUTPUT

  CHR20_FASTA = deepvariant_testdata('ucsc.hg19.chr20.unittest.fasta.gz')
  CHR20_BAM = deepvariant_testdata('NA12878_S1.chr20.10_10p1mb.bam')
  GOLDEN_TRAINING_EXAMPLES = deepvariant_testdata(
      'golden.training_examples.tfrecord')
  GOLDEN_CALLING_CANDIDATES = deepvariant_testdata(
      'golden.calling_examples.tfrecord')
  GOLDEN_CALLING_EXAMPLES = deepvariant_testdata(
      'golden.calling_examples.tfrecord')
  CONFIDENT_REGIONS_BED = deepvariant_testdata(
      'test_nist.b37_chr20_100kbp_at_10mb.bed')
  TRUTH_VARIANTS_VCF = deepvariant_testdata(
      'test_nist.b37_chr20_100kbp_at_10mb.vcf.gz')
  GOLDEN_POSTPROCESS_INPUT = deepvariant_testdata(
      'golden.postprocess_single_site_input.tfrecord')
  GOLDEN_POSTPROCESS_OUTPUT = deepvariant_testdata(
      'golden.postprocess_single_site_output.vcf')
  GOLDEN_POSTPROCESS_GVCF_INPUT = deepvariant_testdata(
      'golden.postprocess_gvcf_input.tfrecord')
  GOLDEN_POSTPROCESS_GVCF_OUTPUT = deepvariant_testdata(
      'golden.postprocess_gvcf_output.g.vcf')
