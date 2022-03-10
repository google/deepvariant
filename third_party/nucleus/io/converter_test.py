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
"""Tests for third_party.nucleus.examples.convert_genomics_file.

These tests do NOT establish the correctness of conversions---tests of the
fidelity of the Reader and Writer classes exist elsewhere in Nucleus.  Rather,
these tests simply exercise that the conversion *runs* for each input/output
file type.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


import os
import unittest

from absl import logging
from absl.testing import absltest
from absl.testing import parameterized
from third_party.nucleus.io import converter
from third_party.nucleus.testing import test_utils

basename = os.path.basename

# Initial (native) input files we will use to begin conversions.
ORIGINAL_TEST_FILES = [
    "test.bed", "test_reads.fastq", "test_features.gff", "test.sam",
    "test_sites.vcf"
]

# These formats require a header, so conversion from tfrecord to a native file
# format cannot be done faithfully.
FORMATS_REQUIRING_HEADER = [".bam", ".gff", ".sam", ".vcf"]


class ConvertGenomicsFileTest(parameterized.TestCase):

  def _convert(self, src, dest):
    logging.info("#### Converting: %s --> %s ... ", basename(src),
                 basename(dest))
    converter.convert(src, dest)

  @parameterized.parameters(*ORIGINAL_TEST_FILES)
  def test_conversion_to_tfrecord_and_back(self, original_input_file):
    """Test conversion from a native file format to tfrecord.gz, then back."""
    input_path = test_utils.genomics_core_testdata(original_input_file)
    tfrecord_output_path = test_utils.test_tmpfile(original_input_file +
                                                   ".tfrecord.gz")
    native_output_path = test_utils.test_tmpfile(original_input_file)

    # Test conversion from native format to tfrecord.
    self._convert(input_path, tfrecord_output_path)

    # TODO: remove this when SAM writer is implemented.
    if native_output_path.endswith(".sam"):
      raise unittest.SkipTest("SAM writing not yet supported")

    # Test conversion from tfrecord format back to native format.  Ensure that
    # conversions where we would need a header, but don't have one from the
    # input, trigger an error message.
    if any(
        native_output_path.endswith(ext) for ext in FORMATS_REQUIRING_HEADER):
      with self.assertRaisesRegexp(
          converter.ConversionError,
          "Input file does not have a header, which is needed to construct "
          "output file"):
        self._convert(tfrecord_output_path, native_output_path)

    else:
      self._convert(tfrecord_output_path, native_output_path)


if __name__ == "__main__":
  absltest.main()
 
