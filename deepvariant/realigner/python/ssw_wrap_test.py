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



from absl.testing import absltest
from deepvariant.realigner.python import ssw

REF = 'CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA'
QUERY = 'CTGAGCCGGTAAATC'

# Expected alignment:
#   CAGCCTTTCTGACCCGG-AAATCAAAATAGGCACAACAAA
#           ||||*|||| |||||
#           CTGAGCCGGTAAATC


class SswWrapTest(absltest.TestCase):
  """Basic tests for the wrapped SSW aligner."""

  def test_Align(self):
    """Tests the Align method."""
    aligner = ssw.Aligner()
    filter_ = ssw.Filter()
    length = aligner.set_reference_sequence(REF)
    self.assertLen(REF, length)
    alignment = aligner.align(QUERY, filter_, 16)
    self.assertEqual(21, alignment.sw_score)
    self.assertEqual(8, alignment.sw_score_next_best)
    self.assertEqual(8, alignment.ref_begin)
    self.assertEqual(21, alignment.ref_end)
    self.assertEqual(0, alignment.query_begin)
    self.assertEqual(14, alignment.query_end)
    self.assertEqual(4, alignment.ref_end_next_best)
    self.assertEqual(2, alignment.mismatches)
    self.assertEqual(b'4=1X4=1I5=', alignment.cigar_string)

  def test_Align2_reversed(self):
    """Tests the Align method, reversing query and ref from above."""
    aligner = ssw.Aligner()
    filter_ = ssw.Filter()
    aligner.set_reference_sequence(QUERY)
    alignment = aligner.align(REF, filter_, 16)
    self.assertEqual(21, alignment.sw_score)
    self.assertEqual(8, alignment.query_begin)
    self.assertEqual(21, alignment.query_end)
    self.assertEqual(0, alignment.ref_begin)
    self.assertEqual(14, alignment.ref_end)
    self.assertEqual(2, alignment.mismatches)
    self.assertEqual(b'8S4=1X4=1D5=17S', alignment.cigar_string)


if __name__ == '__main__':
  absltest.main()
