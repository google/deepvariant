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


def p(obj):
  for x in dir(obj):
    if not x.startswith('_'):
      print(x + ':' + repr(getattr(obj, x, '')))


class SswGccTest(absltest.TestCase):
  """Tests for the wrapped SSW aligner in a way that fails with gcc5.4."""

  def test_short(self):
    """Test very short strings."""
    ref = 'tttt'
    query = 'ttAtt'
    match = 4
    mismatch = 2
    gap_extend_penalty = 2
    gap_open_penalty = 4

    aligner = ssw.Aligner.construct(
        match_score=match,
        mismatch_penalty=mismatch,
        gap_opening_penalty=gap_open_penalty,
        gap_extending_penalty=gap_extend_penalty)
    filter_ = ssw.Filter()
    length = aligner.set_reference_sequence(ref)
    self.assertLen(ref, length)
    alignment = aligner.align(query, filter_, 16)
    p(alignment)
    self.assertEqual(b'2=1I2=', alignment.cigar_string)

  def test_longer(self):
    """Test longer strings, so the second-best alignment is considered."""
    ref = 'TTTTGGGGGGGGGGGGG'
    query = 'TTATTGGGGGGGGGGGGG'
    match = 4
    mismatch = 2
    gap_extend_penalty = 2
    gap_open_penalty = 4

    aligner = ssw.Aligner.construct(
        match_score=match,
        mismatch_penalty=mismatch,
        gap_opening_penalty=gap_open_penalty,
        gap_extending_penalty=gap_extend_penalty)
    filter_ = ssw.Filter()
    length = aligner.set_reference_sequence(ref)
    self.assertLen(ref, length)
    alignment = aligner.align(query, filter_, 16)
    p(alignment)
    self.assertEqual(b'2=1I15=', alignment.cigar_string)


if __name__ == '__main__':
  absltest.main()
