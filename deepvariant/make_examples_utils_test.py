# Copyright 2020 Google LLC.
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
"""Tests for deepvariant .make_examples_utils."""

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import make_examples_utils


class MakeExamplesUtilsTest(parameterized.TestCase):

  def test_init_validation(self):
    with self.assertRaisesRegex(ValueError, 'integer'):
      make_examples_utils.Sample(pileup_height='100')

  def test_print(self):
    sample = make_examples_utils.Sample(
        name='test',
        sam_readers='my sam reader',
        in_memory_sam_reader='my in_memory_sam_reader',
        pileup_height=200)
    # print(sample) should yield:
    # <Sample {'sam_readers': 'my sam reader',
    #     'in_memory_sam_reader': 'my in_memory_sam_reader',
    #     'name': 'test',
    #     'pileup_height': 200}>
    self.assertIn('test', sample.__repr__())
    self.assertIn('my sam reader', sample.__repr__())
    self.assertIn('my in_memory_sam_reader', sample.__repr__())
    self.assertIn('pileup_height', sample.__repr__())
    self.assertIn('200', sample.__repr__())


if __name__ == '__main__':
  absltest.main()
