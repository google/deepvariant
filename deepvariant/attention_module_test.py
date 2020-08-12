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
"""Tests for deepvariant .attention_module."""

from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant import attention_module
from deepvariant import dv_constants


class AttentionModuleTest(parameterized.TestCase):
  """Test to check functions from attention_module.py."""

  @parameterized.parameters([
      (4096, dv_constants.PILEUP_DEFAULT_HEIGHT,
       dv_constants.PILEUP_DEFAULT_WIDTH, dv_constants.PILEUP_NUM_CHANNELS),
      (4096, 8, 8, 2048),
  ])
  def test_se_block(self, batch_size, height, width, channels):
    """Test to check the shape of output tensor after passing SE block."""

    input_shape = (batch_size, height, width, channels)
    input_feature = tf.zeros(shape=input_shape)
    output = attention_module.se_block(
        input_feature=input_feature, name='se_block', ratio=8)

    self.assertEqual(input_shape, output.get_shape())


if __name__ == '__main__':
  absltest.main()
