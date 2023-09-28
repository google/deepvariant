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
"""Tests that call_variants can run in a environment with an accelerator."""



from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant import call_variants_slim
from deepvariant import modeling
from deepvariant import testdata
from third_party.nucleus.testing import test_utils


def setUpModule():
  testdata.init()


class CallVariantsAcceleratorTests(
    tf.test.TestCase, metaclass=parameterized.TestGeneratorMetaclass
):

  @parameterized.parameters(modeling.production_models())
  def test_call_variants_runs_on_gpus(self, model):
    call_variants_slim.call_variants(
        examples_filename=testdata.GOLDEN_CALLING_EXAMPLES,
        checkpoint_path=None,
        model=model,
        execution_hardware='accelerator',
        output_file=test_utils.test_tmpfile('zzz.tfrecord'),
    )


if __name__ == '__main__':
  absltest.main()
