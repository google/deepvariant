# Copyright 2026 Google LLC.
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

import os
import shutil

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant import bagz_index
from deepvariant import testdata
from absl import app
import multiprocessing


def setUpModule():
  testdata.init()


class BagzIndexTest(parameterized.TestCase):

  def test_bagz_index(self):
    temp_dir = self.create_tempdir().full_path
    original_bagz_dir = os.path.dirname(testdata.GOLDEN_TRAINING_EXAMPLES_BAGZ)
    shutil.copytree(original_bagz_dir, temp_dir, dirs_exist_ok=True)
    os.chmod(temp_dir, 0o775)

    # Run the indexer.
    bagz_index.main(['bagz_index', temp_dir])

    # Check that the index files were created.
    self.assertTrue(os.path.exists(os.path.join(temp_dir, 'index.tsv.gz')))
    self.assertTrue(os.path.exists(os.path.join(temp_dir, 'lazy_bag_info.pb')))


if __name__ == '__main__':
  g3_multiprocessing.handle_test_main(absltest.main)
