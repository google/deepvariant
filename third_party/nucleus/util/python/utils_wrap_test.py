# Copyright 2024 Google LLC.
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

from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.util.python import utils as cpp_utils


class UtilsTest(absltest.TestCase):

  def test_read_end_minimal(self):
    read_pb = reads_pb2.Read()
    self.assertEqual(cpp_utils.read_end(read_pb), 0)

  def test_read_range_minimal(self):
    read_pb = reads_pb2.Read()
    range_pb = range_pb2.Range()
    self.assertIsNone(cpp_utils.read_range(read_pb, range_pb))

  def test_read_overlaps_region_minimal(self):
    read_pb = reads_pb2.Read()
    range_pb = range_pb2.Range()
    self.assertFalse(cpp_utils.read_overlaps_region(read_pb, range_pb))


if __name__ == '__main__':
  absltest.main()
