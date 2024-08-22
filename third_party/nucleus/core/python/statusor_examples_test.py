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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from absl.testing import absltest
from third_party.nucleus.core.python import statusor_examples

USING_PYBIND = hasattr(statusor_examples, 'USING_PYBIND')


class StatusorClifWrapTest(absltest.TestCase):

  def test_make_int_ok(self):
    self.assertEqual(statusor_examples.MakeIntOK(), 42)

  def test_make_int_fail(self):
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(ValueError, 'MakeIntFail'):
      statusor_examples.MakeIntFail()

  def test_make_str_ok(self):
    self.assertEqual(statusor_examples.MakeStrOK(), 'hello')

  # See CLIF wrapper for a discussion of why this is commented out.
  # def test_make_str_ok_stripped_type(self):
  #   self.assertEqual(statusor_examples.MakeStrOKStrippedType(), 'hello')

  def test_make_str_fail(self):
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(ValueError, 'MakeStrFail'):
      statusor_examples.MakeStrFail()

  @absltest.skipIf(USING_PYBIND, 'Disabled for now.')
  def test_make_int_unique_ptr_ok(self):
    self.assertEqual(statusor_examples.MakeIntUniquePtrOK(), 421)

  @absltest.skipIf(USING_PYBIND, 'Disabled for now.')
  def test_make_int_unique_ptr_fail(self):
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(ValueError, 'MakeIntUniquePtrFail'):
      statusor_examples.MakeIntUniquePtrFail()

  @absltest.skipIf(USING_PYBIND, 'Disabled for now.')
  def test_make_int_vector_ok(self):
    self.assertEqual(statusor_examples.MakeIntVectorOK(), [1, 2, 42])

  @absltest.skipIf(USING_PYBIND, 'Disabled for now.')
  def test_make_int_vector_fail(self):
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(ValueError, 'MakeIntVectorFail'):
      statusor_examples.MakeIntVectorFail()

  def test_returning_status_ok_returns_none(self):
    self.assertEqual(statusor_examples.FuncReturningStatusOK(), None)

  def test_returning_status_fail_raises(self):
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(ValueError, 'FuncReturningStatusFail'):
      statusor_examples.FuncReturningStatusFail()

  def test_string_owner(self):
    obj = statusor_examples.StringOwner.Factory()
    self.assertEqual(obj.GetText(), 'Factory')


if __name__ == '__main__':
  absltest.main()
