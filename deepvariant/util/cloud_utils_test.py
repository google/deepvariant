# Copyright 2018 Google Inc.
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

"""Tests for deepvariant.util.cloud_utils."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import httplib
import socket


from absl.testing import absltest
import mock

from deepvariant.util import cloud_utils


class CloudUtilsTest(absltest.TestCase):
  """Trivial tests of functionality for fetching the OAuth token.

  These tests currently only test that we cleanly fail to to get the OAuth
  in the event that there is no network access.  We don't test that access is
  definitively available, because this is variable across our testing
  environments.
  """

  def test_with_disabled_network(self):
    # Disable network access, at least via socket.socket.
    with mock.patch('socket.socket'):
      socket.socket.side_effect = httplib.HTTPException(
          'Network disabled by test.')
      token = cloud_utils.oauth2_token()
      self.assertEqual('', token)


if __name__ == '__main__':
  absltest.main()
