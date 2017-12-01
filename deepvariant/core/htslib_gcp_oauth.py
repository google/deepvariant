# Copyright 2017 Google Inc.
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
"""Utility module for making GCP OAuth2 credentials available to htslib.

htslib looks for an OAuth2 token to be available in the env var
GCS_OAUTH_TOKEN.  This module provides helper functions to get the
token and put it there.

This module can be used in two ways:

  1.  The htslib_gcp_oauth2.init() function can be called at your
      program's startup.

  2.  You can run this module as a program, passing your command as
      arguments, and it will set up the environment variable and then
      invoke your command.  Example:

      $ htslib_gcp_oauth.py -- my_exe my_arg1 my_arg2

Note: we do not split hairs here about what scopes are assigned to the
OAuth2 token.  If the token does not have scope access to Google Cloud
Storage, attempt to access gs:// files will fail at file access time
with a permissions message.
"""

# redacted
# front here--but I can't figure out how to use the
# credentials.retrieve_scopes() function.  See:
# https://github.com/google/oauth2client/issues/365

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import logging
import os


import tensorflow as tf

from deepvariant.core import cloud_utils


def init():
  """Gets the GCP OAuth2 token and places it in the ENV var GCS_OAUTH_TOKEN.

  Logs a message about whether or not we succeeded in finding and installing the
  credentials appropriately.

  Returns:
    True if auth was set up, or False otherwise.
  """
  token = cloud_utils.oauth2_token()
  if token:
    logging.info(
        'GCP credentials found; will be able to access non-public gs:// '
        'URIs from htslib')
    os.putenv('GCS_OAUTH_TOKEN', token)
    return True
  else:
    logging.warn(
        'GCP credentials not found; only local files and public gs:// URIs will'
        ' be accessible from htslib')
    return False


def run_command_with_auth(cmd_args):
  cmd = ' '.join(cmd_args)
  logging.info('Running cmd with GCP authentication: %s', cmd)
  init()
  os.system(cmd)


def main(argv):
  if len(argv) < 2:
    raise ValueError('Expecting command as argument.')
  run_command_with_auth(argv[1:])


if __name__ == '__main__':
  tf.app.run()
