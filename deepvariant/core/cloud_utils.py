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
"""Utilities for interoperating with Google Cloud Platform."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from oauth2client import client as oauth2_client
import requests

from google.protobuf import json_format
from google.protobuf import message

# The value of this global is:
#  - None, if we haven't called oauth2_token() yet;
#  - "", if we did call oauth2_token() but cloud credentials were unavailable;
#  - a string encoding the OAuth2 token, otherwise.
_OAUTH_TOKEN = None


def oauth2_token():
  """Return the GCP OAuth2 token as a string.

  Returns:
    GCP OAuth2 token as a string, or "" if cloud credentials are unavailable.
  """
  global _OAUTH_TOKEN
  if _OAUTH_TOKEN is None:
    try:
      credentials = oauth2_client.GoogleCredentials.get_application_default()
      credentials.get_access_token()
      _OAUTH_TOKEN = credentials.access_token
    except oauth2_client.ApplicationDefaultCredentialsError:
      _OAUTH_TOKEN = ''
  return _OAUTH_TOKEN


# Note: "{{" is used to escape "{" characters in this string.
_LOG_JSON_TEMPLATE = """\
{{
  "logName": "projects/{project_name}/logs/{log_name}",
  "resource":  {{
    "type": "global",
    "labels": {{
      "project_id": "{project_name}" }}
  }},
  "entries": [
    {{
       {payload_kv}
    }}
  ],
}}
"""


def log_to_stackdriver(project_name, log_name, msg):
  """Log a text or proto message to a global StackDriver log.

  Args:
    project_name: the GCP project value.
    log_name: the name of the global StackDriver log.
    msg: a string or proto message to log.

  Raises:
    IOError: failure to log to StackDriver.
    ValueError: message has invalid type.
  """
  endpoint = 'https://logging.googleapis.com/v2/entries:write'
  headers = {
      'Authorization': 'Bearer {}'.format(oauth2_token()),
      'Content-Type': 'application/json'
  }
  if isinstance(msg, basestring):
    payload_kv = '"textPayload" : {}'.format(msg)
  elif isinstance(msg, message.Message):
    json_string = json_format.MessageToJson(
        msg,
        preserving_proto_field_name=True,
        including_default_value_fields=True)
    payload_kv = '"jsonPayload" : {}'.format(json_string)
  else:
    raise ValueError('`msg` argument is of invalid type.')

  rest_message = _LOG_JSON_TEMPLATE.format(
      project_name=project_name, log_name=log_name, payload_kv=payload_kv)

  response = requests.post(endpoint, headers=headers, data=rest_message)
  if not response.ok:
    raise IOError('Could not log to StackDriver, error is: %s' % response.text)
