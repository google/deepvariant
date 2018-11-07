#!/usr/bin/python
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

"""Adds a dynamic library load to the top of a Python program.

The only tricky part about this is that Python programs aren't allowed to
have any code before `from __future__` imports, so we have to insert our
load after all of those.

Also, this script is used by nucleus.bzl, so it can't be compiled or
otherwise be the output of a bazel build.

Usage:
  add_cdll_load.py input_program.py some_library.so output_program.py
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import sys

if len(sys.argv) != 4:
  print('Usage: %s <input_program> <so_file> <output_program>' % sys.argv[0])
  sys.exit(-1)

in_py = sys.argv[1]
so_file = sys.argv[2]
out_py = sys.argv[3]

load_str = """

import ctypes
import os
ctypes.CDLL(os.path.dirname(os.path.realpath(__file__)) + \"/%s\",
            ctypes.RTLD_GLOBAL)

""" % so_file

with open(in_py, 'r') as in_file:
  in_str = in_file.read()
  last_future_re = re.compile(r'(.*)(^from __future__ import .+?$)(.*)',
                              re.MULTILINE | re.DOTALL)
  m = re.match(last_future_re, in_str)
  if m is not None:
    out_str = m.group(1) + m.group(2) + load_str + m.group(3)
  else:
    out_str = load_str + in_str

  with open(out_py, 'w') as out_file:
    out_file.write(out_str)
