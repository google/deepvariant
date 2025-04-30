# Copyright 2025 Google LLC.
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

"""Defines a custom build command that generates Python files from .proto files."""

import glob
import os
import shutil
import subprocess
import sys

from setuptools import Command
from setuptools import find_packages
from setuptools import setup

if 'PROTOC' in os.environ and os.path.exists(os.environ['PROTOC']):
  protoc = os.environ['PROTOC']
else:
  protoc = shutil.which('protoc')


def generate_proto(source):
  """Generates Python files from .proto files."""
  if not protoc or not os.path.exists(source):
    return
  if os.path.isdir(source):
    files = glob.glob(os.path.join(source, '*.proto'))
  else:
    files = [source]

  for f in files:
    protoc_command = [protoc, '-I.', '--python_out=.', f]
    print(protoc_command)
    if subprocess.call(protoc_command) != 0:
      sys.exit(-1)


# Define a custom build command that generates Python files from .proto files
class BuildProtoCommand(Command):
  user_options = []

  def initialize_options(self):
    pass

  def finalize_options(self):
    pass

  def run(self):
    generate_proto('deepvariant/protos/')
    generate_proto('third_party/nucleus/protos/')


setup(
    name='deepvariant',
    version='1.9.0',
    packages=find_packages(),
    package_dir={'deepvariant': 'deepvariant', 'third_party': 'third_party'},
    install_requires=[],
    cmdclass={
        'build_proto': BuildProtoCommand,
    },
)
