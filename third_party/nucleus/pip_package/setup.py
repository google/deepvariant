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
"""Setup module for turning Nucleus into a pip package.

Based on
https://github.com/pypa/sampleproject/blob/master/setup.py

This should be invoked through build_pip_package.sh, rather than run
directly.
"""

import fnmatch
import os

from setuptools import find_packages
from setuptools import setup

def find_files(pattern, root):
  """Return all the files matching pattern below root dir."""
  for dirpath, _, files in os.walk(root):
    for filename in fnmatch.filter(files, pattern):
      yield os.path.join(dirpath, filename)


def is_python_file(fn):
  return fn.endswith('.py') or fn.endswith('.pyc')


headers = list(find_files('*.h', 'nucleus'))

matches = ['../' + x for x in find_files('*', 'external')
           if not is_python_file(x)]

so_lib_paths = [
    i for i in os.listdir('.')
    if os.path.isdir(i) and fnmatch.fnmatch(i, '_solib_*')
]

for path in so_lib_paths:
  matches.extend(
      ['../' + x for x in find_files('*', path) if '.py' not in x]
  )

setup(
    name='Nucleus',
    version='0.1.0',
    description='A library for reading and writing genomics data.',
    long_description=
"""
Nucleus is a library of Python and C++ code designed to make it easy to
read, write and analyze data in common genomics file formats like SAM and VCF.
In addition, Nucleus enables painless integration with the TensorFlow machine
learning framework, as anywhere a genomics file is consumed or produced, a
TensorFlow tfrecords file may be substituted.
""",
    url='https://github.com/google/nucleus',
    author='The Genomics team in Google Brain',
    author_email='opensource@google.com',
    license='Apache 2.0',

    # Taken from list of valid classifiers at
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],

    keywords='genomics tensorflow bioinformatics',

    packages=find_packages(exclude=['examples', 'g3doc', 'testdata']),

    # redacted
    # these install_requires.
    install_requires=['contextlib2', 'intervaltree', 'absl-py',
                      'mock', 'numpy', 'scipy', 'six',
                      'tensorflow>=1.7.0'],

    headers=headers,

    include_package_data=True,
    package_data={'nucleus': matches},

    data_files=[],

    entry_points={},

    project_urls={
        'Source': 'https://github.com/google/nucleus',
        'Bug Reports': 'https://github.com/google/nucleus/issues',
    },

    zip_safe=False,
)


