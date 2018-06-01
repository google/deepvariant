#!/bin/bash
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

# Usage:  ./nucleus/pip_package/build_pip_package.sh [optional_dir]
#
# If [optional_dir] is supplied, the created wheel file is placed there.
#
# Important:  You must run
#   source install.sh
# before running this script.

set -e

bazel build -c opt $COPT_FLAGS nucleus/...

function cp_external() {
  local src_dir=$1
  local dest_dir=$2
  for f in `find "$src_dir" -maxdepth 1 -mindepth 1 ! -name '*nucleus*'`; do
    cp -R "$f" "$dest_dir"
  done
}

TMPDIR=$(mktemp -d -t tmp.XXXXXXXXXXX)

RUNFILES=bazel-bin/nucleus/pip_package/build_pip_package.runfiles/nucleus

# $RUNFILES has four subdirectories, each of which gets treated a bit
# differently.

# Subdirectory #1:  Copy /nucleus to top level.
cp -R "${RUNFILES}/nucleus" "${TMPDIR}"

# Subdirectory #2:  Copy /external to /external, but leave out
# /external/nucleus.
mkdir "${TMPDIR}/external"
cp_external "${RUNFILES}" "${TMPDIR}/external"

# Subdirectory #3: Copy /_solib_k8 (or whatever the binary files directory
# is called) to top level.
so_lib_dir=$(ls "$RUNFILES" | grep solib)
if [ -n "${so_lib_dir}" ]; then
  cp -R "${RUNFILES}/${so_lib_dir}" "${TMPDIR}"
fi

# Subdirectory #4: Copy /third_party to /third_party.
mkdir "${TMPDIR}/third_party"
cp -R "${RUNFILES}"/third_party/* "${TMPDIR}/third_party"

cp LICENSE "${TMPDIR}"
cp README.md "${TMPDIR}"
cp nucleus/pip_package/MANIFEST.in "${TMPDIR}"
cp nucleus/pip_package/setup.py "${TMPDIR}"

pushd "${TMPDIR}"
rm -f MANIFEST
echo $(date) : "=== Building wheel in ${TMPDIR}"
python setup.py bdist_wheel
popd

if [ $# -gt 0 ]; then
  DEST=$1
  mkdir -p "${DEST}"
  cp "${TMPDIR}/dist"/* "${DEST}"
else
  DEST="${TMPDIR}/dist"
fi

echo "Output wheel is in ${DEST}"
