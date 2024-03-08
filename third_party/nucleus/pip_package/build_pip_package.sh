#!/bin/bash
# Copyright 2020 Google LLC.
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
# before running this script.  In addition, if you make any changes to the
# source code, you should rebuild:
#   bazel build -c opt $BAZEL_FLAGS nucleus/pip_package:build_pip_package

set -e
set -x

# When changing NUCLEUS_VERSION, be sure to also change it in
# egg_files/PKG-INFO.
NUCLEUS_VERSION="0.7.0"
PACKAGE_NAME="google_nucleus-${NUCLEUS_VERSION}"
PYTHON_VERSION="3.10"

TMPDIR=$(mktemp -d -t tmp.XXXXXXXXXXX)
TOPDIR="${TMPDIR}/${PACKAGE_NAME}"
mkdir -p "${TOPDIR}"

echo $(date) : "=== Copying files to ${TOPDIR}"

RUNFILES=bazel-bin/nucleus/pip_package/build_pip_package.runfiles/nucleus

# $RUNFILES has three subdirectories, each of which gets treated a bit
# differently.

# Subdirectory #1:  Copy /nucleus to top level.
cp -L -R "${RUNFILES}/nucleus" "${TOPDIR}"

# Subdirectory #2:  Copy /third_party to /nucleus/third_party.
mkdir -p "${TOPDIR}/nucleus/third_party"
cp -L -R "${RUNFILES}"/third_party/* "${TOPDIR}/nucleus/third_party"

# Subdirectory #3:  /external.  The only thing we need from it is our
# version of protobuf, which we need to import as google.protobuf.
# See also the top level __init__.py that sets the sys.path to make this
# version of protobuf have precedence over any other installed versions while
# running Nucleus code.
cp -L -R "${RUNFILES}/external/com_google_protobuf/python/google" "${TOPDIR}/google"

# Copy top level files to /nucleus.
cp LICENSE "${TOPDIR}/nucleus"

# Copy setup.py to top level.
cp "${RUNFILES}/nucleus/pip_package/setup.py" "${TOPDIR}"
cp "${RUNFILES}/nucleus/pip_package/setup.cfg" "${TOPDIR}"

# Create egg-info directory.
EGG_DIR="${TOPDIR}/${PACKAGE_NAME}-py${PYTHON_VERSION}.egg-info"
mkdir -p "${EGG_DIR}"
cp "${RUNFILES}/nucleus/pip_package/egg_files"/* "${EGG_DIR}"
pushd "${TOPDIR}"
find . -type f -print > "${EGG_DIR}/SOURCES.txt"
popd

# Fix symbolic links -- any .so file in Nucleus should point to
# google/protobuf/pyext/_message.so with a relative link.
pushd "${TOPDIR}"
find "nucleus" -name '*.so' -exec ln -f -s -r "google/protobuf/pyext/_message.so" {} \;
popd

# Some versions of protobuf have _message.so files named
# _message.cpython-34m.so or _message.cpython-35m-x86_64-linux-gnu.so
# so we create a symbolic link at those filenames so that we overwrite them.
pushd "${TOPDIR}/google/protobuf/pyext"
ln -f -s "_message.so" "_message.cpython-34m.so"
ln -f -s "_message.so" "_message.cpython-35m-x86_64-linux-gnu.so"
popd

# Create tar file
TAR_NAME="${PACKAGE_NAME}.tar.gz"
pushd "${TMPDIR}"
echo $(date) : "=== Building tar file ${TAR_NAME}"
tar cvzf "${TAR_NAME}" "${PACKAGE_NAME}"

# ls the tarfile to see how large it is.
# TODO: Set up monitoring for the pip package size.  It isn't
# allowed to be over 100M, and it's usually a sign of other trouble when
# it is over that size.
ls -lh "${TAR_NAME}"

if [ $# -gt 0 ]; then
  DEST=$1
  mkdir -p "${DEST}"
  cp "${TAR_NAME}" "${DEST}"
else
  DEST="${TAR_NAME}"
fi

popd

echo "Output tar file is in ${DEST}"
