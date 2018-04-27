#!/bin/bash
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
#
# Builds OSS CLIF binary for DeepVariant.
#
# This script should be run on a cloud VM. Known to work on some versions of
# Linux OS.
#
# OSS CLIF takes a very long time to build (10+ minutes) since it needs to
# compile parts of clang and LLVM. To save this build time, we use this script
# to build CLIF, install it in /usr/local/clif, and then packages up
# /usr/local/clif and shared protobuf libraries from /usr/local/lib into a tgz
# called oss_clif.latest.tgz.
#
# This oss_clif.latest.tgz is used by build-prereq.sh to build DeepVariant.
# Various versions that we built and released can be found under:
# https://console.cloud.google.com/storage/browser/deepvariant/packages/oss_clif
#
# We do recognize that this should be temporary, and will update when there is
# an official solution from CLIF.
# GitHub issues such as https://github.com/google/deepvariant/issues/29 has
# some relevant pointers.

set -eux -o pipefail

# Figure out which linux installation we are on to fetch an appropriate version
# of CLIF binary. Note that we only support now Ubuntu (14 and 16), and Debian.
if [[ $(python -mplatform) == *"Ubuntu-16"* ]]; then
  export DV_PLATFORM="ubuntu-16"
  # For ubuntu 16 we install cmake
  sudo -H apt-get -y install cmake
elif [[ $(python -mplatform) == *"Ubuntu-14"* ]]; then
  export DV_PLATFORM="ubuntu-14"
  # For ubuntu 14 we install cmake3
  sudo -H apt-get -y install cmake3
elif [[ $(python -mplatform | grep '[Dd]ebian-\(rodete\|9.*\)') ]]; then
  export DV_PLATFORM="debian"
   # For recent debian, we install cmake.
   sudo -H apt-get -y install cmake
else
  export DV_PLATFORM="unknown"
  exit "unsupported platform"
fi

CLIF_DIR=/usr/local/clif
CLIF_PACKAGE="oss_clif.${DV_PLATFORM}.latest.tgz"

# Install prereqs.
sudo -H apt-get -y install ninja-build subversion
sudo -H apt-get -y install virtualenv python-pip pkg-config
sudo -H pip install 'pyparsing>=2.2.0'
sudo -H pip install 'protobuf>=3.4'

echo === building protobufs

sudo -H apt-get install -y autoconf automake libtool curl make g++ unzip
wget https://github.com/google/protobuf/releases/download/v3.4.1/protobuf-cpp-3.4.1.tar.gz
tar xvzf protobuf-cpp-3.4.1.tar.gz
(cd protobuf-3.4.1 &&
  ./autogen.sh &&
  ./configure &&
  make -j 32 &&
  make -j 32 check &&
  sudo make -j 32 install &&
  sudo ldconfig)

echo === building CLIF

git clone https://github.com/google/clif.git
sed -i 's/\$HOME\/opt/\/usr\/local/g' clif/INSTALL.sh
sed -i 's/-j 2//g' clif/INSTALL.sh
(cd clif && sudo ./INSTALL.sh)

echo === creating package tgz

sudo find ${CLIF_DIR} -type d -exec chmod a+rx {} \;
sudo find ${CLIF_DIR} -type f -exec chmod a+r {} \;
tar czf "${CLIF_PACKAGE}" /usr/local/lib/libproto* "${CLIF_DIR}"

echo === SUCCESS: package is "${CLIF_PACKAGE}"
