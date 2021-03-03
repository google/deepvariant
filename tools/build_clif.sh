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
#
#
# This script is only maintained for Ubuntu 18.04.

set -eux -o pipefail

echo ========== This script is only maintained for Ubuntu 18.04.
echo ========== See https://github.com/google/clif for how to build on different Unix distributions.

UBUNTU_VERSION=18.04
ABSL_VERSION=20200923
PROTOBUF_VERSION=3.13.0
PYTHON_VERSION=3.6

APT_ARGS=(
"-qq"
"-y"
)


sudo apt-get update  "${APT_ARGS[@]}"
sudo apt-get install "${APT_ARGS[@]}" --no-install-recommends \
    autoconf \
    automake \
    cmake \
    curl \
    gpg-agent \
    g++ \
    libtool \
    make \
    pkg-config \
    software-properties-common \
    wget \
    unzip

# Configure LLVM 11 apt repository
wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | sudo apt-key add - && \
 sudo add-apt-repository "deb http://apt.llvm.org/$(lsb_release -sc)/ llvm-toolchain-$(lsb_release -sc)-11 main"

# Install CLIF dependencies
sudo apt-get update "${APT_ARGS[@]}"
sudo apt-get install  "${APT_ARGS[@]}" \
    clang-11 \
    libclang-11-dev \
    libgoogle-glog-dev \
    libgtest-dev \
    libllvm11 \
    llvm-11-dev \
    python3-dev \
    python3-pyparsing \
    zlib1g-dev

# Configure deadsnakes PPA with the more recent versions of python packaged for
# Ubuntu. See https://launchpad.net/~deadsnakes/+archive/ubuntu/ppa
sudo apt-get update "${APT_ARGS[@]}" && \
    sudo apt-get install "${APT_ARGS[@]}" \
    "python$PYTHON_VERSION-dev" \
    "python$PYTHON_VERSION-distutils"

# Install latest version of pip since the version on ubuntu could be outdated
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    "python$PYTHON_VERSION" get-pip.py && \
    rm get-pip.py

# Compile and install absl-cpp from source
wget "https://github.com/abseil/abseil-cpp/archive/$ABSL_VERSION.tar.gz" && \
    tar -xf "$ABSL_VERSION.tar.gz" && \
    mkdir "abseil-cpp-$ABSL_VERSION/build" && \
    cd "abseil-cpp-$ABSL_VERSION/build" && \
    cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE=true && \
    sudo make install && \
    rm -rf "/abseil-cpp-$ABSL_VERSION" "/$ABSL_VERSION.tar.gz"

# Compile and install protobuf from source
wget "https://github.com/protocolbuffers/protobuf/releases/download/v$PROTOBUF_VERSION/protobuf-cpp-$PROTOBUF_VERSION.tar.gz" && \
    tar -xf "protobuf-cpp-$PROTOBUF_VERSION.tar.gz" && \
    cd "protobuf-$PROTOBUF_VERSION" && \
    # Configure and install C++ libraries
    ./autogen.sh && \
    ./configure && \
    make -j"$(nproc)" && \
    sudo make install && \
    sudo ldconfig && \
    rm -rf "/protobuf-$PROTOBUF_VERSION" "/protobuf-cpp-$PROTOBUF_VERSION.tar.gz"

# Install googletest
cd /usr/src/gtest && \
    sudo cmake . && \
    sudo make install

# Install python runtime and test dependencies
"python$PYTHON_VERSION" -m pip install \
    absl-py \
    parameterized \
    protobuf=="$PROTOBUF_VERSION"

sudo "python$PYTHON_VERSION" -m pip uninstall -y pyparsing && \
  "python$PYTHON_VERSION" -m pip install -Iv 'pyparsing==2.2.0'

DV_PLATFORM="ubuntu-${UBUNTU_VERSION}"

sudo ln -sf /usr/bin/python\$PYTHON_VERSION /usr/local/bin/python3

cd && rm -rf clif && git clone https://github.com/google/clif.git && \
  cd clif && \
  sudo ./INSTALL.sh
