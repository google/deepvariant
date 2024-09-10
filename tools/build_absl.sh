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

set -eux -o pipefail

echo ========== This script has been tested on Ubuntu22.04.
echo ========== Run this script in root mode.

ABSL_PIN="${ABSL_PIN-29bf8085f3bf17b84d30e34b3d7ff8248fda404e}"

APT_ARGS=(
"-y"
)


apt-get update  "${APT_ARGS[@]}"
NEEDRESTART_MODE=a apt-get install "${APT_ARGS[@]}" --no-install-recommends \
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

# Install dependencies
apt-get update "${APT_ARGS[@]}"
NEEDRESTART_MODE=a apt-get install "${APT_ARGS[@]}" \
    clang-11 \
    libclang-11-dev \
    libgoogle-glog-dev \
    libgtest-dev \
    libllvm11 \
    llvm-11 \
    llvm-11-dev \
    llvm-11-linker-tools \
    python3-dev \
    zlib1g-dev

# Compile and install absl-cpp from source
git clone https://github.com/abseil/abseil-cpp.git
cd abseil-cpp
if [[ ! -z ${ABSL_PIN} ]]; then
  git checkout "${ABSL_PIN}"
fi
mkdir build && cd build
cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE=true
make install
cd ../..
rm -rf abseil-cpp

curl -o get-pip.py https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --force-reinstall --user
rm -f get-pip.py

export PATH="$HOME/.local/bin":$PATH
echo "$(pip3 --version)"

# Install python runtime and test dependencies
pip3 install \
    absl-py \
    parameterized

# On GPU machines, this might be necessary because of the reason mentioned in:
# https://stackoverflow.com/a/74605488
NEEDRESTART_MODE=a apt-get install "${APT_ARGS[@]}" libstdc++-12-dev
