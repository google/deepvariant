#!/bin/bash
set -euo pipefail

# Copyright 2017 Google LLC.
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

echo ========== This script is only maintained for Ubuntu 18.04.
echo ========== Load config settings.

source settings.sh

################################################################################
# Misc. setup
################################################################################

note_build_stage "Install the runtime packages"

./run-prereq.sh

note_build_stage "Update package list"

sudo -H apt-get -qq -y update

note_build_stage "build-prereq.sh: Install development packages"

# redacted
retry=0
until [[ $retry -ge 3 ]]
do
  sudo -H apt-get -qq -y install pkg-config zip g++ zlib1g-dev unzip curl git wget > /dev/null && break

  echo "apt-get failed. Retrying."
  retry=$(($retry+1))
  sleep 10
done


################################################################################
# bazel
################################################################################

note_build_stage "Install bazel"

function ensure_wanted_bazel_version {
  local wanted_bazel_version=$1
  rm -rf ~/bazel
  mkdir ~/bazel

  if
    v=$(bazel --bazelrc=/dev/null --nomaster_bazelrc version) &&
    echo "$v" | awk -v b="$wanted_bazel_version" '/Build label/ { exit ($3 != b)}'
  then
    echo "Bazel ${wanted_bazel_version} already installed on the machine, not reinstalling"
  else
    pushd ~/bazel
    curl -L -O https://github.com/bazelbuild/bazel/releases/download/"${wanted_bazel_version}"/bazel-"${wanted_bazel_version}"-installer-linux-x86_64.sh
    chmod +x bazel-*.sh
    ./bazel-"${wanted_bazel_version}"-installer-linux-x86_64.sh --user > /dev/null
    rm bazel-"${wanted_bazel_version}"-installer-linux-x86_64.sh
    popd
  fi
}

ensure_wanted_bazel_version "${DV_BAZEL_VERSION}"

################################################################################
# CLIF
################################################################################

note_build_stage "Install CLIF binary"

if [[ -e /usr/local/bin/pyclif ]];
then
  echo "CLIF already installed."
else
  # Build clif binary from scratch. Might not be ideal because it installs a
  # bunch of dependencies, but this works fine when we used this in a Dockerfile
  # because we don't do build-prereq.sh in the final image.
  time ./tools/build_clif.sh
  # redacted
  #                    we can do this better.
  sudo mkdir -p /usr/clang/bin/
  sudo ln -sf /usr/local/bin/clif-matcher /usr/clang/bin/clif-matcher
  sudo mkdir -p /usr/local/clif/bin
  sudo ln -sf /usr/local/bin/pyclif* /usr/local/clif/bin/
  DIST_PACKAGES_DIR=$(python3 -c "import site; print(site.getsitepackages()[0])")
  sudo ln -sf "${DIST_PACKAGES_DIR}"/clif/python /usr/local/clif/
fi

################################################################################
# TensorFlow
################################################################################

note_build_stage "Download and configure TensorFlow sources"

if [[ ! -d ../tensorflow ]]; then
  note_build_stage "Cloning TensorFlow from github as ../tensorflow doesn't exist"
  (cd .. && git clone https://github.com/tensorflow/tensorflow)
fi

export PYTHON_BIN_PATH=$(which python3.6)
export PYTHON_LIB_PATH='/usr/local/lib/python3.6/dist-packages'
(cd ../tensorflow &&
 git checkout "${DV_CPP_TENSORFLOW_TAG}" &&
 echo | ./configure)

# We use TensorFlow's .bazelrc as part of DeepVariant's. In it they use a java
# toolchain flag based on a definition in a BUILD file in the TF repo. This
# causes that flag's usage to raise build errors when building DeepVariant
# unless we also include that BUILD file.
mkdir -p third_party/toolchains/java
cp ../tensorflow/third_party/toolchains/java/BUILD third_party/toolchains/java/

note_build_stage "build-prereq.sh complete"
