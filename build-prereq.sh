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

echo ========== This script is only maintained for Ubuntu 20.04.
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

# Need to wait for dpkg lock (see internal)
wait_for_dpkg_lock
sudo -H apt-get -qq -y install pkg-config zip g++ zlib1g-dev unzip curl git wget > /dev/null


################################################################################
# bazel
################################################################################

note_build_stage "Install bazel"

function ensure_wanted_bazel_version {
  local wanted_bazel_version=$1
  rm -rf ~/bazel
  mkdir ~/bazel

  if
    v=$(bazel --bazelrc=/dev/null --ignore_all_rc_files version) &&
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
  note_build_stage "Build CLIF."
  time sudo ./tools/build_clif.sh
  # TODO:
  # Figure out why these symbolic links are needed and see if
  # we can do this better.
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

# PYTHON_BIN_PATH and PYTHON_LIB_PATH are set in settings.sh.
# I had to remove this line in tensorflow v2.5.0 because I got an ERROR:
# rule() got unexpected keyword argument 'incompatible_use_toolchain_transition'.
# I changed the llvm path to zip to avoid flakiness.
(cd ../tensorflow &&
 git checkout "${DV_CPP_TENSORFLOW_TAG}" &&
 echo | ./configure)

# We want to use a newer absl version. So I grabbed the one from TensorFlow
# r2.13. Eventually we'll want to update to TF 2.13. But for now this works.
# TODO: After updating to v2.13, we can remove this.
wget https://raw.githubusercontent.com/tensorflow/tensorflow/r2.13/third_party/absl/workspace.bzl -O ../tensorflow/third_party/absl/workspace.bzl
wget https://raw.githubusercontent.com/tensorflow/tensorflow/r2.13/third_party/absl/absl_designated_initializers.patch -O ../tensorflow/third_party/absl/absl_design\
ated_initializers.patch

note_build_stage "Set pyparsing to 2.2.0 for CLIF."
export PATH="$HOME/.local/bin":$PATH
pip3 uninstall -y pyparsing && pip3 install -Iv 'pyparsing==2.2.0'

note_build_stage "build-prereq.sh complete"
