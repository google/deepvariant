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

echo ========== This script is only maintained for Ubuntu 22.04.
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
sudo -H NEEDRESTART_MODE=a apt-get -qq -y install pkg-config zip g++ zlib1g-dev unzip curl git wget > /dev/null


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

# This is used for building examples_from_stream.so later.
time sudo ./tools/build_absl.sh

################################################################################
# TensorFlow
################################################################################

note_build_stage "Download and configure TensorFlow sources"

# Getting the directory before switching out.
DV_DIR=$(pwd)

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
rm -f ../tensorflow/third_party/absl/absl_designated_initializers.patch
# To get the @com_google_absl//absl/strings:string_view target:
sed -i -e 's|b971ac5250ea8de900eae9f95e06548d14cd95fe|29bf8085f3bf17b84d30e34b3d7ff8248fda404e|g' ../tensorflow/third_party/absl/workspace.bzl
sed -i -e 's|8eeec9382fc0338ef5c60053f3a4b0e0708361375fe51c9e65d0ce46ccfe55a7|affb64f374b16877e47009df966d0a9403dbf7fe613fe1f18e49802c84f6421e|g' ../tensorflow/third_party/absl/workspace.bzl
sed -i -e 's|patch_file = \["//third_party/absl:absl_designated_initializers.patch"\],||g' ../tensorflow/third_party/absl/workspace.bzl

# Update tensorflow.bzl. This updates the `pybind_extension` rule to use the
# _message.so file.
patch ../tensorflow/tensorflow/tensorflow.bzl "${DV_DIR}"/third_party/tensorflow.bzl.patch

# I want to replace this part in ../tensorflow/tensorflow/workspace2.bzl
# From:
# tf_http_archive(
#     name = "pybind11",
#     urls = tf_mirror_urls("https://github.com/pybind/pybind11/archive/v2.10.0.tar.gz"),
#     sha256 = "eacf582fa8f696227988d08cfc46121770823839fe9e301a20fbce67e7cd70ec",
#     strip_prefix = "pybind11-2.10.0",
#     build_file = "//third_party:pybind11.BUILD",
#     system_build_file = "//third_party/systemlibs:pybind11.BUILD",
# )
# To:
# tf_http_archive(
#     name = "pybind11",
#     urls = tf_mirror_urls("https://github.com/pybind/pybind11/archive/a7b91e33269ab6f3f90167291af2c4179fc878f5.zip"),
#     sha256 = "09d2ab67e91457c966eb335b361bdc4d27ece2d4dea681d22e5d8307e0e0c023",
#     strip_prefix = "pybind11-a7b91e33269ab6f3f90167291af2c4179fc878f5",
#     build_file = "//third_party:pybind11.BUILD",
#     system_build_file = "//third_party/systemlibs:pybind11.BUILD",
# )
sed -i -e 's|v2.10.0.tar.gz|a7b91e33269ab6f3f90167291af2c4179fc878f5.zip|g' ../tensorflow/tensorflow/workspace2.bzl
sed -i -e 's|eacf582fa8f696227988d08cfc46121770823839fe9e301a20fbce67e7cd70ec|09d2ab67e91457c966eb335b361bdc4d27ece2d4dea681d22e5d8307e0e0c023|g' ../tensorflow/tensorflow/workspace2.bzl
sed -i -e 's|pybind11-2.10.0|pybind11-a7b91e33269ab6f3f90167291af2c4179fc878f5|g' ../tensorflow/tensorflow/workspace2.bzl

# TODO: Test removing this version pinning.
note_build_stage "Set pyparsing to 2.2.2 for CLIF."
export PATH="$HOME/.local/bin":$PATH
pip3 uninstall -y pyparsing && pip3 install -Iv 'pyparsing==2.2.2'

note_build_stage "build-prereq.sh complete"
