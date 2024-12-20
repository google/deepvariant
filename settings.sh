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

# Source this file---these options are needed for TF config and for
# successive bazel runs.

# Set this to 1 if the system image already has TensorFlow preinstalled.  This
# will skip the installation of TensorFlow.
export DV_USE_PREINSTALLED_TF="${DV_USE_PREINSTALLED_TF:-0}"

export TF_NEED_GCP=1

export CUDNN_INSTALL_PATH="/usr/lib/x86_64-linux-gnu"

# The version of bazel we want to build DeepVariant.
# https://www.tensorflow.org/install/source#tested_build_configurations
DV_BAZEL_VERSION="5.3.0"

# We need to make sure that $HOME/bin is first in the binary search path so that
# `bazel` will find the latest version of bazel installed in the user's home
# directory. This is set in setting.sh as all DeepVariant scripts source
# settings.sh and assume that `bazel` will find the right version.
export PATH="$HOME/bin:$PATH"

# Path to the public bucket containing DeepVariant-related artifacts.
export DEEPVARIANT_BUCKET="gs://deepvariant"
export DV_PACKAGE_BUCKET_PATH="${DEEPVARIANT_BUCKET}/packages"
export DV_PACKAGE_CURL_PATH="https://storage.googleapis.com/deepvariant/packages"

# Set this to 1 to use the nightly (latest) build of TensorFlow instead of a
# named release version. Set it to an already existing value in the environment
# (allowing command line control of the build), defaulting to 0 (release build).
# Note that setting this to 1 implies that the C++ code in DeepVariant will be
# build using the master branch and not the pinned version to avoid
# incompatibilities between TensorFlow C++ used to build DeepVariant and the
# tf-nightly wheel.
export DV_TF_NIGHTLY_BUILD="${DV_TF_NIGHTLY_BUILD:-0}"

# The branch/tag we checkout to build our C++ dependencies against. This is not
# the same as the python version of TensorFlow we use, but should be similar or
# we risk having version incompatibilities between our C++ code and the Python
# code we use at runtime.
if [[ "${DV_TF_NIGHTLY_BUILD}" = "1" ]]; then
  export DV_CPP_TENSORFLOW_TAG="master"
else
  export DV_CPP_TENSORFLOW_TAG="v2.13.1"
fi
# These WHL_VERSIONs determine the Python version of TensorFlow we use.
export DV_GCP_OPTIMIZED_TF_WHL_VERSION="2.13.1"
export DV_TENSORFLOW_STANDARD_GPU_WHL_VERSION="2.13.1"
export DV_TENSORFLOW_STANDARD_CPU_WHL_VERSION="2.13.1"

# Set this to 1 to use DeepVariant with GPUs. Set it to an already existing
# value in the environment (allowing command line control of the build),
# defaulting to 0 (CPU only build).
export DV_GPU_BUILD="${DV_GPU_BUILD:-0}"

# NOTE: CPU TensorFlow has a TF_ENABLE_ONEDNN_OPTS option that can be used to
# enable Intel-specific optimization.
export GCP_OPTIMIZED_TF_WHL_FILENAME="tensorflow-${DV_GCP_OPTIMIZED_TF_WHL_VERSION}.deepvariant_gcp-cp27-none-linux_x86_64.whl"
export GCP_OPTIMIZED_TF_WHL_PATH="${DV_PACKAGE_BUCKET_PATH}/tensorflow"
export GCP_OPTIMIZED_TF_WHL_CURL_PATH="${DV_PACKAGE_CURL_PATH}/tensorflow"

# Set this to 1 to make our prereq scripts install the CUDA libraries.
# If you already have CUDA installed, such as on a properly provisioned
# Docker image, it shouldn't be necessary.
export DV_INSTALL_GPU_DRIVERS="${DV_INSTALL_GPU_DRIVERS:-0}"

export PYTHON_VERSION=3.10
# shellcheck disable=SC2155
export PYTHON_BIN_PATH="$(which python${PYTHON_VERSION})"
export PYTHON_LIB_PATH="/usr/local/lib/python${PYTHON_VERSION}/dist-packages"
export USE_DEFAULT_PYTHON_LIB_PATH=1
# N.B. The --experimental_build_setting_api had to be added on protobuf
# upgrade to 3.9.2 to avoid error in bazel_skylib:
#   "parameter 'build_setting' is experimental and thus unavailable with the
#    current flags. It may be enabled by setting
#    --experimental_build_setting_api"
# Presumably it won't be needed at some later point when bazel_skylib is
# upgraded again.
export DV_COPT_FLAGS="--copt=-march=corei7 --copt=-Wno-sign-compare --copt=-Wno-write-strings --experimental_build_setting_api --java_runtime_version=remotejdk_11"

function note_build_stage {
  echo "========== [$(date)] Stage '${1}' starting"
}

function wait_for_dpkg_lock {
  # Wait for at most 5 minutes.
  echo "Calling wait_for_dpkg_lock."
  max_wait=300
  i=0
  while sudo fuser /var/lib/dpkg/{lock,lock-frontend} >/dev/null 2>&1 ; do
    echo "Waiting to obtain dpkg lock.."
    sleep 10
    ((i=i+10))
    if (( i > max_wait )); then
      echo "ERROR: Waited for dpkg lock for 5 minutes."
      exit 1
    fi
  done
}
