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

# Source this file---these options are needed for TF config and for
# successive bazel runs.

# Set this to 1 if the system image already has TensorFlow preinstalled.  This
# will skip the installation of TensorFlow.
export DV_USE_PREINSTALLED_TF="${DV_USE_PREINSTALLED_TF:-0}"

export TF_CUDA_CLANG=0
export TF_ENABLE_XLA=0
export TF_NEED_CUDA=0
export TF_NEED_GCP=1
export TF_NEED_GDR=0
export TF_NEED_HDFS=0
export TF_NEED_JEMALLOC=0
export TF_NEED_MKL=0
export TF_NEED_MPI=0
export TF_NEED_OPENCL=0
export TF_NEED_OPENCL_SYCL=0
export TF_NEED_S3=0
export TF_NEED_VERBS=0

# Used if TF_NEED_CUDA=1
export TF_CUDA_VERSION="8.0"
export CUDA_TOOLKIT_PATH="/usr/local/cuda"
export TF_CUDNN_VERSION="6"
export CUDNN_INSTALL_PATH="/usr/lib/x86_64-linux-gnu"

# Path to the public bucket containing DeepVariant-related artifacts.
export DEEPVARIANT_BUCKET="gs://deepvariant"
export DV_PACKAGE_BUCKET_PATH="${DEEPVARIANT_BUCKET}/packages"
export DV_PACKAGE_CURL_PATH="https://storage.googleapis.com/deepvariant/packages"

# Set this to 1 to use DeepVariant with GPUs. Set it to an already existing
# value in the environment (allowing command line control of the build),
# defaulting to 0 (CPU only build).
export DV_GPU_BUILD="${DV_GPU_BUILD:-0}"

# If this variable is set to 1, DeepVariant will use a TensorFlow wheel file
# compiled to use AVX and SSE instructions. This instructions require Sandy
# Bridge or better chipsets on the host machine. The default TensorFlow wheel
# files don't contain these instructions (and thereby run on a broader set of
# CPUs). Using this optimized wheel reduces the runtime of DeepVariant's
# call_variants step by ~20%. This is called the GCP (Google Cloud Platform)
# optimized wheel because all GCP instances have at least Sandy Bridge or better
# chipsets, so this wheel should run anywhere on GCP.
export DV_USE_GCP_OPTIMIZED_TF_WHL="${DV_USE_GCP_OPTIMIZED_TF_WHL:-1}"
export GCP_OPTIMIZED_TF_WHL_FILENAME="tensorflow-1.4.1.deepvariant_gcp-cp27-none-linux_x86_64.whl"
export GCP_OPTIMIZED_TF_WHL_PATH="${DV_PACKAGE_BUCKET_PATH}/tensorflow"
export GCP_OPTIMIZED_TF_WHL_CURL_PATH="${DV_PACKAGE_CURL_PATH}/tensorflow"

# Set this to 1 to use the nightly (latest) build of TensorFlow instead of a
# named release version. Set it to an already existing value in the environment
# (allowing command line control of the build), defaulting to 0 (release build).
export DV_TF_NIGHTLY_BUILD="${DV_TF_NIGHTLY_BUILD:-0}"

# Set this to 1 to make our prereq scripts install the CUDA libraries.
# If you already have CUDA installed, such as on a properly provisioned
# Docker image, it shouldn't be necessary.
export DV_INSTALL_GPU_DRIVERS="${DV_INSTALL_GPU_DRIVERS:-0}"

export PYTHON_BIN_PATH=$(which python)
export USE_DEFAULT_PYTHON_LIB_PATH=1
export DV_COPT_FLAGS="--copt=-msse4.1 --copt=-msse4.2 --copt=-mavx --copt=-O3"

function note_build_stage {
  echo "========== [$(date)] Stage '${1}' starting"
}
