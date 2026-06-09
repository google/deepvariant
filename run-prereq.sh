#!/bin/bash

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

# This installs all the libraries (python, dso, etc) that are needed
# by DeepVariant at runtime (except for tensorflow, which is special).
# Some extra stuff may also be included.

set -euo pipefail

echo ========== This script is only maintained for Ubuntu 22.04.
echo ========== Load config settings.

source settings.sh

################################################################################
# misc setup
################################################################################

note_build_stage "Misc setup"

APT_ARGS=(
"-qq"
"-y"
)

UV_ARGS=()
if [[ "$EUID" = "0" ]]; then
  # Just in case:
  # https://github.com/NVIDIA/nvidia-docker/issues/1632#issuecomment-1112667716
  rm -f /etc/apt/sources.list.d/cuda.list
  rm -f /etc/apt/sources.list.d/nvidia-ml.list
  # Ensure sudo exists, even if we don't need it.
  apt-get update "${APT_ARGS[@]}" > /dev/null
  apt-get install "${APT_ARGS[@]}" sudo > /dev/null
  UV_ARGS+=("--system")
fi

note_build_stage "Update package list"

sudo -H apt-get update "${APT_ARGS[@]}" > /dev/null

note_build_stage "run-prereq.sh: Install development packages"

# Need to wait for dpkg lock (see internal)
wait_for_dpkg_lock

# See https://askubuntu.com/questions/909277.
sudo -H DEBIAN_FRONTEND=noninteractive apt-get install "${APT_ARGS[@]}" \
  gcc \
  git \
  curl \
  pkg-config \
  python3-distutils \
  python3-testresources \
  unzip \
  wget \
  zip \
  zlib1g-dev > /dev/null

# Install uv
note_build_stage "Install uv package manager"

curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"

# Verify uv is installed
uv --version

echo "$(python3 --version)"

################################################################################
# python packages
################################################################################

note_build_stage "Install python3 packages (batch with uv)"

# NOTE: Some packages are excluded from this batch because they have
# transitive dependency conflicts that uv's strict resolver catches.
# They are installed separately after TensorFlow.
cat > /tmp/requirements-dv.txt << 'EOF'
contextlib2
etils
importlib_resources
enum34==1.1.8
sortedcontainers==2.1.0
intervaltree==3.1.0
mock>=2.0.0
ml_collections
PyYAML
clu==0.0.9
protobuf==4.21.9
argparse==1.4.0
pyasn1<0.5.0,>=0.4.6
requests>=2.18
oauth2client>=4.0.0
crcmod>=1.7
six>=1.11.0
joblib
psutil
google-api-python-client==2.187.0
google-auth==2.47.0
google-auth-httplib2==0.3.0
httplib2==0.31.0
pandas==1.3.4
altair==5.5.0
jsonschema==4.17.3
Pillow==9.5.0
ipython==8.22.2
pysam==0.20.0
scikit-learn==1.0.2
setuptools==61.0.0
packaging==25.0
pyparsing>=3.0.0,<4.0.0
EOF

uv pip install "${UV_ARGS[@]}" -r /tmp/requirements-dv.txt
rm -f /tmp/requirements-dv.txt

# tf-models-official has transitive deps that conflict with TF 2.16.1's
# typing-extensions requirements. Install it separately without deps.
uv pip install "${UV_ARGS[@]}" --no-deps "tf-models-official==2.13.1"

################################################################################
# TensorFlow
################################################################################

note_build_stage "Install TensorFlow pip package"

if [[ "${DV_USE_PREINSTALLED_TF}" = "1" ]]; then
  echo "Skipping TensorFlow installation at user request; will use pre-installed TensorFlow."
else
  # Also pip install the latest TensorFlow with cpu support. We don't build the
  # full TF from source, but instead using prebuilt version. However, we still
  # need the full source version to build DeepVariant.

  # Gets the nightly TF build: https://pypi.python.org/pypi/tf-nightly which is
  # necessary right now if we aren't pinning the TF source. We have observed
  # runtime failures if there's too much skew between the released TF package and
  # the source.
  if [[ "${DV_TF_NIGHTLY_BUILD}" = "1" ]]; then
    if [[ "${DV_GPU_BUILD}" = "1" ]]; then
      echo "Installing GPU-enabled TensorFlow nightly wheel"
      uv pip install "${UV_ARGS[@]}" --upgrade tf_nightly_gpu
    else
      echo "Installing CPU-only TensorFlow nightly wheel"
      uv pip install "${UV_ARGS[@]}" --upgrade tf_nightly
    fi
  else
    # Use the official TF release pip package.
    if [[ "${DV_GPU_BUILD}" = "1" ]]; then
      echo "Installing GPU-enabled TensorFlow ${DV_TENSORFLOW_STANDARD_GPU_WHL_VERSION} wheel"
      uv pip install "${UV_ARGS[@]}" --upgrade "tensorflow[and-cuda]==${DV_TENSORFLOW_STANDARD_GPU_WHL_VERSION}"
    else
      echo "Installing CPU TensorFlow ${DV_TENSORFLOW_STANDARD_CPU_WHL_VERSION} wheel"
      uv pip install "${UV_ARGS[@]}" --upgrade "tensorflow==${DV_TENSORFLOW_STANDARD_CPU_WHL_VERSION}"
    fi
  fi
fi

# Post-TF install fixups.
# These packages have mutually incompatible transitive dependencies with TF,
# so they must be installed after TF in controlled groups.
uv pip install "${UV_ARGS[@]}" 'markupsafe==2.0.1' 'tf_keras==2.16.0' 'protobuf==4.21.9'
# jax must come after TF because it needs ml-dtypes>=0.4.0 which conflicts with
# TF's ml-dtypes<0.4.0 pin. Install jax without deps, then jaxlib+ml-dtypes.
uv pip install "${UV_ARGS[@]}" --no-deps 'jax==0.4.35'
uv pip install "${UV_ARGS[@]}" 'jaxlib==0.4.35' 'ml-dtypes>=0.4.0'

################################################################################
# CUDA
################################################################################

note_build_stage "Install CUDA"

# TF 2.16.1 requires CUDA 12.3 and CuDNN 8.9.
# See https://www.tensorflow.org/install/source#gpu for versions required.
if [[ "${DV_GPU_BUILD}" = "1" ]]; then
  if [[ "${DV_INSTALL_GPU_DRIVERS}" = "1" ]]; then
    # This script is only maintained for Ubuntu 22.04.
    echo "Checking for CUDA..."
    if ! dpkg-query -W cuda-12-3; then
      echo "Installing CUDA..."
      UBUNTU_VERSION="2204"
      curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VERSION}/x86_64/cuda-ubuntu${UBUNTU_VERSION}.pin
      sudo mv cuda-ubuntu${UBUNTU_VERSION}.pin /etc/apt/preferences.d/cuda-repository-pin-600

      curl https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub | gpg --dearmor | sudo tee /usr/share/keyrings/nvidia-cuda-archive-keyring.gpg > /dev/null
      echo \
        "deb [signed-by=/usr/share/keyrings/nvidia-cuda-archive-keyring.gpg] https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /" | \
        sudo tee /etc/apt/sources.list.d/cuda.list > /dev/null
      sudo -H NEEDRESTART_MODE=a apt-get update "${APT_ARGS[@]}"
      sudo -H DEBIAN_FRONTEND=noninteractive NEEDRESTART_MODE=a apt-get full-upgrade "${APT_ARGS[@]}"
      sudo -H DEBIAN_FRONTEND=noninteractive NEEDRESTART_MODE=a apt-get install "${APT_ARGS[@]}" cuda-12-3
    fi
    echo "Checking for CUDNN..."
    if [[ ! -e /usr/local/cuda-12/include/cudnn.h ]]; then
      echo "Installing CUDNN..."
      CUDNN_TAR_FILE="cudnn-linux-x86_64-8.9.0.131_cuda12-archive.tar.xz"
      wget -q https://developer.download.nvidia.com/compute/cudnn/redist/cudnn/linux-x86_64/${CUDNN_TAR_FILE}
      tar -xvf ${CUDNN_TAR_FILE}
      sudo cp -P cudnn-linux-x86_64-8.9.0.131_cuda12-archive/include/cudnn.h /usr/local/cuda-12/include
      sudo cp -P cudnn-linux-x86_64-8.9.0.131_cuda12-archive/lib/libcudnn* /usr/local/cuda-12/lib64/
      sudo chmod a+r /usr/local/cuda-12/lib64/libcudnn*
      sudo ldconfig
    fi
    sudo -H NEEDRESTART_MODE=a apt-get install "${APT_ARGS[@]}" libcupti-dev > /dev/null
  fi

  nvidia-smi || :
fi

################################################################################
# TensorRT
################################################################################

note_build_stage "Install TensorRT"

# Address the issue:
# 'dlerror: libnvinfer.so.7: cannot open shared object file: No such file or directory'
# It's unclear whether we need this or not. Setting up to get rid of the errors.

if [[ "${DV_GPU_BUILD}" = "1" ]]; then
  uv pip install "${UV_ARGS[@]}" tensorrt==8.5.3.1
  echo "For debugging:"
  uv pip show tensorrt
  TENSORRT_PATH=$(python3 -c 'import tensorrt; print(tensorrt.__path__[0])')
  sudo ln -sf "${TENSORRT_PATH}/libnvinfer.so.8" "${TENSORRT_PATH}/libnvinfer.so.7"
  sudo ln -sf "${TENSORRT_PATH}/libnvinfer_plugin.so.8" "${TENSORRT_PATH}/libnvinfer_plugin.so.7"
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH-}:${TENSORRT_PATH}"
  sudo ldconfig
  # Just in case this still doesn't work, we link them.
  # This is a workaround that we might want to get rid of, if we can make sure
  # setting LD_LIBRARY_PATH and `sudo ldconfig`` works.
  if [[ ! -e /usr/local/nvidia/lib ]]; then
    sudo mkdir -p /usr/local/nvidia/lib
    sudo ln -sf "${TENSORRT_PATH}/libnvinfer.so.7" /usr/local/nvidia/lib/libnvinfer.so.7
    sudo ln -sf "${TENSORRT_PATH}/libnvinfer_plugin.so.7" /usr/local/nvidia/lib/libnvinfer_plugin.so.7
  fi
fi

################################################################################
# Misc dependencies
################################################################################

note_build_stage "Install other packages"

sudo -H NEEDRESTART_MODE=a apt-get install "${APT_ARGS[@]}" \
  libboost-graph-dev \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblz-dev \
  liblzma-dev \
  libssl-dev > /dev/null

note_build_stage "Linking Cuda and TF shared libraries"
if [[ "${DV_GPU_BUILD}" = "1" ]]; then
  echo "####### For debugging linking of shared libraries in GPU build:"
  # Adding symbolic links to NVIDIA shared libraries. https://www.tensorflow.org/install/pip
  TENSORFLOW_DIR=$(dirname "$(python3 -c 'print(__import__("tensorflow").__file__)')")
  (cd "${TENSORFLOW_DIR}" && sudo ln -svf ../nvidia/*/lib/*.so* .)
fi

note_build_stage "run-prereq.sh complete"
