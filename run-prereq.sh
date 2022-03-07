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

echo ========== This script is only maintained for Ubuntu 20.04.
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

if [[ "$EUID" = "0" ]]; then
  # Ensure sudo exists, even if we don't need it.
  apt-get update "${APT_ARGS[@]}" > /dev/null
  apt-get install "${APT_ARGS[@]}" sudo > /dev/null
  PIP_ARGS=(
    "-qq")
else
  PIP_ARGS=(
    "--user"
    "-qq")
fi

note_build_stage "Update package list"

sudo -H apt-get update "${APT_ARGS[@]}" > /dev/null

note_build_stage "run-prereq.sh: Install development packages"

# Need to wait for dpkg lock (see internal)
wait_for_dpkg_lock

# See https://askubuntu.com/questions/909277.
sudo -H DEBIAN_FRONTEND=noninteractive apt-get install "${APT_ARGS[@]}" pkg-config zip zlib1g-dev unzip curl git wget > /dev/null
sudo -H apt-get install "${APT_ARGS[@]}" python3-distutils > /dev/null

note_build_stage "Install python3 packaging infrastructure"

sudo -H apt-get install "${APT_ARGS[@]}" "python${PYTHON_VERSION}-dev"
sudo update-alternatives --install /usr/bin/python3 python3 "/usr/bin/python${PYTHON_VERSION}" 0
sudo update-alternatives --install /usr/bin/python python "/usr/bin/python${PYTHON_VERSION}" 0
# If we install python3-pip directly, the pip3 version points to:
#   pip 8.1.1 from /usr/lib/python3/dist-packages (python 3.5)
# Use the following lines to ensure correct Python version.
curl -o get-pip.py https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --force-reinstall --user
rm -f get-pip.py

echo "$(python3 --version)"

export PATH="$HOME/.local/bin":$PATH
echo "$(pip3 --version)"

################################################################################
# python packages
################################################################################

note_build_stage "Install python3 packages"

pip3 install "${PIP_ARGS[@]}" contextlib2
pip3 install "${PIP_ARGS[@]}" 'enum34==1.1.8'
pip3 install "${PIP_ARGS[@]}" 'sortedcontainers==2.1.0'
pip3 install "${PIP_ARGS[@]}" 'intervaltree==3.0.2'
pip3 install "${PIP_ARGS[@]}" 'mock>=2.0.0'
# Note that protobuf installed with pip needs to be 3.13 because of the pyclif
# version we're using. This is currently inconsistent with C++ protobuf version
# in WORKSPACE and protobuf.BUILD, but we can't update those, because those
# files need to be consistent with what TensorFlow needs, which is currently
# still 3.9.2.
# Ideally we want to make these protobuf versions all match, eventually.
pip3 install "${PIP_ARGS[@]}" 'protobuf==3.13.0'
pip3 install "${PIP_ARGS[@]}" 'argparse==1.4.0'
pip3 install "${PIP_ARGS[@]}" git+https://github.com/google-research/tf-slim.git

pip3 install "${PIP_ARGS[@]}" "numpy==${DV_TF_NUMPY_VERSION}"

# Reason:
# ========== [Wed Dec 11 19:57:32 UTC 2019] Stage 'Install python3 packages' starting
# ERROR: pyasn1-modules 0.2.7 has requirement pyasn1<0.5.0,>=0.4.6, but you'll have pyasn1 0.1.9 which is incompatible.
pip3 install "${PIP_ARGS[@]}" 'pyasn1<0.5.0,>=0.4.6'
pip3 install "${PIP_ARGS[@]}" 'requests>=2.18'
pip3 install "${PIP_ARGS[@]}" --ignore-installed 'oauth2client>=4.0.0'
pip3 install "${PIP_ARGS[@]}" 'crcmod>=1.7'
pip3 install "${PIP_ARGS[@]}" 'six>=1.11.0'
pip3 install "${PIP_ARGS[@]}" joblib
pip3 install "${PIP_ARGS[@]}" psutil
pip3 install "${PIP_ARGS[@]}" --upgrade google-api-python-client
pip3 install "${PIP_ARGS[@]}" 'pandas==1.2.4'
# We manually install jsonschema here to pin it to v3.2.0, since
# the latest v4.0.1 has issues with Altair v4.1.0.
# See https://github.com/altair-viz/altair/issues/2496
# If Altair version is updated below, the jsonschema version
# should also be updated accordingly.
pip3 install "${PIP_ARGS[@]}" 'jsonschema==3.2.0'
pip3 install "${PIP_ARGS[@]}" 'altair==4.1.0'
pip3 install "${PIP_ARGS[@]}" 'Pillow>=5.4.1'
pip3 install "${PIP_ARGS[@]}" 'ipython>=7.9.0'


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
      pip3 install "${PIP_ARGS[@]}" --upgrade tf_nightly_gpu
    else
      echo "Installing CPU-only TensorFlow nightly wheel"
      pip3 install "${PIP_ARGS[@]}" --upgrade tf_nightly
    fi
  else
    # Use the official TF release pip package.
    if [[ "${DV_GPU_BUILD}" = "1" ]]; then
      echo "Installing GPU-enabled TensorFlow ${DV_TENSORFLOW_STANDARD_GPU_WHL_VERSION} wheel"
      pip3 install "${PIP_ARGS[@]}" --upgrade "tensorflow-gpu==${DV_TENSORFLOW_STANDARD_GPU_WHL_VERSION}"
    elif [[ "${DV_USE_GCP_OPTIMIZED_TF_WHL}" = "1" ]]; then
      echo "Installing Intel's CPU-only MKL TensorFlow ${DV_GCP_OPTIMIZED_TF_WHL_VERSION} wheel"
      pip3 install "${PIP_ARGS[@]}" --upgrade "intel-tensorflow==${DV_GCP_OPTIMIZED_TF_WHL_VERSION}"
    else
      echo "Installing standard CPU-only TensorFlow ${DV_TENSORFLOW_STANDARD_CPU_WHL_VERSION} wheel"
      pip3 install "${PIP_ARGS[@]}" --upgrade "tensorflow==${DV_TENSORFLOW_STANDARD_CPU_WHL_VERSION}"
    fi
  fi
fi


################################################################################
# CUDA
################################################################################

# See https://www.tensorflow.org/install/source#gpu for versions required.
if [[ "${DV_GPU_BUILD}" = "1" ]]; then
  if [[ "${DV_INSTALL_GPU_DRIVERS}" = "1" ]]; then
    # This script is only maintained for Ubuntu 20.04.
    UBUNTU_VERSION="2004"
    # https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=20.04&target_type=deb_local
    echo "Checking for CUDA..."
    if ! dpkg-query -W cuda-11-3; then
      echo "Installing CUDA..."
      CUDA_DEB="cuda-repo-ubuntu${UBUNTU_VERSION}-11-3-local_11.3.0-465.19.01-1_amd64.deb"
      curl -O "https://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VERSION}/x86_64/cuda-ubuntu2004.pin"
      sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
      curl -O https://developer.download.nvidia.com/compute/cuda/11.3.0/local_installers/"${CUDA_DEB}"
      sudo -H apt-key adv --fetch-keys "http://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VERSION}/x86_64/7fa2af80.pub"
      sudo -H dpkg -i "./${CUDA_DEB}"
      sudo apt-key add /var/cuda-repo-ubuntu2004-11-3-local/7fa2af80.pub
      sudo -H apt-get update "${APT_ARGS[@]}" > /dev/null
      sudo -H apt-get install "${APT_ARGS[@]}" cuda
    fi
    echo "Checking for CUDNN..."
    if [[ ! -e /usr/local/cuda-11/include/cudnn.h ]]; then
      echo "Installing CUDNN..."
      CUDNN_TAR_FILE="cudnn-11.3-linux-x64-v8.2.0.53.tgz"
      wget -q https://developer.download.nvidia.com/compute/redist/cudnn/v8.2.0/${CUDNN_TAR_FILE}
      tar -xzvf ${CUDNN_TAR_FILE}
      sudo cp -P cuda/include/cudnn.h /usr/local/cuda-11/include
      sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-11/lib64/
      sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-11/lib64/
      sudo chmod a+r /usr/local/cuda-11/lib64/libcudnn*
      sudo ldconfig
    fi
    # Tensorflow says to do this.
    sudo -H apt-get install "${APT_ARGS[@]}" libcupti-dev > /dev/null
  fi

  # If we are doing a gpu-build, nvidia-smi should be install. Run it so we
  # can see what gpu is installed.
  nvidia-smi || :
fi

################################################################################
# OpenVINO
################################################################################
if [[ "${DV_OPENVINO_BUILD}" = "1" ]]; then
  pip3 install "${PIP_ARGS[@]}" --no-deps networkx defusedxml openvino-dev~=2022.1.0.dev
  pip3 install "${PIP_ARGS[@]}" openvino~=2022.1.0.dev
fi

################################################################################
# Misc dependencies
################################################################################

note_build_stage "Install other packages"

# for htslib
sudo -H apt-get install "${APT_ARGS[@]}" libssl-dev libcurl4-openssl-dev liblz-dev libbz2-dev liblzma-dev > /dev/null

# for the debruijn graph
sudo -H apt-get install "${APT_ARGS[@]}" libboost-graph-dev > /dev/null

note_build_stage "run-prereq.sh complete"
