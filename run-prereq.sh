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

note_build_stage "Install development packages"

sudo -H apt-get install "${APT_ARGS[@]}" pkg-config zip zlib1g-dev unzip curl git lsb-release wget > /dev/null

note_build_stage "Install python3 packaging infrastructure"

# For altair to work, we need Python to be >= 3.5.3.
# https://github.com/altair-viz/altair/issues/972
# On Ubuntu 16.04, default is 3.5.2.
# So, install Python 3.6.
sudo -H apt-get install "${APT_ARGS[@]}" software-properties-common

# Install Python 3.6.
# Reference: https://askubuntu.com/a/1069303
sudo -H -E add-apt-repository -y ppa:deadsnakes/ppa
sudo -H apt update "${APT_ARGS[@]}"
sudo -H apt-get install "${APT_ARGS[@]}" python3.6
sudo -H apt-get install "${APT_ARGS[@]}" python3.6-dev
sudo -H apt-get install "${APT_ARGS[@]}" python3.6-venv
sudo ln -sf /usr/bin/python3.6 /usr/local/bin/python3
sudo ln -sf /usr/bin/python3.6 /usr/bin/python
# If we install python3-pip directly, the pip3 version points to:
#   pip 8.1.1 from /usr/lib/python3/dist-packages (python 3.5)
# Use the following lines to ensure 3.6.
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
pip3 install "${PIP_ARGS[@]}" 'protobuf==3.12.0'
pip3 install "${PIP_ARGS[@]}" 'argparse==1.4.0'
pip3 install "${PIP_ARGS[@]}" git+https://github.com/google-research/tf-slim.git

# Because of an issue with pypi's numpy on Ubuntu 14.04. we need to compile from
# source. But we know that on 16.04 we don't need to compile from source
# See https://github.com/tensorflow/tensorflow/issues/6968#issuecomment-279061085
if [[ "$(lsb_release -d)" == *Ubuntu*16.04.* ]]; then
  pip3 install "${PIP_ARGS[@]}" "numpy==${DV_TF_NUMPY_VERSION}"
else
  echo "Installing numpy with -no-binary=:all:. This will take a bit longer."
  pip3 install "${PIP_ARGS[@]}" --no-binary=:all: "numpy==${DV_TF_NUMPY_VERSION}"
fi

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
pip3 install "${PIP_ARGS[@]}" 'pandas==0.24.1'
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

if [[ "${DV_GPU_BUILD}" = "1" ]]; then
  if [[ "${DV_INSTALL_GPU_DRIVERS}" = "1" ]]; then
    if [[ "$(lsb_release -d)" != *Ubuntu*16.*.* ]]; then
      echo "CUDA installation only configured for Ubuntu 16"
      exit 1
    fi

    # https://www.tensorflow.org/install/gpu?hl=en#ubuntu_1604_cuda_101
    echo "Checking for CUDA..."
    if ! dpkg-query -W cuda-10-1; then
      echo "Installing CUDA..."
      CUDA_DEB="cuda-repo-ubuntu1604_10.1.243-1_amd64.deb"
      curl -O http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/${CUDA_DEB}
      sudo -H apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
      sudo -H dpkg -i "./${CUDA_DEB}"
      sudo -H apt-get update "${APT_ARGS[@]}" > /dev/null
      sudo -H apt-get install "${APT_ARGS[@]}" cuda-10-1 > /dev/null
    fi
    echo "Checking for CUDNN..."
    if [[ ! -e /usr/local/cuda-10.1/include/cudnn.h ]]; then
      echo "Installing CUDNN..."
      CUDNN_TAR_FILE="cudnn-10.1-linux-x64-v7.6.5.32.tgz"
      wget -q https://developer.download.nvidia.com/compute/redist/cudnn/v7.6.5/${CUDNN_TAR_FILE}
      tar -xzvf ${CUDNN_TAR_FILE}
      sudo cp -P cuda/include/cudnn.h /usr/local/cuda-10.1/include
      sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-10.1/lib64/
      sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-10.1/lib64/
      sudo chmod a+r /usr/local/cuda-10.1/lib64/libcudnn*
      # https://stackoverflow.com/questions/55224016/importerror-libcublas-so-10-0-cannot-open-shared-object-file-no-such-file-or/64472380#64472380
      export LD_LIBRARY_PATH=/usr/local/cuda-10.2/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
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
  pip3 install "${PIP_ARGS[@]}" git+https://github.com/openvinotoolkit/openvino.git@releases/2021/3#subdirectory=model-optimizer
  pip3 install "${PIP_ARGS[@]}" openvino==2021.3
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
