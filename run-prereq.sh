#!/bin/bash

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

if [[ "$EUID" = "0" ]]; then
  # Ensure sudo exists, even if we don't need it.
  apt-get update
  apt-get -y install sudo
fi

note_build_stage "Update package list"

sudo -H apt-get -qq -y update

note_build_stage "Install development packages"

sudo -H apt-get -y install pkg-config zip zlib1g-dev unzip curl

note_build_stage "Install python packaging infrastructure"

sudo -H apt-get -y install python-dev python-pip python-wheel
sudo -H pip install --upgrade pip

################################################################################
# python packages
################################################################################

note_build_stage "Install python packages"

sudo -H pip install contextlib2
sudo -H pip install enum34
# sortedcontainers>=2.0.0 breaks intervaltree=2.1.0
# Remove this when https://github.com/chaimleib/intervaltree/pull/69
# is resolved.  Tracked internally at b/80085543.
sudo -H pip install 'sortedcontainers==1.5.3'
sudo -H pip install intervaltree
sudo -H pip install 'mock>=2.0.0'
sudo -H pip install 'numpy==1.14'  # To match GCP_OPTIMIZED_TF_WHL_FILENAME
sudo -H pip install 'requests>=2.18'
sudo -H pip install 'scipy==1.0'
sudo -H pip install 'oauth2client>=4.0.0'
sudo -H pip install 'crcmod>=1.7'
sudo -H pip install six
sudo -H pip install sklearn
sudo -H pip install pandas
sudo -H pip install psutil
sudo -H pip install --upgrade google-api-python-client


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
  # First we uninstall any and all tensorflow's already present to ensure our
  # subsequent installs work.
  sudo -H pip uninstall -y tensorflow tf_nightly_gpu tf_nightly tensorflow-gpu
  if [[ "${DV_TF_NIGHTLY_BUILD}" = "1" ]]; then
    if [[ "${DV_GPU_BUILD}" = "1" ]]; then
      echo "Installing GPU-enabled TensorFlow nightly wheel"
      sudo -H pip install --upgrade tf_nightly_gpu
    else
      echo "Installing CPU-only TensorFlow nightly wheel"
      sudo -H pip install --upgrade tf_nightly
    fi
  else
    # Use the official TF release pip package.
    if [[ "${DV_GPU_BUILD}" = "1" ]]; then
      echo "Installing GPU-enabled TensorFlow wheel"
      sudo -H pip install --upgrade "tensorflow-gpu==${DV_TENSORFLOW_STANDARD_GPU_WHL_VERSION}"
    elif [[ "${DV_USE_GCP_OPTIMIZED_TF_WHL}" = "1" ]]; then
      echo "Installing Google Cloud Platform optimized CPU-only TensorFlow wheel"
      curl "${GCP_OPTIMIZED_TF_WHL_CURL_PATH}/${GCP_OPTIMIZED_TF_WHL_FILENAME}" \
        > "/tmp/${GCP_OPTIMIZED_TF_WHL_FILENAME}"
      sudo -H pip install --upgrade "/tmp/${GCP_OPTIMIZED_TF_WHL_FILENAME}"
    else
      echo "Installing standard CPU-only TensorFlow wheel"
      sudo -H pip install --upgrade "tensorflow==${DV_TENSORFLOW_STANDARD_CPU_WHL_VERSION}"
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

    # from https://cloud.google.com/compute/docs/gpus/add-gpus
    echo "Checking for CUDA..."
    if ! dpkg-query -W cuda-9-0; then
      echo "Installing CUDA..."
      CUDA_DEB="cuda-repo-ubuntu1604_9.0.176-1_amd64.deb"
      curl -O http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/${CUDA_DEB}
      sudo -H apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
      sudo -H dpkg -i ./cuda-repo-ubuntu1604_9.0.176-1_amd64.deb
      sudo -H apt-get update
      sudo -H apt-get -y install cuda-9-0
    fi

    echo "Checking for CUDNN..."
    if [[ ! -e /usr/local/cuda-9.0/include/cudnn.h ]]; then
      echo "Installing CUDNN..."
      CUDNN_TAR_FILE="cudnn-9.0-linux-x64-v7.tgz"
      wget http://developer.download.nvidia.com/compute/redist/cudnn/v7.0.5/${CUDNN_TAR_FILE}
      tar -xzvf ${CUDNN_TAR_FILE}
      sudo cp -P cuda/include/cudnn.h /usr/local/cuda-9.0/include
      sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-9.0/lib64/
      sudo chmod a+r /usr/local/cuda-9.0/lib64/libcudnn*
      sudo ldconfig
    fi

    # Tensorflow says to do this.
    sudo -H apt-get -y install libcupti-dev
  fi

  # If we are doing a gpu-build, nvidia-smi should be install. Run it so we
  # can see what gpu is installed.
  nvidia-smi || :
fi


################################################################################
# Misc dependencies
################################################################################

note_build_stage "Install other packages"

# for htslib
sudo -H apt-get -y install libssl-dev libcurl4-openssl-dev liblz-dev libbz2-dev liblzma-dev

# for the debruijn graph
sudo -H apt-get -y install libboost-graph-dev

note_build_stage "run-prereq.sh complete"
