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
# This script is used to install nvidia docker on Ubuntu 22.04.
# For different Linux distributions and versions, modifications might be needed.

set -euo pipefail

APT_ARGS=(
"-qq"
"-y"
)

# Installing nvidia docker to use deepvariant_gpu Docker image.
# (1) Install nvidia driver:
# https://linuxhint.com/install-cuda-ubuntu/
sudo apt-get "${APT_ARGS[@]}" update
sudo apt-get "${APT_ARGS[@]}" install \
  build-essential \
  curl \
  "linux-headers-$(uname -r)" \
  nvidia-cuda-toolkit

# See https://www.tensorflow.org/install/source#gpu for versions required.
if ! dpkg-query -W cuda-11-8; then
  echo "Installing CUDA..."
  UBUNTU_VERSION="2204"
  curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VERSION}/x86_64/cuda-ubuntu${UBUNTU_VERSION}.pin
  sudo mv cuda-ubuntu${UBUNTU_VERSION}.pin /etc/apt/preferences.d/cuda-repository-pin-600

  curl https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub | gpg --dearmor | sudo tee /usr/share/keyrings/nvidia-cuda-archive-keyring.gpg > /dev/null
  echo \
    "deb [signed-by=/usr/share/keyrings/nvidia-cuda-archive-keyring.gpg] https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /" | \
    sudo tee /etc/apt/sources.list.d/cuda.list > /dev/null
  sudo -H apt-get update "${APT_ARGS[@]}"
  sudo -H apt-get full-upgrade "${APT_ARGS[@]}"
  sudo -H apt-get install "${APT_ARGS[@]}" cuda-11-8
fi

echo "Checking for CUDNN..."
if [[ ! -e /usr/local/cuda-11/include/cudnn.h ]]; then
  echo "Installing CUDNN..."
  CUDNN_TAR_FILE="cudnn-linux-x86_64-8.6.0.163_cuda11-archive.tar.xz"
  wget -q https://developer.download.nvidia.com/compute/redist/cudnn/v8.6.0/local_installers/11.8/${CUDNN_TAR_FILE}
  tar -xvf ${CUDNN_TAR_FILE}
  sudo cp -P cudnn-linux-x86_64-8.6.0.163_cuda11-archive/include/cudnn.h /usr/local/cuda-11/include
  sudo cp -P cudnn-linux-x86_64-8.6.0.163_cuda11-archive/lib/libcudnn* /usr/local/cuda-11/lib64/
  sudo chmod a+r /usr/local/cuda-11/lib64/libcudnn*
  sudo ldconfig
fi

# (2) Install Docker CE:
# https://docs.docker.com/engine/install/ubuntu/
sudo apt-get "${APT_ARGS[@]}" update
sudo apt-get "${APT_ARGS[@]}" install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get "${APT_ARGS[@]}" update
sudo apt-get "${APT_ARGS[@]}" install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# (3) Install nvidia docker:
# https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt-get "${APT_ARGS[@]}" update
sudo apt-get "${APT_ARGS[@]}" install nvidia-container-toolkit

sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

#### Test nvidia-smi with the latest official CUDA image
sudo docker run --gpus 1 nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04 nvidia-smi
