#!/bin/bash
# Copyright 2019 Google LLC.
# This script is used to install `nvidia-docker` on Ubutun 16.04.
# For different Linux distributions and versions, modifications might be needed.

set -euo pipefail

# Installing `nvidia-docker` to use deepvariant_gpu Docker image.
# (1) Install nvidia driver:
# https://github.com/NVIDIA/nvidia-docker/wiki/Frequently-Asked-Questions#how-do-i-install-the-nvidia-driver
sudo apt-get -qq -y update
# From: https://docs.docker.com/install/linux/docker-ce/ubuntu/#set-up-the-repository
sudo apt-get -qq -y install \
  apt-transport-https \
  ca-certificates \
  curl \
  gnupg-agent \
  software-properties-common

echo "Installing CUDA..."
if ! dpkg-query -W cuda-10-0; then
  echo "Installing CUDA..."
  CUDA_DEB="cuda-repo-ubuntu1604_10.0.130-1_amd64.deb"
  curl -O http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/${CUDA_DEB}
  sudo -H apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
  sudo -H dpkg -i "./${CUDA_DEB}"
  sudo -H apt-get -qq -y update > /dev/null
  sudo -H apt-get -qq -y install cuda-10-0 > /dev/null
fi

# (2) Install Docker CE:
# https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo apt-key fingerprint 0EBFCD88
sudo add-apt-repository \
 "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
 $(lsb_release -cs) \
 stable"
sudo apt-get -qq -y update
sudo apt-get -qq -y install docker-ce

# (3) Install nvidia-docker:
# https://github.com/NVIDIA/nvidia-docker#ubuntu-140416041804-debian-jessiestretch
# "If you have nvidia-docker 1.0 installed: we need to remove it and all existing GPU containers"
sudo docker volume ls -q -f driver=nvidia-docker | xargs -r -I{} -n1 sudo docker ps -q -a -f volume={} | xargs -r sudo docker rm -f
sudo apt-get purge -y nvidia-docker || true

# Add the package repositories
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | \
  sudo apt-key add -
distribution=$(. /etc/os-release;echo "${ID}${VERSION_ID}")
curl -s -L "https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list" | \
  sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get -qq -y update

# Install nvidia-docker2 and reload the Docker daemon configuration
sudo apt-get install -qq -y nvidia-docker2
sudo pkill -SIGHUP dockerd

# Test nvidia-smi with the latest official CUDA image
sudo docker run --runtime=nvidia --rm nvidia/cuda nvidia-smi
