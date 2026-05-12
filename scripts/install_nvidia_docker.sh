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

echo "====================================================================="
echo "Waiting for background APT processes (unattended-upgrades) to finish..."
echo "====================================================================="
while sudo fuser /var/lib/dpkg/lock-frontend >/dev/null 2>&1; do
    echo "Waiting for other apt-get processes to finish..."
    sleep 5
done

# (1) Install Nvidia driver NATIVELY via Ubuntu
echo "====================================================================="
echo "(1) Installing NVIDIA driver natively..."
echo "====================================================================="
sudo apt-get "${APT_ARGS[@]}" update
sudo DEBIAN_FRONTEND=noninteractive apt-get "${APT_ARGS[@]}" install ubuntu-drivers-common linux-headers-"$(uname -r)"

# Automatically detect the attached GPU and install the recommended driver
sudo DEBIAN_FRONTEND=noninteractive ubuntu-drivers autoinstall

# Explicitly install nvidia-modprobe (GCP VMs often skip this due to minimal image settings)
sudo DEBIAN_FRONTEND=noninteractive apt-get "${APT_ARGS[@]}" install nvidia-modprobe


# (1.5) Hot-load the NVIDIA drivers (Avoids needing a reboot)
echo "====================================================================="
echo "(1.5) Hot-loading the NVIDIA drivers into the kernel..."
echo "====================================================================="

# STOP the daemon that locks older drivers in memory
sudo systemctl stop nvidia-persistenced || true

# FORCEFULLY UNLOAD any older NVIDIA modules that GCP or Ubuntu pre-loaded
sudo modprobe -r nvidia_drm nvidia_modeset nvidia_uvm nvidia || true
sudo modprobe -r nouveau || true

# LOAD the NVIDIA kernel modules we just installed
sudo modprobe nvidia
sudo modprobe nvidia_uvm
sudo modprobe nvidia_drm

# Explicitly create the /dev/nvidia* device nodes that Docker needs
sudo nvidia-modprobe
sudo nvidia-modprobe -u -c=0

# Restart the daemon with the new drivers
sudo systemctl start nvidia-persistenced || true

# Verify driver on the host
echo "====================================================================="
echo "Testing nvidia-smi on Host..."
echo "====================================================================="
sudo nvidia-smi


# (2) Install Docker CE:
echo "====================================================================="
echo "(2) Installing Docker CE..."
echo "====================================================================="
sudo apt-get "${APT_ARGS[@]}" update
sudo apt-get "${APT_ARGS[@]}" install ca-certificates curl gnupg
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get "${APT_ARGS[@]}" update
sudo DEBIAN_FRONTEND=noninteractive NEEDRESTART_MODE=a \
  apt-get "${APT_ARGS[@]}" install \
  docker-ce docker-ce-cli \
  containerd.io \
  docker-buildx-plugin \
  docker-compose-plugin


# (3) Install nvidia docker:
echo "====================================================================="
echo "(3) Installing NVIDIA Container Toolkit..."
echo "====================================================================="
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor --yes -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list > /dev/null

sudo apt-get "${APT_ARGS[@]}" update
sudo DEBIAN_FRONTEND=noninteractive NEEDRESTART_MODE=a \
  apt-get "${APT_ARGS[@]}" install nvidia-container-toolkit

sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker


# (4) Test
echo "====================================================================="
echo "(4) Testing nvidia-smi with the latest official CUDA image..."
echo "====================================================================="
sudo docker run --rm --gpus 1 nvidia/cuda:12.3.2-cudnn9-devel-ubuntu22.04 nvidia-smi

echo "====================================================================="
echo "✅ SUCCESS! Your GPU is fully functional inside Docker."
echo "====================================================================="
