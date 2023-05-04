#!/bin/bash
# Copyright 2019 Google LLC.
# This script is used to install `singularity` 3.7.0.
# It has been tested on Ubutun 20.04.
# For different Linux distributions and versions, modifications might be needed.
# Installation instructions are from: https://sylabs.io/docs/

set -euo pipefail

sudo apt-get update && sudo apt-get install -y \
  build-essential \
  libssl-dev \
  uuid-dev \
  libgpgme11-dev \
  squashfs-tools \
  libseccomp-dev \
  wget \
  pkg-config \
  git

export VERSION=1.15.6 OS=linux ARCH=amd64
# Downloads the required Go package
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
# Extracts the archive
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz
# Deletes the ``tar`` file
rm go$VERSION.$OS-$ARCH.tar.gz

export VERSION=3.7.0
wget https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
tar -xzf singularity-${VERSION}.tar.gz
pushd singularity
export PATH=/usr/local/go/bin:$PATH

./mconfig
make -C builddir
sudo make -C builddir install

# Returns to the original directory.
popd
