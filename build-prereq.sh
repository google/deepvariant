#!/bin/bash
set -euo pipefail

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

echo ========== Load config settings.

source settings.sh

################################################################################
# Misc. setup
################################################################################

note_build_stage "Install the runtime packages"

./run-prereq.sh

note_build_stage "Update package list"

sudo -H apt-get -qq -y update

note_build_stage "Install development packages"

sudo -H apt-get -y install pkg-config zip g++ zlib1g-dev unzip curl git

################################################################################
# Java
################################################################################

note_build_stage "Install Java and friends"

# Java is available on Kokoro, so we add this cutout.
if ! java -version 2>&1 | fgrep "1.8"; then
  echo "No Java 8, will install."
  sudo -H apt-get install -y software-properties-common debconf-utils
  # Debian needs authentication.
  # (http://www.webupd8.org/2014/03/how-to-install-oracle-java-8-in-debian.html)
  [[ $(lsb_release -d | grep 'Debian') ]] && \
    sudo -H apt-get install -y gnupg dirmngr && \
    sudo -H apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886
  sudo add-apt-repository -y ppa:webupd8team/java
  sudo -H apt-get -qq -y update
  echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections
  sudo -H apt-get install -y oracle-java8-installer
  sudo -H apt-get -y install ca-certificates-java
  sudo update-ca-certificates -f
else
  echo "Java 8 found, will not reinstall."
fi

################################################################################
# bazel
################################################################################

note_build_stage "Install bazel"

function update_bazel_linux {
  BAZEL_VERSION=$1
  rm -rf ~/bazel
  mkdir ~/bazel

  pushd ~/bazel
  curl -L -O https://github.com/bazelbuild/bazel/releases/download/"${BAZEL_VERSION}"/bazel-"${BAZEL_VERSION}"-installer-linux-x86_64.sh
  chmod +x bazel-*.sh
  ./bazel-"${BAZEL_VERSION}"-installer-linux-x86_64.sh --user
  rm bazel-"${BAZEL_VERSION}"-installer-linux-x86_64.sh
  popd

  PATH="$HOME/bin:$PATH"
}

bazel_ver="0.11.0"
if
  v=$(bazel --bazelrc=/dev/null --nomaster_bazelrc version) &&
  echo "$v" | awk -v b="$bazel_ver" '/Build label/ { exit ($3 != b)}'
then
  echo "Bazel $bazel_ver already installed on the machine, not reinstalling"
else
  update_bazel_linux "$bazel_ver"
fi

################################################################################
# CLIF
################################################################################

note_build_stage "Install CLIF binary"

if [[ -e /usr/local/clif/bin/pyclif ]];
then
  echo "CLIF already installed."
else
  # Figure out which linux installation we are on to fetch an appropriate
  # version of the pre-built CLIF binary. Note that we only support now Ubuntu
  # 14 and 16.
  case "$(lsb_release -d)" in
    *Ubuntu*16.*.*) export DV_PLATFORM="ubuntu-16" ;;
    *Ubuntu*14.*.*) export DV_PLATFORM="ubuntu-14" ;;
    *Debian*9.*)    export DV_PLATFORM="debian" ;;
    *Debian*rodete) export DV_PLATFORM="debian" ;;
    *) echo "CLIF is not installed on this machine and a prebuilt binary is not
unavailable for this platform. Please install CLIF at
https://github.com/google/clif before continuing."
    exit 1
  esac

  OSS_CLIF_CURL_ROOT="${DV_PACKAGE_CURL_PATH}/oss_clif"
  OSS_CLIF_PKG="oss_clif.${DV_PLATFORM}.latest.tgz"

  if [[ ! -f "/tmp/${OSS_CLIF_PKG}" ]]; then
    curl "${OSS_CLIF_CURL_ROOT}/${OSS_CLIF_PKG}" > /tmp/${OSS_CLIF_PKG}
  fi

  (cd / && sudo tar xzf "/tmp/${OSS_CLIF_PKG}")
  sudo ldconfig  # Reload shared libraries.
fi

################################################################################
# TensorFlow
################################################################################

note_build_stage "Download and configure TensorFlow sources"

(cd .. &&
 git clone https://github.com/tensorflow/tensorflow &&
 cd tensorflow &&
 git checkout v1.7.0 &&
 echo | ./configure)

note_build_stage "build-prereq.sh complete"

