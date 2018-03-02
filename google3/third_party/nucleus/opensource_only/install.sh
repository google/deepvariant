#!/bin/bash

# Usage:  ./install.sh
#
# This script installs all the packages required to build Nucleus, and
# then builds Nucleus.
#
# This script will run as-is on Ubuntu 14, Ubuntu 16, and Debian 9 systems.
# On all other systems, you will need to first install CLIF by following the
# instructions at https://github.com/google/clif#installation
#
# We also assume that apt-get is already installed and available.

function note_build_stage {
  echo "========== [$(date)] Stage '${1}' starting"
}

# Update package list
################################################################################
note_build_stage "Update package list"
sudo -H apt-get -qq -y update

# Install generic dependencies
################################################################################
note_build_stage "Update misc. dependencies"
sudo -H apt-get -y install pkg-config zip g++ zlib1g-dev unzip curl git

# Install htslib dependencies
################################################################################
note_build_stage "Install htslib dependencies"
sudo -H apt-get -y install libssl-dev libcurl4-openssl-dev liblz-dev libbz2-dev liblzma-dev

# Install pip
################################################################################
note_build_stage "Update pip"
sudo -H apt-get -y install python-dev python-pip python-wheel
sudo -H install --upgrade pip

# Install python packages used by Nucleus
################################################################################
sudo -H pip install contextlib2
sudo -H pip install intervaltree
sudo -H pip install 'mock>=2.0.0'
sudo -H pip install 'numpy==1.12'
sudo -H pip install 'scipy==1.0'

# redacted
sudo -H pip install 'requests>=2.18'
sudo -H pip install 'oauth2client>=4.0.0'

# Install Java (required for Bazel)
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

# Install Bazel
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

bazel_ver="0.9.0"
if
  v=$(bazel --bazelrc=/dev/null --nomaster_bazelrc version) &&
  echo "$v" | awk -v b="$bazel_ver" '/Build label/ { exit ($3 != b)}'
then
  echo "Bazel $bazel_ver already installed on the machine, not reinstalling"
else
  update_bazel_linux "$bazel_ver"
fi

# Install CLIF
################################################################################

note_build_stage "Install CLIF binary"

if [[ -e /usr/local/clif/bin/pyclif ]];
then
  echo "CLIF already installed."
else
  # Figure out which linux installation we are on to fetch an appropriate
  # version of the pre-built CLIF binary. Note that we only support now Ubuntu
  # 14, Ubuntu 16, and Debian 9.
  case "$(lsb_release -d)" in
    *Ubuntu*16.*.*)  PLATFORM="ubuntu-16" ;;
    *Ubuntu*14.*.*)  PLATFORM="ubuntu-14" ;;
    *Debian*9.*)     PLATFORM="debian" ;;
    *Debian*rodete*) PLATFORM="debian" ;;
    *) echo "CLIF is not installed on this machine and a prebuilt binary is not
unavailable for this platform. Please install CLIF at
https://github.com/google/clif before continuing."
    exit 1
  esac

  PACKAGE_CURL_PATH="https://storage.googleapis.com/deepvariant/packages"
  OSS_CLIF_CURL_ROOT="${PACKAGE_CURL_PATH}/oss_clif"
  OSS_CLIF_PKG="oss_clif.${PLATFORM}.latest.tgz"

  if [[ ! -f "/tmp/${OSS_CLIF_PKG}" ]]; then
    curl "${OSS_CLIF_CURL_ROOT}/${OSS_CLIF_PKG}" > /tmp/${OSS_CLIF_PKG}
  fi

  (cd / && sudo tar xzf "/tmp/${OSS_CLIF_PKG}")
  sudo ldconfig  # Reload shared libraries.
fi

# Download and build TensorFlow
################################################################################
note_build_stage "Download and build TensorFlow"

(cd .. &&
 git clone https://github.com/tensorflow/tensorflow &&
 cd tensorflow &&
 echo | ./configure)

# Build Nucleus
################################################################################
note_build_stage "Building Nucleus"

COPT_FLAGS="--copt=-msse4.1 --copt=-msse4.2 --copt=-mavx --copt=-O3"
bazel build -c opt ${COPT_FLAGS} nucleus:all

bazel build :licenses_zip

# Done!
################################################################################

echo "Installation complete at $(date)!"
