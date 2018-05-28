#!/bin/bash
set -euo pipefail
set -x

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

# A helper script for building python binaries required for the deepvariant
# docker image.

source settings.sh  # Make sure we define DV_COPT_FLAGS used below.
./build-prereq.sh   # implies run-prereq.sh

# For bazel.
PATH="$HOME/bin:$PATH"

# Build all required binaries as python zip files. Note that par executables
# that are using subpar library does not yet work due to c-extensions not being
# supported. Also symlinks do not work in Dockerfile, so copy them explicitly
# to //deepvariant/docker directory.
# shellcheck disable=SC2086
# because DV_COPT_FLAGS contains (potentially) multiple optimization flags but
# is set as an environment variable in setting.sh via build-prereq.sh and
# therefore cannot be an array. So in order to expand correctly we do not
# surround the variable with quotes.
bazel build --build_python_zip -c opt ${DV_COPT_FLAGS} \
    //deepvariant:make_examples \
    //deepvariant:call_variants \
    //deepvariant:postprocess_variants \
    //deepvariant:model_train \
    //deepvariant:model_eval
cp bazel-bin/deepvariant/*.zip deepvariant/docker/
cp run-prereq.sh settings.sh deepvariant/docker/
cp LICENSE AUTHORS deepvariant/docker/
