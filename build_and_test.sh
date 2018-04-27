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

# NOLINT
set -eux -o pipefail

source settings.sh

# bazel should have been installed in build-prereq.sh, but the PATH might
# need to be added in this script.
if ! bazel; then
  PATH="$HOME/bin:$PATH"
fi

# Run all deepvariant tests.  Take bazel options from args, if any.
# Note: If running with GPU, tests must be executed serially due to a GPU
# contention issue.
if [[ "${DV_GPU_BUILD:-0}" = "1" ]]; then
  bazel test -c opt --local_test_jobs=1 ${DV_COPT_FLAGS} "$@" \
    deepvariant/...
  # Also run the GPU-specific tests.
  bazel test -c opt --local_test_jobs=1 ${DV_COPT_FLAGS} "$@" \
    deepvariant:gpu_tests
else
  # Running parallel tests on CPU.
  bazel test -c opt ${DV_COPT_FLAGS} "$@" deepvariant/...
fi

# Create some build artifacts.
bazel build -c opt ${DV_COPT_FLAGS} "$@" deepvariant:binaries

# Bundle the licenses.
bazel build :licenses_zip
