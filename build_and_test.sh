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

# NOLINT
set -eux -o pipefail

source settings.sh

# bazel should have been installed in build-prereq.sh, but the PATH might
# need to be added in this script.
if ! bazel; then
  PATH="$HOME/bin:$PATH"
fi

# Building examples_from_stream.so C++ library. It cannot be built correctly
# with the default bazel setup, so we build it manually.
# examples_from_stream.so is used by call_variants target therefore it has to
# be built before :binaries.
TF_CFLAGS=( $(python3 -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_compile_flags()))') )
TF_LFLAGS=( $(python3 -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_link_flags()))') )

# shellcheck disable=SC2068
g++ -std=c++14 -shared \
        deepvariant/stream_examples_kernel.cc  \
        deepvariant/stream_examples_ops.cc \
        -o deepvariant/examples_from_stream.so \
        -fPIC \
        -l:libtensorflow_framework.so.2  \
        -I. \
        ${TF_CFLAGS[@]} \
        ${TF_LFLAGS[@]} \
        -D_GLIBCXX_USE_CXX11_ABI=1 \
        --std=c++17 \
        -DEIGEN_MAX_ALIGN_BYTES=64 \
        -O2

# Run all deepvariant tests.  Take bazel options from args, if any.
# Note: If running with GPU, tests must be executed serially due to a GPU
# contention issue.
if [[ "${DV_GPU_BUILD:-0}" = "1" ]]; then
  bazel test -c opt --local_test_jobs=1 ${DV_COPT_FLAGS} "$@" \
    deepvariant/...
  # GPU tests are commented out for now.
  # Because they seem to be all filtered out, and as a result causing an error.
  # See internal#comment5.
  # TODO: Uncomment this once it's resolved.
  # bazel test -c opt --local_test_jobs=1 ${DV_COPT_FLAGS} "$@" \
  #   deepvariant:gpu_tests
else
  # Running parallel tests on CPU.
  bazel test -c opt ${DV_COPT_FLAGS} "$@" deepvariant/...
fi

# Build the binary.
./build_release_binaries.sh

echo 'Expect a usage message:'
(python3 bazel-out/k8-opt/bin/deepvariant/call_variants.zip --help || : ) | grep '/call_variants.py:'

