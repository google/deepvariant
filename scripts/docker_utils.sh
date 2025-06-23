#!/bin/bash
# Copyright 2025 Google LLC.
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

# Utility functions for Docker operations.

# Constructs the full Docker image name based on source, version, and GPU use.
#
# Arguments:
#   DOCKER_SOURCE: The base name of the repository (e.g., google/deepsomatic).
#   BIN_VERSION: The version tag (e.g., 1.9.0).
#   USE_GPU: Boolean ('true' or 'false') indicating if a GPU image is needed.
# Returns:
#   Prints the constructed Docker image name to stdout.
function get_deepsomatic_docker_image_name() {
  local DOCKER_SOURCE="$1"
  local BIN_VERSION="$2"
  local USE_GPU_ARG="$3"
  local IMAGE_NAME
  local -i use_gpu=0 # Default to false

  if [[ "${DOCKER_SOURCE}" = "gcr.io/google.com/brain-genomics/deepvariant" ]]; then
    IMAGE_NAME="${DOCKER_SOURCE}:deepsomatic-${BIN_VERSION}"
  else
    IMAGE_NAME="${DOCKER_SOURCE}:${BIN_VERSION}"
  fi

  # echo "Debug in get_deepsomatic_docker_image_name: USE_GPU_ARG='${USE_GPU_ARG}'" >&2

  # Convert common true values to 1
  if [[ "${USE_GPU_ARG}" == "1" || "${USE_GPU_ARG}" == "true" || "${USE_GPU_ARG}" == "True" || "${USE_GPU_ARG}" == "TRUE" ]]; then
    use_gpu=1
  fi

  if (( use_gpu )); then
    IMAGE_NAME="${IMAGE_NAME}-gpu"
  fi

  echo "${IMAGE_NAME}"
}

function get_deepvariant_docker_image_name() {
  local DOCKER_SOURCE="$1"
  local BIN_VERSION="$2"
  local USE_GPU_ARG="$3"
  local IMAGE_NAME
  local -i use_gpu=0 # Default to false

  IMAGE_NAME="${DOCKER_SOURCE}:${BIN_VERSION}"

  # echo "Debug in get_deepvariant_docker_image_name: USE_GPU_ARG='${USE_GPU_ARG}'" >&2

  # Convert common true values to 1
  if [[ "${USE_GPU_ARG}" == "1" || "${USE_GPU_ARG}" == "true" || "${USE_GPU_ARG}" == "True" || "${USE_GPU_ARG}" == "TRUE" ]]; then
    use_gpu=1
  fi

  if (( use_gpu )); then
    IMAGE_NAME="${IMAGE_NAME}-gpu"
  fi

  echo "${IMAGE_NAME}"
}

function get_deeptrio_docker_image_name() {
  local DOCKER_SOURCE="$1"
  local BIN_VERSION="$2"
  local USE_GPU_ARG="$3"
  local IMAGE_NAME
  local -i use_gpu=0 # Default to false

  IMAGE_NAME="${DOCKER_SOURCE}:deeptrio-${BIN_VERSION}"

  # echo "Debug in get_deepvariant_docker_image_name: USE_GPU_ARG='${USE_GPU_ARG}'" >&2

  # Convert common true values to 1
  if [[ "${USE_GPU_ARG}" == "1" || "${USE_GPU_ARG}" == "true" || "${USE_GPU_ARG}" == "True" || "${USE_GPU_ARG}" == "TRUE" ]]; then
    use_gpu=1
  fi

  if (( use_gpu )); then
    IMAGE_NAME="${IMAGE_NAME}-gpu"
  fi

  echo "${IMAGE_NAME}"
}
