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

# Sets up GPU drivers and pulls the specified Docker image.
#
# Arguments:
#   vm_project: The Google Cloud project ID for the virtual machine.
#   host: The name of the virtual machine instance.
#   zone: The Google Cloud zone where the VM is located.
#   test_script: The name of the test script being run.
#   docker_source: The base name of the repository (e.g., google/deepvariant).
#   bin_version: The version tag (e.g., 1.6.0).
#   use_gpu: Boolean ('true' or 'false') indicating if a GPU is being used.
function setup_gpu_and_pull_docker_image() {
  local vm_project="$1"
  local host="$2"
  local zone="$3"
  local test_script="$4"
  local docker_source="$5"
  local bin_version="$6"
  local use_gpu="$7"

  echo "Waiting for NVIDIA driver to be installed on '${host}'..."
  echo "This can take 2-3 minutes. The script will check every 10 seconds."

  while ! gcloud compute ssh --project="${vm_project}" "${host}" --zone "${zone}" --command="nvidia-smi" > /dev/null 2>&1; do
    sleep 10
  done

  echo "✅ GPU drivers are active! The instance is ready to use."

  local IMAGE
  if [[ "${test_script}" == "inference_deepvariant.sh" ]]; then
    IMAGE=$(get_deepvariant_docker_image_name "${docker_source}" "${bin_version}" "${use_gpu}")
  elif [[ "${test_script}" == "inference_deepsomatic.sh" ]]; then
    IMAGE=$(get_deepsomatic_docker_image_name "${docker_source}" "${bin_version}" "${use_gpu}")
  elif [[ "${test_script}" == "inference_deetrio.sh" ]]; then
    IMAGE=$(get_deeptrio_docker_image_name "${docker_source}" "${bin_version}" "${use_gpu}")
  else
    echo "ERROR: Unknown test script '${test_script}'." >&2
    exit 1
  fi
  local -r DONE_FILE="/tmp/docker_pull_is_done"
  local -r FAIL_FILE="/tmp/docker_pull_is_failed"
  local -r LOG_FILE="/tmp/docker_pull.log"
  local -a GCLOUD_SSH_ARGS_SYNC=("--ssh-flag=-o ServerAliveInterval=600")

  echo "Kicking off background Docker image pull on ${host}..."
  echo "Logs for the pull will be available at ${host}:${LOG_FILE}"

  gcloud compute ssh --project="${vm_project}" "${host}" --zone="${zone}" --command="rm -f ${DONE_FILE} ${FAIL_FILE} ${LOG_FILE}"

  # NOTE: I was hoping that "ServerAliveInterval" would prevent the SSH connection from timing out, but it does not work.
  # So, I still had to use the nohup command.
  local REMOTE_COMMAND="nohup sh -c '(gcloud auth print-access-token | sudo docker login -u oauth2accesstoken --password-stdin https://gcr.io && sudo docker pull ${IMAGE} && touch ${DONE_FILE}) || touch ${FAIL_FILE}' > ${LOG_FILE} 2>&1 &"

  gcloud compute ssh --project="${vm_project}" "${host}" --zone="${zone}" "${GCLOUD_SSH_ARGS_SYNC[@]}" --command="${REMOTE_COMMAND}"

  echo "Waiting for remote Docker pull to complete. This can take several minutes."
  echo "You can check progress with: gcloud compute ssh ${host} --zone ${zone} --command=\"tail -f ${LOG_FILE}\""

  local counter=0
  while ! gcloud compute ssh --project="${vm_project}" "${host}" --zone="${zone}" "${GCLOUD_SSH_ARGS_SYNC[@]}" --command="test -f ${DONE_FILE} || test -f ${FAIL_FILE}" 2>/dev/null; do
    echo "Waiting for Docker to install on ${host}. Sleep for 30 seconds."
    sleep 30
    ((counter+=30))
    if [[ $counter -gt 900 ]]; then
      echo "You've waited more than 900 seconds. You might want to ctrl-c, and directly run 'gcloud compute ssh --project=${vm_project} --zone=${zone} ${host}' to check."
      exit 1
    fi
  done

  if gcloud compute ssh --project="${vm_project}" "${host}" --zone="${zone}" "${GCLOUD_SSH_ARGS_SYNC[@]}" --command="test -f ${FAIL_FILE}" 2>/dev/null; then
    echo "ERROR: Docker pull failed on remote host '${host}'." >&2
    echo "Displaying the contents of the log file (${LOG_FILE}):" >&2
    gcloud compute ssh --project="${vm_project}" "${host}" --zone="${zone}" --command="cat ${LOG_FILE}"
    exit 1
  else
    echo "✅ Docker image pull completed successfully."
  fi

  gcloud compute ssh --project="${vm_project}" "${host}" --zone="${zone}" "${GCLOUD_SSH_ARGS_SYNC[@]}" --command="rm -f ${DONE_FILE}"
}
