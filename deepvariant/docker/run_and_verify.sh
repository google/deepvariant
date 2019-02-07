#!/bin/bash
set -euo pipefail

# Copyright 2018 Google LLC.
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
# Used by Cloud container build, runs DeepVariant dockers on GCP, and checks if
# expected files are generated.

function init() {
  local project_id="$1"
  local image_tag="$2"
  local model="$3"
  local platform="$4"
  local extra_flags=""
  if [[ "${platform}" == "gpu" ]]; then
    extra_flags="--gpu"
  elif [[ "${platform}" == "tpu" ]]; then
    local extra_flags="--tpu --gke_cluster_version 1.10 --gke_cluster_zone us-central1-c"
  fi

  readonly DATE=$(date '+%Y-%m-%d-%H-%M-%S')
  readonly PROJECT_ID="${project_id}"
  readonly WORKING_BUCKET="gs://${PROJECT_ID}"
  readonly ZONES="us-west1-b,us-east1-c"
  readonly OUTDIR="${WORKING_BUCKET}/deepvariant-cbi-${DATE}"
  readonly OUTFILE="${OUTDIR}/output.vcf"
  readonly OUTFILE_GVCF="${OUTDIR}/output.gvcf"
  readonly MODEL="${model}"
  readonly BAM="gs://deepvariant/quickstart-testdata/NA12878_S1.chr20.10_10p1mb.bam"
  readonly REF="gs://deepvariant/quickstart-testdata/ucsc.hg19.chr20.unittest.fasta.gz"
  readonly REGION="chr20:10,000,000-10,010,000"
  readonly RUNNER_IMAGE="gcr.io/${PROJECT_ID}/deepvariant_runner:${image_tag}"
  readonly DEEPVARIANT_IMAGE="gcr.io/${PROJECT_ID}/deepvariant:${image_tag}"
  readonly DEEPVARIANT_IMAGE_GPU="gcr.io/${PROJECT_ID}/deepvariant_gpu:${image_tag}"
  readonly COMMAND="
    ./opt/deepvariant_runner/bin/gcp_deepvariant_runner \
        --project ${PROJECT_ID} \
        --zones ${ZONES} \
        --docker_image ${DEEPVARIANT_IMAGE} \
        --docker_image_gpu ${DEEPVARIANT_IMAGE_GPU} \
        --outfile ${OUTFILE} \
        --gvcf_outfile ${OUTFILE_GVCF} \
        --staging ${OUTDIR} \
        --model ${MODEL} \
        --bam ${BAM} \
        --ref ${REF} \
        --regions chr20:10,000,000-10,010,000 \
        --gcsfuse \
        ${extra_flags}
  "
}

err() {
  echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: ERROR: $*" >&2
  exit 1
}

info() {
  echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: INFO: $*" >&2
}

function print_useful_info() {
  info "Logs and staging files can be found in ${OUTDIR}."
  info "Output VCF file will be generated in ${OUTFILE}."
  info "Output gVCF file will be generated in ${OUTFILE_GVCF}."
}

########################################
# Runs DeepVariant on GCP using pipelines tool.
# Globals:
#   COMMAND
#   ZONES
#   RUNNER_IMAGE
#   OUTDIR
# Arguments:
#   None
# Returns:
#   None
########################################
function run_deepvariant() {
  info "Running DeepVariant ..."
  pipelines \
    --project "${PROJECT_ID}" \
    run \
    --scopes="https://www.googleapis.com/auth/cloud-platform" \
    --zones "${ZONES}" \
    --image "${RUNNER_IMAGE}" \
    --output "${OUTDIR}/runner.log" \
    --attempts 1 \
    --pvm-attempts 0 \
    --command "${COMMAND}"
}

########################################
# Checks if the given GCS file exists and is non-empty.
# Arguments:
#   GCS file path.
# Returns:
#   None
########################################
function check_file() {
  local file_path="$1"
  local file_size="$(gsutil stat "${file_path}" | awk '$1 == "Content-Length:" {print $2}')"
  if [[ -z "${file_size}" ]]; then
    err "${file_path} does not exist."
  elif [[ "${file_size}" == 0 ]]; then
    err "${file_path} is empty. Expected a non-empty file."
  fi
}

########################################
# Checks if expected files are generated.
# Globals:
#   OUTFILE
#   OUTFILE_GVCF
#   EXPECTED_OUTFILE
#   EXPECTED_OUTFILE_GVCF
# Arguments:
#   None
# Returns:
#   None
########################################
function verify() {
  info "Verifying outputs ..."
  check_file "${OUTFILE}"
  check_file "${OUTFILE_GVCF}"
}

function main() {
  if [[ $# -ne 4 ]]; then
    err "Usage: $0 <project_id> <image-tag> <model> <cpu|gpu|tpu>"
  fi
  local project_id="$1"
  local image_tag="$2"
  local model="$3"
  local platform="$4"

  init "${project_id}" "${image_tag}" "${model}" "${platform}"
  print_useful_info
  run_deepvariant
  verify
  info "SUCCESS: All tests passed."
}

main "$@"
