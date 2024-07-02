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

set -euo pipefail


USAGE=$'
Example usage:
inference_deepvariant.sh --model_preset "WGS" --docker_build true --use_gpu true

Flags:
--docker_build (true|false)  Whether to build docker image. (default: false)
--dry_run (true|false)  If true, print out the main commands instead of running. (default: false)
--use_gpu (true|false)   Whether to use GPU when running case study. Make sure to specify vm_zone that is equipped with GPUs. (default: false)
--docker_source Where to pull the Docker image from. Default: google/deepvariant.
--bin_version Version of DeepVariant model to use
--customized_model Path to checkpoint directory containing model checkpoint.
--customized_small_model Path to checkpoint directory containing small model checkpoint.
--regions Regions passed into both variant calling and hap.py.
--make_examples_extra_args Flags for make_examples, specified as "flag1=param1,flag2=param2".
--call_variants_extra_args Flags for call_variants, specified as "flag1=param1,flag2=param2".
--postprocess_variants_extra_args Flags for postprocess_variants, specified as "flag1=param1,flag2=param2".
--model_preset Preset case study to run: WGS, WES, PACBIO, ONT_R104, ONT_R104_DUPLEX_CHR20, or HYBRID_PACBIO_ILLUMINA. ONT_R104_DUPLEX_CHR20 will use the ONT_R104 model_type.
--par_regions_bed Path to BED containing Human Pseudoautosomal Region (PAR) regions. This is used in postprocess_variants. We separate it out as a flag because we need to copy data from gs://.
--population_vcfs Path to VCFs containing population allele frequencies. Use wildcard pattern.
--proposed_variants Path to VCF containing proposed variants. In make_examples_extra_args, you must also specify variant_caller=vcf_candidate_importer but not proposed_variants.
--save_intermediate_results (true|false) If True, keep intermediate outputs from make_examples and call_variants.
--skip_happy (true|false) If True, skip the hap.py evaluation.
--report_title Optional title for reports (VCF stats report and make_examples runtime report).
# TODO: It might be worth considering having a separate inference_pangenome_aware_deepvariant.sh script.
# Because there are specific optional flags, such as --pangenome that could be added.
# However, right now it seems most convenient for me to re-use inference_deepvariant.sh
# for my experiment.
--main_binary_name (run_deepvariant|run_pangenome_aware_deepvariant)  Default is run_deepvariant. If using the pangenome-aware DeepVariant Docker image, use run_pangenome_aware_deepvariant.

If model_preset is not specified, the below flags are required:
--model_type Type of DeepVariant model to run (WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA)
--ref Path to GCP bucket containing ref file (.fa)
--bam Path to GCP bucket containing BAM
--truth_vcf Path to GCP bucket containing truth VCF
--truth_bed Path to GCP bucket containing truth BED
--capture_bed Path to GCP bucket containing captured file (only needed for WES model_type)

Note: All paths to dataset must be of the form "gs://..."
'

# Specify default values.
# Booleans; sorted alphabetically.
BUILD_DOCKER=false
DRY_RUN=false
USE_GPU=false
SAVE_INTERMEDIATE_RESULTS=false
SKIP_HAPPY=false
# Strings; sorted alphabetically.
BAM=""
DOCKER_SOURCE="google/deepvariant"
BIN_VERSION="1.6.1"
CALL_VARIANTS_ARGS=""
CAPTURE_BED=""
CUSTOMIZED_MODEL=""
CUSTOMIZED_SMALL_MODEL=""
MAIN_BINARY_NAME="run_deepvariant"
MAKE_EXAMPLES_ARGS=""
MODEL_PRESET=""
MODEL_TYPE=""
PAR_REGIONS_BED=""
POPULATION_VCFS=""
POSTPROCESS_VARIANTS_ARGS=""
PROPOSED_VARIANTS=""
REF=""
REGIONS=""
REPORT_TITLE=""
TRUTH_BED=""
TRUTH_VCF=""

while (( "$#" )); do
  case "$1" in
    --docker_build)
      BUILD_DOCKER="$2"
      if [[ ${BUILD_DOCKER} != "true" ]] && [[ ${BUILD_DOCKER} != "false" ]]; then
        echo "Error: --docker_build needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --dry_run)
      DRY_RUN="$2"
      if [[ ${DRY_RUN} != "true" ]] && [[ ${DRY_RUN} != "false" ]]; then
        echo "Error: --dry_run needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --use_gpu)
      USE_GPU="$2"
      if [[ ${USE_GPU} != "true" ]] && [[ ${USE_GPU} != "false" ]]; then
        echo "Error: --use_gpu needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --save_intermediate_results)
      SAVE_INTERMEDIATE_RESULTS="$2"
      if [[ "${SAVE_INTERMEDIATE_RESULTS}" != "true" ]] && [[ "${SAVE_INTERMEDIATE_RESULTS}" != "false" ]]; then
        echo "Error: --save_intermediate_results needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --skip_happy)
      SKIP_HAPPY="$2"
      if [[ "${SKIP_HAPPY}" != "true" ]] && [[ "${SKIP_HAPPY}" != "false" ]]; then
        echo "Error: --SKIP_HAPPY needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --regions)
      REGIONS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --customized_model)
      CUSTOMIZED_MODEL="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --main_binary_name)
      MAIN_BINARY_NAME="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --make_examples_extra_args)
      MAKE_EXAMPLES_ARGS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --call_variants_extra_args)
      CALL_VARIANTS_ARGS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --postprocess_variants_extra_args)
      POSTPROCESS_VARIANTS_ARGS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --report_title)
      REPORT_TITLE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --model_preset)
      MODEL_PRESET="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --model_type)
      MODEL_TYPE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --docker_source)
      DOCKER_SOURCE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bin_version)
      BIN_VERSION="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --ref)
      REF="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bam)
      BAM="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --truth_vcf)
      TRUTH_VCF="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --truth_bed)
      TRUTH_BED="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --capture_bed)
      CAPTURE_BED="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --par_regions_bed)
      PAR_REGIONS_BED="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --population_vcfs)
      POPULATION_VCFS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --proposed_variants)
      PROPOSED_VARIANTS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --customized_small_model)
      CUSTOMIZED_SMALL_MODEL="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --help)
      echo "$USAGE" >&2
      exit 1
      ;;
    -*|--*=) # other flags not supported
      echo "Error: unrecognized flag $1" >&2
      echo "Run with --help to see usage." >&2
      exit 1
      ;;
    *)
      echo "Error: unrecognized extra args $1" >&2
      echo "Run with --help to see usage." >&2
      exit 1
      ;;
  esac
done

## Presets
# These settings specify the commonly run case studies
GCS_DATA_DIR="gs://deepvariant"  # Or you can use "https://storage.googleapis.com/deepvariant"
BASE="${HOME}/custom-case-study"

declare -a extra_args
declare -a happy_args
declare -a docker_args

if [[ "${MODEL_PRESET}" = "PACBIO" ]]; then
  MODEL_TYPE="PACBIO"
  BASE="${HOME}/pacbio-case-study"

  REF="${REF:=${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
  BAM="${BAM:=${GCS_DATA_DIR}/pacbio-case-study-testdata/HG003.pfda_challenge.35x.grch38.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
elif [[ "${MODEL_PRESET}" = "ONT_R104" ]]; then
  MODEL_TYPE="ONT_R104"
  BASE="${HOME}/ont-case-study"

  REF="${REF:=${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
  BAM="${BAM:=${GCS_DATA_DIR}/ont-case-study-testdata/HG003_R104_sup_merged.80x.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
elif [[ "${MODEL_PRESET}" = "ONT_R104_DUPLEX_CHR20" ]]; then
  MODEL_TYPE="ONT_R104"
  BASE="${HOME}/ont-duplex-case-study"

  REF="${REF:=${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
  BAM="${BAM:=${GCS_DATA_DIR}/ont-case-study-testdata/HG002_R1041_Duplex_all_Dorado_v0.1.1_400bps_pass_2_GRCh38.chr20.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/case-study-testdata/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/case-study-testdata/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
  REGIONS="chr20"
elif [[ "${MODEL_PRESET}" = "WGS" ]]; then
  MODEL_TYPE="WGS"
  BASE="${HOME}/wgs-case-study"

  REF="${REF:=${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
  if [[ "${MAIN_BINARY_NAME}" = "run_pangenome_aware_deepvariant" ]]; then
    echo "Use VG BAM for pangenome-aware DeepVariant."
    BAM="${BAM:=gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam}"
    echo "Add VG default make_examples args for pangenome-aware DeepVariant."
    MAKE_EXAMPLES_ARGS="${MAKE_EXAMPLES_ARGS:+$MAKE_EXAMPLES_ARGS,}min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true"
  else
    BAM="${BAM:=${GCS_DATA_DIR}/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam}"
  fi
elif [[ "${MODEL_PRESET}" = "WES" ]]; then
  MODEL_TYPE="WES"
  BASE="${HOME}/wes-case-study"

  REF="${REF:=${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
  BAM="${BAM:=${GCS_DATA_DIR}/exome-case-study-testdata/HG003.novaseq.wes_idt.100x.dedup.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
  CAPTURE_BED="${CAPTURE_BED:=${GCS_DATA_DIR}/exome-case-study-testdata/idt_capture_novogene.grch38.bed}"
elif [[ "${MODEL_PRESET}" = "HYBRID_PACBIO_ILLUMINA" ]]; then
  MODEL_TYPE="HYBRID_PACBIO_ILLUMINA"
  BASE="${HOME}/hybrid-pacbio-illumina-case-study"

  REF="${REF:=${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
  BAM="${BAM:=${GCS_DATA_DIR}/hybrid-case-study-testdata/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
else
  if [[ -n "${MODEL_PRESET}" ]]; then
    echo "Error: --model_preset must be one of WGS, WES, PACBIO, HYBRID_PACBIO_ILLUMINA." >&2
    exit 1
  fi

fi

INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="deepvariant.output.vcf.gz"
OUTPUT_GVCF="deepvariant.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

if [[ "${MODEL_TYPE}" = "WES" ]]; then
  if [[ -n "${REGIONS}" ]]; then
    echo "Error: --regions is not used with model_type WES. Please use --capture_bed." >&2
    exit 1
  fi
  extra_args+=( --regions "/input/$(basename $CAPTURE_BED)")
fi

# CAPTURE_BED can be specified even when it's not WES.
# But for now, we only use this in hap.py if it's not WES.
if [[ -n "${CAPTURE_BED}" ]]; then
  happy_args+=( -T "${INPUT_DIR}/$(basename $CAPTURE_BED)")
fi

if [[ "${SAVE_INTERMEDIATE_RESULTS}" == "true" ]]; then
  extra_args+=( --intermediate_results_dir "/output/intermediate_results_dir")
fi

echo "========================="
echo "# Booleans; sorted alphabetically."
echo "BUILD_DOCKER: ${BUILD_DOCKER}"
echo "DRY_RUN: ${DRY_RUN}"
echo "USE_GPU: ${USE_GPU}"
echo "SAVE_INTERMEDIATE_RESULTS: ${SAVE_INTERMEDIATE_RESULTS}"
echo "SKIP_HAPPY: ${SKIP_HAPPY}"
echo "# Strings; sorted alphabetically."
echo "BAM: ${BAM}"
echo "DOCKER_SOURCE: ${DOCKER_SOURCE}"
echo "BIN_VERSION: ${BIN_VERSION}"
echo "CALL_VARIANTS_ARGS: ${CALL_VARIANTS_ARGS}"
echo "CAPTURE_BED: ${CAPTURE_BED}"
echo "CUSTOMIZED_MODEL: ${CUSTOMIZED_MODEL}"
echo "CUSTOMIZED_SMALL_MODEL: ${CUSTOMIZED_SMALL_MODEL}"
echo "MAIN_BINARY_NAME: ${MAIN_BINARY_NAME}"
echo "MAKE_EXAMPLES_ARGS: ${MAKE_EXAMPLES_ARGS}"
echo "MODEL_PRESET: ${MODEL_PRESET}"
echo "MODEL_TYPE: ${MODEL_TYPE}"
echo "PAR_REGIONS_BED: ${PAR_REGIONS_BED}"
echo "POPULATION_VCFS: ${POPULATION_VCFS}"
echo "POSTPROCESS_VARIANTS_ARGS: ${POSTPROCESS_VARIANTS_ARGS}"
echo "PROPOSED_VARIANTS: ${PROPOSED_VARIANTS}"
echo "REF: ${REF}"
echo "REGIONS: ${REGIONS}"
echo "REPORT_TITLE: ${REPORT_TITLE}"
echo "TRUTH_BED: ${TRUTH_BED}"
echo "TRUTH_VCF: ${TRUTH_VCF}"
echo "========================="

function run() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    # Prints out command to stdout and a [DRY RUN] tag to stderr.
    # This allows the users to use the dry_run mode to get a list of
    # executable commands by redirecting stdout to a file.
    1>&2 printf "[DRY RUN] " && echo "$*"
  else
    echo "$*"
    eval "$@"
  fi
}

function copy_gs_or_http_file() {
  if [[ "$1" == http* ]]; then
    if curl --output /dev/null --silent --head --fail "$1"; then
      run echo "Copying from \"$1\" to \"$2\""
      run aria2c -c -x10 -s10 "$1" -d "$2"
    else
      run echo "File $1 does not exist. Skip copying."
    fi
  elif [[ "$1" == gs://* ]]; then
    status=0
    gsutil -q stat "$1" || status=1
    if [[ $status == 0 ]]; then
      run echo "Copying from \"$1\" to \"$2\""
      run gcloud storage cp "$1" "$2"
    else
      run echo "File $1 does not exist. Skip copying."
    fi
  else
    echo "Unrecognized file format: $1" >&2
    exit 1
  fi
}

function copy_correct_index_file() {
  BAM="$1"
  INPUT_DIR="$2"
  # Index files have two acceptable naming patterns. We explicitly check for
  # both since we cannot use wildcard paths with http files.
  if [[ "${BAM}" == *".cram" ]]; then
    copy_gs_or_http_file "${BAM%.cram}.crai" "${INPUT_DIR}"
    copy_gs_or_http_file "${BAM}.crai" "${INPUT_DIR}"
  else
    copy_gs_or_http_file "${BAM%.bam}.bai" "${INPUT_DIR}"
    copy_gs_or_http_file "${BAM}.bai" "${INPUT_DIR}"
  fi
}

function copy_data() {
  copy_gs_or_http_file "${TRUTH_BED}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF}.tbi" "${INPUT_DIR}"
  copy_gs_or_http_file "${BAM}" "${INPUT_DIR}"
  copy_correct_index_file "${BAM}" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz.fai" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz.gzi" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gzi" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.fai" "${INPUT_DIR}"
  if [[ -n "${CAPTURE_BED}" ]]; then
    copy_gs_or_http_file "${CAPTURE_BED}" "${INPUT_DIR}"
  fi
  if [[ -n "${PROPOSED_VARIANTS}" ]]; then
    copy_gs_or_http_file "${PROPOSED_VARIANTS}" "${INPUT_DIR}"
    copy_gs_or_http_file "${PROPOSED_VARIANTS}.tbi" "${INPUT_DIR}"
  fi
  if [[ -n "${REGIONS}" ]]; then
    if [[ "${REGIONS}" = http* ]] || [[ "${REGIONS}" = gs://* ]]; then
      copy_gs_or_http_file "${REGIONS}" "${INPUT_DIR}"
    fi
  fi
  if [[ -n "${PAR_REGIONS_BED}" ]]; then
    copy_gs_or_http_file "${PAR_REGIONS_BED}" "${INPUT_DIR}"
  fi
  if [[ -n "${POPULATION_VCFS}" ]]; then
    copy_gs_or_http_file "${POPULATION_VCFS}" "${INPUT_DIR}"
    copy_gs_or_http_file "${POPULATION_VCFS}.tbi" "${INPUT_DIR}"
  fi
}

function setup_test() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    return
  fi

  ## Create local directory structure
  mkdir -p "${OUTPUT_DIR}"
  mkdir -p "${INPUT_DIR}"
  mkdir -p "${LOG_DIR}"

  ## Download extra packages
  # Install aria2 to download data files.
  sudo apt-get -qq -y update
  sudo apt-get -qq -y install aria2

  if ! hash docker 2>/dev/null; then
    echo "'docker' was not found in PATH. Installing docker..."
    # Install docker using instructions on:
    # https://docs.docker.com/install/linux/docker-ce/ubuntu/
    sudo apt-get -qq -y install \
      apt-transport-https \
      ca-certificates \
      curl \
      gnupg-agent \
      software-properties-common
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    sudo add-apt-repository \
      "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
      $(lsb_release -cs) \
      stable"
    sudo apt-get -qq -y update
    sudo apt-get -qq -y install docker-ce
  fi
}

function get_docker_image() {
  DOCKERFILE_NAME="Dockerfile"
  if [[ "${MAIN_BINARY_NAME}" = "run_pangenome_aware_deepvariant" ]]; then
    DOCKERFILE_NAME="Dockerfile.pangenome_aware_deepvariant"
  fi
  if [[ "${BUILD_DOCKER}" = true ]]; then
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="deepvariant_gpu:latest"
      run "sudo docker build \
        -f ${DOCKERFILE_NAME} \
        --build-arg=FROM_IMAGE=nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04 \
        --build-arg=DV_GPU_BUILD=1 -t deepvariant_gpu ."
      run echo "Done building GPU Docker image ${IMAGE}."
      docker_args+=( --gpus 1 )
    else
      IMAGE="deepvariant:latest"
      # Building twice in case the first one times out.
      run "sudo docker build -f ${DOCKERFILE_NAME} -t deepvariant . || \
        (sleep 5 ; sudo docker build -f ${DOCKERFILE_NAME} -t deepvariant .)"
      run echo "Done building Docker image ${IMAGE}."
    fi

  else
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="${DOCKER_SOURCE}:${BIN_VERSION}-gpu"
      # shellcheck disable=SC2027
      # shellcheck disable=SC2086
      run "sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")"
      docker_args+=( --gpus 1 )
    else
      IMAGE="${DOCKER_SOURCE}:${BIN_VERSION}"
      # shellcheck disable=SC2027
      # shellcheck disable=SC2086
      run "sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")"
    fi
  fi
  if [[ "${USE_GPU}" = true ]]; then
    # shellcheck disable=SC2027
    # shellcheck disable=SC2086
    run "sudo docker run --gpus 1 "${IMAGE}" \
      python3 -c 'import tensorflow as tf; \
      print(\"is_gpu_available=\" + str(tf.test.is_gpu_available()))'"
    # shellcheck disable=SC2027
    # shellcheck disable=SC2086
    run "sudo docker run --gpus 1 "${IMAGE}" \
      python3 -c 'import tensorflow as tf; \
      tf.test.is_gpu_available() or exit(1)' \
      2> /dev/null || exit 1"
  fi
}

function setup_args() {
  if [[ -n "${CUSTOMIZED_MODEL}" ]]; then
    run echo "Copy from gs:// path ${CUSTOMIZED_MODEL} to ${INPUT_DIR}/"
    # Check if it's saved Model
    saved_modelpath=${CUSTOMIZED_MODEL}/saved_model.pb
    using_saved_model=$(gsutil -q stat "$saved_modelpath" || echo 1)
    if [[ $using_saved_model != 1 ]]; then
      echo "Using saved model"
      run mkdir -p "${INPUT_DIR}/savedmodel"
      run gcloud storage cp -R "${CUSTOMIZED_MODEL}"/'*' "${INPUT_DIR}"/savedmodel/
      run gcloud storage cp "${CUSTOMIZED_MODEL}"/example_info.json "${INPUT_DIR}"/savedmodel/example_info.json
      extra_args+=( --customized_model "/input/savedmodel")
    else
      echo "Using checkpoint"
      run gcloud storage cp "${CUSTOMIZED_MODEL}".data-00000-of-00001 "${INPUT_DIR}/model.ckpt.data-00000-of-00001"
      run gcloud storage cp "${CUSTOMIZED_MODEL}".index "${INPUT_DIR}/model.ckpt.index"
      # Starting from v1.7.0, example_info.json is required.
      CUSTOMIZED_MODEL_DIR="$(dirname "${CUSTOMIZED_MODEL}")"
      run "gcloud storage cp ${CUSTOMIZED_MODEL_DIR}/example_info.json ${INPUT_DIR}/example_info.json"
      extra_args+=( --customized_model "/input/model.ckpt")
    fi
  else
    run echo "No custom model specified."
  fi
  if [[ -n "${CUSTOMIZED_SMALL_MODEL}" ]]; then
    echo "Using customized small model"
    run mkdir -p "${INPUT_DIR}/smallmodel"
    run echo "Copy from gs:// path ${CUSTOMIZED_SMALL_MODEL} to ${INPUT_DIR}/smallmodel"
    run gcloud storage cp -R "${CUSTOMIZED_SMALL_MODEL}"/'*' "${INPUT_DIR}"/smallmodel/
    extra_args+=( --customized_small_model "/input/smallmodel")
  fi
  if [[ -n "${POPULATION_VCFS}" ]]; then
    MAKE_EXAMPLES_ARGS="${MAKE_EXAMPLES_ARGS:+${MAKE_EXAMPLES_ARGS},}population_vcfs=/input/$(basename "$POPULATION_VCFS")"
  fi
  if [[ -n "${MAKE_EXAMPLES_ARGS}" ]]; then
    # In order to use proposed variants, we have to pass vcf_candidate_importer
    # to make_examples_extra_args, so we know that we will enter this if
    # statement.
    if [[ -n "${PROPOSED_VARIANTS}" ]]; then
      MAKE_EXAMPLES_ARGS="${MAKE_EXAMPLES_ARGS},proposed_variants=/input/$(basename "$PROPOSED_VARIANTS")"
    fi
    extra_args+=( --make_examples_extra_args "\"${MAKE_EXAMPLES_ARGS}\"")
  fi
  if [[ -n "${CALL_VARIANTS_ARGS}" ]]; then
    extra_args+=( --call_variants_extra_args "\"${CALL_VARIANTS_ARGS}\"")
  fi
  if [[ -n "${POSTPROCESS_VARIANTS_ARGS}" ]]; then
    if [[ -n "${PAR_REGIONS_BED}" ]]; then
      # Note: currently this following line won't be run if
      # POSTPROCESS_VARIANTS_ARGS is empty, which means `par_regions_bed` won't
      # be set if POSTPROCESS_VARIANTS_ARGS is empty. That is the intended
      # behavior because --par_regions_bed isn't useful unless --haploid_contigs
      # is set.
      POSTPROCESS_VARIANTS_ARGS="${POSTPROCESS_VARIANTS_ARGS},par_regions_bed=/input/$(basename "$PAR_REGIONS_BED")"
    fi
    extra_args+=( --postprocess_variants_extra_args "\"${POSTPROCESS_VARIANTS_ARGS}"\")
  fi
  if [[ -n "${REPORT_TITLE}" ]]; then
    extra_args+=( --report_title "${REPORT_TITLE}")
  fi
  if [[ -n "${REGIONS}" ]]; then
    if [[ "${REGIONS}" = http* ]] || [[ "${REGIONS}" = gs://* ]]; then
      extra_args+=( --regions "/input/$(basename $REGIONS)")
      happy_args+=( -T "${INPUT_DIR}/$(basename $REGIONS)")
    else
      extra_args+=( --regions "${REGIONS}")
      happy_args+=( -l "${REGIONS}")
    fi
  fi
  # If you're running an older version (before 1.2) that doesn't have this flag,
  # you'll need to comment out this line.
  extra_args+=( --runtime_report )
}

function run_deepvariant_with_docker() {
  run echo "Run DeepVariant..."
  run echo "using IMAGE=${IMAGE}"
  # shellcheck disable=SC2027
  # shellcheck disable=SC2046
  # shellcheck disable=SC2068
  # shellcheck disable=SC2086
  # shellcheck disable=SC2145
  run "(time (sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    ${docker_args[@]-} \
    "${IMAGE}" \
    /opt/deepvariant/bin/${MAIN_BINARY_NAME} \
    --model_type="${MODEL_TYPE}" \
    --ref="/input/$(basename $REF).gz" \
    --reads="/input/$(basename $BAM)" \
    --output_vcf="/output/${OUTPUT_VCF}" \
    --output_gvcf="/output/${OUTPUT_GVCF}" \
    --num_shards "$(nproc)" \
    --logging_dir="/output/logs" \
    "${extra_args[@]-}" && \
  echo "Done.")) 2>&1 | tee "${LOG_DIR}/deepvariant_runtime.log""
  echo
}


function run_happy() {
  ## Evaluation: run hap.py
  run echo "Start evaluation with hap.py..."
  UNCOMPRESSED_REF="${INPUT_DIR}/$(basename $REF)"
  # hap.py cannot read the compressed fa, so uncompress
  # into a writable directory. Index file was downloaded earlier.
  # shellcheck disable=SC2027
  # shellcheck disable=SC2046
  # shellcheck disable=SC2086
  run "zcat <"${INPUT_DIR}/$(basename $REF).gz" >"${UNCOMPRESSED_REF}""

  HAPPY_VERSION="v0.3.12"
  # Pulling twice in case the first one times out.
  run "sudo docker pull jmcdani20/hap.py:${HAPPY_VERSION} || \
    (sleep 5 ; sudo docker pull jmcdani20/hap.py:${HAPPY_VERSION})"
  # shellcheck disable=SC2027
  # shellcheck disable=SC2046
  # shellcheck disable=SC2086
  # shellcheck disable=SC2145
  run "( sudo docker run -i \
  -v "${INPUT_DIR}:${INPUT_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  jmcdani20/hap.py:${HAPPY_VERSION} /opt/hap.py/bin/hap.py \
    "${INPUT_DIR}/$(basename $TRUTH_VCF)" \
    "${OUTPUT_DIR}/${OUTPUT_VCF}" \
    -f "${INPUT_DIR}/$(basename $TRUTH_BED)" \
    -r "${UNCOMPRESSED_REF}" \
    -o "${OUTPUT_DIR}/happy.output" \
    --engine=vcfeval \
    --pass-only \
    ${happy_args[@]-} \
  ) 2>&1 | tee "${LOG_DIR}/happy.log""
  if [[ "${DRY_RUN}" != "true" ]]; then
    echo "${HAPPY_VERSION}" > "${LOG_DIR}/happy_version.log"
  fi
  run echo "Done."
}

function main() {
  run echo 'Starting the test...'

  setup_test
  copy_data
  if [[ ${DOCKER_SOURCE} =~ ^gcr.io ]]; then
    run "gcloud auth print-access-token | sudo docker login -u oauth2accesstoken --password-stdin https://gcr.io"
  fi
  get_docker_image
  setup_args
  run_deepvariant_with_docker
  if [[ "${SKIP_HAPPY}" == "false" ]]; then
    if [[ "${DRY_RUN}" == "true" ]]; then
      run_happy
    else
      run_happy 2>&1 | tee "${LOG_DIR}/happy.log"
    fi
  fi
}

main "$@"
