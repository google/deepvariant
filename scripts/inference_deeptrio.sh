#!/bin/bash
# Copyright 2020 Google LLC.

set -euo pipefail


USAGE=$'
Example usage:
inference_deeptrio.sh --model_preset "WGS" --docker_build true --use_gpu true

The only trio supported is: child=HG002, parents=HG003,HG004.

This is because otherwise there are too many flags needed to set VCFs and BEDs for 3 samples.

Flags:
--docker_build (true|false)  Whether to build docker image. (default: false)
--dry_run (true|false)  If true, print out the main commands instead of running. (default: false)
--use_gpu (true|false)   Whether to use GPU when running case study. Make sure to specify vm_zone that is equipped with GPUs. (default: false)
--use_hp_information (true|false) Use to set --use_hp_information. Only set this for PACBIO*.
--bin_version Version of DeepTrio model to use.
--customized_model Path to checkpoint directory containing model checkpoint.
--regions Regions passed into both variant calling and hap.py.
--make_examples_extra_args Flags for make_examples, specified as "flag1=param1,flag2=param2".
--call_variants_extra_args Flags for call_variants, specified as "flag1=param1,flag2=param2".
--postprocess_variants_extra_args Flags for postprocess_variants, specified as "flag1=param1,flag2=param2".
--model_preset Preset case study to run: WGS, WGS_CHR20, WES, PACBIO, or PACBIO_CHR20.
--proposed_variants Path to VCF containing proposed variants. In make_examples_extra_args, you must also specify variant_caller=vcf_candidate_importer but not proposed_variants.
--save_intermediate_results (true|false) If True, keep intermediate outputs from make_examples and call_variants.


If model_preset is not specified, the below flags are required:
--model_type Type of DeepTrio model to run (WGS, WES, PACBIO)
--bam_child Path to bam for HG002 on GCP.
--bam_parent1 Path to bam for HG003 on GCP.
--bam_parent2 Path to bam for HG004 on GCP.

--capture_bed Path to GCP bucket containing captured file (only needed for WES model_type)

Note: All paths to dataset must be of the form "gs://..."
'

# Specify default values.
# Booleans; sorted alphabetically.
BUILD_DOCKER=false
DRY_RUN=false
USE_GPU=false
USE_HP_INFORMATION="unset"  # To distinguish whether this flag is set explicitly or not.
SAVE_INTERMEDIATE_RESULTS=false
# Strings; sorted alphabetically.
BAM=""
BAM_PARENT1=""
BAM_PARENT2=""
BIN_VERSION="1.3.0"
CALL_VARIANTS_ARGS=""
CAPTURE_BED=""
CUSTOMIZED_MODEL=""
MAKE_EXAMPLES_ARGS=""
MODEL_PRESET=""
MODEL_TYPE=""
POSTPROCESS_VARIANTS_ARGS=""
PROPOSED_VARIANTS=""
REGIONS=""

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
    --use_hp_information)
      USE_HP_INFORMATION="$2"
      if [[ "${USE_HP_INFORMATION}" != "true" ]] && [[ "${USE_HP_INFORMATION}" != "false" ]]; then
        echo "Error: --use_hp_information needs to have value (true|false)." >&2
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
    --bin_version)
      BIN_VERSION="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bam_child)
      BAM_CHILD="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bam_parent1)
      BAM_PARENT1="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bam_parent2)
      BAM_PARENT2="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --capture_bed)
      CAPTURE_BED="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --proposed_variants)
      PROPOSED_VARIANTS="$2"
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
GCS_DATA_DIR="https://storage.googleapis.com/deepvariant"
BASE="${HOME}/custom-case-study"

declare -a extra_args
declare -a happy_args
declare -a docker_args

if [[ "${MODEL_PRESET}" = "PACBIO" ]]; then
  MODEL_TYPE="PACBIO"
  BASE="${HOME}/pacbio-case-study"

  BAM_CHILD="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG002.pfda_challenge.grch38.phased.bam"
  BAM_PARENT1="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG003.pfda_challenge.grch38.phased.bam"
  BAM_PARENT2="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG004.pfda_challenge.grch38.phased.bam"
elif [[ "${MODEL_PRESET}" = "PACBIO_CHR20" ]]; then
  MODEL_TYPE="PACBIO"
  BASE="${HOME}/pacbio-case-study"

  BAM_CHILD="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG002.pfda_challenge.grch38.phased.chr20.bam"
  BAM_PARENT1="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG003.pfda_challenge.grch38.phased.chr20.bam"
  BAM_PARENT2="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG004.pfda_challenge.grch38.phased.chr20.bam"

  if [[ -n "${REGIONS}" ]]; then
    echo "For --model_preset=${MODEL_PRESET}, regions will be set to chr20."
  fi
  REGIONS=chr20
elif [[ "${MODEL_PRESET}" = "WGS" ]]; then
  MODEL_TYPE="WGS"
  BASE="${HOME}/wgs-case-study"

  BAM_CHILD="${GCS_DATA_DIR}/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
  BAM_PARENT1="${GCS_DATA_DIR}/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
  BAM_PARENT2="${GCS_DATA_DIR}/case-study-testdata/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
elif [[ "${MODEL_PRESET}" = "WGS_CHR20" ]]; then
  MODEL_TYPE="WGS"
  BASE="${HOME}/wgs-case-study"

  BAM_CHILD="${GCS_DATA_DIR}/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
  BAM_PARENT1="${GCS_DATA_DIR}/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
  BAM_PARENT2="${GCS_DATA_DIR}/case-study-testdata/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"

  if [[ -n "${REGIONS}" ]]; then
    echo "For --model_preset=${MODEL_PRESET}, regions will be set to chr20."
  fi
  REGIONS=chr20
elif [[ "${MODEL_PRESET}" = "WES" ]]; then
  MODEL_TYPE="WES"
  BASE="${HOME}/exome-case-study"

  BAM_CHILD="${GCS_DATA_DIR}/exome-case-study-testdata/HG002.novaseq.wes_idt.100x.dedup.bam"
  BAM_PARENT1="${GCS_DATA_DIR}/exome-case-study-testdata/HG003.novaseq.wes_idt.100x.dedup.bam"
  BAM_PARENT2="${GCS_DATA_DIR}/exome-case-study-testdata/HG004.novaseq.wes_idt.100x.dedup.bam"
  CAPTURE_BED="${GCS_DATA_DIR}/exome-case-study-testdata/idt_capture_novogene.grch38.bed"
else
  if [[ -n "${MODEL_PRESET}" ]]; then
    echo "Error: --model_preset must be one of WGS, WGS_CHR20, WES, PACBIO, PACBIO_CHR20." >&2
    exit 1
  fi
fi

## Flag consistency sanity checks.

## The user should have not set --use_hp_information if it's not PACBIO.
if [[ "${USE_HP_INFORMATION}" != "unset" ]] && [[ "${MODEL_TYPE}" != "PACBIO" ]]; then
  echo "Error: Only set --use_hp_information for PACBIO." >&2
  exit 1
fi
if [[ "${USE_HP_INFORMATION}" == "unset" ]] && [[ "${MODEL_TYPE}" == "PACBIO" ]]; then
  # This is for backward-compatibility.
  echo "For PACBIO, set use_hp_information=true if it's not explicitly set."
  USE_HP_INFORMATION="true"
fi

N_SHARDS="64"
INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"
MODELS_DIR="${INPUT_DIR}/models"
CHILD_MODEL_DIR="${MODELS_DIR}/child"
PARENT_MODEL_DIR="${MODELS_DIR}/parent"
CHILD_MODEL="${CHILD_MODEL_DIR}/model.ckpt"
PARENT_MODEL="${PARENT_MODEL_DIR}/model.ckpt"
EXAMPLES_PREFIX="${OUTPUT_DIR}/make_examples.tfrecord@${N_SHARDS}.gz"
EXAMPLES_CHILD="${OUTPUT_DIR}/make_examples_child.tfrecord@${N_SHARDS}.gz"
EXAMPLES_PARENT1="${OUTPUT_DIR}/make_examples_parent1.tfrecord@${N_SHARDS}.gz"
EXAMPLES_PARENT2="${OUTPUT_DIR}/make_examples_parent2.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS_PREFIX="${OUTPUT_DIR}/gvcf.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS_CHILD="${OUTPUT_DIR}/gvcf_child.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS_PARENT1="${OUTPUT_DIR}/gvcf_parent1.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS_PARENT2="${OUTPUT_DIR}/gvcf_parent2.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT_CHILD="${OUTPUT_DIR}/call_variants_output_child.tfrecord.gz"
CALL_VARIANTS_OUTPUT_PARENT1="${OUTPUT_DIR}/call_variants_output_parent1.tfrecord.gz"
CALL_VARIANTS_OUTPUT_PARENT2="${OUTPUT_DIR}/call_variants_output_parent2.tfrecord.gz"
OUTPUT_VCF_CHILD="HG002.output.vcf.gz"
OUTPUT_VCF_PARENT1="HG003.output.vcf.gz"
OUTPUT_VCF_PARENT2="HG004.output.vcf.gz"
OUTPUT_VCF_MERGED="HG002-3-4.output.vcf.gz"
OUTPUT_GVCF_CHILD="HG002.output.g.vcf.gz"
OUTPUT_GVCF_PARENT1="HG003.output.g.vcf.gz"
OUTPUT_GVCF_PARENT2="HG004.output.g.vcf.gz"
OUTPUT_STATS_CHILD="HG002.output.vcf_stats"
OUTPUT_STATS_PARENT1="HG003.output.vcf_stats"
OUTPUT_STATS_PARENT2="HG004.output.vcf_stats"
REF="${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
SAMPLE_NAME_CHILD="HG002"
SAMPLE_NAME_PARENT1="HG003"
SAMPLE_NAME_PARENT2="HG004"
TRUTH_VCF_CHILD="${GCS_DATA_DIR}/case-study-testdata/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_CHILD="${GCS_DATA_DIR}/case-study-testdata/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF_PARENT1="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_PARENT1="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF_PARENT2="${GCS_DATA_DIR}/case-study-testdata/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_PARENT2="${GCS_DATA_DIR}/case-study-testdata/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
LOG_DIR="${OUTPUT_DIR}/logs"

if [[ "${MODEL_TYPE}" = "WES" ]]; then
  if [[ -n "${REGIONS}" ]]; then
    echo "Error: --regions is not used with model_type WES. Please use --capture_bed." >&2
    exit 1
  fi
  extra_args+=( --regions "/input/$(basename $CAPTURE_BED)")
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
echo "USE_HP_INFORMATION: ${USE_HP_INFORMATION}"
echo "SAVE_INTERMEDIATE_RESULTS: ${SAVE_INTERMEDIATE_RESULTS}"
echo "# Strings; sorted alphabetically."
echo "BAM_CHILD: ${BAM_CHILD}"
echo "BAM_PARENT1: ${BAM_PARENT1}"
echo "BAM_PARENT2: ${BAM_PARENT2}"
echo "BIN_VERSION: ${BIN_VERSION}"
echo "CALL_VARIANTS_ARGS: ${CALL_VARIANTS_ARGS}"
echo "CAPTURE_BED: ${CAPTURE_BED}"
echo "CUSTOMIZED_MODEL: ${CUSTOMIZED_MODEL}"
echo "MAKE_EXAMPLES_ARGS: ${MAKE_EXAMPLES_ARGS}"
echo "MODEL_PRESET: ${MODEL_PRESET}"
echo "MODEL_TYPE: ${MODEL_TYPE}"
echo "POSTPROCESS_VARIANTS_ARGS: ${POSTPROCESS_VARIANTS_ARGS}"
echo "PROPOSED_VARIANTS: ${PROPOSED_VARIANTS}"
echo "REF: ${REF}"
echo "REGIONS: ${REGIONS}"
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
      run gsutil -m cp "$1" "$2"
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
  copy_gs_or_http_file "${TRUTH_BED_CHILD}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF_CHILD}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF_CHILD}.tbi" "${INPUT_DIR}"
  # Copy HG003 truth
  copy_gs_or_http_file "${TRUTH_BED_PARENT1}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF_PARENT1}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF_PARENT1}.tbi" "${INPUT_DIR}"
  # Copy HG004 truth
  copy_gs_or_http_file "${TRUTH_BED_PARENT2}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF_PARENT2}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF_PARENT2}.tbi" "${INPUT_DIR}"

  copy_gs_or_http_file "${BAM_CHILD}" "${INPUT_DIR}"
  copy_correct_index_file "${BAM_CHILD}" "${INPUT_DIR}"
  copy_gs_or_http_file "${BAM_PARENT1}" "${INPUT_DIR}"
  copy_correct_index_file "${BAM_PARENT1}" "${INPUT_DIR}"
  copy_gs_or_http_file "${BAM_PARENT2}" "${INPUT_DIR}"
  copy_correct_index_file "${BAM_PARENT2}" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz.fai" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz.gzi" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gzi" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.fai" "${INPUT_DIR}"
  if [[ "${MODEL_TYPE}" = "WES" ]]; then
    copy_gs_or_http_file "${CAPTURE_BED}" "${INPUT_DIR}"
  fi
  if [[ -n "${PROPOSED_VARIANTS}" ]]; then
    copy_gs_or_http_file "${PROPOSED_VARIANTS}" "${INPUT_DIR}"
    copy_gs_or_http_file "${PROPOSED_VARIANTS}.tbi" "${INPUT_DIR}"
  fi
}

function setup_test() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    return
  fi

  ## Create local directory structure
  mkdir -p "${OUTPUT_DIR}"
  mkdir -p "${INPUT_DIR}"
  mkdir -p "${MODELS_DIR}"
  mkdir -p "${CHILD_MODEL_DIR}"
  mkdir -p "${PARENT_MODEL_DIR}"
  mkdir -p "${LOG_DIR}"

  ## Download extra packages
  # Install aria2 to download data files.
  sudo apt-get -qq -y update
  sudo apt-get -qq -y install aria2
  sudo apt-get -qq -y install bcftools
  sudo apt-get -qq -y install tabix

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
  if [[ "${BUILD_DOCKER}" = true ]]; then
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="deeptrio_gpu:latest"
      run "sudo docker build \
        -f Dockerfile.deeptrio \
        --build-arg=FROM_IMAGE=nvidia/cuda:11.3.0-cudnn8-devel-ubuntu20.04 \
        --build-arg=DV_GPU_BUILD=1 -t deeptrio_gpu ."
      run echo "Done building GPU Docker image ${IMAGE}."
      docker_args+=( --gpus 1 )
    else
      IMAGE="deeptrio:latest"
      # Building twice in case the first one times out.
      run "sudo docker build -f Dockerfile.deeptrio -t deeptrio . || \
        (sleep 5 ; sudo docker build -f Dockerfile.deeptrio -t deeptrio . )"
      run echo "Done building Docker image ${IMAGE}."
    fi
  else
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="google/deepvariant:deeptrio-${BIN_VERSION}-gpu"
      # shellcheck disable=SC2027
      # shellcheck disable=SC2086
      run "sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")"
      docker_args+=( --gpus 1 )
    else
      IMAGE="google/deepvariant:deeptrio-${BIN_VERSION}"
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
    # redacted
    run echo "Copy from gs:// path ${CUSTOMIZED_MODEL} to ${INPUT_DIR}/"
    run gsutil cp "${CUSTOMIZED_MODEL}".data-00000-of-00001 "${INPUT_DIR}/model.ckpt.data-00000-of-00001"
    run gsutil cp "${CUSTOMIZED_MODEL}".index "${INPUT_DIR}/model.ckpt.index"
    run gsutil cp "${CUSTOMIZED_MODEL}".meta "${INPUT_DIR}/model.ckpt.meta"
    run "gsutil cp ${CUSTOMIZED_MODEL}.input_shape ${INPUT_DIR}/model.ckpt.input_shape || echo 'skip input_shape'"
    extra_args+=( --customized_model "/input/model.ckpt")
    # redacted
    # echo "Copy from gs:// path $CUSTOMIZED_MODEL to ${INPUT_DIR}/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/child/model.ckpt.data-00000-of-00001 "${INPUT_DIR}/child/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/child/model.ckpt.index "${INPUT_DIR}child/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/child/model.ckpt.meta "${INPUT_DIR}child/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/parent/model.ckpt.data-00000-of-00001 "${INPUT_DIR}/parent/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/parent/model.ckpt.index "${INPUT_DIR}parent/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/parent/model.ckpt.meta "${INPUT_DIR}parent/"
    # extra_args+=( --customized_model "/input/model.ckpt")
  else
    run echo "No custom model specified."
  fi
  if [[ -n "${MAKE_EXAMPLES_ARGS}" ]]; then
    # In order to use proposed variants, we have to pass vcf_candidate_importer
    # to make_examples_extra_args, so we know that we will enter this if
    # statement.
    if [[ -n "${PROPOSED_VARIANTS}" ]]; then
      MAKE_EXAMPLES_ARGS="${MAKE_EXAMPLES_ARGS},proposed_variants=/input/$(basename "$PROPOSED_VARIANTS")"
    fi
    extra_args+=( --make_examples_extra_args "${MAKE_EXAMPLES_ARGS}")
  fi
  if [[ -n "${CALL_VARIANTS_ARGS}" ]]; then
    extra_args+=( --call_variants_extra_args "${CALL_VARIANTS_ARGS}")
  fi
  if [[ "${USE_HP_INFORMATION}" == "true" ]] || [[ "${USE_HP_INFORMATION}" == "false" ]]; then
    # Note that because --use_hp_information is a binary flag, I need to use
    # the --flag=(true|false) format, or use --[no]flag format.
    extra_args+=( "--use_hp_information=${USE_HP_INFORMATION}" )
  fi
  if [[ -n "${POSTPROCESS_VARIANTS_ARGS}" ]]; then
    extra_args+=( --postprocess_variants_extra_args "${POSTPROCESS_VARIANTS_ARGS}")
  fi
  if [[ -n "${REGIONS}" ]]; then
    extra_args+=( --regions "${REGIONS}")
    happy_args+=( -l "${REGIONS}")
  fi
  if [[ "${BUILD_DOCKER}" = true ]] || [[ "${BIN_VERSION}" =~ ^1\.[2-9]\.0$ ]]; then
    extra_args+=( --runtime_report )
  fi
}

function run_deeptrio() {
  run echo "Run DeepTrio..."
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
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
    --model_type "${MODEL_TYPE}" \
    --ref="/input/$(basename $REF).gz" \
    --reads_child "/input/$(basename $BAM_CHILD)" \
    --reads_parent1 "/input/$(basename $BAM_PARENT1)" \
    --reads_parent2 "/input/$(basename $BAM_PARENT2)" \
    --output_vcf_child "/output/${OUTPUT_VCF_CHILD}" \
    --output_vcf_parent1 "/output/${OUTPUT_VCF_PARENT1}" \
    --output_vcf_parent2 "/output/${OUTPUT_VCF_PARENT2}" \
    --sample_name_child "${SAMPLE_NAME_CHILD}" \
    --sample_name_parent1 "${SAMPLE_NAME_PARENT1}" \
    --sample_name_parent2 "${SAMPLE_NAME_PARENT2}" \
    --num_shards "$(nproc)" \
    --intermediate_results_dir /output/intermediate_results_dir \
    --output_gvcf_child "/output/${OUTPUT_GVCF_CHILD}" \
    --output_gvcf_parent1 "/output/${OUTPUT_GVCF_PARENT1}" \
    --output_gvcf_parent2 "/output/${OUTPUT_GVCF_PARENT2}" \
    --logging_dir="/output/logs" \
    "${extra_args[@]-}" && \
  echo "Done.")) 2>&1 | tee "${LOG_DIR}/deeptrio_runtime.log""
  echo
}


function run_happy() {
  ## Evaluation: run hap.py
  local -r truth_vcf=$(basename "${1}")
  local -r truth_bed=$(basename "${2}")
  local -r vcf_output=$(basename "${3}")
  local -r optional_chr="${4-}"

  run echo "Start evaluation with hap.py..."
  declare -a optional_chr_args
  if [[ -n ${optional_chr} ]]; then
    optional_chr_args+=(-l "${optional_chr}")
    run echo "... evaluating on ${optional_chr}"
  fi
  UNCOMPRESSED_REF=${INPUT_DIR}/$(basename $REF)
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
  # shellcheck disable=SC2086
  # shellcheck disable=SC2145
  run "( sudo docker run -i \
  -v "${INPUT_DIR}:${INPUT_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  jmcdani20/hap.py:${HAPPY_VERSION} /opt/hap.py/bin/hap.py \
    "${INPUT_DIR}/${truth_vcf}" \
    "${OUTPUT_DIR}/${vcf_output}" \
    -f "${INPUT_DIR}/${truth_bed}" \
    -r ${UNCOMPRESSED_REF} \
    -o "${OUTPUT_DIR}/happy.output${optional_chr}-${vcf_output}" \
    --engine=vcfeval \
    --pass-only \
    ${happy_args[@]-} \
    ${optional_chr_args[@]-} \
  ) 2>&1 | tee "${LOG_DIR}/happy.output${optional_chr}-${vcf_output}.log""
  if [[ "${DRY_RUN}" != "true" ]]; then
    echo "${HAPPY_VERSION}" > "${LOG_DIR}/happy_version.log"
  fi
  run echo "Done."
}

function run_happy_reports() {
  local -r optional_chr="${1-}"

  if [[ ${DRY_RUN} == "true" ]]; then
    run_happy "${TRUTH_VCF_CHILD}" "${TRUTH_BED_CHILD}" "${OUTPUT_VCF_CHILD}" "${optional_chr}"
    run_happy "${TRUTH_VCF_PARENT1}" "${TRUTH_BED_PARENT1}" "${OUTPUT_VCF_PARENT1}" "${optional_chr}"
    run_happy "${TRUTH_VCF_PARENT2}" "${TRUTH_BED_PARENT2}" "${OUTPUT_VCF_PARENT2}" "${optional_chr}"
    run_happy "${TRUTH_VCF_CHILD}" "${TRUTH_BED_CHILD}" "HG002_split.vcf" "${optional_chr}"
    run_happy "${TRUTH_VCF_PARENT1}" "${TRUTH_BED_PARENT1}" "HG003_split.vcf" "${optional_chr}"
    run_happy "${TRUTH_VCF_PARENT2}" "${TRUTH_BED_PARENT2}" "HG004_split.vcf" "${optional_chr}"
  else
    run_happy "${TRUTH_VCF_CHILD}" "${TRUTH_BED_CHILD}" "${OUTPUT_VCF_CHILD}" "${optional_chr}" 2>&1 | tee "${LOG_DIR}/happy_child.log"
    run_happy "${TRUTH_VCF_PARENT1}" "${TRUTH_BED_PARENT1}" "${OUTPUT_VCF_PARENT1}" "${optional_chr}" 2>&1 | tee "${LOG_DIR}/happy_parent1.log"
    run_happy "${TRUTH_VCF_PARENT2}" "${TRUTH_BED_PARENT2}" "${OUTPUT_VCF_PARENT2}" "${optional_chr}" 2>&1 | tee "${LOG_DIR}/happy_parent2.log"
    run_happy "${TRUTH_VCF_CHILD}" "${TRUTH_BED_CHILD}" "HG002_split.vcf" "${optional_chr}" 2>&1 | tee "${LOG_DIR}/happy_child_split.log"
    run_happy "${TRUTH_VCF_PARENT1}" "${TRUTH_BED_PARENT1}" "HG003_split.vcf" "${optional_chr}" 2>&1 | tee "${LOG_DIR}/happy_parent1_split.log"
    run_happy "${TRUTH_VCF_PARENT2}" "${TRUTH_BED_PARENT2}" "HG004_split.vcf" "${optional_chr}" 2>&1 | tee "${LOG_DIR}/happy_parent2_split.log"
  fi
}

function run_glnexus() {
  if [[ ${DRY_RUN} == "true" ]]; then
    return
  fi

  sudo docker pull quay.io/mlin/glnexus:v1.2.7 || \
    (sleep 5 ; sudo docker pull quay.io/mlin/glnexus:v1.2.7)

  time sudo docker run \
    -v "${OUTPUT_DIR}":"/output" \
    quay.io/mlin/glnexus:v1.2.7 \
    /usr/local/bin/glnexus_cli \
    --config DeepVariant_unfiltered \
    "/output/${OUTPUT_GVCF_PARENT2}" "/output/${OUTPUT_GVCF_PARENT1}" "/output/${OUTPUT_GVCF_CHILD}" \
    | bcftools view - | bgzip -c > "${OUTPUT_DIR}/${OUTPUT_VCF_MERGED}"
}

function extract_samples() {
  if [[ ${DRY_RUN} == "true" ]]; then
    return
  fi

  bcftools view -s HG002 "${OUTPUT_DIR}/${OUTPUT_VCF_MERGED}" > "${OUTPUT_DIR}/HG002_split.vcf"
  bcftools view -s HG003 "${OUTPUT_DIR}/${OUTPUT_VCF_MERGED}" > "${OUTPUT_DIR}/HG003_split.vcf"
  bcftools view -s HG004 "${OUTPUT_DIR}/${OUTPUT_VCF_MERGED}" > "${OUTPUT_DIR}/HG004_split.vcf"
}

function main() {
  run echo 'Starting the test...'

  setup_test
  copy_data
  get_docker_image
  setup_args
  run_deeptrio
  run_glnexus
  extract_samples
  run_happy_reports
  if [[ -z "${REGIONS}" ]]; then
    run_happy_reports "chr20"
  fi
}

main "$@"
