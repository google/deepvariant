#!/bin/bash
# Copyright 2020 Google LLC.

set -euo pipefail


USAGE=$'
Example usage:
inference_deepvariant.sh --model_preset "WGS" --docker_build true --use_gpu true

Flags:
--docker_build (true|false)  Whether to build docker image. (default: false)
--use_gpu (true|false)   Whether to use GPU when running case study. Make sure to specify vm_zone that is equipped with GPUs. (default: false)
--customized_model Path to checkpoint directory containing model checkpoint.
--regions Regions passed into both variant calling and hap.py.
--make_examples_extra_args Flags for make_examples, specified as "flag1=param1,flag2=param2".
--call_variants_extra_args Flags for call_variants, specified as "flag1=param1,flag2=param2".
--postprocess_variants_extra_args Flags for postprocess_variants, specified as "flag1=param1,flag2=param2".
--model_preset Preset case study to run: WGS, WES, PACBIO, or HYBRID_PACBIO_ILLUMINA.

If model_preset is not specified, the below flags are required:
--model_type Type of DeepVariant model to run (WGS, WES, PACBIO, HYBRID_PACBIO_ILLUMINA)
--bin_version Version of DeepVariant model to use
--ref Path to GCP bucket containing ref file (.fa)
--bam Path to GCP bucket containing BAM
--truth_vcf Path to GCP bucket containing truth VCF
--truth_bed Path to GCP bucket containing truth BED
--capture_bed Path to GCP bucket containing captured file (only needed for WES model_type)

Note: All paths to dataset must be of the form "gs://..."
'

# Specify default values.
BUILD_DOCKER=false
USE_GPU=false
REGIONS=""
CUSTOMIZED_MODEL=""
MAKE_EXAMPLES_ARGS=""
CALL_VARIANTS_ARGS=""
POSTPROCESS_VARIANTS_ARGS=""
MODEL_PRESET=""
MODEL_TYPE=""
BIN_VERSION=""
REF=""
BAM=""
TRUTH_VCF=""
TRUTH_BED=""
CAPTURE_BED=""

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
    -*|--*=) # other flags not supported
      echo "Error: unrecognized flag $1" >&2
      echo "$USAGE" >&2
      exit 1
      ;;
    *)
      echo "Error: unrecognized extra args $1" >&2
      echo "$USAGE" >&2
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

if [[ -n $MODEL_PRESET ]]; then
  BIN_VERSION="1.1.0"
fi

if [[ "${MODEL_PRESET}" = "PACBIO" ]]; then
  MODEL_TYPE="PACBIO"
  BASE="${HOME}/pacbio-case-study"

  REF="${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  BAM="${GCS_DATA_DIR}/pacbio-case-study-testdata/HG003.pfda_challenge.grch38.phased.bam"
  TRUTH_VCF="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
  TRUTH_BED="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.bed"
elif [[ "${MODEL_PRESET}" = "WGS" ]]; then
  MODEL_TYPE="WGS"
  BASE="${HOME}/wgs-case-study"

  REF="${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  BAM="${GCS_DATA_DIR}/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
  TRUTH_VCF="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
  TRUTH_BED="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.bed"
elif [[ "${MODEL_PRESET}" = "WES" ]]; then
  MODEL_TYPE="WES"
  BASE="${HOME}/exome-case-study"

  REF="${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  BAM="${GCS_DATA_DIR}/exome-case-study-testdata/HG003.novaseq.wes_idt.100x.dedup.bam"
  TRUTH_VCF="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
  TRUTH_BED="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.bed"
  CAPTURE_BED="${GCS_DATA_DIR}/exome-case-study-testdata/idt_capture_novogene.grch38.bed"
elif [[ "${MODEL_PRESET}" = "HYBRID_PACBIO_ILLUMINA" ]]; then
  MODEL_TYPE="HYBRID_PACBIO_ILLUMINA"
  BASE="${HOME}/hybrid-pacbio-illumina-case-study"

  REF="${GCS_DATA_DIR}/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  BAM="${GCS_DATA_DIR}/hybrid-case-study-testdata/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.bam"
  TRUTH_VCF="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
  TRUTH_BED="${GCS_DATA_DIR}/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.bed"
else
  if [[ -n $MODEL_PRESET ]]; then
    echo "Error: --model_preset must be one of WGS, WES, PACBIO, HYBRID_PACBIO_ILLUMINA." >&2
    exit 1
  fi

fi

N_SHARDS="64"
INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="deepvariant.output.vcf.gz"
OUTPUT_GVCF="deepvariant.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

if [[ "${MODEL_TYPE}" = "PACBIO" ]]; then
  extra_args+=( --use_hp_information )
fi
if [[ "${MODEL_TYPE}" = "WES" ]]; then
  if [[ -n $REGIONS ]]; then
    echo "Error: --regions is not used with model_type WES. Please use --capture_bed." >&2
    exit 1
  fi
  extra_args+=( --regions "/input/$(basename $CAPTURE_BED)")
  happy_args+=( -T "${INPUT_DIR}/$(basename $CAPTURE_BED)")
fi

echo "========================="
echo "BUILD_DOCKER: $BUILD_DOCKER"
echo "CUSTOMIZED_MODEL: $CUSTOMIZED_MODEL"
echo "REGIONS: $REGIONS"
echo "MAKE_EXAMPLES_ARGS: $MAKE_EXAMPLES_ARGS"
echo "CALL_VARIANTS_ARGS: $CALL_VARIANTS_ARGS"
echo "POSTPROCESS_VARIANTS_ARGS: $POSTPROCESS_VARIANTS_ARGS"
echo "USE_GPU: $USE_GPU"
echo "MODEL_PRESET: $MODEL_PRESET"
echo "MODEL_TYPE: $MODEL_TYPE"
echo "BIN_VERSION: $BIN_VERSION"
echo "REF: $REF"
echo "BAM: $BAM"
echo "TRUTH_VCF: $TRUTH_VCF"
echo "TRUTH_BED: $TRUTH_BED"
echo "CAPTURE_BED: $CAPTURE_BED"
echo "========================="

function copy_data() {
  # For the presets, we use `aria2c https://storage.googleapis.com/...` since
  # some users have had difficulty installing gsutil in the past.
  # However, in the general use case, we prefer to use `gsutil cp`, so to use
  # custom data with this script gsutil is required.
  if [[ -n "${MODEL_PRESET}" ]]; then
    aria2c -c -x10 -s10 "${TRUTH_BED}" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${TRUTH_VCF}" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${TRUTH_VCF}.tbi" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${BAM}" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${BAM}.bai" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${REF}.gz" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${REF}.gz.fai" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${REF}.gz.gzi" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${REF}.gzi" -d "${INPUT_DIR}"
    aria2c -c -x10 -s10 "${REF}.fai" -d "${INPUT_DIR}"
    if [[ "${MODEL_PRESET}" = "WES" ]]; then
      aria2c -c -x10 -s10 "${CAPTURE_BED}" -d "${INPUT_DIR}"
    fi
  else
    gsutil -m cp "${TRUTH_BED}" "${INPUT_DIR}"
    gsutil -m cp "${TRUTH_VCF}" "${INPUT_DIR}"
    gsutil -m cp "${TRUTH_VCF}.tbi" "${INPUT_DIR}"
    gsutil -m cp "${BAM}" "${INPUT_DIR}"
    gsutil -m cp "${BAM}.bai" "${INPUT_DIR}"
    gsutil -m cp "${REF}.gz" "${INPUT_DIR}"
    gsutil -m cp "${REF}.gz.fai" "${INPUT_DIR}"
    gsutil -m cp "${REF}.gz.gzi" "${INPUT_DIR}"
    gsutil -m cp "${REF}.gzi" "${INPUT_DIR}"
    gsutil -m cp "${REF}.fai" "${INPUT_DIR}"
    if [[ "${MODEL_TYPE}" = "WES" ]]; then
      gsutil -m cp "${CAPTURE_BED}" "${INPUT_DIR}"
    fi
  fi
}

function setup_test() {

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

  copy_data

  if [[ "${BUILD_DOCKER}" = true ]]; then
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="deepvariant_gpu:latest"
      sudo docker build \
        --build-arg=FROM_IMAGE=nvidia/cuda:10.0-cudnn7-devel-ubuntu18.04 \
        --build-arg=DV_GPU_BUILD=1 -t deepvariant_gpu .
      echo "Done building GPU Docker image ${IMAGE}."
      docker_args+=( --gpus 1 )
    else
      IMAGE="deepvariant:latest"
      # Building twice in case the first one times out.
      sudo docker build -t deepvariant . --build-arg DV_OPENVINO_BUILD=1 || \
        (sleep 5 ; sudo docker build -t deepvariant . --build-arg DV_OPENVINO_BUILD=1)
      echo "Done building Docker image ${IMAGE}."
    fi
  else
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="google/deepvariant:${BIN_VERSION}-gpu"
      sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")
      docker_args+=( --gpus 1 )
    else
      IMAGE="google/deepvariant:${BIN_VERSION}"
      sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")
    fi
  fi
}

function run_deepvariant_with_docker() {
  echo "Run DeepVariant..."
  echo "using IMAGE=$IMAGE"

  if [[ -n $CUSTOMIZED_MODEL ]]; then
    echo "Copy from gs:// path $CUSTOMIZED_MODEL to ${INPUT_DIR}/"
    gsutil cp "${CUSTOMIZED_MODEL}"/model.ckpt.data-00000-of-00001 "${INPUT_DIR}"
    gsutil cp "${CUSTOMIZED_MODEL}"/model.ckpt.index "${INPUT_DIR}"
    gsutil cp "${CUSTOMIZED_MODEL}"/model.ckpt.meta "${INPUT_DIR}"
    extra_args+=( --customized_model "/input/model.ckpt")
  else
      echo "No custom model specified."
  fi
  if [[ -n $MAKE_EXAMPLES_ARGS ]]; then
    extra_args+=( --make_examples_extra_args "${MAKE_EXAMPLES_ARGS}")
  fi
  if [[ -n $CALL_VARIANTS_ARGS ]]; then
    extra_args+=( --call_variants_extra_args "${CALL_VARIANTS_ARGS}")
  fi
  if [[ -n $POSTPROCESS_VARIANTS_ARGS ]]; then
    extra_args+=( --postprocess_variants_extra_args "${POSTPROCESS_VARIANTS_ARGS}")
  fi
  if [[ -n $REGIONS ]]; then
    extra_args+=( --regions "${REGIONS}")
    happy_args+=( -l "${REGIONS}")
  fi

  # shellcheck disable=SC2068
  (time (sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    ${docker_args[@]-} \
    "${IMAGE}" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type="${MODEL_TYPE}" \
    --ref="/input/$(basename $REF).gz" \
    --reads="/input/$(basename $BAM)" \
    --output_vcf="/output/${OUTPUT_VCF}" \
    --output_gvcf="/output/${OUTPUT_GVCF}" \
    --num_shards=${N_SHARDS} \
    --logging_dir="/output/logs" \
    "${extra_args[@]-}"
  echo "Done.")) 2>&1 | tee "${LOG_DIR}/deepvariant_runtime.log"
  echo
}

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${INPUT_DIR}/$(basename $REF)"

function run_happy() {
  # hap.py cannot read the compressed fa, so uncompress
  # into a writable directory. Index file was downloaded earlier.
  zcat <"${INPUT_DIR}/$(basename $REF).gz" >"${UNCOMPRESSED_REF}"

  sudo docker pull pkrusche/hap.py
  # shellcheck disable=SC2068
  ( sudo docker run -i \
  -v "${INPUT_DIR}:${INPUT_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
    "${INPUT_DIR}/$(basename $TRUTH_VCF)" \
    "${OUTPUT_DIR}/${OUTPUT_VCF}" \
    -f "${INPUT_DIR}/$(basename $TRUTH_BED)" \
    -r "${UNCOMPRESSED_REF}" \
    -o "${OUTPUT_DIR}/happy.output" \
    --engine=vcfeval \
    ${happy_args[@]-}
  ) 2>&1 | tee "${LOG_DIR}/happy.log"
  echo "Done."
}

function main() {
  echo 'Starting the test...'

  setup_test
  run_deepvariant_with_docker
  run_happy 2>&1 | tee "${LOG_DIR}/happy.log"
}

main "$@"
