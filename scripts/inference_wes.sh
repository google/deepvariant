#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/exome-case-study"
BIN_VERSION="1.0.0"

INPUT_DIR="${BASE}/input/data"
REF="hs37d5.fa"
BAM="151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
TRUTH_VCF="HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz"
TRUTH_BED="HG002_GRCh37_1_22_v4.1_draft_benchmark.bed"

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/HG002.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="HG002.output.vcf.gz"
OUTPUT_GVCF="HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

CAPTURE_BED="agilent_sureselect_human_all_exon_v5_b37_targets.bed"

# Whether to build docker image.
BUILD_DOCKER="${1:-false}"

# Optional extra flags for DeepVariant.
CUSTOMIZED_MODEL="${2:-}"
MAKE_EXAMPLES_ARGS="${3:-}"
CALL_VARIANTS_ARGS="${4:-}"
POSTPROCESS_VARIANTS_ARGS="${5:-}"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"
mkdir -p "${LOG_DIR}"

declare -a extra_args

if [[ -n $CUSTOMIZED_MODEL ]]
then
  echo "Copy from gs:// path $CUSTOMIZED_MODEL to ${INPUT_DIR}/"
  gsutil cp "${CUSTOMIZED_MODEL}"/model.ckpt.data-00000-of-00001 "${INPUT_DIR}"
  gsutil cp "${CUSTOMIZED_MODEL}"/model.ckpt.index "${INPUT_DIR}"
  gsutil cp "${CUSTOMIZED_MODEL}"/model.ckpt.meta "${INPUT_DIR}"
  extra_args+=( --customized_model "/input/model.ckpt")
else
    echo "No custom model specified."
fi
if [[ -n $MAKE_EXAMPLES_ARGS ]]
then
  extra_args+=( --make_examples_extra_args "${MAKE_EXAMPLES_ARGS}")
fi
if [[ -n $CALL_VARIANTS_ARGS ]]
then
  extra_args+=( --call_variants_extra_args "${CALL_VARIANTS_ARGS}")
fi
if [[ -n $POSTPROCESS_VARIANTS_ARGS ]]
then
  extra_args+=( --postprocess_variants_extra_args "${POSTPROCESS_VARIANTS_ARGS}")
fi

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

# Copy the data
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${BAM}"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${TRUTH_BED}"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${TRUTH_VCF}"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${TRUTH_VCF}.tbi"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${CAPTURE_BED}"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${REF}.gz"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${REF}.gz.fai"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${REF}.gz.gzi"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${REF}.gzi"
aria2c -c -x10 -s10 -d "${INPUT_DIR}" "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/${REF}.fai"

if [[ "${BUILD_DOCKER}" = true ]]
then
  IMAGE="deepvariant:latest"
  # Pulling twice in case the first one times out.
  sudo docker build -t deepvariant . || \
    (sleep 5 ; sudo docker build -t deepvariant .)
  echo "Done building Docker image ${IMAGE}."
else
  IMAGE="google/deepvariant:${BIN_VERSION}"
  sudo docker pull "${IMAGE}"
fi

echo "Run DeepVariant..."
echo "using IMAGE=$IMAGE"
(time ( sudo docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  "${IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="/input/${REF}.gz" \
  --reads="/input/${BAM}" \
  --regions="/input/${CAPTURE_BED}" \
  --output_vcf="/output/${OUTPUT_VCF}" \
  --output_gvcf="/output/${OUTPUT_GVCF}" \
  --num_shards=${N_SHARDS} \
  "${extra_args[@]-}"
echo "Done.")) 2>&1 | tee "${LOG_DIR}/deepvariant_runtime.log"
echo

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${INPUT_DIR}/${REF}"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory. Index file was downloaded earlier.
zcat <"${INPUT_DIR}/${REF}.gz" >"${UNCOMPRESSED_REF}"

sudo docker pull pkrusche/hap.py
( sudo docker run -i \
-v "${INPUT_DIR}:${INPUT_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${INPUT_DIR}/${TRUTH_VCF}" \
  "${OUTPUT_DIR}/${OUTPUT_VCF}" \
  -f "${INPUT_DIR}/${TRUTH_BED}" \
  -T "${INPUT_DIR}/${CAPTURE_BED}" \
  -r "${UNCOMPRESSED_REF}" \
  -o "${OUTPUT_DIR}/happy.output" \
  --engine=vcfeval
) 2>&1 | tee "${LOG_DIR}/happy.log"
echo "Done."
