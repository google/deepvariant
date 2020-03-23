#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/pacbio-case-study"
BIN_VERSION="0.10.0"

INPUT_DIR="${BASE}/input/data"
REF="hs37d5.fa.gz"
BAM="pacbio.8M.30x.bam"
TRUTH_VCF="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="HG002.output.vcf.gz"
OUTPUT_GVCF="HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

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

  # Copy the data
  aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${TRUTH_BED}" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${TRUTH_VCF}" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${TRUTH_VCF}.tbi" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/pacbio-case-study-testdata/${BAM}" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/pacbio-case-study-testdata/${BAM}.bai" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.fai -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.gzi -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gzi -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.fai -d "${INPUT_DIR}"

  ## Pull the docker image.
  sudo docker pull google/deepvariant:"${BIN_VERSION}"
}

function run_deepvariant_with_docker() {
  echo "Run DeepVariant..."

  # If a customized model is specified, copy the model files to your local disk.
  GS_ADDRESS="^gs:\/\/.*"
  declare -a extra_args
  extra_args=( )

  if [[ -z $model_http_dir ]];
  then
    echo 'Use default model path.'
  elif [[ $model_http_dir =~ $GS_ADDRESS ]];
  then
    echo "Copy from gs:// path $model_http_dir to ${INPUT_DIR}/"
    gsutil cp "${model_http_dir}"/model.ckpt.data-00000-of-00001 "${INPUT_DIR}"
    gsutil cp "${model_http_dir}"/model.ckpt.index "${INPUT_DIR}"
    gsutil cp "${model_http_dir}"/model.ckpt.meta "${INPUT_DIR}"
    extra_args+=( --customized_model "/input/model.ckpt")
  else
    echo "Copy from local path $model_http_dir to ${INPUT_DIR}/"
    cp -f "${model_http_dir}"/model.ckpt.data-00000-of-00001 "${INPUT_DIR}"
    cp -f "${model_http_dir}"/model.ckpt.index "${INPUT_DIR}"
    cp -f "${model_http_dir}"/model.ckpt.meta "${INPUT_DIR}"
    extra_args+=( --customized_model "/input/model.ckpt")
  fi

  sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    google/deepvariant:"${BIN_VERSION}" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref="/input/${REF}" \
    --reads="/input/${BAM}" \
    --output_vcf=/output/${OUTPUT_VCF} \
    --output_gvcf=/output/${OUTPUT_GVCF} \
    --num_shards=${N_SHARDS} \
    "${extra_args[@]-}"
  echo "Done."
  echo
}

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${INPUT_DIR}/hs37d5.fa"

function run_happy() {
  # hap.py cannot read the compressed fa, so uncompress
  # into a writable directory. Index file was downloaded earlier.
  zcat <"${INPUT_DIR}/${REF}" >"${UNCOMPRESSED_REF}"

  sudo docker pull pkrusche/hap.py
  ( sudo docker run -i \
  -v "${INPUT_DIR}:${INPUT_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
    "${INPUT_DIR}/${TRUTH_VCF}" \
    "${OUTPUT_DIR}/${OUTPUT_VCF}" \
    -f "${INPUT_DIR}/${TRUTH_BED}" \
    -r "${UNCOMPRESSED_REF}" \
    -l 20 \
    -o "${OUTPUT_DIR}/happy.output" \
    --engine=vcfeval
  ) 2>&1 | tee "${LOG_DIR}/happy.log"
  echo "Done."
}

function main() {
  echo 'Starting the test...'

  local -r model_http_dir="${1:-}"
  if [[ ! -z $model_http_dir ]]; then
    echo "Using model from: ${model_http_dir}"
  fi

  setup_test
  run_deepvariant_with_docker
  run_happy 2>&1 | tee "${LOG_DIR}/happy.log"
}

main "$@"
