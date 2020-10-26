#!/bin/bash
# Copyright 2020 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/case-study"

INPUT_DIR="${BASE}/input"
MODELS_DIR="${INPUT_DIR}/models"
CHILD_MODEL_DIR="${MODELS_DIR}/child"
PARENT_MODEL_DIR="${MODELS_DIR}/parent"
CHILD_MODEL="${CHILD_MODEL_DIR}/model.ckpt"
PARENT_MODEL="${PARENT_MODEL_DIR}/model.ckpt"
DATA_DIR="${INPUT_DIR}/data"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
BAM_CHILD="HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
BAM_PARENT1="HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
BAM_PARENT2="HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
SAMPLE_NAME_CHILD="HG002"
SAMPLE_NAME_PARENT1="HG003"
SAMPLE_NAME_PARENT2="HG004"
TRUTH_VCF_CHILD="HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz"
TRUTH_BED_CHILD="HG002_GRCh38_1_22_v4.2_benchmark.bed"
TRUTH_VCF_PARENT1="HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
TRUTH_BED_PARENT1="HG003_GRCh38_1_22_v4.2_benchmark.bed"
TRUTH_VCF_PARENT2="HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz"
TRUTH_BED_PARENT2="HG004_GRCh38_1_22_v4.2_benchmark.bed"

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"
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
LOG_DIR="${OUTPUT_DIR}/logs"

# Whether to build docker image.
BUILD_DOCKER="${1:-false}"

# Optional extra flags for DeepVariant.
CUSTOMIZED_MODEL="${2:-}"
MAKE_EXAMPLES_ARGS="${3:-}"
CALL_VARIANTS_ARGS="${4:-}"
POSTPROCESS_VARIANTS_ARGS="${5:-}"
REGIONS="${6:-}"

declare -a extra_args
declare -a happy_args

function setup_test() {
  ## Create local directory structure
  mkdir -p "${OUTPUT_DIR}"
  mkdir -p "${DATA_DIR}"
  mkdir -p "${MODELS_DIR}"
  mkdir -p "${CHILD_MODEL_DIR}"
  mkdir -p "${PARENT_MODEL_DIR}"
  mkdir -p "${LOG_DIR}"

  ## Download extra packages
  # There are some extra programs we will need.
  # We are going to use [GNU Parallel](https://www.gnu.org/software/parallel/) to
  # run `make_examples`.
  sudo apt-get -y update
  sudo apt-get -y install aria2
  sudo apt-get -y install bcftools
  sudo apt-get -y install tabix

  # Install docker if needed
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

  # Copy the data, using http:// because of https://github.com/aria2/aria2/issues/1012.
  # Copy HG002 truth
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh38_1_22_v4.2_benchmark.bed -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi -d "${INPUT_DIR}"
  # Copy HG003 truth
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.bed -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi -d "${INPUT_DIR}"
  # Copy HG004 truth
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004_GRCh38_1_22_v4.2_benchmark.bed -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi -d "${INPUT_DIR}"

  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "http://storage.googleapis.com/deepvariant/case-study-testdata/${REF}.gz" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "http://storage.googleapis.com/deepvariant/case-study-testdata/${REF}.gz.fai" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "http://storage.googleapis.com/deepvariant/case-study-testdata/${REF}.gz.gzi" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "http://storage.googleapis.com/deepvariant/case-study-testdata/${REF}.gzi" -d "${INPUT_DIR}"
  aria2c -c -x10 -s10 "http://storage.googleapis.com/deepvariant/case-study-testdata/${REF}.fai" -d "${INPUT_DIR}"
}

function setup_args() {
  if [[ -n $CUSTOMIZED_MODEL ]]
  then
    # redacted
    echo "No custom model specified."
    # echo "Copy from gs:// path $CUSTOMIZED_MODEL to ${INPUT_DIR}/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/child/model.ckpt.data-00000-of-00001 "${INPUT_DIR}/child/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/child/model.ckpt.index "${INPUT_DIR}child/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/child/model.ckpt.meta "${INPUT_DIR}child/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/parent/model.ckpt.data-00000-of-00001 "${INPUT_DIR}/parent/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/parent/model.ckpt.index "${INPUT_DIR}parent/"
    # gsutil cp "${CUSTOMIZED_MODEL}"/parent/model.ckpt.meta "${INPUT_DIR}parent/"
    # extra_args+=( --customized_model "/input/model.ckpt")
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

  if [[ -n $REGIONS ]]
  then
    extra_args+=( --regions "${REGIONS}")
    happy_args+=( -l "${REGIONS}")
  fi
}

function get_docker_image() {
  if [[ "${BUILD_DOCKER}" = true ]]
  then
    IMAGE="deeptrio:latest"
    # Pulling twice in case the first one times out.
    sudo docker build -f Dockerfile.deeptrio -t deeptrio . || \
      (sleep 5 ; sudo docker -f Dockerfile.deeptrio build -t deeptrio .)
    echo "Done building Docker image ${IMAGE}."
  else
    IMAGE="google/deeptrio:${BIN_VERSION}"
    sudo docker pull "${IMAGE}"
  fi
}

function run_deeptrio() {
  echo "Run DeepTrio..."
  echo "using IMAGE=$IMAGE"
  (time ( sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    "${IMAGE}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
    --model_type WGS \
    --ref="/input/${REF}.gz" \
    --reads_child "/input/${BAM_CHILD}" \
    --reads_parent1 "/input/${BAM_PARENT1}" \
    --reads_parent2 "/input/${BAM_PARENT2}" \
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
    --regions chr20 \
    --logging_dir="/output/logs" \
    "${extra_args[@]-}"
  echo "Done.")) 2>&1 | tee "${LOG_DIR}/deeptrio_runtime.log"
  echo
}

function run_glnexus() {
  sudo docker pull quay.io/mlin/glnexus:v1.2.7 || \
    (sleep 5 ; sudo docker pull quay.io/mlin/glnexus:v1.2.7)

  time sudo docker run \
    -v "${OUTPUT_DIR}":"/output" \
    quay.io/mlin/glnexus:v1.2.7 \
    /usr/local/bin/glnexus_cli \
    --config DeepVariantWGS \
    "/output/${OUTPUT_GVCF_PARENT2}" "/output/${OUTPUT_GVCF_PARENT1}" "/output/${OUTPUT_GVCF_CHILD}" \
    | bcftools view - | bgzip -c > "${OUTPUT_DIR}/${OUTPUT_VCF_MERGED}"
}

## Evaluation: run hap.py
function run_happy() {
  local -r truth_vcf="${1}"
  local -r truth_bed="${2}"
  local -r vcf_output="${3}"

  echo "Start evaluation with hap.py..."

  # hap.py cannot read the compressed fa, so uncompress
  # into a writable directory. Index file was downloaded earlier.
  zcat <"${INPUT_DIR}/${REF}.gz" >"${INPUT_DIR}/${REF}"

  # Pulling twice in case the first one times out.
  sudo docker pull pkrusche/hap.py || \
    (sleep 5 ; sudo docker pull pkrusche/hap.py)
  sudo docker run -i \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
    "/input/${truth_vcf}" \
    "/output/${vcf_output}" \
    -f "/input/${truth_bed}" \
    -r "/input/${REF}" \
    -o "/output/happy.output" \
    --engine=vcfeval \
    -l chr20
  echo "Done."
}

function run_all_happy_reports() {
  run_happy "${TRUTH_VCF_CHILD}" "${TRUTH_BED_CHILD}" "${OUTPUT_VCF_CHILD}" 2>&1 | tee "${LOG_DIR}/happy_child.log"
  run_happy "${TRUTH_VCF_PARENT1}" "${TRUTH_BED_PARENT1}" "${OUTPUT_VCF_PARENT1}" 2>&1 | tee "${LOG_DIR}/happy_parent1.log"
  run_happy "${TRUTH_VCF_PARENT2}" "${TRUTH_BED_PARENT2}" "${OUTPUT_VCF_PARENT2}" 2>&1 | tee "${LOG_DIR}/happy_parent2.log"
}

## Run `vcf_stats_report`
# The report is also created by postprocess_variants, but this separate runner
# script works on existing VCF files.
function run_vcf_stats_report() {
  local -r vcf_output="${1}"
  local -r stats_output="${2}"
  echo "Start running vcf_stats_report...Log will be in the terminal and also to ${LOG_DIR}/vcf_stats_report.log."
  sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    "${IMAGE}" \
    python /opt/deepvariant/vcf_stats_report.zip \
      --input_vcf "/output/${vcf_output}" \
      --outfile_base "/output/${stats_output}"
  echo "Done."
  echo
}

function run_all_vcf_stats_report() {
  (time run_vcf_stats_report "${OUTPUT_VCF_CHILD}" "${OUTPUT_STATS_CHILD}") 2>&1 | tee "${LOG_DIR}/vcf_stats_report_child.log"
  (time run_vcf_stats_report "${OUTPUT_VCF_PARENT1}" "${OUTPUT_STATS_PARENT1}") 2>&1 | tee "${LOG_DIR}/vcf_stats_report_parent1.log"
  (time run_vcf_stats_report "${OUTPUT_VCF_PARENT2}" "${OUTPUT_STATS_PARENT2}") 2>&1 | tee "${LOG_DIR}/vcf_stats_report_parent2.log"
}

function main() {
  echo 'Starting the test...'

  setup_test
  setup_args
  get_docker_image
  run_deeptrio
  run_glnexus
  # redacted
  #run_all_vcf_stats_report
  run_all_happy_reports
}

main "$@"
