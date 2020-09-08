#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/case-study"
MODEL_VERSION=redacted
# redacted
DEFAULT_CHILD_MODEL_HTTP_DIR=redacted
DEFAULT_PARENT_MODEL_HTTP_DIR=redacted

INPUT_DIR="${BASE}/input"
MODELS_DIR="${INPUT_DIR}/models"
CHILD_MODEL_DIR="${MODELS_DIR}/child"
PARENT_MODEL_DIR="${MODELS_DIR}/parent"
# redacted
CHILD_MODEL="${CHILD_MODEL_DIR}/model"
PARENT_MODEL="${PARENT_MODEL_DIR}/model"
DATA_DIR="${INPUT_DIR}/data"
REF="${DATA_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
BAM_CHILD="${DATA_DIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
BAM_PARENT1="${DATA_DIR}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
BAM_PARENT2="${DATA_DIR}/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
SAMPLE_NAME_CHILD="HG002"
SAMPLE_NAME_PARENT1="HG003"
SAMPLE_NAME_PARENT2="HG004"
TRUTH_VCF_CHILD="${DATA_DIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz"
TRUTH_BED_CHILD="${DATA_DIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed"
TRUTH_VCF_PARENT1="${DATA_DIR}/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
TRUTH_BED_PARENT1="${DATA_DIR}/HG003_GRCh38_1_22_v4.2_benchmark.bed"
TRUTH_VCF_PARENT2="${DATA_DIR}/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz"
TRUTH_BED_PARENT2="${DATA_DIR}/HG004_GRCh38_1_22_v4.2_benchmark.bed"

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
OUTPUT_VCF_CHILD="${OUTPUT_DIR}/HG002.output.vcf.gz"
OUTPUT_VCF_PARENT1="${OUTPUT_DIR}/HG003.output.vcf.gz"
OUTPUT_VCF_PARENT2="${OUTPUT_DIR}/HG004.output.vcf.gz"
OUTPUT_VCF_MERGED="${OUTPUT_DIR}/HG002-3-4.output.vcf.gz"
OUTPUT_GVCF_CHILD="${OUTPUT_DIR}/HG002.output.g.vcf.gz"
OUTPUT_GVCF_PARENT1="${OUTPUT_DIR}/HG003.output.g.vcf.gz"
OUTPUT_GVCF_PARENT2="${OUTPUT_DIR}/HG004.output.g.vcf.gz"
OUTPUT_STATS_CHILD="${OUTPUT_DIR}/HG002.output.vcf_stats"
OUTPUT_STATS_PARENT1="${OUTPUT_DIR}/HG003.output.vcf_stats"
OUTPUT_STATS_PARENT2="${OUTPUT_DIR}/HG004.output.vcf_stats"
LOG_DIR="${OUTPUT_DIR}/logs"

# Build binaries.
# If you're using the pre-built binaries, you can skip these and just run
# ./run-prereq.sh instead. And update the script to point to your *zip binaries.
function build_binaries() {
  ./build-prereq.sh
  ./build_release_binaries.sh
}

function copy_model() {
  # Expected params: model source directory, models destination directory

  model_http_dir="${1}"
  MODELS_DIR="${2}"

  ## Download models, and test data
  # Copy the model files to your local disk.
  HTTPS_ADDRESS="^https:\/\/.*"
  GS_ADDRESS="^gs:\/\/.*"
  if [[ $model_http_dir =~ $HTTPS_ADDRESS ]];
  then
    aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${model_http_dir}"/model.data-00000-of-00001
    aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${model_http_dir}"/model.index
    aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${model_http_dir}"/model.meta
  elif [[ $model_http_dir =~ $GS_ADDRESS ]];
  then
    gsutil cp "${model_http_dir}"/model.data-00000-of-00001 "${MODELS_DIR}"
    gsutil cp "${model_http_dir}"/model.index "${MODELS_DIR}"
    gsutil cp "${model_http_dir}"/model.meta "${MODELS_DIR}"
  else
    echo "Could not copy model. Unknown address prefix: ${model_http_dir}"
    exit 1
  fi
}

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
  sudo apt-get -y install parallel
  sudo apt-get -y install docker.io
  sudo apt-get -y install aria2
  sudo apt-get -y install bcftools
  sudo apt-get -y install tabix

  copy_model "${child_model_http_dir}" "${CHILD_MODEL_DIR}"
  copy_model "${parent_model_http_dir}" "${PARENT_MODEL_DIR}"

  # Copy the data, using http:// because of https://github.com/aria2/aria2/issues/1012.
  # Copy HG002 truth
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh38_1_22_v4.2_benchmark.bed -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi -d "${DATA_DIR}"
  # Copy HG003 truth
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.bed -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi -d "${DATA_DIR}"
  # Copy HG004 truth
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004_GRCh38_1_22_v4.2_benchmark.bed -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi -d "${DATA_DIR}"

  # Add a checksum for the biggest file:
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.fai -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.gzi -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gzi -d "${DATA_DIR}"
  aria2c -c -x10 -s10 http://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai -d "${DATA_DIR}"
}

## Run `make_examples`
function run_make_examples() {
  echo "Start running make_examples...Log will be in the terminal and also to ${LOG_DIR}/make_examples.log."
  seq 0 $((N_SHARDS-1)) | \
    parallel --halt 2 --line-buffer \
      python ./bazel-bin/deeptrio/make_examples.zip \
        --mode calling \
        --ref "${REF}" \
        --reads "${BAM_CHILD}" \
        --reads_parent1 "${BAM_PARENT1}" \
        --reads_parent2 "${BAM_PARENT2}" \
        --sample_name "${SAMPLE_NAME_CHILD}" \
        --sample_name_parent1 "${SAMPLE_NAME_PARENT1}" \
        --sample_name_parent2 "${SAMPLE_NAME_PARENT2}" \
        --examples "${EXAMPLES_PREFIX}" \
        --pileup_image_height_child 60 \
        --pileup_image_height_parent 40 \
        --gvcf "${GVCF_TFRECORDS_PREFIX}" \
        --task {}
  echo "Done."
  echo
}

## Run `call_variants`
function run_call_variants() {
  echo "Start running call_variants...Log will be in the terminal and also to ${LOG_DIR}/call_variants.log."
  local -r examples="${1}"
  local -r model="${2}"
  local -r cvo_output="${3}"
  python ./bazel-bin/deepvariant/call_variants.zip \
      --outfile "${cvo_output}" \
      --examples "${examples}" \
      --checkpoint "${model}"
  echo "Done."
  echo
}

## Run `postprocess_variants`, without gVCFs.
function run_postprocess_variants() {
  echo "Start running postprocess_variants (without gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.log."
  local -r cvo_output="${1}"
  local -r vcf_output="${2}"
  python ./bazel-bin/deepvariant/postprocess_variants.zip \
      --ref "${REF}" \
      --infile "${cvo_output}" \
      --outfile "${vcf_output}"
  echo "Done."
  echo
}

## Run `postprocess_variants`, with gVCFs.
function run_postprocess_variants_gVCF() {
  echo "Start running postprocess_variants (with gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.withGVCF.log."
  local -r cvo_output="${1}"
  local -r vcf_output="${2}"
  local -r gvcf_output="${3}"
  local -r gvcf_tfrecords="${4}"

  python ./bazel-bin/deepvariant/postprocess_variants.zip \
      --ref "${REF}" \
      --infile "${cvo_output}" \
      --outfile "${vcf_output}" \
      --nonvariant_site_tfrecord_path "${gvcf_tfrecords}" \
      --gvcf_outfile "${gvcf_output}"
  echo "Done."
  echo
}

## Run `vcf_stats_report`
# The report is also created by postprocess_variants, but this separate runner
# script works on existing VCF files.
function run_vcf_stats_report() {
  local -r vcf_output="${1}"
  local -r stats_output="${2}"
  echo "Start running vcf_stats_report...Log will be in the terminal and also to ${LOG_DIR}/vcf_stats_report.log."
  python ./bazel-bin/deepvariant/vcf_stats_report.zip \
      --input_vcf "${vcf_output}" \
      --outfile_base "${stats_output}"
  echo "Done."
  echo
}

## Evaluation: run hap.py
function run_happy() {
  local -r truth_vcf="${1}"
  local -r truth_bed="${2}"
  local -r vcf_output="${3}"

  echo "Start evaluation with hap.py..."
  UNCOMPRESSED_REF="${DATA_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

  # hap.py cannot read the compressed fa, so uncompress
  # into a writable directory. Index file was downloaded earlier.
  zcat <"${REF}" >"${UNCOMPRESSED_REF}"

  # Pulling twice in case the first one times out.
  sudo docker pull pkrusche/hap.py || \
    (sleep 5 ; sudo docker pull pkrusche/hap.py)
  sudo docker run -i \
  -v "${DATA_DIR}:${DATA_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
    "${truth_vcf}" \
    "${vcf_output}" \
    -f "${truth_bed}" \
    -r "${UNCOMPRESSED_REF}" \
    -o "${OUTPUT_DIR}/happy.output" \
    --engine=vcfeval
  echo "Done."
}

function run_glnexus() {
  sudo docker pull quay.io/mlin/glnexus:v1.2.7 || \
    (sleep 5 ; sudo docker pull quay.io/mlin/glnexus:v1.2.7)

  time sudo docker run \
    -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
    quay.io/mlin/glnexus:v1.2.7 \
    /usr/local/bin/glnexus_cli \
    --config DeepVariantWGS \
    "${OUTPUT_GVCF_PARENT2}" "${OUTPUT_GVCF_PARENT1}" "${OUTPUT_GVCF_CHILD}" \
    | bcftools view - | bgzip -c > "${OUTPUT_VCF_MERGED}"
}

function run_all_call_variants() {
  (time run_call_variants "${EXAMPLES_CHILD}" "${CHILD_MODEL}" "${CALL_VARIANTS_OUTPUT_CHILD}") >> "${LOG_DIR}/call_variants.log" 2>&1
  (time run_call_variants "${EXAMPLES_PARENT1}" "${PARENT_MODEL}" "${CALL_VARIANTS_OUTPUT_PARENT1}") >> "${LOG_DIR}/call_variants.log" 2>&1
  (time run_call_variants "${EXAMPLES_PARENT2}" "${PARENT_MODEL}" "${CALL_VARIANTS_OUTPUT_PARENT2}") >> "${LOG_DIR}/call_variants.log" 2>&1
}

function run_all_postprocess_variants() {
  (time run_postprocess_variants "${CALL_VARIANTS_OUTPUT_CHILD}" "${OUTPUT_VCF_CHILD}") >> "${LOG_DIR}/postprocess_variants.log" 2>&1
  (time run_postprocess_variants "${CALL_VARIANTS_OUTPUT_PARENT1}" "${OUTPUT_VCF_PARENT1}") >> "${LOG_DIR}/postprocess_variants.log" 2>&1
  (time run_postprocess_variants "${CALL_VARIANTS_OUTPUT_PARENT2}" "${OUTPUT_VCF_PARENT2}") >> "${LOG_DIR}/postprocess_variants.log" 2>&1
}

function run_all_postprocess_variants_gVCF() {
  (time run_postprocess_variants_gVCF "${CALL_VARIANTS_OUTPUT_CHILD}" "${OUTPUT_VCF_CHILD}" "${OUTPUT_GVCF_CHILD}" "${GVCF_TFRECORDS_CHILD}") >> "${LOG_DIR}/postprocess_variants.withGVCF.log" 2>&1
  (time run_postprocess_variants_gVCF "${CALL_VARIANTS_OUTPUT_PARENT1}" "${OUTPUT_VCF_PARENT1}" "${OUTPUT_GVCF_PARENT1}" "${GVCF_TFRECORDS_PARENT1}") >> "${LOG_DIR}/postprocess_variants.withGVCF.log" 2>&1
  (time run_postprocess_variants_gVCF "${CALL_VARIANTS_OUTPUT_PARENT2}" "${OUTPUT_VCF_PARENT2}" "${OUTPUT_GVCF_PARENT2}" "${GVCF_TFRECORDS_PARENT2}") >> "${LOG_DIR}/postprocess_variants.withGVCF.log" 2>&1
}

function run_all_vcf_stats_report() {
  (time run_vcf_stats_report "${OUTPUT_VCF_CHILD}" "${OUTPUT_STATS_CHILD}") 2>&1 | tee "${LOG_DIR}/vcf_stats_report_child.log"
  (time run_vcf_stats_report "${OUTPUT_VCF_PARENT1}" "${OUTPUT_STATS_PARENT1}") 2>&1 | tee "${LOG_DIR}/vcf_stats_report_parent1.log"
  (time run_vcf_stats_report "${OUTPUT_VCF_PARENT2}" "${OUTPUT_STATS_PARENT2}") 2>&1 | tee "${LOG_DIR}/vcf_stats_report_parent2.log"
}

function run_all_happy_reports() {
  run_happy "${TRUTH_VCF_CHILD}" "${TRUTH_BED_CHILD}" "${OUTPUT_VCF_CHILD}" 2>&1 | tee "${LOG_DIR}/happy_child.log"
  run_happy "${TRUTH_VCF_PARENT1}" "${TRUTH_BED_PARENT1}" "${OUTPUT_VCF_PARENT1}" 2>&1 | tee "${LOG_DIR}/happy_parent1.log"
  run_happy "${TRUTH_VCF_PARENT2}" "${TRUTH_BED_PARENT2}" "${OUTPUT_VCF_PARENT2}" 2>&1 | tee "${LOG_DIR}/happy_parent2.log"
}

function main() {
  echo 'Starting the test...'
  local -r child_model_http_dir="${1:-$DEFAULT_CHILD_MODEL_HTTP_DIR}"
  echo "Using child model from: ${child_model_http_dir}"
  local -r parent_model_http_dir="${1:-$DEFAULT_PARENT_MODEL_HTTP_DIR}"
  echo "Using parent model from: ${parent_model_http_dir}"

  build_binaries
  setup_test
  (time run_make_examples) 2>&1 | tee "${LOG_DIR}/make_examples.log"
  (time run_all_call_variants) >> "${LOG_DIR}/call_variants.log" 2>&1
  (time run_all_postprocess_variants) >> "${LOG_DIR}/postprocess_variants.log" 2>&1
  (time run_all_postprocess_variants_gVCF) >> "${LOG_DIR}/postprocess_variants.withGVCF.log" 2>&1
  run_glnexus
  run_all_vcf_stats_report
  run_all_happy_reports
}

main "$@"
