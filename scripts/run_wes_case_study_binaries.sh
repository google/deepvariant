#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/exome-case-study"
MODEL_VERSION="0.7.2"
MODEL_NAME="DeepVariant-inception_v3-${MODEL_VERSION}+data-wes_standard"
DEFAULT_MODEL_HTTP_DIR="https://storage.googleapis.com/deepvariant/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}"

INPUT_DIR="${BASE}/input"
MODELS_DIR="${INPUT_DIR}/models"
MODEL="${MODELS_DIR}/model.ckpt"
DATA_DIR="${INPUT_DIR}/data"
REF="${DATA_DIR}/hs37d5.fa.gz"
BAM="${DATA_DIR}/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
TRUTH_VCF="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/HG002.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="${OUTPUT_DIR}/HG002.output.vcf.gz"
OUTPUT_GVCF="${OUTPUT_DIR}/HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

CAPTURE_BED="${DATA_DIR}/agilent_sureselect_human_all_exon_v5_b37_targets.bed"

# Build binaries.
# If you're using the pre-built binaries, you can skip these and just run
# ./run-prereq.sh instead. And update the script to point to your *zip binaries.
function build_binaries() {
  ./build-prereq.sh
  ./build_release_binaries.sh
}

function setup_test() {
  ## Create local directory structure
  mkdir -p "${OUTPUT_DIR}"
  mkdir -p "${DATA_DIR}"
  mkdir -p "${MODELS_DIR}"
  mkdir -p "${LOG_DIR}"

  ## Download extra packages
  # There are some extra programs we will need.
  # We are going to use [GNU Parallel](https://www.gnu.org/software/parallel/) to
  # run `make_examples`.
  sudo apt-get -y update
  sudo apt-get -y install parallel
  sudo apt-get -y install docker.io
  sudo apt-get -y install aria2

  ## Download models, and test data
  # Copy the model files to your local disk.
  HTTPS_ADDRESS="^https:\/\/.*"
  GS_ADDRESS="^gs:\/\/.*"
  if [[ $model_http_dir =~ $HTTPS_ADDRESS ]];
  then
    aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${model_http_dir}"/model.ckpt.data-00000-of-00001
    aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${model_http_dir}"/model.ckpt.index
    aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${model_http_dir}"/model.ckpt.meta
  elif [[ $model_http_dir =~ $GS_ADDRESS ]];
  then
    gsutil cp "${model_http_dir}"/model.ckpt.data-00000-of-00001 "${MODELS_DIR}"
    gsutil cp "${model_http_dir}"/model.ckpt.index "${MODELS_DIR}"
    gsutil cp "${model_http_dir}"/model.ckpt.meta "${MODELS_DIR}"
  else
    echo 'Could not copy model. Unknown address prefix.'
  fi

  # Copy the data
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/agilent_sureselect_human_all_exon_v5_b37_targets.bed
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz.fai
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz.gzi
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gzi
  aria2c -c -x10 -s10 -d "${DATA_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.fai
}

## Run `make_examples`
# In this step, we used the `--regions` flag to constrain the regions we processed
# to the capture region BED file:
function run_make_examples() {
  echo "Start running make_examples...Log will be in the terminal and also to ${LOG_DIR}/make_examples.log."
  seq 0 $((N_SHARDS-1)) | \
    parallel -k --line-buffer \
      python ./bazel-bin/deepvariant/make_examples.zip \
        --mode calling \
        --ref "${REF}" \
        --reads "${BAM}" \
        --examples "${EXAMPLES}" \
        --regions "${CAPTURE_BED}" \
        --gvcf "${GVCF_TFRECORDS}" \
        --task {}
  echo "Done."
  echo
}

## Run `call_variants`
function run_call_variants() {
  echo "Start running call_variants...Log will be in the terminal and also to ${LOG_DIR}/call_variants.log."
  python ./bazel-bin/deepvariant/call_variants.zip \
      --outfile "${CALL_VARIANTS_OUTPUT}" \
      --examples "${EXAMPLES}" \
      --checkpoint "${MODEL}"
  echo "Done."
  echo
}

## Run `postprocess_variants`, without gVCFs.
function run_postprocess_variants() {
  echo "Start running postprocess_variants (without gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.log."
  python ./bazel-bin/deepvariant/postprocess_variants.zip \
      --ref "${REF}" \
      --infile "${CALL_VARIANTS_OUTPUT}" \
      --outfile "${OUTPUT_VCF}"
  echo "Done."
  echo
}

## Run `postprocess_variants`, with gVCFs.
function run_postprocess_variants_gVCF() {
  echo "Start running postprocess_variants (with gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.withGVCF.log."
  python ./bazel-bin/deepvariant/postprocess_variants.zip \
      --ref "${REF}" \
      --infile "${CALL_VARIANTS_OUTPUT}" \
      --outfile "${OUTPUT_VCF}" \
      --nonvariant_site_tfrecord_path "${GVCF_TFRECORDS}" \
      --gvcf_outfile "${OUTPUT_GVCF}"
  echo "Done."
  echo
}

## Evaluation: run hap.py
function run_happy() {
  echo "Start evaluation with hap.py..."
  UNCOMPRESSED_REF="${DATA_DIR}/hs37d5.fa"

  # hap.py cannot read the compressed fa, so uncompress
  # into a writable directory. Index file was downloaded earlier.
  zcat <"${REF}" >"${UNCOMPRESSED_REF}"

  sudo docker pull pkrusche/hap.py
  sudo docker run -i \
  -v "${DATA_DIR}:${DATA_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
    "${TRUTH_VCF}" \
    "${OUTPUT_VCF}" \
    -f "${TRUTH_BED}" \
    -T "${CAPTURE_BED}" \
    -r "${UNCOMPRESSED_REF}" \
    -o "${OUTPUT_DIR}/happy.output" \
    --engine=vcfeval
  echo "Done."
}

function main() {
  echo 'Starting the test...'
  local -r model_http_dir="${1:-$DEFAULT_MODEL_HTTP_DIR}"
  echo "Using model from: ${model_http_dir}"

  build_binaries
  setup_test
  (time run_make_examples) 2>&1 | tee "${LOG_DIR}/make_examples.log"
  (time run_call_variants) 2>&1 | tee "${LOG_DIR}/call_variants.log"
  (time run_postprocess_variants) 2>&1 | tee "${LOG_DIR}/postprocess_variants.log"
  (time run_postprocess_variants_gVCF) 2>&1 | tee "${LOG_DIR}/postprocess_variants.withGVCF.log"
  run_happy 2>&1 | tee "${LOG_DIR}/happy.log"
}

main "$@"
