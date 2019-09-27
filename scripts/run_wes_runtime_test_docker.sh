#!/bin/bash
# Copyright 2019 Google LLC.
# This script builds a Docker image and runs DeepVariant for whole exome.
# Main purpose of this script is to evaluate the total runtime of DeepVariant on
# different computer (cloud instance) types.
# Runtime measurements do not include the time for building the Docker image and
# localizing test data.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/exome-case-study"

INPUT_DIR="${BASE}/input"
DATA_DIR="${INPUT_DIR}/data"
REF="${DATA_DIR}/hs37d5.fa.gz"
BAM="${DATA_DIR}/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
TRUTH_VCF="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS=$(nproc)

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/make_examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/call_variants_output.tfrecord.gz"
OUTPUT_VCF="${OUTPUT_DIR}/HG002.output.vcf.gz"
OUTPUT_GVCF="${OUTPUT_DIR}/HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

CAPTURE_BED="${DATA_DIR}/agilent_sureselect_human_all_exon_v5_b37_targets.bed"

# Build Docker image.
function build_docker_image() {
  echo "Start building Docker image."
  sudo docker build -t deepvariant .
  echo "Done building Docker image."
}

function setup_test() {
  ## Create local directory structure
  mkdir -p "${OUTPUT_DIR}"
  mkdir -p "${DATA_DIR}"
  mkdir -p "${LOG_DIR}"

  ## Download extra packages
  sudo apt-get -y update
  sudo apt-get -y install docker.io
  sudo apt-get -y install aria2

  # Copy the data
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/agilent_sureselect_human_all_exon_v5_b37_targets.bed
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz.fai
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz.gzi
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gzi
  aria2c -c -x10 -s10 -d "${DATA_DIR}" http://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.fai
}

function run_deepvariant() {
  echo "Start running run_deepvariant...Log will be in the terminal and also to ${LOG_DIR}/deepvariant_runtime.log."
  sudo docker run \
    -v "${INPUT_DIR}":"${INPUT_DIR}" \
    -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
    deepvariant:latest \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type "WES" \
    --num_shards "${N_SHARDS}" \
    --output_gvcf "${OUTPUT_GVCF}" \
    --output_vcf "${OUTPUT_VCF}" \
    --reads "${BAM}" \
    --ref "${REF}" \
    --regions "${CAPTURE_BED}"
  echo "Done."
  echo
}

function main() {
  echo 'Starting the test...'

  setup_test
  build_docker_image

  (time run_deepvariant) 2>&1 | tee "${LOG_DIR}/deepvariant_runtime.log"
}

main "$@"
