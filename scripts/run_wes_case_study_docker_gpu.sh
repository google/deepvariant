#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/exome-case-study"
BIN_VERSION="0.9.0"

INPUT_DIR="${BASE}/input/data"
REF="hs37d5.fa.gz"
BAM="151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
TRUTH_VCF="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS="16"

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/HG002.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="HG002.output.vcf.gz"
OUTPUT_GVCF="HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

CAPTURE_BED="agilent_sureselect_human_all_exon_v5_b37_targets.bed"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"
mkdir -p "${LOG_DIR}"

## Download extra packages
# Install aria2 to download data files.
sudo apt-get -qq -y update
sudo apt-get -qq -y install aria2

if ! hash nvidia-docker 2>/dev/null; then
  echo "'nvidia-docker' was not found in PATH. Installing nvidia-docker..."
  ./scripts/install_nvidia_docker.sh
fi

# Copy the data
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/agilent_sureselect_human_all_exon_v5_b37_targets.bed
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz.fai
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz.gzi
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gzi
aria2c -c -x10 -s10 -d "${INPUT_DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.fai

## Pull the docker image.
sudo nvidia-docker pull google/deepvariant:"${BIN_VERSION}-gpu"

echo "Run DeepVariant..."
sudo nvidia-docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --regions="/input/${CAPTURE_BED}" \
  --output_vcf="/output/${OUTPUT_VCF}" \
  --output_gvcf="/output/${OUTPUT_GVCF}" \
  --num_shards=${N_SHARDS}
echo "Done."
echo

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${INPUT_DIR}/hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory. Index file was downloaded earlier.
zcat <"${INPUT_DIR}/${REF}" >"${UNCOMPRESSED_REF}"

sudo nvidia-docker pull pkrusche/hap.py
( sudo nvidia-docker run -i \
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
