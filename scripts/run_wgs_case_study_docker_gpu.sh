#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/case-study"
BIN_VERSION="0.8.0"

INPUT_DIR="${BASE}/input/data"
REF="hs37d5.fa.gz"
BAM="HG002_NIST_150bp_50x.bam"
TRUTH_VCF="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS="16"

OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="HG002.output.vcf.gz"
OUTPUT_GVCF="HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

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
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam.bai -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.fai -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.gzi -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gzi -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.fai -d "${INPUT_DIR}"

## Pull the docker image.
sudo nvidia-docker pull gcr.io/deepvariant-docker/deepvariant_gpu:"${BIN_VERSION}"

echo "Run DeepVariant..."
sudo nvidia-docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  gcr.io/deepvariant-docker/deepvariant_gpu:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/${OUTPUT_VCF} \
  --output_gvcf=/output/${OUTPUT_GVCF} \
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
  -r "${UNCOMPRESSED_REF}" \
  -o "${OUTPUT_DIR}/happy.output" \
  --engine=vcfeval
) 2>&1 | tee "${LOG_DIR}/happy.log"
echo "Done."
