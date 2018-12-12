#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BASE="${HOME}/case-study"
BIN_VERSION="0.7.2"
MODEL_VERSION="0.7.2"
MODEL_NAME="DeepVariant-inception_v3-${MODEL_VERSION}+data-wgs_standard"
MODEL_HTTP_DIR="https://storage.googleapis.com/deepvariant/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}"

INPUT_DIR="${BASE}/input"
MODELS_DIR="${INPUT_DIR}/models"
MODEL="${MODELS_DIR}/model.ckpt"
DATA_DIR="${INPUT_DIR}/data"
REF="${DATA_DIR}/hs37d5.fa.gz"
BAM="${DATA_DIR}/HG002_NIST_150bp_50x.bam"
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

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${MODELS_DIR}"
mkdir -p "${LOG_DIR}"


## Download extra packages
# There are some extra programs we will need.
# We are going to use [GNU Parallel](https://www.gnu.org/software/parallel/) to
# run `make_examples`. We are going to install `samtools` and `docker.io` to help
# do some analysis at the end.
sudo apt-get -y update
sudo apt-get -y install parallel
sudo apt-get -y install samtools
sudo apt-get -y install docker.io
sudo apt-get -y install aria2


## Download models, and test data
# Copy the model files to your local disk.
aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${MODEL_HTTP_DIR}"/model.ckpt.data-00000-of-00001
aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${MODEL_HTTP_DIR}"/model.ckpt.index
aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${MODEL_HTTP_DIR}"/model.ckpt.meta

# Copy the data
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam.bai -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.fai -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.gzi -d "${DATA_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gzi -d "${DATA_DIR}"

## Pull the docker image.
sudo docker pull gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"

## Run `make_examples`
echo "Start running make_examples...Log will be in the terminal and also to ${LOG_DIR}/make_examples.log."
( time seq 0 $((N_SHARDS-1)) | \
  parallel -k --line-buffer \
    sudo docker run \
      -v "${BASE}":"${BASE}" \
      gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
      /opt/deepvariant/bin/make_examples \
      --mode calling \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${EXAMPLES}" \
      --gvcf "${GVCF_TFRECORDS}" \
      --task {} \
) 2>&1 | tee "${LOG_DIR}/make_examples.log"
echo "Done."
echo

## Run `call_variants`
echo "Start running call_variants...Log will be in the terminal and also to ${LOG_DIR}/call_variants.log."
( time sudo docker run \
    -v "${BASE}":"${BASE}" \
    gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
    /opt/deepvariant/bin/call_variants \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --examples "${EXAMPLES}" \
    --checkpoint "${MODEL}"
) 2>&1 | tee "${LOG_DIR}/call_variants.log"
echo "Done."
echo

## Run `postprocess_variants`, without gVCFs.
echo "Start running postprocess_variants (without gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.log."
( time sudo docker run \
    -v "${BASE}":"${BASE}" \
    gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
    /opt/deepvariant/bin/postprocess_variants \
    --ref "${REF}" \
    --infile "${CALL_VARIANTS_OUTPUT}" \
    --outfile "${OUTPUT_VCF}"
) 2>&1 | tee "${LOG_DIR}/postprocess_variants.log"
echo "Done."
echo

## Run `postprocess_variants`, with gVCFs.
echo "Start running postprocess_variants (with gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.withGVCF.log."
( time sudo docker run \
    -v "${BASE}":"${BASE}" \
    gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
    /opt/deepvariant/bin/postprocess_variants \
    --ref "${REF}" \
    --infile "${CALL_VARIANTS_OUTPUT}" \
    --outfile "${OUTPUT_VCF}" \
    --nonvariant_site_tfrecord_path "${GVCF_TFRECORDS}" \
    --gvcf_outfile "${OUTPUT_GVCF}"
) 2>&1 | tee "${LOG_DIR}/postprocess_variants.withGVCF.log"
echo "Done."
echo

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${OUTPUT_DIR}/hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory and index it.
zcat <"${REF}" >"${UNCOMPRESSED_REF}"
samtools faidx "${UNCOMPRESSED_REF}"

sudo docker pull pkrusche/hap.py
( sudo docker run -i \
-v "${DATA_DIR}:${DATA_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${OUTPUT_VCF}" \
  -f "${TRUTH_BED}" \
  -r "${UNCOMPRESSED_REF}" \
  -o "${OUTPUT_DIR}/happy.output" \
  --engine=vcfeval
) 2>&1 | tee "${LOG_DIR}/happy.log"
echo "Done."
