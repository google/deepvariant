#!/bin/bash
model=${1}

INPUT_DIR="${PWD}/input"
REF_OUTPUT_DIR="${PWD}/reference"
docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${REF_OUTPUT_DIR}:/output" \
    deepvariant \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${model} \
    --ref=/input/GRCh38_no_alt_analysis_set.fasta \
    --reads=/input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
    --regions "chr20:10,000,000-10,010,000" \
    --output_vcf=/output/HG002.output.vcf.gz \
    --output_gvcf=/output/HG002.output.g.vcf.gz \
    --call_variants_extra_args="use_openvino=False" \
    --num_shards=$(nproc)

OUTPUT_DIR="${PWD}/quickstart-output-ovino"
docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    deepvariant \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${model} \
    --ref=/input/GRCh38_no_alt_analysis_set.fasta \
    --reads=/input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
    --regions "chr20:10,000,000-10,010,000" \
    --output_vcf=/output/HG002.output.vcf.gz \
    --output_gvcf=/output/HG002.output.g.vcf.gz \
    --call_variants_extra_args="use_openvino=True" \
    --num_shards=$(nproc)

files=$(ls $REF_OUTPUT_DIR)
for f in $files; do
    if [[ $(cmp $REF_OUTPUT_DIR/$f $OUTPUT_DIR/$f) ]]; then
        echo "$REF_OUTPUT_DIR/$f and $OUTPUT_DIR/$f are different!"
        exit 1
    fi
done
sudo rm -r $REF_OUTPUT_DIR $OUTPUT_DIR

