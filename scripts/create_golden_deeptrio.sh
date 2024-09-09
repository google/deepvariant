#!/bin/bash
# Copyright 2024 Google LLC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Usage:
# ./scripts/create_golden_deeptrio.sh
# to generate data under deeptrio/testdata.
#
# See scripts/create_golden.sh for more details.
#
# NOTE: Please run until you see this line:
#
# ================ create_golden_deeptrio.sh completed! ================
#
# Otherwise, some commands might have failed in between.

set -euo pipefail
set -x

TESTDATA_DIR=$(pwd)/deeptrio/testdata

# First, clean up all files on the top directory.

find "${TESTDATA_DIR}" -maxdepth 1 -type f -delete

# Inputs
REF="${TESTDATA_DIR}"/input/hs37d5.chr20.fa.gz
READS="${TESTDATA_DIR}"/input/HG001.chr20.10_10p1mb_sorted.bam
READS_PARENT1="${TESTDATA_DIR}"/input/NA12891.chr20.10_10p1mb_sorted.bam
READS_PARENT2="${TESTDATA_DIR}"/input/NA12892.chr20.10_10p1mb_sorted.bam
VARIANTS="${TESTDATA_DIR}"/input/test_hg001_giab_grch37_chr20_100kbp_at_10mb.vcf.gz
# This is a simplified version, manually modified from
# test_nist.b37_chr20_100kbp_at_10mb.vcf.gz. "type=1" and "type=2" are manually
# added to the INFO field.
REGIONS="${TESTDATA_DIR}"/input/test_giab.b37_chr20_100kbp_at_10mb.bed
WITH_TYPE_VARIANTS="${TESTDATA_DIR}"/input/with_types.test_nist.b37_chr20_4kbp_at_10mb.vcf.gz

# Outputs
GOLDEN_TRAINING_EXAMPLES="${TESTDATA_DIR}"/golden.training_examples.tfrecord.gz
GOLDEN_TRAINING_EXAMPLES_VCF="${TESTDATA_DIR}"/golden.training_examples.vcf
GOLDEN_CALLING_CANDIDATES="${TESTDATA_DIR}"/golden.calling_candidates.tfrecord.gz
# This file doesn't get generated as is.
# Instead, golden_{child,parent1,parent2} are generated.
GOLDEN_CALLING_EXAMPLES="${TESTDATA_DIR}"/golden.calling_examples.tfrecord.gz
CALLING_EXAMPLES_TEMP=${TESTDATA_DIR}/tmp.examples.tfrecord.gz
GOLDEN_CALLING_EXAMPLES_CHILD="${TESTDATA_DIR}"/golden_child.calling_examples.tfrecord.gz
GOLDEN_CANDIDATE_POSITIONS="${TESTDATA_DIR}"/golden_child.candidate_positions
GOLDEN_POSTPROCESS_INPUT="${TESTDATA_DIR}"/golden.postprocess_single_site_input.tfrecord.gz
GOLDEN_POSTPROCESS_OUTPUT="${TESTDATA_DIR}"/golden.postprocess_single_site_output.vcf
# This file doesn't get generated as is.
# Instead, golden_{child,parent1,parent2} are generated.
GOLDEN_POSTPROCESS_GVCF_INPUT="${TESTDATA_DIR}"/golden.postprocess_gvcf_input.tfrecord.gz
GOLDEN_POSTPROCESS_GVCF_INPUT_CHILD="${TESTDATA_DIR}"/golden_child.postprocess_gvcf_input.tfrecord.gz
GOLDEN_POSTPROCESS_GVCF_OUTPUT="${TESTDATA_DIR}"/golden.postprocess_gvcf_output.g.vcf

MODEL=gs://deepvariant/models/DeepTrio/1.6.0/savedmodels/deeptrio.wgs_child.savedmodel
# Speed up by copying to /tmp/
rm -rf /tmp/deeptrio.wgs_child.savedmodel
gsutil -m cp -R ${MODEL} /tmp/
MODEL=/tmp/deeptrio.wgs_child.savedmodel

source settings.sh

# Need to run this first, otherwise the bazel build command below will complain
# about deepvariant/examples_from_stream.so.
./build_release_binaries.sh

# shellcheck disable=SC2086
bazel build -c opt ${DV_COPT_FLAGS} //deeptrio:make_examples \
  //deepvariant:make_examples \
  //deepvariant:call_variants \
  //deepvariant:postprocess_variants \
  //deepvariant/labeler:labeled_examples_to_vcf

# Makes the training output.
time ./bazel-bin/deeptrio/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --sample_name_to_train "child" \
  --regions "20:10,000,000-10,010,000" \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${GOLDEN_TRAINING_EXAMPLES}" \
  --channel_list='BASE_CHANNELS,insert_size' \
  --write_run_info \
  --deterministic_serialization

# Makes the training output; 3 shards.
for i in $(seq 0 2);
do
time ./bazel-bin/deeptrio/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --sample_name_to_train "child" \
  --regions "20:10,000,000-10,010,000" \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${GOLDEN_TRAINING_EXAMPLES}"@3 \
  --channel_list='BASE_CHANNELS,insert_size' \
  --nowrite_run_info \
  --task "${i}" \
  --deterministic_serialization
done

# Now make the calling output.
time ./bazel-bin/deeptrio/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --regions "20:10,000,000-10,010,000" \
  --confident_regions "${REGIONS}" \
  --candidates "${GOLDEN_CALLING_CANDIDATES}" \
  --examples "${GOLDEN_CALLING_EXAMPLES}" \
  --channel_list='BASE_CHANNELS,insert_size' \
  --gvcf "${GOLDEN_POSTPROCESS_GVCF_INPUT}" \
  --nowrite_run_info \
  --deterministic_serialization

time ./bazel-bin/deeptrio/make_examples \
  --mode candidate_sweep \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --regions "20:10,000,000-10,010,000" \
  --confident_regions "${REGIONS}" \
  --examples "${CALLING_EXAMPLES_TEMP}" \
  --channel_list='BASE_CHANNELS,insert_size' \
  --candidate_positions "${GOLDEN_CANDIDATE_POSITIONS}" \
  --nowrite_run_info \
  --deterministic_serialization

# Now make the calling output; 3 shards.
for i in $(seq 0 2);
do
time ./bazel-bin/deeptrio/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --regions "20:10,000,000-10,010,000" \
  --confident_regions "${REGIONS}" \
  --examples "${GOLDEN_CALLING_EXAMPLES}"@3 \
  --channel_list='BASE_CHANNELS,insert_size' \
  --gvcf "${GOLDEN_POSTPROCESS_GVCF_INPUT}"@3 \
  --nowrite_run_info \
  --task "${i}" \
  --deterministic_serialization

# Generate golden set for candidate_sweep mode
time ./bazel-bin/deeptrio/make_examples \
  --mode candidate_sweep \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --regions "20:10,000,000-10,010,000" \
  --confident_regions "${REGIONS}" \
  --examples "${CALLING_EXAMPLES_TEMP}"@3 \
  --channel_list='BASE_CHANNELS,insert_size' \
  --candidate_positions "${GOLDEN_CANDIDATE_POSITIONS}"@3 \
  --nowrite_run_info \
  --task "${i}" \
  --deterministic_serialization
done

# This calls `call_variants` and generate the input for `postprocess_variants`.
time ./bazel-bin/deepvariant/call_variants \
  --outfile "${GOLDEN_POSTPROCESS_INPUT}" \
  --examples "${GOLDEN_CALLING_EXAMPLES_CHILD}" \
  --checkpoint "${MODEL}"

# Finally, we run `postprocess_variants`.
time ./bazel-bin/deepvariant/postprocess_variants \
  --infile "${GOLDEN_POSTPROCESS_INPUT}" \
  --outfile "${GOLDEN_POSTPROCESS_OUTPUT}.gz" \
  --ref "${REF}" \
  --nonvariant_site_tfrecord_path "${GOLDEN_POSTPROCESS_GVCF_INPUT_CHILD}" \
  --gvcf_outfile "${GOLDEN_POSTPROCESS_GVCF_OUTPUT}" \
  --cpus 0

# The VCF version is also used in our unit test.
zcat "${GOLDEN_POSTPROCESS_OUTPUT}.gz" > "${GOLDEN_POSTPROCESS_OUTPUT}"

time ./bazel-bin/deepvariant/labeler/labeled_examples_to_vcf \
  --ref "${REF}" \
  --examples "${GOLDEN_TRAINING_EXAMPLES}" \
  --output_vcf "${GOLDEN_TRAINING_EXAMPLES_VCF}"


####### Data for CustomizedClassesVariantLabler #######
CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES="${TESTDATA_DIR}"/customized_classes.golden.training_examples.tfrecord.gz
time ./bazel-bin/deeptrio/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --sample_name_to_train "child" \
  --regions "20:10,000,000-10,004,000" \
  --truth_variants "${WITH_TYPE_VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES}" \
  --channel_list='BASE_CHANNELS,insert_size' \
  --labeler_algorithm "customized_classes_labeler" \
  --customized_classes_labeler_classes_list "ref,class1,class2" \
  --customized_classes_labeler_info_field_name "type" \
  --nowrite_run_info \
  --deterministic_serialization

####### Data for running with alt_aligned_pileup #######
ALT_ALIGNED_PILEUP_IMAGE_GOLDEN_TRAINING_EXAMPLES="${TESTDATA_DIR}"/alt_aligned_pileup.golden.training_examples.tfrecord.gz
time ./bazel-bin/deeptrio/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --sample_name_to_train "child" \
  --regions "20:10,000,000-10,010,000" \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${ALT_ALIGNED_PILEUP_IMAGE_GOLDEN_TRAINING_EXAMPLES}" \
  --channel_list='BASE_CHANNELS' \
  --alt_aligned_pileup "diff_channels" \
  --pileup_image_height_child=60 \
  --pileup_image_height_parent=40 \
  --pileup_image_width=199 \
  --nowrite_run_info \
  --deterministic_serialization


GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES="${TESTDATA_DIR}"/golden.vcf_candidate_importer.training_examples.tfrecord.gz
time ./bazel-bin/deeptrio/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --sample_name_to_train "child" \
  --regions "20:10,000,000-10,010,000" \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --variant_caller="vcf_candidate_importer" \
  --examples "${GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES}" \
  --channel_list='BASE_CHANNELS,insert_size' \
  --pileup_image_height_child=60 \
  --pileup_image_height_parent=40 \
  --nowrite_run_info \
  --deterministic_serialization

# Just for the sake of creating test data, the --proposed_variants_parent{1,2}
# flags were set to the same VCF file. Depending on the use case, they are
# likely to be different files in real scenarios.
GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES="${TESTDATA_DIR}"/golden.vcf_candidate_importer.calling_examples.tfrecord.gz
time ./bazel-bin/deeptrio/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --reads_parent1 "${READS_PARENT1}" \
  --reads_parent2 "${READS_PARENT2}" \
  --sample_name "child" \
  --sample_name_parent1 "parent1" \
  --sample_name_parent2 "parent2" \
  --regions "20:10,000,000-10,010,000" \
  --proposed_variants_child "${VARIANTS}" \
  --proposed_variants_parent1 "${VARIANTS}" \
  --proposed_variants_parent2 "${VARIANTS}" \
  --variant_caller="vcf_candidate_importer" \
  --examples "${GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES}" \
  --channel_list='BASE_CHANNELS,insert_size' \
  --pileup_image_height_child=60 \
  --pileup_image_height_parent=40 \
  --nowrite_run_info \
  --deterministic_serialization

# ONT make examples golden test with phasing enabled
ONT_HG002_READS=${TESTDATA_DIR}/input/HG002_R10_chr20_5050000_5075000.bam
ONT_HG003_READS=${TESTDATA_DIR}/input/HG003_R10_chr20_5050000_5075000.bam
ONT_HG004_READS=${TESTDATA_DIR}/input/HG004_R10_chr20_5050000_5075000.bam
HG002_CONFIDENT_BED=${TESTDATA_DIR}/input/HG002_GRCh38_1_22_v4.2.1_benchmark.chr20.bed
HG002_CONFIDENT_VCF=${TESTDATA_DIR}/input/HG002_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz
ONT_REF=${TESTDATA_DIR}/input/grch38.chr20_5050000_5075000.masked.fa.gz
ONT_EXAMPLES_OUT=${TESTDATA_DIR}/HG002_ONT_deeptrio.examples.tfrecord.gz

time ./bazel-bin/deeptrio/make_examples \
--examples="${ONT_EXAMPLES_OUT}" \
--regions="chr20:5050000-5075000" \
--channel_list='BASE_CHANNELS,haplotype' \
--alt_aligned_pileup="diff_channels" \
--confident_regions="${HG002_CONFIDENT_BED}" \
--min_mapping_quality=1 \
--mode="training" \
--parse_sam_aux_fields=true \
--partition_size=25000 \
--phase_reads=true \
--pileup_image_height_child=100 \
--pileup_image_height_parent=100 \
--pileup_image_width=199 \
--reads="${ONT_HG002_READS}" \
--reads_parent1="${ONT_HG003_READS}"  \
--reads_parent2="${ONT_HG004_READS}" \
--realign_reads=false \
--ref="${ONT_REF}" \
--sample_name="HG002" \
--sample_name_parent1="HG003" \
--sample_name_parent2="HG004" \
--sample_name_to_train="HG002" \
--skip_parent_calling=true \
--sort_by_haplotypes=true \
--track_ref_reads=true \
--truth_variants="${HG002_CONFIDENT_VCF}" \
--vsc_min_fraction_indels=0.12 \
--vsc_min_fraction_snps=0.1 \
--nowrite_run_info \
--deterministic_serialization

# ONT make_examples with denovo regions enabled
HG002_DENOVO_REGIONS=${TESTDATA_DIR}/input/HG002_GRCh38_1_22_v4.2.1_benchmark.chr20.denovo_regions.bed
ONT_EXAMPLES_OUT=${TESTDATA_DIR}/HG002_ONT_deeptrio.denovo.examples.tfrecord.gz

time ./bazel-bin/deeptrio/make_examples \
--examples="${ONT_EXAMPLES_OUT}" \
--regions="chr20:5050000-5075000" \
--channel_list='BASE_CHANNELS,haplotype' \
--alt_aligned_pileup="diff_channels" \
--confident_regions="${HG002_CONFIDENT_BED}" \
--min_mapping_quality=1 \
--mode="training" \
--parse_sam_aux_fields=true \
--partition_size=25000 \
--phase_reads=true \
--pileup_image_height_child=100 \
--pileup_image_height_parent=100 \
--pileup_image_width=199 \
--reads="${ONT_HG002_READS}" \
--reads_parent1="${ONT_HG003_READS}"  \
--reads_parent2="${ONT_HG004_READS}" \
--realign_reads=false \
--ref="${ONT_REF}" \
--sample_name="HG002" \
--sample_name_parent1="HG003" \
--sample_name_parent2="HG004" \
--sample_name_to_train="HG002" \
--skip_parent_calling=true \
--sort_by_haplotypes=true \
--track_ref_reads=true \
--truth_variants="${HG002_CONFIDENT_VCF}" \
--denovo_regions="${HG002_DENOVO_REGIONS}" \
--vsc_min_fraction_indels=0.12 \
--vsc_min_fraction_snps=0.1 \
--write_run_info \
--deterministic_serialization

# For now, remove any golden_parent1 and golden_parent2 files because they're
# not used in testing.
rm -f "${TESTDATA_DIR}"/golden_parent*

# Remove any tmp files because they're not used.
rm -f "${TESTDATA_DIR}"/tmp.*
rm -f "${TESTDATA_DIR}"/tmp_*

# Clean up run_info to remove internal info.
find "${TESTDATA_DIR}" -name "*run_info.pbtxt" \
 -exec sed -i \
 -e 's|/home/.*/testdata|deeptrio/testdata|g' \
 -e 's/host_name: ".*"/host_name: "host"/g' {} \;

echo "================ create_golden_deeptrio.sh completed! ================"
