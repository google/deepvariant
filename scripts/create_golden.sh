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
# ./scripts/create_golden.sh
# to generate data under deepvariant/testdata.
#
# This script is used to generate the files in testdata/ which are useful for
# the e2e-stype unit tests.
#
# The reason why we don't run it everytime is because the tfrecord.gz files
# won't be the same everytime we generate it (even without code changes).
# Which is why this script is usually only run when we know our commits would
# change these files.
#
# Recommendation: run this script with every commit and make sure the
# behavior is as expected. Can revert the binary files if you know there is not
# actual content change.
#
# Before running this script, set up as usual:
#   sudo su; ./build-preqreq.sh
#
# NOTE: Please run until you see this line:
#
# ================ create_golden.sh completed! ================
#
# Otherwise, some commands might have failed in between.

set -euo pipefail
set -x

TESTDATA_DIR=$(pwd)/deepvariant/testdata

# First, clean up all files on the top directory.

find "${TESTDATA_DIR}" -maxdepth 1 -type f -delete

# Inputs
REF=${TESTDATA_DIR}/input/ucsc.hg19.chr20.unittest.fasta.gz
READS=${TESTDATA_DIR}/input/NA12878_S1.chr20.10_10p1mb.bam
# CRAM_READS=${TESTDATA_DIR}/NA12878_S1.chr20.10_10p1mb.cram
# Here is how I prepared the CRAM file:
# samtools view -Ch -T /brain-genomics/biotf/pfda/ucsc_hg19.fa ${READS} > ${CRAM_READS}
# samtools index ${CRAM_READS}
VARIANTS=${TESTDATA_DIR}/input/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz
# This is a simplified version, manually modified from
# test_nist.b37_chr20_100kbp_at_10mb.vcf.gz. "type=1" and "type=2" are manually
# added to the INFO field.
WITH_TYPE_VARIANTS=${TESTDATA_DIR}/input/with_types.test_nist.b37_chr20_4kbp_at_10mb.vcf.gz
REGIONS=${TESTDATA_DIR}/input/test_nist.b37_chr20_100kbp_at_10mb.bed

# Files below were created using commands in internal#comment2.
GRCH38_FASTA=${TESTDATA_DIR}/input/grch38.chr20_and_21_10M.fa.gz
AF_VCF_CHR20=${TESTDATA_DIR}/input/cohort-chr20_100k.vcf.gz
AF_VCF_CHR21=${TESTDATA_DIR}/input/cohort-chr21_100k.vcf.gz
AF_VCF_CHR20_AND_21=${TESTDATA_DIR}/input/cohort-chr20_and_chr21_100k.vcf.gz
GRCH38_CHR20_AND_21_BAM=${TESTDATA_DIR}/input/grch38_1k_subset_chr20_and_chr21.bam


# Outputs
GOLDEN_TRAINING_EXAMPLES=${TESTDATA_DIR}/golden.training_examples.tfrecord.gz
GOLDEN_TRAINING_EXAMPLES_VCF=${TESTDATA_DIR}/golden.training_examples.vcf
GOLDEN_CALLING_CANDIDATES=${TESTDATA_DIR}/golden.calling_candidates.tfrecord.gz
GOLDEN_CALLING_EXAMPLES=${TESTDATA_DIR}/golden.calling_examples.tfrecord.gz
GOLDEN_CALLING_EMPTY_EXAMPLES=${TESTDATA_DIR}/golden.calling_examples_empty.tfrecord.gz
CALLING_EXAMPLES_TEMP=${TESTDATA_DIR}/tmp.examples.tfrecord.gz
GOLDEN_CANDIDATE_POSITIONS="${TESTDATA_DIR}"/golden.candidate_positions
GOLDEN_POSTPROCESS_INPUT=${TESTDATA_DIR}/golden.postprocess_single_site_input.tfrecord.gz
GOLDEN_POSTPROCESS_OUTPUT=${TESTDATA_DIR}/golden.postprocess_single_site_output.vcf
GOLDEN_POSTPROCESS_OUTPUT_PASS_ONLY=${TESTDATA_DIR}/golden.postprocess_single_site_output.pass_only.vcf
GOLDEN_POSTPROCESS_GVCF_INPUT=${TESTDATA_DIR}/golden.postprocess_gvcf_input.tfrecord.gz
GOLDEN_POSTPROCESS_GVCF_OUTPUT=${TESTDATA_DIR}/golden.postprocess_gvcf_output.g.vcf
GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES=${TESTDATA_DIR}/golden.vcf_candidate_importer.training_examples.tfrecord.gz

MODEL=gs://deepvariant/models/DeepVariant/1.6.0/savedmodels/deepvariant.wgs.savedmodel
# Speed up by copying to /tmp/
rm -rf /tmp/deepvariant.wgs.savedmodel
gsutil -m cp -R ${MODEL} /tmp/
MODEL=/tmp/deepvariant.wgs.savedmodel

source settings.sh

# Need to run this first, otherwise the bazel build command below will complain
# about deepvariant/examples_from_stream.so.
./build_release_binaries.sh

# shellcheck disable=SC2086
bazel build -c opt ${DV_COPT_FLAGS} //deepvariant:make_examples \
  //deepvariant:call_variants \
  //deepvariant:postprocess_variants \
  //deepvariant/labeler:labeled_examples_to_vcf

# Makes the training output.
time ./bazel-bin/deepvariant/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${GOLDEN_TRAINING_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --write_run_info \
  --deterministic_serialization

# Currently shuffle_tfrecords_test.cc uses nucleus::ReadProtosFromTFRecord,
# which doesn't read from gzipped tfrecord. So we unzip this file for now.
# Remove this line once we can read gzipped tfrecord in Nucleus.
gunzip -k -f "${GOLDEN_TRAINING_EXAMPLES}"

# Makes the training output; 3 shards.
for i in $(seq 0 2);
do
time ./bazel-bin/deepvariant/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${GOLDEN_TRAINING_EXAMPLES}"@3 \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --task "${i}" \
  --deterministic_serialization
done

# Now make the calling output.
# NOTE:
# --proposed_variants is not used unless --variant_caller=vcf_candidate_importer
# is set.
time ./bazel-bin/deepvariant/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --proposed_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --candidates "${GOLDEN_CALLING_CANDIDATES}" \
  --examples "${GOLDEN_CALLING_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --gvcf "${GOLDEN_POSTPROCESS_GVCF_INPUT}" \
  --deterministic_serialization

# Generate make_examples output for candidate_sweep mode
time ./bazel-bin/deepvariant/make_examples \
  --mode candidate_sweep \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --confident_regions "${REGIONS}" \
  --examples "${CALLING_EXAMPLES_TEMP}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --candidate_positions "${GOLDEN_CANDIDATE_POSITIONS}" \
  --deterministic_serialization

# Now make the calling output; 3 shards.
for i in $(seq 0 2);
do
time ./bazel-bin/deepvariant/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --proposed_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${GOLDEN_CALLING_EXAMPLES}"@3 \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --gvcf "${GOLDEN_POSTPROCESS_GVCF_INPUT}"@3 \
  --task "${i}" \
  --deterministic_serialization

# Generate make_examples output for candidate_sweep mode
time ./bazel-bin/deepvariant/make_examples \
  --mode candidate_sweep \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --confident_regions "${REGIONS}" \
  --examples "${CALLING_EXAMPLES_TEMP}"@3 \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref' \
  --candidate_positions "${GOLDEN_CANDIDATE_POSITIONS}"@3 \
  --task "${i}" \
  --deterministic_serialization
done

# Create empty examples file.
touch "${GOLDEN_CALLING_EMPTY_EXAMPLES}"

# This calls `call_variants` and generate the input for `postprocess_variants`.
time ./bazel-bin/deepvariant/call_variants \
  --outfile "${GOLDEN_POSTPROCESS_INPUT}" \
  --examples "${GOLDEN_CALLING_EXAMPLES}" \
  --checkpoint "${MODEL}"

# Finally, we run `postprocess_variants`.
time ./bazel-bin/deepvariant/postprocess_variants \
  --infile "${GOLDEN_POSTPROCESS_INPUT}" \
  --outfile "${GOLDEN_POSTPROCESS_OUTPUT}.gz" \
  --ref "${REF}" \
  --nonvariant_site_tfrecord_path "${GOLDEN_POSTPROCESS_GVCF_INPUT}" \
  --gvcf_outfile "${GOLDEN_POSTPROCESS_GVCF_OUTPUT}" \
  --novcf_stats_report \
  --cpus 0

# Make a version where we make all chr20 a haploid contig.
GOLDEN_POSTPROCESS_OUTPUT_HAPLOID=${TESTDATA_DIR}/golden.haploid_chr20.postprocess_single_site_output.vcf
GOLDEN_POSTPROCESS_GVCF_OUTPUT_HAPLOID=${TESTDATA_DIR}/golden.haploid_chr20.postprocess_gvcf_output.g.vcf

time ./bazel-bin/deepvariant/postprocess_variants \
  --infile "${GOLDEN_POSTPROCESS_INPUT}" \
  --outfile "${GOLDEN_POSTPROCESS_OUTPUT_HAPLOID}" \
  --ref "${REF}" \
  --nonvariant_site_tfrecord_path "${GOLDEN_POSTPROCESS_GVCF_INPUT}" \
  --gvcf_outfile "${GOLDEN_POSTPROCESS_GVCF_OUTPUT_HAPLOID}" \
  --haploid_contigs chr20 \
  --novcf_stats_report \
  --cpus 0

# The VCF version is also used in our unit test.
zcat "${GOLDEN_POSTPROCESS_OUTPUT}.gz" > "${GOLDEN_POSTPROCESS_OUTPUT}"

# Also run `postprocess_variants` with --only_keep_pass
time ./bazel-bin/deepvariant/postprocess_variants \
  --infile "${GOLDEN_POSTPROCESS_INPUT}" \
  --outfile "${GOLDEN_POSTPROCESS_OUTPUT_PASS_ONLY}" \
  --ref "${REF}" \
  --novcf_stats_report \
  --only_keep_pass \
  --cpus 0

time ./bazel-bin/deepvariant/labeler/labeled_examples_to_vcf \
  --ref "${REF}" \
  --examples "${GOLDEN_TRAINING_EXAMPLES}" \
  --output_vcf "${GOLDEN_TRAINING_EXAMPLES_VCF}"


####### Data for CustomizedClassesVariantLabler #######
CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES=${TESTDATA_DIR}/customized_classes.golden.training_examples.tfrecord.gz
time ./bazel-bin/deepvariant/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions "chr20:10,000,000-10,004,000" \
  --truth_variants "${WITH_TYPE_VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --labeler_algorithm "customized_classes_labeler" \
  --customized_classes_labeler_classes_list "ref,class1,class2" \
  --customized_classes_labeler_info_field_name "type" \
  --deterministic_serialization

####### Data for VcfCaller #######
# Input:
VCF_CANDIDATE_IMPORTER_VARIANTS=${TESTDATA_DIR}/input/vcf_candidate_importer.indels.chr20.vcf.gz

GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES=${TESTDATA_DIR}/golden.vcf_candidate_importer_calling_examples.tfrecord.gz
GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT=${TESTDATA_DIR}/golden.vcf_candidate_importer_postprocess_single_site_input.tfrecord.gz
GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_OUTPUT=${TESTDATA_DIR}/golden.vcf_candidate_importer_postprocess_single_site_output.vcf

time ./bazel-bin/deepvariant/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --variant_caller "vcf_candidate_importer" \
  --truth_variants "${VARIANTS}" \
  --examples "${GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --deterministic_serialization

time ./bazel-bin/deepvariant/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --variant_caller "vcf_candidate_importer" \
  --proposed_variants "${VCF_CANDIDATE_IMPORTER_VARIANTS}" \
  --examples "${GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --regions "chr20:59,777,000-60,000,000" \
  --norealign_reads \
  --deterministic_serialization

# This calls `call_variants` and generate the input for `postprocess_variants`.
time ./bazel-bin/deepvariant/call_variants \
  --outfile "${GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT}" \
  --examples "${GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES}" \
  --checkpoint "${MODEL}"

# Finally, we run `postprocess_variants`.
time ./bazel-bin/deepvariant/postprocess_variants \
  --infile "${GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT}" \
  --outfile "${GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_OUTPUT}" \
  --nogroup_variants \
  --ref "${REF}" \
  --novcf_stats_report \
  --cpus 0

# Run with --alt_aligned_pileup "rows"
ALT_ALIGNED_ROWS_EXAMPLES=${TESTDATA_DIR}/golden.alt_aligned_pileup_rows_examples.tfrecord.gz
time ./bazel-bin/deepvariant/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions "chr20:10,000,000-10,010,000" \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${ALT_ALIGNED_ROWS_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref' \
  --alt_aligned_pileup "rows" \
  --deterministic_serialization

# Run with --alt_aligned_pileup "diff_channels"
ALT_ALIGNED_DIFF_CHANNELS_EXAMPLES=${TESTDATA_DIR}/golden.alt_aligned_pileup_diff_channels_examples.tfrecord.gz
time ./bazel-bin/deepvariant/make_examples \
  --mode training \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions "chr20:10,000,000-10,010,000" \
  --truth_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --examples "${ALT_ALIGNED_DIFF_CHANNELS_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref' \
  --alt_aligned_pileup "diff_channels" \
  --deterministic_serialization

# Run with allele frequency.
# Regions supported by these files: chr20:61001-62000,chr21:5114000-5114999.
GOLDEN_ALLELE_FREQUENCY_EXAMPLES=${TESTDATA_DIR}/golden.allele_frequency_examples.tfrecord.gz
time ./bazel-bin/deepvariant/make_examples \
  --mode "calling" \
  --ref "${GRCH38_FASTA}" \
  --regions 'chr20:61001-62000' \
  --population_vcfs "${AF_VCF_CHR20_AND_21}" \
  --reads "${GRCH38_CHR20_AND_21_BAM}" \
  --examples "${GOLDEN_ALLELE_FREQUENCY_EXAMPLES}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size,allele_frequency' \
  --deterministic_serialization

## Testing the 3 steps with include_med_dp:
GOLDEN_CALLING_CANDIDATES_WITH_MED_DP=${TESTDATA_DIR}/DELETME.golden.calling_candidates.med_dp.tfrecord.gz
GOLDEN_CALLING_EXAMPLES_WITH_MED_DP=${TESTDATA_DIR}/DELETME.golden.calling_examples.med_dp.tfrecord.gz
GOLDEN_POSTPROCESS_GVCF_INPUT_WITH_MED_DP=${TESTDATA_DIR}/DELETME.golden.postprocess_gvcf_input.med_dp.tfrecord.gz
GOLDEN_POSTPROCESS_INPUT_WITH_MED_DP=${TESTDATA_DIR}/DELETME.golden.postprocess_single_site_input.med_dp.tfrecord.gz
GOLDEN_POSTPROCESS_OUTPUT_WITH_MED_DP=${TESTDATA_DIR}/DELETME.golden.postprocess_single_site_output.med_dp.vcf
# This file is created just for the sake of diff during CL review. It is not
# used in any unit tests.
GOLDEN_POSTPROCESS_GVCF_OUTPUT_WITH_MED_DP=${TESTDATA_DIR}/golden.postprocess_gvcf_output.med_dp.g.vcf

time ./bazel-bin/deepvariant/make_examples \
  --mode calling \
  --ref "${REF}" \
  --reads "${READS}" \
  --regions chr20:10,000,000-10,010,000 \
  --proposed_variants "${VARIANTS}" \
  --confident_regions "${REGIONS}" \
  --candidates "${GOLDEN_CALLING_CANDIDATES_WITH_MED_DP}" \
  --examples "${GOLDEN_CALLING_EXAMPLES_WITH_MED_DP}" \
  --channel_list='read_base,base_quality,mapping_quality,strand,read_supports_variant,base_differs_from_ref,insert_size' \
  --gvcf "${GOLDEN_POSTPROCESS_GVCF_INPUT_WITH_MED_DP}" \
  --include_med_dp \
  --deterministic_serialization

time ./bazel-bin/deepvariant/call_variants \
  --outfile "${GOLDEN_POSTPROCESS_INPUT_WITH_MED_DP}" \
  --examples "${GOLDEN_CALLING_EXAMPLES_WITH_MED_DP}" \
  --checkpoint "${MODEL}"

time ./bazel-bin/deepvariant/postprocess_variants \
  --infile "${GOLDEN_POSTPROCESS_INPUT_WITH_MED_DP}" \
  --outfile "${GOLDEN_POSTPROCESS_OUTPUT_WITH_MED_DP}" \
  --ref "${REF}" \
  --nonvariant_site_tfrecord_path "${GOLDEN_POSTPROCESS_GVCF_INPUT_WITH_MED_DP}" \
  --gvcf_outfile "${GOLDEN_POSTPROCESS_GVCF_OUTPUT_WITH_MED_DP}" \
  --novcf_stats_report \
  --cpus 0

# Remove anything that started with DELETEME.
rm -f "${TESTDATA_DIR}"/DELETME.*

# Remove any tmp.* files because they're not used.
rm -f "${TESTDATA_DIR}"/tmp.*

# Clean up run_info to remove internal info.
find "${TESTDATA_DIR}" -name "*run_info.pbtxt" \
 -exec sed -i \
 -e 's|/home/.*/testdata|deepvariant/testdata|g' \
 -e 's/host_name: ".*"/host_name: "host"/g' {} \;
echo "================ create_golden.sh completed! ================"
