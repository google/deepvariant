#!/bin/bash
# Copyright 2023 Google LLC.
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

set -euo pipefail


USAGE=$'
Example usage:
inference_deepsomatic.sh --model_preset "WGS" --docker_build true --use_gpu true

Flags:
--docker_build (true|false)  Whether to build docker image. (default: false)
--dry_run (true|false)  If true, print out the main commands instead of running. (default: false)
--use_gpu (true|false)   Whether to use GPU when running case study. Make sure to specify vm_zone that is equipped with GPUs. (default: false)
--docker_source Where to pull the Docker image from. Default: google/deepsomatic.
--bin_version Version of DeepSomatic docker to use
--customized_model Path to checkpoint directory containing model checkpoint.
--regions Regions passed into both variant calling and som.py.
--sample_name_normal Sample name to use for normal bam.
--sample_name_tumor Sample name to use for tumor bam.
--make_examples_extra_args Flags for make_examples, specified as "flag1=param1,flag2=param2".
--call_variants_extra_args Flags for call_variants, specified as "flag1=param1,flag2=param2".
--postprocess_variants_extra_args Flags for postprocess_variants, specified as "flag1=param1,flag2=param2".
--model_preset Preset case study to run: WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES.
--pon_filtering Path to PON (panel of normal) VCF. If set, pass --pon_filtering to postprocess_variants.
--population_vcfs Path to VCFs containing population allele frequencies. Use wildcard pattern.
--proposed_variants Path to VCF containing proposed variants. In make_examples_extra_args, you must also specify variant_caller=vcf_candidate_importer but not proposed_variants.
--save_intermediate_results (true|false) If True, keep intermediate outputs from make_examples and call_variants.
--sompy_docker_source Where to pull the som.py Docker image from. Default: pkrusche/hap.py.
--skip_sompy (true|false) If True, skip the som.py evaluation.
--report_title Optional title for reports (VCF stats report and make_examples runtime report).

If model_preset is not specified, the below flags are required:
--model_type Type of DeepSomatic model to run (WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES)
--ref Path to GCP bucket containing ref file (.fa)
--bam_normal Path to GCP bucket containing BAM_NORMAL
--bam_tumor Path to GCP bucket containing BAM_TUMOR
--truth_vcf Path to GCP bucket containing truth VCF
--truth_bed Path to GCP bucket containing truth BED
--capture_bed Path to GCP bucket containing captured file (only needed for WES model_type)

Note: All paths to dataset must be of the form "gs://..."
'

# Specify default values.
# Booleans; sorted alphabetically.
BUILD_DOCKER=false
DRY_RUN=false
USE_CANDIDATE_PARTITION=false
USE_DEFAULT_PON_FILTERING=false
USE_GPU=false
SAVE_INTERMEDIATE_RESULTS=false
SKIP_SOMPY=false
# Strings; sorted alphabetically.
unset BAM_NORMAL
BAM_TUMOR=""
DOCKER_SOURCE="google/deepsomatic"
BIN_VERSION="1.6.1"
CALL_VARIANTS_ARGS=""
CAPTURE_BED=""
CUSTOMIZED_MODEL=""
MAKE_EXAMPLES_ARGS=""
MODEL_PRESET=""
MODEL_TYPE=""
PON_FILTERING=""
POPULATION_VCFS=""
POSTPROCESS_VARIANTS_ARGS=""
PROPOSED_VARIANTS=""
REF=""
REGIONS=""
REPORT_TITLE=""
SAMPLE_NAME_NORMAL=""
SAMPLE_NAME_TUMOR=""
SOMPY_DOCKER_SOURCE="pkrusche/hap.py"
TRUTH_BED=""
TRUTH_VCF=""

while (( "$#" )); do
  case "$1" in
    --docker_build)
      BUILD_DOCKER="$2"
      if [[ ${BUILD_DOCKER} != "true" ]] && [[ ${BUILD_DOCKER} != "false" ]]; then
        echo "Error: --docker_build needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --dry_run)
      DRY_RUN="$2"
      if [[ ${DRY_RUN} != "true" ]] && [[ ${DRY_RUN} != "false" ]]; then
        echo "Error: --dry_run needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --use_gpu)
      USE_GPU="$2"
      if [[ ${USE_GPU} != "true" ]] && [[ ${USE_GPU} != "false" ]]; then
        echo "Error: --use_gpu needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --save_intermediate_results)
      SAVE_INTERMEDIATE_RESULTS="$2"
      if [[ "${SAVE_INTERMEDIATE_RESULTS}" != "true" ]] && [[ "${SAVE_INTERMEDIATE_RESULTS}" != "false" ]]; then
        echo "Error: --save_intermediate_results needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --skip_sompy)
      SKIP_SOMPY="$2"
      if [[ "${SKIP_SOMPY}" != "true" ]] && [[ "${SKIP_SOMPY}" != "false" ]]; then
        echo "Error: --SKIP_SOMPY needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --use_candidate_partition)
      USE_CANDIDATE_PARTITION="$2"
      if [[ ${USE_CANDIDATE_PARTITION} != "true" ]] && [[ ${USE_CANDIDATE_PARTITION} != "false" ]]; then
        echo "Error: --use_candidate_partition needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --use_default_pon_filtering)
      USE_DEFAULT_PON_FILTERING="$2"
      if [[ ${USE_DEFAULT_PON_FILTERING} != "true" ]] && [[ ${USE_DEFAULT_PON_FILTERING} != "false" ]]; then
        echo "Error: --use_default_pon_filtering needs to have value (true|false)." >&2
        echo "$USAGE" >&2
        exit 1
      fi
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --regions)
      REGIONS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --customized_model)
      CUSTOMIZED_MODEL="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --make_examples_extra_args)
      MAKE_EXAMPLES_ARGS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --call_variants_extra_args)
      CALL_VARIANTS_ARGS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --postprocess_variants_extra_args)
      POSTPROCESS_VARIANTS_ARGS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --report_title)
      REPORT_TITLE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --model_preset)
      MODEL_PRESET="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --model_type)
      MODEL_TYPE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --docker_source)
      DOCKER_SOURCE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bin_version)
      BIN_VERSION="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --ref)
      REF="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bam_normal)
      BAM_NORMAL="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --bam_tumor)
      BAM_TUMOR="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --truth_vcf)
      TRUTH_VCF="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --truth_bed)
      TRUTH_BED="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --capture_bed)
      CAPTURE_BED="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --pon_filtering)
      PON_FILTERING="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --population_vcfs)
      POPULATION_VCFS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --proposed_variants)
      PROPOSED_VARIANTS="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --sompy_docker_source)
      SOMPY_DOCKER_SOURCE="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --sample_name_normal)
      SAMPLE_NAME_NORMAL="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --sample_name_tumor)
      SAMPLE_NAME_TUMOR="$2"
      shift # Remove argument name from processing
      shift # Remove argument value from processing
      ;;
    --help)
      echo "$USAGE" >&2
      exit 1
      ;;
    -*|--*=) # other flags not supported
      echo "Error: unrecognized flag $1" >&2
      echo "Run with --help to see usage." >&2
      exit 1
      ;;
    *)
      echo "Error: unrecognized extra args $1" >&2
      echo "Run with --help to see usage." >&2
      exit 1
      ;;
  esac
done

## Presets
# These settings specify the commonly run case studies
GCS_DATA_DIR="gs://deepvariant"  # Or you can use "https://storage.googleapis.com/deepvariant"
BASE="${HOME}/custom-case-study"

declare -a extra_args
declare -a sompy_args
declare -a docker_args

if [[ "${MODEL_PRESET}" = "WGS" ]]; then
  # Only set to default if MODEL_TYPE is not set.
  # This will allow MODEL_TYPE to be set to WGS_TUMOR_ONLY, and allow usage of
  # the tumor-only model
  if [[ -z "${MODEL_TYPE}" ]]; then
    MODEL_TYPE="WGS"
  fi
  BASE="${HOME}/deepsomatic-case-studies"

  REF="${REF:=${GCS_DATA_DIR}/deepsomatic-case-studies/GRCh38_no_alt_analysis_set.fasta}"
  # Only use the default if BAM_NORMAL is unset.
  # This will allow BAM_NORMAL to be set to an empty string, in order to enable
  # tumor-only model
  if [[ "${BAM_NORMAL+set}" != set ]]; then
    BAM_NORMAL="${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-wgs-case-study/S1395_WGS_NS_N_1.bwa.dedup.bam"
  fi
  BAM_TUMOR="${BAM_TUMOR:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-wgs-case-study/S1395_WGS_NS_T_1.bwa.dedup.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/High-Confidence_Regions_v1.2.bed}"
elif [[ "${MODEL_PRESET}" = "WES" ]]; then
  # Only set to default if MODEL_TYPE is not set.
  # This will allow MODEL_TYPE to be set to WES_TUMOR_ONLY, and allow usage of
  # the tumor-only model
  if [[ -z "${MODEL_TYPE}" ]]; then
    MODEL_TYPE="WES"
  fi
  BASE="${HOME}/deepsomatic-case-studies"

  REF="${REF:=${GCS_DATA_DIR}/deepsomatic-case-studies/GRCh38_no_alt_analysis_set.fasta}"
  # Only use the default if BAM_NORMAL is unset.
  # This will allow BAM_NORMAL to be set to an empty string, in order to enable
  # tumor-only model
  if [[ "${BAM_NORMAL+set}" != set ]]; then
    BAM_NORMAL="${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-wes-case-study/WES_IL_N_1.bwa.dedup.bam"
  fi
  BAM_TUMOR="${BAM_TUMOR:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-wes-case-study/WES_IL_T_1.bwa.dedup.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/High-Confidence_Regions_v1.2.bed}"
  CAPTURE_BED="${CAPTURE_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-wes-case-study/seqc2_hg38.exome_regions.bed}"
elif [[ "${MODEL_PRESET}" = "PACBIO" ]]; then
  # Only set to default if MODEL_TYPE is not set.
  # This will allow MODEL_TYPE to be set to WES_TUMOR_ONLY, and allow usage of
  # the tumor-only model
  if [[ -z "${MODEL_TYPE}" ]]; then
    MODEL_TYPE="PACBIO"
  fi
  BASE="${HOME}/deepsomatic-case-studies"

  REF="${REF:=${GCS_DATA_DIR}/deepsomatic-case-studies/GRCh38_no_alt_analysis_set.fasta}"
  # Only use the default if BAM_NORMAL is unset.
  # This will allow BAM_NORMAL to be set to an empty string, in order to enable
  # tumor-only model
  if [[ "${BAM_NORMAL+set}" != set ]]; then
    BAM_NORMAL="${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-pacbio-case-study/HCC1395-BL.pacbio.normal.GRCh38.bam"
  fi
  BAM_TUMOR="${BAM_TUMOR:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-pacbio-case-study/HCC1395.pacbio.tumor.GRCh38.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/High-Confidence_Regions_v1.2.bed}"
elif [[ "${MODEL_PRESET}" = "ONT" ]]; then
  # Only set to default if MODEL_TYPE is not set.
  # This will allow MODEL_TYPE to be set to WES_TUMOR_ONLY, and allow usage of
  # the tumor-only model
  if [[ -z "${MODEL_TYPE}" ]]; then
    MODEL_TYPE="ONT"
  fi
  BASE="${HOME}/deepsomatic-case-studies"

  REF="${REF:=${GCS_DATA_DIR}/deepsomatic-case-studies/GRCh38_no_alt_analysis_set.fasta}"
  # Only use the default if BAM_NORMAL is unset.
  # This will allow BAM_NORMAL to be set to an empty string, in order to enable
  # tumor-only model
  # TODO: Update to a path in gs://deepvariant.
  if [[ "${BAM_NORMAL+set}" != set ]]; then
    BAM_NORMAL="${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-ont-case-study/1395_Normal_ONT.GRCh38.sorted.bam"
  fi

  SAMPLE_NAME_NORMAL="1395_normal_ont"
  SAMPLE_NAME_TUMOR="1395_tumor_ont"
  BAM_TUMOR="${BAM_TUMOR:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-ont-case-study/1395_Tumor_ONT.50x.GRCh38.sorted.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/High-Confidence_Regions_v1.2.bed}"

elif [[ "${MODEL_PRESET}" = "FFPE_WGS" ]]; then
  MODEL_TYPE="FFPE_WGS"
  BASE="${HOME}/deepsomatic-case-studies"

  REF="${REF:=${GCS_DATA_DIR}/deepsomatic-case-studies/GRCh38_no_alt_analysis_set.fasta}"
  # Only use the default if BAM_NORMAL is unset.
  # This will allow BAM_NORMAL to be set to an empty string, in order to enable
  # tumor-only model
  # TODO: Update to a path in gs://deepvariant.
  if [[ "${BAM_NORMAL+set}" != set ]]; then
    BAM_NORMAL="${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-ffpe-wgs-case-study/FFG_IL_N_6h.bwa.dedup.bam"
  fi

  SAMPLE_NAME_NORMAL="1395_normal_ffpe_wgs"
  SAMPLE_NAME_TUMOR="1395_tumor_ffpe_wgs"
  BAM_TUMOR="${BAM_TUMOR:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-ffpe-wgs-case-study/FFG_IL_T_6h.bwa.dedup.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/High-Confidence_Regions_v1.2.bed}"

elif [[ "${MODEL_PRESET}" = "FFPE_WES" ]]; then
  MODEL_TYPE="FFPE_WES"
  BASE="${HOME}/deepsomatic-case-studies"

  REF="${REF:=${GCS_DATA_DIR}/deepsomatic-case-studies/GRCh38_no_alt_analysis_set.fasta}"
  # Only use the default if BAM_NORMAL is unset.
  # This will allow BAM_NORMAL to be set to an empty string, in order to enable
  # tumor-only model
  # TODO: Update to a path in gs://deepvariant.
  if [[ "${BAM_NORMAL+set}" != set ]]; then
    BAM_NORMAL="${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-ffpe-wes-case-study/FFX_IL_N_6h_2.bwa.dedup.bam"
  fi

  SAMPLE_NAME_NORMAL="1395_normal_ffpe_wes"
  SAMPLE_NAME_TUMOR="1395_tumor_ffpe_wes"
  CAPTURE_BED="${CAPTURE_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-wes-case-study/seqc2_hg38.exome_regions.bed}"
  BAM_TUMOR="${BAM_TUMOR:=${GCS_DATA_DIR}/deepsomatic-case-studies/deepsomatic-ffpe-wes-case-study/FFX_IL_T_6h_1.bwa.dedup.bam}"
  TRUTH_VCF="${TRUTH_VCF:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/high-confidence_sINDEL_sSNV_in_HC_regions_v1.2.1.merged.vcf.gz}"
  TRUTH_BED="${TRUTH_BED:=${GCS_DATA_DIR}/deepsomatic-case-studies/SEQC2-S1395-truth/High-Confidence_Regions_v1.2.bed}"

else
  if [[ -n "${MODEL_PRESET}" ]]; then
    echo "Error: --model_preset must be one of WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES." >&2
    exit 1
  fi
fi

# Sanity check: Tumor BAM is required.
if [[ -z "${BAM_TUMOR}" ]]; then
  echo "Error: Need to set --bam_tumor" >&2
  exit 1
fi

# Print a warning if BAM_NORMAL is unset or empty.
if [[ "${BAM_NORMAL+set}" != set ]] || [[ -z "${BAM_NORMAL}" ]]; then
  echo "--bam_normal is not specified. Please make sure you use a model that was trained with corresponding tumor-only mode."
  # Set BAM_NORMAL to empty string here, in case it's unset.
  BAM_NORMAL=""
fi

INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="deepsomatic.output.vcf.gz"
OUTPUT_GVCF="deepsomatic.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

if [[ "${MODEL_TYPE}" = "WES" || "${MODEL_TYPE}" = "FFPE_WES" ]]; then
  if [[ -n "${REGIONS}" ]]; then
    echo "Error: --regions is not used with model_type WES. Please use --capture_bed." >&2
    exit 1
  fi
  extra_args+=( --regions "/input/$(basename "$CAPTURE_BED")")
  sompy_args+=( -T "${INPUT_DIR}/$(basename $CAPTURE_BED)")
fi

if [[ "${SAVE_INTERMEDIATE_RESULTS}" == "true" ]]; then
  extra_args+=( --intermediate_results_dir "/output/intermediate_results_dir")
fi

# Because the default is False, we only set it if
# set to true.
if [[ "${USE_CANDIDATE_PARTITION}" == "true" ]]; then
  extra_args+=( --use_candidate_partition )
fi
if [[ "${USE_DEFAULT_PON_FILTERING}" == "true" ]]; then
  extra_args+=( --use_default_pon_filtering )
fi

echo "========================="
echo "# Booleans; sorted alphabetically."
echo "BUILD_DOCKER: ${BUILD_DOCKER}"
echo "DRY_RUN: ${DRY_RUN}"
echo "USE_GPU: ${USE_GPU}"
echo "SAVE_INTERMEDIATE_RESULTS: ${SAVE_INTERMEDIATE_RESULTS}"
echo "SKIP_SOMPY: ${SKIP_SOMPY}"
echo "USE_CANDIDATE_PARTITION: ${USE_CANDIDATE_PARTITION}"
echo "USE_DEFAULT_PON_FILTERING: ${USE_DEFAULT_PON_FILTERING}"
echo "# Strings; sorted alphabetically."
echo "BAM_NORMAL: ${BAM_NORMAL}"
echo "BAM_TUMOR: ${BAM_TUMOR}"
echo "DOCKER_SOURCE: ${DOCKER_SOURCE}"
echo "BIN_VERSION: ${BIN_VERSION}"
echo "CALL_VARIANTS_ARGS: ${CALL_VARIANTS_ARGS}"
echo "CAPTURE_BED: ${CAPTURE_BED}"
echo "CUSTOMIZED_MODEL: ${CUSTOMIZED_MODEL}"
echo "MAKE_EXAMPLES_ARGS: ${MAKE_EXAMPLES_ARGS}"
echo "MODEL_PRESET: ${MODEL_PRESET}"
echo "MODEL_TYPE: ${MODEL_TYPE}"
echo "PON_FILTERING: ${PON_FILTERING}"
echo "POPULATION_VCFS: ${POPULATION_VCFS}"
echo "POSTPROCESS_VARIANTS_ARGS: ${POSTPROCESS_VARIANTS_ARGS}"
echo "PROPOSED_VARIANTS: ${PROPOSED_VARIANTS}"
echo "REF: ${REF}"
echo "REGIONS: ${REGIONS}"
echo "REPORT_TITLE: ${REPORT_TITLE}"
echo "SOMPY_DOCKER_SOURCE: ${SOMPY_DOCKER_SOURCE}"
echo "TRUTH_BED: ${TRUTH_BED}"
echo "TRUTH_VCF: ${TRUTH_VCF}"
echo "========================="

function run() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    # Prints out command to stdout and a [DRY RUN] tag to stderr.
    # This allows the users to use the dry_run mode to get a list of
    # executable commands by redirecting stdout to a file.
    1>&2 printf "[DRY RUN] " && echo "$*"
  else
    echo "$*"
    eval "$@"
  fi
}

function copy_gs_or_http_file() {
  if [[ "$1" == http* ]]; then
    if curl --output /dev/null --silent --head --fail "$1"; then
      run echo "Copying from \"$1\" to \"$2\""
      run aria2c -c -x10 -s10 "$1" -d "$2"
    else
      run echo "File $1 does not exist. Skip copying."
    fi
  elif [[ "$1" == gs://* ]]; then
    status=0
    gsutil -q stat "$1" || status=1
    if [[ $status == 0 ]]; then
      run echo "Copying from \"$1\" to \"$2\""
      # Skip the file if it exists.
      run gcloud storage cp -n "$1" "$2"
    else
      run echo "File $1 does not exist. Skip copying."
    fi
  else
    echo "Unrecognized file format: $1" >&2
    exit 1
  fi
}

function copy_correct_index_file() {
  BAM="$1"
  INPUT_DIR="$2"
  # Index files have two acceptable naming patterns. We explicitly check for
  # both since we cannot use wildcard paths with http files.
  if [[ "${BAM}" == *".cram" ]]; then
    copy_gs_or_http_file "${BAM%.cram}.crai" "${INPUT_DIR}"
    copy_gs_or_http_file "${BAM}.crai" "${INPUT_DIR}"
  else
    copy_gs_or_http_file "${BAM%.bam}.bai" "${INPUT_DIR}"
    copy_gs_or_http_file "${BAM}.bai" "${INPUT_DIR}"
  fi
}

function copy_data() {
  copy_gs_or_http_file "${TRUTH_BED}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF}" "${INPUT_DIR}"
  copy_gs_or_http_file "${TRUTH_VCF}.tbi" "${INPUT_DIR}"
  # Only copy the BAM_NORMAL files if BAM_NORMAL file is not empty.
  if [[ ! -z "${BAM_NORMAL}" ]]; then
    copy_gs_or_http_file "${BAM_NORMAL}" "${INPUT_DIR}"
    copy_correct_index_file "${BAM_NORMAL}" "${INPUT_DIR}"
  fi
  copy_gs_or_http_file "${BAM_TUMOR}" "${INPUT_DIR}"
  copy_correct_index_file "${BAM_TUMOR}" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz.fai" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gz.gzi" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.gzi" "${INPUT_DIR}"
  copy_gs_or_http_file "${REF}.fai" "${INPUT_DIR}"
  if [[ "${MODEL_TYPE}" = "WES" || "${MODEL_TYPE}" = "FFPE_WES" ]]; then
    copy_gs_or_http_file "${CAPTURE_BED}" "${INPUT_DIR}"
  fi
  if [[ -n "${PROPOSED_VARIANTS}" ]]; then
    copy_gs_or_http_file "${PROPOSED_VARIANTS}" "${INPUT_DIR}"
    copy_gs_or_http_file "${PROPOSED_VARIANTS}.tbi" "${INPUT_DIR}"
  fi
  if [[ -n "${REGIONS}" ]]; then
    if [[ "${REGIONS}" = http* ]] || [[ "${REGIONS}" = gs://* ]]; then
      copy_gs_or_http_file "${REGIONS}" "${INPUT_DIR}"
    fi
  fi
  if [[ -n "${PON_FILTERING}" ]]; then
    copy_gs_or_http_file "${PON_FILTERING}" "${INPUT_DIR}"
    copy_gs_or_http_file "${PON_FILTERING}.tbi" "${INPUT_DIR}"
  fi
  if [[ -n "${POPULATION_VCFS}" ]]; then
    copy_gs_or_http_file "${POPULATION_VCFS}" "${INPUT_DIR}"
    copy_gs_or_http_file "${POPULATION_VCFS}.tbi" "${INPUT_DIR}"
  fi
}

function setup_test() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    return
  fi

  ## Create local directory structure
  mkdir -p "${OUTPUT_DIR}"
  mkdir -p "${INPUT_DIR}"
  mkdir -p "${LOG_DIR}"

  ## Download extra packages
  # Install aria2 to download data files.
  sudo apt-get -qq -y update
  sudo apt-get -qq -y install aria2

  if ! hash docker 2>/dev/null; then
    echo "'docker' was not found in PATH. Installing docker..."
    # This installs an older version, which currently can work as a workaround
    # for the issue in internal#comment23.
    sudo apt update -y && sudo apt install -y docker.io
  fi
}

function get_docker_image() {
  if [[ "${BUILD_DOCKER}" = true ]]; then
    if [[ "${USE_GPU}" = true ]]; then
      IMAGE="deepsomatic_gpu:latest"
      run "sudo docker build \
        -f Dockerfile.deepsomatic \
        --build-arg=FROM_IMAGE=nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04 \
        --build-arg=DV_GPU_BUILD=1 -t deepsomatic_gpu ."
      run echo "Done building GPU Docker image ${IMAGE}."
      docker_args+=( --gpus 1 )
    else
      IMAGE="deepsomatic:latest"
      # Building twice in case the first one times out.
      run "sudo docker build -f Dockerfile.deepsomatic -t deepsomatic . || \
        (sleep 5 ; sudo docker build -f Dockerfile.deepsomatic -t deepsomatic .)"
      run echo "Done building Docker image ${IMAGE}."
    fi

  else
    if [[ "${USE_GPU}" = true ]]; then
      if [[ "${DOCKER_SOURCE}" = "google/deepsomatic" ]]; then
        IMAGE="${DOCKER_SOURCE}:${BIN_VERSION}-gpu"
      else
        IMAGE="${DOCKER_SOURCE}:deepsomatic-${BIN_VERSION}-gpu"
      fi
      # shellcheck disable=SC2027
      # shellcheck disable=SC2086
      run "sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")"
      docker_args+=( --gpus 1 )
    else
      if [[ "${DOCKER_SOURCE}" = "google/deepsomatic" ]]; then
        IMAGE="${DOCKER_SOURCE}:${BIN_VERSION}"
      else
        IMAGE="${DOCKER_SOURCE}:deepsomatic-${BIN_VERSION}"
      fi
      # shellcheck disable=SC2027
      # shellcheck disable=SC2086
      run "sudo docker pull "${IMAGE}" || \
        (sleep 5 ; sudo docker pull "${IMAGE}")"
    fi
  fi
  if [[ "${USE_GPU}" = true ]]; then
    # shellcheck disable=SC2027
    # shellcheck disable=SC2086
    run "sudo docker run --gpus 1 "${IMAGE}" \
      python3 -c 'import tensorflow as tf; \
      print(\"is_gpu_available=\" + str(tf.test.is_gpu_available()))'"
    # shellcheck disable=SC2027
    # shellcheck disable=SC2086
    run "sudo docker run --gpus 1 "${IMAGE}" \
      python3 -c 'import tensorflow as tf; \
      tf.test.is_gpu_available() or exit(1)' \
      2> /dev/null || exit 1"
  fi
  # Pull som.py image earlier too. If there's any issues we can fail early.
  if [[ "${SKIP_SOMPY}" == "false" ]]; then
    SOMPY_VERSION="v0.3.9"
    # Pulling twice in case the first one times out.
    run "sudo docker pull ${SOMPY_DOCKER_SOURCE}:${SOMPY_VERSION} || \
      (sleep 5 ; sudo docker pull ${SOMPY_DOCKER_SOURCE}:${SOMPY_VERSION})"
  fi
}

function setup_args() {
  if [[ -n "${CUSTOMIZED_MODEL}" ]]; then
    run echo "Copy from gs:// path ${CUSTOMIZED_MODEL} to ${INPUT_DIR}/"
    # Check if it's saved Model
    saved_modelpath=${CUSTOMIZED_MODEL}/saved_model.pb
    using_saved_model=$(gsutil -q stat "$saved_modelpath" || echo 1)
    if [[ $using_saved_model != 1 ]]; then
      echo "Using saved model"
      run mkdir -p "${INPUT_DIR}/savedmodel"
      run gcloud storage cp -R "${CUSTOMIZED_MODEL}"/'*' "${INPUT_DIR}"/savedmodel/
      run gcloud storage cp "${CUSTOMIZED_MODEL}"/example_info.json "${INPUT_DIR}"/savedmodel/example_info.json
      extra_args+=( --customized_model "/input/savedmodel")
    else
      echo "Using checkpoint"
      run gcloud storage cp "${CUSTOMIZED_MODEL}".data-00000-of-00001 "${INPUT_DIR}/model.ckpt.data-00000-of-00001"
      run gcloud storage cp "${CUSTOMIZED_MODEL}".index "${INPUT_DIR}/model.ckpt.index"
      # Starting from v1.7.0, example_info.json is required.
      CUSTOMIZED_MODEL_DIR="$(dirname "${CUSTOMIZED_MODEL}")"
      run "gcloud storage cp ${CUSTOMIZED_MODEL_DIR}/example_info.json ${INPUT_DIR}/example_info.json"
      extra_args+=( --customized_model "/input/model.ckpt")
    fi
  else
    run echo "No custom model specified."
  fi
  if [[ -n "${POPULATION_VCFS}" ]]; then
    MAKE_EXAMPLES_ARGS="${MAKE_EXAMPLES_ARGS:+${MAKE_EXAMPLES_ARGS},}population_vcfs=/input/$(basename "$POPULATION_VCFS")"
  fi
  if [[ -n "${MAKE_EXAMPLES_ARGS}" ]]; then
    # In order to use proposed variants, we have to pass vcf_candidate_importer
    # to make_examples_extra_args, so we know that we will enter this if
    # statement.
    if [[ -n "${PROPOSED_VARIANTS}" ]]; then
      MAKE_EXAMPLES_ARGS="${MAKE_EXAMPLES_ARGS},proposed_variants=/input/$(basename "$PROPOSED_VARIANTS")"
    fi
    extra_args+=( --make_examples_extra_args "\"${MAKE_EXAMPLES_ARGS}\"")
  fi
  if [[ -n "${CALL_VARIANTS_ARGS}" ]]; then
    extra_args+=( --call_variants_extra_args "\"${CALL_VARIANTS_ARGS}\"")
  fi
  if [[ -n "${POSTPROCESS_VARIANTS_ARGS}" ]]; then
    extra_args+=( --postprocess_variants_extra_args "\"${POSTPROCESS_VARIANTS_ARGS}\"")
  fi
  if [[ -n "${REPORT_TITLE}" ]]; then
    extra_args+=( --report_title "${REPORT_TITLE}")
  fi
  if [[ -n "${REGIONS}" ]]; then
    if [[ "${REGIONS}" = http* ]] || [[ "${REGIONS}" = gs://* ]]; then
      extra_args+=( --regions "/input/$(basename "$REGIONS")")
      sompy_args+=( -T "${INPUT_DIR}/$(basename "$REGIONS")")
    else
      extra_args+=( --regions "${REGIONS}")
      sompy_args+=( -l "${REGIONS}")
    fi
  fi
  # If you're running an older version (before 1.2) that doesn't have this flag,
  # you'll need to comment out this line.
  extra_args+=( --runtime_report )
}

function run_deepsomatic_with_docker() {
  run echo "Run DeepSomatic..."
  run echo "using IMAGE=${IMAGE}"
  if [[ ! -z "${BAM_NORMAL}" ]]; then
    extra_args+=( --reads_normal "/input/$(basename "$BAM_NORMAL")" )
  fi
  if [[ ! -z "${PON_FILTERING}" ]]; then
    extra_args+=( --pon_filtering "/input/$(basename "$PON_FILTERING")" )
  fi
  if [[ ! -z "${SAMPLE_NAME_NORMAL}" ]]; then
    extra_args+=( --sample_name_normal "${SAMPLE_NAME_NORMAL}" )
  fi
  if [[ ! -z "${SAMPLE_NAME_TUMOR}" ]]; then
    extra_args+=( --sample_name_tumor "${SAMPLE_NAME_TUMOR}" )
  fi
  # shellcheck disable=SC2027
  # shellcheck disable=SC2046
  # shellcheck disable=SC2068
  # shellcheck disable=SC2086
  # shellcheck disable=SC2145
  run "(time (sudo docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    ${docker_args[@]-} \
    "${IMAGE}" \
    run_deepsomatic \
    --model_type="${MODEL_TYPE}" \
    --ref="/input/$(basename $REF).gz" \
    --reads_tumor="/input/$(basename $BAM_TUMOR)" \
    --output_vcf="/output/${OUTPUT_VCF}" \
    --output_gvcf="/output/${OUTPUT_GVCF}" \
    --num_shards "$(nproc)" \
    --logging_dir="/output/logs" \
    "${extra_args[@]-}" && \
  echo "Done.")) 2>&1 | tee "${LOG_DIR}/deepsomatic_runtime.log""
  echo
}

function run_sompy() {
  ## Evaluation: run som.py
  run echo "Start evaluation with som.py..."
  UNCOMPRESSED_REF="${INPUT_DIR}/$(basename $REF)"
  # Uncompress fa into  a writable directory. Index file was downloaded earlier.
  # shellcheck disable=SC2027
  # shellcheck disable=SC2046
  # shellcheck disable=SC2086
  run "zcat <"${INPUT_DIR}/$(basename $REF).gz" >"${UNCOMPRESSED_REF}""

  SOMPY_VERSION="v0.3.9"
  # shellcheck disable=SC2027
  # shellcheck disable=SC2046
  # shellcheck disable=SC2086
  # shellcheck disable=SC2145
  run "( sudo docker run -i \
  -v "${INPUT_DIR}:${INPUT_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  ${SOMPY_DOCKER_SOURCE}:${SOMPY_VERSION} /opt/hap.py/bin/som.py \
    "${INPUT_DIR}/$(basename $TRUTH_VCF)" \
    "${OUTPUT_DIR}/$1" \
    --restrict-regions "${INPUT_DIR}/$(basename $TRUTH_BED)" \
    -r "${UNCOMPRESSED_REF}" \
    -o "${OUTPUT_DIR}/$2.output" \
    --feature-table generic \
    ${sompy_args[@]-} \
  ) 2>&1 | tee "${LOG_DIR}/$2.log""
  if [[ "${DRY_RUN}" != "true" ]]; then
    echo "${SOMPY_VERSION}" > "${LOG_DIR}/sompy_version.log"
  fi
  run echo "Done."
}

function main() {
  run echo 'Starting the test...'

  setup_test
  # Get or build Docker image first, before downloading data.
  # They're independent steps, but might be nice to know if the Docker soource
  # doesn't exist.
  if [[ ${DOCKER_SOURCE} =~ ^gcr.io ]] || [[ ${SOMPY_DOCKER_SOURCE} =~ ^gcr.io ]]; then
    run "gcloud auth print-access-token | sudo docker login -u oauth2accesstoken --password-stdin https://gcr.io"
  fi
  get_docker_image
  copy_data
  setup_args
  run_deepsomatic_with_docker
  if [[ "${SKIP_SOMPY}" == "false" ]]; then
    if [[ "${DRY_RUN}" == "true" ]]; then
      run_sompy "${OUTPUT_VCF}" "sompy"
    else
      run_sompy "${OUTPUT_VCF}" "sompy" 2>&1 | tee "${LOG_DIR}/sompy.log"
    fi
  fi
}

main "$@"
