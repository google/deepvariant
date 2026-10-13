#!/bin/bash
# Shuffle DeepVariant examples on an EMR on EC2 cluster (YARN).
#
# 1) Creates a transient cluster that installs deps via emr_bootstrap.sh, runs
#    the shuffle as a single Spark step, then terminates itself.
# 2) Set CLUSTER_ID to submit to an existing cluster instead (it must have been
#    created with emr_bootstrap.sh as a bootstrap action).
#
# Disk, not RAM, is what a large shuffle needs: Spark spills the shuffle to the
# EBS volumes of core/task nodes. Give the cluster in total roughly 3x the
# compressed size of your input examples.
set -euo pipefail

BUCKET="${BUCKET:?set BUCKET, e.g. s3://my-bucket}"
INPUT="${INPUT:-${BUCKET}/training-case-study/output/training_set.with_label.tfrecord-?????-of-?????.gz}"
OUT_DIR="${OUT_DIR:-${BUCKET}/training-case-study/output/shuffled}"
NAME="${NAME:-training_set}"
DATASET_NAME="${DATASET_NAME:-HG001}"
REGION="${REGION:-us-east-1}"
PROFILE="${PROFILE:-}"  # optional AWS CLI profile
CLUSTER_ID="${CLUSTER_ID:-}"
# Leave empty to size shards automatically (~256 MB compressed each).
NUM_OUTPUT_SHARDS="${NUM_OUTPUT_SHARDS:-0}"

# Transient-cluster settings (ignored when CLUSTER_ID is set).
# Default roles come from `aws emr create-default-roles`.
LOG_URI="${LOG_URI:-${BUCKET}/emr-logs/}"
SUBNET_ID="${SUBNET_ID:-}"  # required unless CLUSTER_ID is set
INSTANCE_PROFILE="${INSTANCE_PROFILE:-EMR_EC2_DefaultRole}"
SERVICE_ROLE="${SERVICE_ROLE:-EMR_DefaultRole}"
WORKER_TYPE="${WORKER_TYPE:-m5.2xlarge}"
WORKER_COUNT="${WORKER_COUNT:-4}"
WORKER_EBS_GB="${WORKER_EBS_GB:-500}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHUFFLE_PY="${SCRIPT_DIR}/../shuffle_tfrecords_spark.py"
[[ -f "${SHUFFLE_PY}" ]] || SHUFFLE_PY="${SCRIPT_DIR}/shuffle_tfrecords_spark.py"
CODE_PREFIX="${BUCKET}/code/aws_shuffle"
AWS=(aws --region "${REGION}" ${PROFILE:+--profile "${PROFILE}"})

"${AWS[@]}" s3 cp "${SHUFFLE_PY}" "${CODE_PREFIX}/"
"${AWS[@]}" s3 cp "${SCRIPT_DIR}/emr_bootstrap.sh" "${CODE_PREFIX}/"

# JSON (not shorthand) so '?' / '*' / ',' in args survive intact.
STEPS_JSON="$(mktemp)"
cat > "${STEPS_JSON}" <<EOF
[{
  "Type": "Spark",
  "Name": "shuffle-${NAME}",
  "ActionOnFailure": "CONTINUE",
  "Args": [
    "--deploy-mode", "cluster",
    "--conf", "spark.python.worker.memory=1g",
    "--conf", "spark.executor.memoryOverhead=3g",
    "--conf", "spark.yarn.appMasterEnv.AWS_DEFAULT_REGION=${REGION}",
    "--conf", "spark.executorEnv.AWS_DEFAULT_REGION=${REGION}",
    "${CODE_PREFIX}/shuffle_tfrecords_spark.py",
    "--input_pattern_list=${INPUT}",
    "--output_pattern_prefix=${OUT_DIR}/${NAME}.with_label.shuffled",
    "--output_dataset_name=${DATASET_NAME}",
    "--output_dataset_config_pbtxt=${OUT_DIR}/${NAME}.dataset_config.pbtxt",
    "--num_output_shards=${NUM_OUTPUT_SHARDS}"
  ]
}]
EOF

if [[ -n "${CLUSTER_ID}" ]]; then
  "${AWS[@]}" emr add-steps --cluster-id "${CLUSTER_ID}" --steps "file://${STEPS_JSON}"
else
  : "${SUBNET_ID:?set SUBNET_ID (a subnet for the new cluster)}"
  EBS="{\"EbsBlockDeviceConfigs\":[{\"VolumeSpecification\":{\"VolumeType\":\"gp3\",\"SizeInGB\":${WORKER_EBS_GB}},\"VolumesPerInstance\":1}]}"
  "${AWS[@]}" emr create-cluster \
    --name "deepvariant-shuffle-${NAME}" \
    --release-label emr-7.5.0 \
    --applications Name=Hadoop Name=Spark \
    --log-uri "${LOG_URI}" \
    --service-role "${SERVICE_ROLE}" \
    --ec2-attributes "{\"InstanceProfile\":\"${INSTANCE_PROFILE}\",\"SubnetId\":\"${SUBNET_ID}\"}" \
    --instance-groups \
      "[{\"InstanceGroupType\":\"MASTER\",\"InstanceCount\":1,\"InstanceType\":\"m5.xlarge\"},
        {\"InstanceGroupType\":\"CORE\",\"InstanceCount\":${WORKER_COUNT},\"InstanceType\":\"${WORKER_TYPE}\",\"EbsConfiguration\":${EBS}}]" \
    --bootstrap-actions "[{\"Name\":\"shuffle-deps\",\"Path\":\"${CODE_PREFIX}/emr_bootstrap.sh\"}]" \
    --steps "file://${STEPS_JSON}" \
    --auto-terminate
fi
rm -f "${STEPS_JSON}"
