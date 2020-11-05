#!/bin/bash

set -e

if [ $# -eq 0 ]; then
	echo "Usage: $0 <Pattern> <Working directory> <Number of buckets> <Number of threads> [Mount Pattern]"
	exit
fi

IMAGE_NAME=dvorak.maas/oddjobs/flink-base
PATTERN=$1
WORKDIR=$2
NUM_BUCKETS=$3
NUM_THREADS=$4
MOUNT_PATTERN="$HOME/DockerStorage:/root/storage"

if [ $# -eq 5 ]; then
	MOUNT_PATTERN=$5
fi

get_dockerized_path() {
	path=$1
	source_dir=$(echo $MOUNT_PATTERN | cut -f1 -d ":")
	dest_dir=$(echo $MOUNT_PATTERN | cut -f2 -d ":")
	echo $(echo "$path" | sed "s?^$source_dir?$dest_dir?g")
}

if [ -d $WORKDIR ]; then
	echo "$WORKDIR exists. Provide a directory path that doesn't exist."
	exit
fi

mkdir -p $WORKDIR

# First create partitions
echo "Performing partitioning"
input_pattern_list=$(get_dockerized_path "$PATTERN")
workdir=$(get_dockerized_path "$WORKDIR/partitions")

docker run -it --rm -v $MOUNT_PATTERN --entrypoint python $IMAGE_NAME \
	/opt/tools/partition_tfrecords_beam.py \
		--input_pattern_list="$input_pattern_list" \
		--workdir="$workdir" \
		--num_buckets=$NUM_BUCKETS \
		--direct_num_workers=$NUM_THREADS \
		--direct_running_mode="multi_processing" \
		--streaming

# Shuffle each partition
echo "Performing shuffling"
mkdir -p $WORKDIR/shuffled

for i in `seq 0 $((NUM_BUCKETS - 1))`; do
	partition_prefix=$(printf "%05d" "$i")
	input_pattern_list=$(get_dockerized_path "$WORKDIR/partitions/partition${partition_prefix}-*")
	output_pattern_name=$(printf "$WORKDIR/shuffled/shuffled_partition-%05d-" "$i")
	output_pattern_prefix=$(get_dockerized_path "$output_pattern_name")
	output_dataset_name="partition_shuffle${i}"
	echo "Shuffling partition $i"
	docker run --rm -it -v $MOUNT_PATTERN --entrypoint python $IMAGE_NAME \
		/opt/tools/shuffle_tfrecords_beam.py \
			--input_pattern_list="$input_pattern_list" \
			--output_pattern_prefix="$output_pattern_prefix" \
			--output_dataset_name="$output_dataset_name" \
			--direct_num_workers=$NUM_THREADS \
			--direct_running_mode"multi_processing" \
			--streaming
done

# Combine all partitions into the output files
echo "Combining outputs"
mkdir -p $WORKDIR/results
input_pattern_list=$(get_dockerized_path "$WORKDIR/shuffled/shuffled_partition-?????-*")
output_pattern_prefix=$(get_dockerized_path "$WORKDIR/results/combined")
output_dataset_config_pbtxt=$(get_dockerized_path "$WORKDIR/results/combined_pbtxt.txt")

docker run --rm -t -v $MOUNT_PATTERN --entrypoint python $IMAGE_NAME \
	/opt/tools/combine_tfrecords_beam.py \
		--input_pattern_list="$input_pattern_list" \
		--output_pattern_prefix="$output_pattern_prefix" \
		--output_dataset_config_pbtxt="$output_dataset_config_pbtxt" \
		--output_dataset_name="results" \
		--direct_num_workers=$NUM_THREADS \
		--direct_running_mode="multi_processing" \
		--streaming
