# Shuffling DeepVariant training examples without Google Cloud

The existing `tools/shuffle_tfrecords_beam.py` runs on Google Cloud Dataflow or
on Beam's DirectRunner. Without Google Cloud, that leaves DirectRunner, which
holds the whole shuffle in memory: a large training set needs a machine with
very large RAM, or fails with out-of-memory errors.

`tools/shuffle_tfrecords_spark.py` does the same SHA1-keyed shuffle with plain
PySpark and produces the same outputs (sharded `*.tfrecord.gz`,
`dataset_config.pbtxt`, `example_info.json`). Spark spills the shuffle to disk,
so memory stays roughly constant no matter how big the dataset is. **Plan for
disk, not RAM:** about 1.5x the gzipped input size (see measurements below).

There are two ways to run it:

- **One machine** (recommended). Simplest; a large dataset just takes longer.
- **An AWS EMR cluster**, if you don't want to wait on a single instance.

## Option 1: one machine

Any Linux or macOS machine with Java 8, 11 or 17 works, for example a single
EC2 instance. Give it a data volume large enough for the spill (about 1.5x the
gzipped input), and point Spark at that volume with `spark.local.dir`. By
default Spark spills to `/tmp`, which on EC2 is the small root volume.

```bash
pip install pyspark google-crc32c psutil boto3
spark-submit --master 'local[*]' --conf spark.local.dir=/mnt/data/spark-tmp \
  tools/shuffle_tfrecords_spark.py \
  --input_pattern_list='s3://my-bucket/out/training_set.with_label.tfrecord-?????-of-?????.gz' \
  --output_pattern_prefix='s3://my-bucket/out/training_set.with_label.shuffled' \
  --output_dataset_name=HG001 \
  --output_dataset_config_pbtxt='s3://my-bucket/out/training_set.dataset_config.pbtxt'
```

Inputs and outputs can be `s3://` paths or local paths. If the output goes to
the same disk as the spill, add another 1x the input size of disk space.

## Option 2: AWS EMR cluster (faster)

A cluster spreads the same job across several machines, so it finishes sooner.
It also retries failed tasks on other nodes, which helps on long runs.
`run_emr_ec2.sh` creates a temporary cluster, runs the shuffle as one Spark
step, then terminates the cluster:

```bash
BUCKET=s3://my-bucket \
INPUT='s3://my-bucket/out/training_set.with_label.tfrecord-?????-of-?????.gz' \
OUT_DIR=s3://my-bucket/out/shuffled NAME=training_set \
SUBNET_ID=subnet-xxxxxxxx WORKER_COUNT=4 WORKER_EBS_GB=500 \
  tools/aws_shuffle/run_emr_ec2.sh
```

The cluster installs the dependencies (`google-crc32c psutil boto3`) through
`emr_bootstrap.sh`. It uses the standard EMR roles (create them once with
`aws emr create-default-roles`). To use an existing cluster that was created
with that bootstrap action, set `CLUSTER_ID=j-...`. On EMR, Spark already
spills to the EBS volumes under `/mnt`, so `WORKER_COUNT x WORKER_EBS_GB` is
the total spill space.

## Resource needs (measured)

Test: 80,000 DeepVariant-sized examples, 3.0 GB gzipped (~12 GB raw),
`local[4]` with `--driver-memory 1g` and `spark.python.worker.memory=256m`.

| Peak | Value |
|---|---|
| Total memory, all Spark processes | 2.1 GB |
| Largest Python worker | 0.46 GB |
| Spill disk | 4.4 GB (~1.5x gzipped input) |
| Open files per Python worker | 37 |

Memory does not grow with dataset size. Disk does: plan for ~1.5x the gzipped
input for spill, plus 1x more if the output is written to the same disk.

## Tuning

- `--num_output_shards` (default: about 256 MB of compressed input per shard).
  Each shard is one Spark partition, so you need enough shards to keep every
  core busy.
- `spark.python.worker.memory` (Spark's default is 512m; the EMR script sets
  1g) is the threshold at which the Python-side sort spills to disk.
- Input records are CRC-checked by default. Use `--no_verify_crc` to skip that.

## Why not run the Beam script on Spark?

Beam's Python `SparkRunner` is a *portable* runner: the Python process starts a
Java job server (downloaded from Maven at runtime), and every Spark executor
launches a Beam SDK harness, by default in a Docker container. On EMR that
breaks in several ways:

- `spark-submit my_beam_pipeline.py` does not turn a Beam pipeline into a Spark
  job; the script just runs as a PySpark driver that then tries to start a
  second (Beam) runner.
- The Python `SparkRunner` cannot target `--spark_master=yarn` directly. It only
  supports `local[*]` and `spark://host:port`, or building a jar with
  `--output_executable_path` and spark-submitting that jar.
- Executors need Docker (or a hand-installed SDK `boot` binary with
  `--environment_type=PROCESS`), plus internet access for pip/Maven.

The PySpark script avoids all of that.
