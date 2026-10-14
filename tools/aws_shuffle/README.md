# Shuffling DeepVariant training examples without Google Cloud

DeepVariant training needs the examples shuffled globally before training. The
upstream script, `tools/shuffle_tfrecords_beam.py`, supports two ways to run it:

- **Google Cloud Dataflow.** Scales to large datasets, but needs a GCP account.
- **Beam DirectRunner, on one machine.** Holds the whole shuffle in memory, so
  memory grows with the dataset. In our test below, the 1.97 GB (gzipped)
  training set needed 32.6 GB of RAM. Larger training sets need a very
  large-memory machine or fail with out-of-memory errors.

`tools/shuffle_tfrecords_spark.py` does the same SHA1-keyed shuffle with plain
PySpark and writes the same outputs (sharded `*.tfrecord.gz`,
`dataset_config.pbtxt`, `example_info.json`). Spark spills the shuffle to disk,
so memory depends on how many cores you use, not on dataset size. It runs on a
single machine or on an AWS EMR cluster. **Plan for disk, not RAM.**

## Results on real data

We shuffled the training case study's examples (HG001, made with
`make_examples` from DeepVariant 1.6.1) with both scripts on the same 64-core,
187 GB Linux server:

- Validation set: 16 files, 0.34 GB gzipped, 59,401 examples.
- Training set: 7 files, 1.97 GB gzipped, 342,758 examples.

Beam used the v1.10 script with the pinned `apache_beam==2.71.0`. Spark used
PySpark 3.5.3 on Java 17. Peak memory is the total over all processes of the
job.

**Training set:**

| | Beam DirectRunner | Spark `local[4]` | Spark `local[*]` | Spark `local[*]`, 64 shards |
|---|---|---|---|---|
| Wall time | 9:52 | 13:09 | 6:33 | **1:45** |
| Peak memory | 32.6 GB | **3.9 GB** | 6.1 GB | 38.9 GB |
| Average cores busy | ~4 | ~4 | ~8 | ~41 |
| Output shards | 5 | 8 | 8 | 64 |

**Validation set:**

| | Beam DirectRunner | Spark `local[*]` |
|---|---|---|
| Wall time | 2:26 | 4:08 |
| Peak memory | 5.0-6.0 GB | 3.1 GB |
| Output shards | 2 | 2 |

Spark at its default settings is slower on small inputs: the default gives 2
shards for validation and 8 for training, and each shard is sorted and
gzipped by one core (see [Tuning](#tuning)). With `--num_output_shards=64`
the training shuffle took 1:45, with reading done in about 30 seconds.

### Randomness checks

We checked every output, Beam and Spark, against its input:

1. **Same records.** The output holds exactly the input records (compared by
   SHA1), and the label counts and `num_examples` in `dataset_config.pbtxt`
   match the output.
2. **Order check.** For the first records of each output shard, we printed the
   input file and position each record came from. They jump between input
   files and positions, with no runs of neighbouring input records.
3. **Neighbour rate.** We measured the share of consecutive output records that
   were within 100 records of each other in the same input file, and compared
   it with random shuffles of the same records:

   | Output | Neighbour rate | Random shuffles |
   |---|---|---|
   | Spark, validation | 0.31% | 0.34% |
   | Spark, training (8 and 64 shards) | 0.06% | 0.05-0.06% |
   | Beam, validation (3 runs) | 0.30-0.36% | 0.29-0.39% (95% range of 200 shuffles) |
   | Beam, training | 0.06% | 0.06% |

All outputs passed all three checks, so both scripts produce a properly random
order.

### Reproducibility

The Spark output is deterministic. Each record's shard comes from its SHA1, and
each shard is sorted by SHA1, so the same inputs and shard count always give
the same records in the same order:

- Two Spark runs on the validation set gave identical records in identical
  order in every shard.
- Training runs on 4 cores (`local[4]`) and 64 cores (`local[*]`) gave
  identical output in all 8 shards. The core count does not change the result.
- The same shuffle on AWS EMR (2 nodes) and on the local server (64 cores) gave
  identical output.
- The `.gz` files themselves differ by checksum between runs, because gzip
  stores a timestamp in each file header. Compare decompressed records instead.

Beam's output order changes on every run, even with identical settings.

### AWS EMR test

`run_emr_ec2.sh` shuffled the validation set from S3 on a temporary EMR
cluster: emr-7.5.0, 1 m5.xlarge primary node and 2 m5.2xlarge core nodes with
100 GB gp3 each, using the account's existing EMR roles through
`SERVICE_ROLE` and `INSTANCE_PROFILE`. Provisioning and bootstrap took about
3.5 minutes and the shuffle step about 5.7 minutes, for roughly $0.25. The
output passed the same checks as above, and was identical, record for record
and in the same order, to the output of the same shuffle on the local server.

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
the same disk as the spill, add another 1x the input size of disk space. To cap
memory, use fewer cores, for example `--master 'local[4]'`.

## Option 2: AWS EMR cluster

A cluster spreads the same job across several machines, and retries failed
tasks on other nodes, which helps on long runs. `run_emr_ec2.sh` creates a
temporary cluster, runs the shuffle as one Spark step, then terminates the
cluster:

```bash
BUCKET=s3://my-bucket \
INPUT='s3://my-bucket/out/training_set.with_label.tfrecord-?????-of-?????.gz' \
OUT_DIR=s3://my-bucket/out/shuffled NAME=training_set \
SUBNET_ID=subnet-xxxxxxxx WORKER_COUNT=4 WORKER_EBS_GB=500 \
  tools/aws_shuffle/run_emr_ec2.sh
```

The cluster installs the dependencies (`google-crc32c psutil boto3`) through
`emr_bootstrap.sh`. By default it uses the standard EMR roles,
`EMR_DefaultRole` and `EMR_EC2_DefaultRole`, which
`aws emr create-default-roles` creates. Creating them needs IAM permissions,
so ask your AWS admin if you don't have them. If your account uses other EMR
roles, set `SERVICE_ROLE` and `INSTANCE_PROFILE`. To use an existing cluster
that was created with that bootstrap action, set `CLUSTER_ID=j-...`. On EMR,
Spark already spills to the EBS volumes under `/mnt`, so
`WORKER_COUNT x WORKER_EBS_GB` is the total spill space.

## Resource needs

- **Memory** depends on the number of tasks running at once (about 0.6 GB each)
  plus the Spark JVM, not on dataset size. On the training set, 4 cores peaked
  at 3.9 GB and 64 cores with 64 shards at 38.9 GB.
- **Disk** for the spill is about 1.5x the gzipped input. This was measured on
  80,000 synthetic DeepVariant-sized examples (3.0 GB gzipped, about 12 GB
  raw), which spilled 4.4 GB.

## Tuning

- `--num_output_shards` (default: about 256 MB of gzipped input per shard).
  Each shard is sorted and gzipped by one core, so with fewer shards than cores
  the last stage leaves cores idle. Set it to at least the number of cores for
  speed. The shard count changes which records go into which shard, so keep it
  fixed if you need identical output across runs.
- `spark.python.worker.memory` (Spark's default is 512m; the EMR script sets
  1g) is the threshold at which the Python-side sort spills to disk.
- Input records are CRC-checked by default. Use `--no_verify_crc` to skip that.

## Notes on Beam DirectRunner

With `apache_beam` 2.61 and later, DirectRunner runs the shuffle pipeline on
Beam's Prism runner. We found two things to watch for:

- **Long jobs fail** with `StatusCode.DEADLINE_EXCEEDED` once they run past
  `--job_server_timeout` (300 s by default in 2.71). Set it to a large value,
  for example `--job_server_timeout=86400`. The validation set (2:26) passes
  without it; the training set (9:52) failed at 5 minutes.
- **`--direct_num_workers` has no effect.** On the validation set, `=0` (all 64
  cores) and `=1` both took 2:26 and kept about 2-3 cores busy.

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
