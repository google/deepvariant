# DeepVariant usage guide

## Overview

DeepVariant is a set of programs used to transform aligned sequencing reads into
variant calls. At the highest level, a user needs to provide three inputs:

1.  A reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
    format and its corresponding [.fai index
    file](http://www.htslib.org/doc/faidx.html) generated using the `samtools
    faidx` command.

1.  An aligned reads file in [BAM](http://genome.sph.umich.edu/wiki/BAM) format
    and its corresponding index file (.bai). The reads must be aligned to the
    reference genome described above.

1.  A model checkpoint for DeepVariant.

The output of DeepVariant is a list of all variant calls in
[VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format.

DeepVariant is composed of three programs: `make_examples`, `call_variants`, and
`postprocess_variants`. Each program's function and example usage is described
in detail in the [quick start].

## Inputs and outputs

### General notes

*   Sharded files are a single logical collection of files with a common naming
    convention. For example, we talk about `filename@10` as a single 10-way
    sharded file named `filename`. On most filesystems this actually looks like
    10 distinct files `filename-00000-of-00010`, ..., `filename-00009-of-00010`.
    DeepVariant can write sharded files using their `filename@10`-style name and
    can read sharded files using both that style as well as the glob form, such
    as `filename-*` or `filename-*-of-00010`.
*   Files with the `.gz` suffix are interpreted as being compressed with gzip
    and are read/written accordingly.
*   All of the DeepVariant tools can read/write their inputs/outputs from local
    disk as well as directly from a cloud storage bucket. Directly accessing BAM
    and reference files from cloud storage is currently not as efficient as we'd
    like and will certainly improve over time, both in the current
    implementation as well as when [htsget] is widely available. See the
    discussion below on reading directly from cloud buckets for details.

### make_examples

`make_examples` consumes reads and the reference genome to create TensorFlow
examples for evaluation with our deep learning models. The tf.Example protos are
written out in TFRecord format.

`make_examples` is a single-threaded program using 1-2 GB of RAM. Since the
process of generating examples is embarrassingly parallel across the genome,
`make_examples` supports sharding of its input and output via the `--task`
argument with a sharded output specification. For example, if the output is
specified as `--examples examples.tfrecord@10.gz` and `--task 0`, the input to
the program will be 10% of the regions and the output will be written to
`examples.tfrecord-00000-of-00010.gz`. A concrete example of using multiple
processes and sharded data is given in the [quick start].

The reference genome FASTA, passed in using the `--ref` flag, must be indexed
and can either be uncompressed or compressed with bgzip.

### call_variants

`call_variants` consumes TFRecord file(s) of tf.Examples protos created by
`make_examples` and a deep learning model checkpoint and evaluates the model on
each example in the input TFRecord. The output here is a TFRecord of
CallVariantsOutput protos. `call_variants` doesn't directly support sharding its
outputs, but accepts a glob or shard-pattern for its inputs.

`call_variants` uses around 4 GB per process and uses TensorFlow for evaluation.
When evaluating a model in CPU mode, TensorFlow can make use of multiple cores,
but scaling is sub-linear. In other words, `call_variants` on a 64 core machine
is less than 8x faster than running on a 8 core machine.

When using a GPU, `call_variants` is both faster, more efficient, and needs
fewer CPUs. Based on a small number of experiments, currently the most efficient
configuration for `call_variants` on a GPU instance is 4-8 CPUs and 1 GPU.
Compared to our setting in the [whole genome case study], we noticed a 2.5x
speedup on the call_variants step using a single K80 GPU and 8 CPUs. Note that
currently `call_variants` can only use one GPU at most. So it doesn't improve
the speed if you get a multiple-GPU machine.

More speed and cost details when running DeepVariant with GPU can be found on
[Running DeepVariant on Google Cloud Platform]. The following steps show what to
do if you want to manually create an instance with GPU. In order to run with
GPU, here are a few changes you want to make in addition to the instructions in
the [whole genome case study] or the [exome case study].

1.  To use GPUs, please read the [GPUs on Compute
    Engine](https://cloud.google.com/compute/docs/gpus/) page first. For
    example, you'll need to request GPU quota first.

1.  Request a machine with GPU.

    If you want to use command line to start a machine, here is an example we
    used:

    ```
    gcloud beta compute instances create "deepvariant-gpu-run" \
    --scopes "compute-rw,storage-full,cloud-platform" \
    --image-family "ubuntu-1604-lts" --image-project "ubuntu-os-cloud" \
    --machine-type "n1-standard-8" --boot-disk-size "300" \
    --boot-disk-type "pd-ssd" \
    --boot-disk-device-name "deepvariant-gpu-run" \
    --zone "us-west1-b" \
    --accelerator type=nvidia-tesla-k80,count=1 \
    --maintenance-policy "TERMINATE"
    ```

1.  When setting up DeepVariant, make sure you run `run-prereq.sh` with
    `DV_GPU_BUILD=1` and `DV_INSTALL_GPU_DRIVERS=1`. `DV_GPU_BUILD` gets you a
    GPU-enabled TensorFlow build, and `DV_INSTALL_GPU_DRIVERS` installs the GPU
    drivers.

    ```
    DV_GPU_BUILD=1 DV_INSTALL_GPU_DRIVERS=1 bash run-prereq.sh
    ```

    This helps you set up an environment with GPU-enabled TensorFlow on
    Ubuntu 16. You won't need `DV_INSTALL_GPU_DRIVERS` if you already installed
    and configured your machine with GPU capabilities. You can read [NVIDIA
    requirements to run TensorFlow with GPU
    support](https://www.tensorflow.org/install/install_linux#nvidia_requirements_to_run_tensorflow_with_gpu_support)
    for more details.

1.  (Optional) You can install and run `nvidia-smi` to confirm that your GPU
    drivers are correctly installed.

    ```
    sudo apt-get -y install nvidia-smi
    ```

    After you install, run `nvidia-smi` to make sure that your GPU drivers are
    set up properly before you proceed.

### postprocess_variants

`postprocess_variants` reads all of the output TFRecord files from
`call_variants`, sorts them, combines multi-allelic records, and writes out a
VCF file.

Because `postprocess_variants` combines and sorts the output of
`call_variants`, it needs to see all of the outputs from `call_variants` for a
single sample to merge into a final VCF. `postprocess_variants` is
single-threaded and needs a non-trivial amount of memory to run (20-30 GB), so
it is best run on a single/dual core machine with sufficient memory.

## Updates on DeepVariant since precisionFDA truth challenge and bioRxiv preprint

The DeepVariant team has been hard at work since we first presented the method.
Key changes and improvements include:

*   Rearchitected with open source release in mind
*   Built on [TensorFlow]
*   Increased variant calling accuracy, especially for indels
*   Vastly faster with reduced memory usage

We have made a number of improvements to the methodology as well. The biggest
change was to move away from RGB-encoded (3-channel) pileup images and instead
represent the aligned read data using a multi-channel tensor data layout. We
current represent the data as a 7-channel raw tensor in which we encode:

*   The read base (A, C, G, T)
*   The base's quality score
*   The read's mapping quality score
*   The read's strand (positive or negative)
*   Does the read support the allele being evaluated?
*   Does the base match the reference genome at this position?
*   The CIGAR operation length for the current op

These are all readily derived from the information found in the BAM file
encoding of each read.

Additional modeling changes were to move to the inception-v3 architecture and to
train on many more independent sequencing replicates of the ground truth
training samples, including 50% downsampled versions of each of those read sets.
In our testing this allowed the model to better generalize to other data types.

In the end these changes reduced our error rate by more than 50% on the held out
evaluation sample (NA24385 / HG002) as compared to our results in the
[PrecisionFDA Truth Challenge](https://precision.fda.gov/challenges/truth/results/):

DeepVariant April 2016 (HG002, GIAB v3.2.2, b37):

Type  | Recall   | Precision | F1_Score
----- | -------- | --------- | --------
INDEL | 0.987882 | 0.991728  | 0.989802
SNP   | 0.999447 | 0.999728  | 0.999587

DeepVariant December 2017 (HG002, GIAB v3.2.2, b37):

Type  | Recall   | Precision | F1_Score
----- | -------- | --------- | --------
INDEL | 0.993081 | 0.994954  | 0.994017
SNP   | 0.999759 | 0.999881  | 0.999820

See the [whole genome case study], which we update with each release of
DeepVariant, for the latest results.

## Optimizing DeepVariant

For educational purposes, the DeepVariant Case Study uses the simplest, but
least efficient, computing configuration for DeepVariant. At scale, DeepVariant
can be run far more efficiently by following a few simple guidelines:

-   The work performed by `make_examples` can be split across on any number of
    machines, using all NCPU processes on each machine. For example, to run
    efficiently on `M` machines, each with `C` cores, the number of shards in
    the output file should be `N = M*C`. Each process should pass as its output
    file flag `--examples output@N`, with the `--task` flag set uniquely for
    each process (e.g., for core `c` on machine `m`, where `0 <= c < C` and `0
    <= m < M`), the task should be `m*C + c`. The number of jobs / shards should
    be balanced by the time needed to start up a running instance of DeepVariant
    (especially when localizing the inputs to the instance).

-   `make_examples` with `--task i` processes every ith interval on the genome.
    This means that each process will access the entire genome (in pieces). This
    provides excellent load-balancing across tasks, but it means that each
    process needs access to the entire reference and BAM file. Localization
    costs can be overcome by directly accessing files on cloud storage, or
    splitting up the inputs by chromosome and then processing those with
    `--task` and `--regions` set as appropriate. Note that the DeepVariant team
    doesn't recommend going down the per-chromosome BAM route, as direct cloud
    storage access is a more attractive long-term approach.

-   `call_variants` can read any number of TFRecord files from `make_examples`.
    As in the [whole genome case study], we run `make_examples` 64 ways
    parallel, writing an `examples@64` file. We then run `call_variants
    --examples examples@64` to evaluate all of the examples across all 64
    example files. Given the sub-linear scaling with CPUs, the most
    cost-efficient way to run `call_variants` is either on a single GPU instance
    or one with 4 CPUs.

### Optimizing for turnaround time

*   Run `make_examples` on a large number of machines, reading the BAM and
    reference genome from cloud storage directly, to a sharded output file
    `examples@N`.
*   Run multiple `call_variants` jobs on 4-8 CPU / 1 GPU machines, providing
    each separate `call_variants` job a subset of the example outputs. For
    example, if `N` separate `call_variants` jobs are run, each job `J` reads
    `examples-0000J-of-0000N`, and writes output to `calls-0000J-of-0000N`).
*   Run `postprocess_variants` when all `call_variants` jobs have finished,
    merging the `calls@N` into a single VCF on a 30 GB (or so) machine with 1-2
    cores.

### Optimizing for cost

*   Use a job manager that supports preemptable VMs (i.e., can detect a machine
    going down and restart the dead job on a new machine) for `make_examples`.
*   Explore the optimal cost trade-off between preemptable CPU/GPU instances of
    `call_variants`.
*   Run `postprocess_variants` on a preemptable VM as noted above.

### Example: Running the Case Study distributed across machines

To be totally concrete, suppose we wanted to rerun the Case Study but
distributed across many machines.

- Run `make_examples` on 8 64 core machines, each running `make_examples` with
`parallel` or equivalent with 64 processes on each machine. This would complete
in roughly 45 minutes each. The output TFRecords should be written directly to
cloud storage.
- Run `call_variants` on 8 8-CPU/1-GPU instances, each processing 64 shards of
the 8x64 sharded examples. This completes in roughly 75 minutes. The inputs
should be read directly from cloud storage and the outputs should be written
directly to cloud storage. An alternative CPU-only workflow would be to run on
64 4-CPU machines.
- Run a single `postprocess_variants` on a 2 core 30 GB machine, completing in
30-45 minutes.

It's possible to run on more machines for `make_examples` and `call_variants` if
a shorter turn-around-time is needed. At some point scaling limits will be hit,
either due to localizing BAMs, startup costs on the machines, reading too many
sharded inputs, etc. The optimal configuration for minimizing cost or turnaround
time is currently an open question and likely depends on input data (e.g. whole
genome vs exome, sequencing coverage, etc.).

Note: All of the specific machine type recommendations here are based on
empirical results on Google Cloud Platform and are likely to be different on
different machine types. Moreover, the optimal configuration may change over
time as DeepVariant, TensorFlow, and the CPU, GPU, and
[TPU](https://cloud.google.com/tpu/) hardware evolves.

## Reading directly from Google Cloud Storage (GCS) buckets

DeepVariant can read SAM/BAM and fasta sources from GCS buckets, as well as read
and write TFRecord (i.e., TensorFlow examples) to/from GCS buckets. In fact, in
general DeepVariant can read SAM/BAM/FASTA files from any source supported by
[htslib] and/or TensorFlow.

DeepVariant runs slightly slower when processing reads directly from GCS; our
current best estimates are roughly a ~3-5% slowdown for an end-to-end run of
`make_examples` compared to first localizing the entire BAM to a local
persistent SSD on the cloud instance. Consequently, the choice of whether to
first copy the BAM to local disk or directly read from a GCS bucket depends on
how one is running DeepVariant. If running `make_examples` on a single big
multi-core machine as in the Case Study examples, it is likely advantageous to
first localize the BAM as copying the full 50x BAM to local disk only takes
10-20 minutes, a small overhead compared to the end-to-end runtime of
`make_examples`. However, more advanced distributed runs of DeepVariant, such as
running distributed DeepVariant with [dsub] and the [Pipelines API] as
exemplified by the [DeepVariant with docker], can benefit
enormously from this capability. When there are N instances running
`make_examples` in a distributed fashion, the overhead of localization is N
times larger, all the worse since each `make_examples` run is only accessing
roughly 1 / N of the BAM reads. Moreover, distributing `make_examples` over many
machines enables far more parallelizing than even the largest instance, so the
per-instance runtime can be only marginally longer than the cost to localize the
BAM itself. In this situation direct GCS access is likely the best option
despite its slight overhead.

As of 12/4/17, there are some limitations to remote file access that we are
working to overcome:

*   GCS BAM: there's a race condition in [htslib] when running multiple
    `make_examples` jobs `parallel` on a single instance. [htslib] downloads the
    BAI file to local disk, and multiple parallel executions of `make_examples`
    will all try to load this BAI file at the same time, resulting in data
    corruption. In this situation it is safest to simply localize the BAI file
    (but not the BAM) to the directory where you are running `make_examples`.
    This is not necessary if you are just running a single `make_examples`
    process on the machine or run each `make_examples` job in separate working
    directories.
*   GCS fasta: there is a problem in [htslib] where the connection to the FASTA
    file is getting dropped, so long-running `make_examples` can fail using a
    remote FASTA can fail. We are looking at how to fix this issue.

Note these concerns do not apply to persistent mounted remote filesystems like
NFS or Fuse, only to non-POSIX remote systems like GCS or FTP.

As of today, the safest, performant approach to using remote reads and reference
files with DeepVariant `make_examples` is to localize the reference fasta and
associated fai/gzi metadata as well as the BAI, while leaving the BAM itself
remote, which is almost always the largest data source needed by DeepVariant.

The examples output of `make_examples` and both the input examples and output
calls of `call_variants` can be read from / written to GCS.

Currently `postprocess_variants` supports reading its inputs directly from GCS
but the output must go to a local disk.

[htsget]: http://samtools.github.io/hts-specs/htsget.html
[htslib]: http://www.htslib.org/
[dsub]: https://cloud.google.com/genomics/v1alpha2/dsub
[Pipelines API]: https://cloud.google.com/genomics/v1alpha2/pipelines
[DeepVariant with docker]: deepvariant-docker.md
[exome case study]: deepvariant-exome-case-study.md
[whole genome case study]: deepvariant-case-study.md
[quick start]: deepvariant-quick-start.md
[Running DeepVariant on Google Cloud Platform]: https://cloud.google.com/genomics/deepvariant
[TensorFlow]: http://www.tensorflow.org/
