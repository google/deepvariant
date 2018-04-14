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
    disk as well as directly from a cloud storage bucket (through gcsfuse). In
    the [whole genome case study] and [exome case study], we didn't use gcsfuse
    because there is not much benefit on a single machine compared to just
    copying the files first.

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

#### Input assumptions

`make_examples` requires its input files to satisfy a few basic requirements to
be processed correctly.

First, the reference genome FASTA, passed in using the `--ref` flag, must be
indexed and can either be uncompressed or compressed with bgzip.

Second, the BAM file provided to `--reads` should be aligned to a "compatible"
version of the genome reference provided as the `--ref`. By compatible here we
mean the BAM and FASTA share at least a reasonable set of common contigs, as
DeepVariant will only process contigs shared by both the BAM and reference. As
an example, suppose you have a BAM file mapped to b37 + decoy FASTA and you
provide just the vanilla b37 fasta to `make_examples`. DeepVariant will only
process variants on the shared contigs, effectively excluding the hs37d5 contig
present in the BAM but not in the reference.

The BAM file must be also sorted and indexed. It must exist on disk, so you
cannot pipe it into DeepVariant. We currently recommend that the BAM be
duplicate marked, but it's unclear if this is even necessary. Finally, it's not
necessary to recalibrate the base qualities or do any form of indel realignment.

Third, if you are providing `--regions` or other similar arguments these should
refer to contigs present in the reference genome. These arguments accept
space-separated lists, so all of the follow examples are valid arguments for
`--regions` or similar arguments:

*   `--regions chr20` => only process all of chromosome 20
*   `--regions chr20:10,000,000-11,000,000` => only process 10-11mb of chr20
*   `--regions "chr20 chr21"` => only process chromosomes 20 and 21

Fourth and finally, if running in training mode the `truth_vcf` and
`confident_regions` arguments should point to VCF and BED files containing the
true variants and regions where we are confident in our calls (i.e., calls
within these regions and not in the truth_vcf are considered false positives).
These should be bgzipped and tabix indexed and be on a reference consistent with
the one provided with the `--ref` argument.

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
VCF file. When [gVCF output](deepvariant-gvcf-support.md) is requested, it also
outputs a gVCF file which merges the VCF with the non-variant sites.

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
currently represent the data as a 6-channel raw tensor in which we encode:

*   The read base (A, C, G, T)
*   The base's quality score
*   The read's mapping quality score
*   The read's strand (positive or negative)
*   Does the read support the allele being evaluated?
*   Does the base match the reference genome at this position?

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

Type  | # FN | # FP | Recall   | Precision | F1_Score
----- | ---- | ---- | -------- | --------- | --------
INDEL | 4175 | 2839 | 0.987882 | 0.991728  | 0.989802
SNP   | 1689 | 832  | 0.999447 | 0.999728  | 0.999587

DeepVariant December 2017 (HG002, GIAB v3.2.2, b37):

Type  | # FN | # FP | Recall   | Precision | F1_Score
----- | ---- | ---- | -------- | --------- | --------
INDEL | 2384 | 1811 | 0.993081 | 0.994954  | 0.994017
SNP   | 735  | 363  | 0.999759 | 0.999881  | 0.999820

See the [whole genome case study], which we update with each release of
DeepVariant, for the latest results.

You can also see the [Datalab example] to see how you can visualize the pileup
images.

## Training data over time

For the models we've released over time, here are more details of the training
data we used.

### WGS models

version | Replicates                             | #examples
------- | -------------------------------------- | -----------
v0.4    | 9 HG001                                | 85,323,867
v0.5    | 9 HG001<br>2 HG005<br>78 HG001 WES<br>1 HG005 WES<sup>[(1)](#myfootnote1)</sup> | 115,975,740
v0.6    | 9 HG001<br>2 HG005<br>5 HG001 PCR+     | 156,571,227

### WES models

version | Replicates                  | #examples
------- | --------------------------- | ------------------------------
v0.5    | 78 HG001 WES<br>1 HG005 WES | 15,714,062
v0.6    | 78 HG001 WES<br>1 HG005 WES<sup>[(2)](#myfootnote2)</sup> | 15,705,449


<a name="myfootnote1">(1)</a>: In v0.5, we experimented with adding whole exome sequencing data into training data. In v0.6, we took it out because it didn't improve the WGS accuracy.

<a name="myfootnote2">(2)</a>: The training data are from the same replicates as v0.5. The number of examples changed because of the update in [haplotype_labeler](https://github.com/google/deepvariant/tree/r0.6/deepvariant/labeler/haplotype_labeler.py).

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

## Reading from Google Cloud Storage (GCS) using gcsfuse

[gcsfuse](https://github.com/GoogleCloudPlatform/gcsfuse) allows access to files
in the GCS bucket much like local files.

### Install gcsfuse

Install gcsfuse package as explained
[here](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/installing.md)
(recommended).

### Mount GCS Bucket

First ensure your "application default" credential for GCP is available. The
easiest way to set this up is to run gcloud tool.

```
gcloud auth application-default login
```

Then mount the bucket to a dir on your local machine.

```
mkdir -p $local_dir                  # E.g. local_dir=$HOME/input
gcsfuse $gcs_bucket_name $local_dir  # E.g. gcs_bucket_name=deepvariant
```

Confirm it is mounted successfully by doing `ls` on `$local_dir`. If there is a
problem, you may want to unmount and mount again.

```
fusermount -u $local_dir  # Linux
umount $local_dir         # OS X
```

Now you can run DeepVariant as if the input files are stored locally under
`$local_dir`.

### Known Issue

Currently gcsfuse does not work well with DeepVariant in multi-tasking scenarios
, i.e. one gcsfuse daemon and multiple DeepVariant tasks (see the
[issue](https://github.com/GoogleCloudPlatform/gcsfuse/issues/262)). There are
two workarounds for this:

1.  Run multiple gcsfuse; one gcsfuse daemon for every DeepVariant task.

    ```
    for i in `seq 1 $num_task`;do
      mkdir -p $local_dir_task_$i
      gcsfuse $gcs_bucket_name $local_dir_task_$i
    done
    ```

2.  Reduce read sizes in DeepVariant tasks. For example, we observed that the
    issue disappears for `make_examples` if hts block size is reduced to `128KB`
    (from default `128MB`). Use `--hts_block_size` flag for this.

### Useful gcsfuse links

*   [gcsfuse installing
    doc](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/installing.md)
*   [mounting
    doc](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/mounting.md)
    (if you exprienced credentional issue)

[dsub]: https://cloud.google.com/genomics/v1alpha2/dsub
[Pipelines API]: https://cloud.google.com/genomics/v1alpha2/pipelines
[DeepVariant with docker]: deepvariant-docker.md
[exome case study]: deepvariant-exome-case-study.md
[whole genome case study]: deepvariant-case-study.md
[quick start]: deepvariant-quick-start.md
[Running DeepVariant on Google Cloud Platform]: https://cloud.google.com/genomics/deepvariant
[TensorFlow]: http://www.tensorflow.org/
[Datalab example]: visualizing_examples.ipynb
