f# DeepVariant usage guide

## Overview

DeepVariant is a set of programs used to transform aligned sequencing reads into
variant calls. At the highest level, a user needs to provide three inputs:

1.  A reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
    format and its corresponding
    [.fai index file](http://www.htslib.org/doc/faidx.html) generated using the
    `samtools faidx` command.

1.  An aligned reads file in [BAM](http://genome.sph.umich.edu/wiki/BAM) format
    and its corresponding index file (.bai). The reads must be aligned to the
    reference genome described above.

1.  A model checkpoint for DeepVariant.

The output of DeepVariant is a list of all variant calls in
[VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format.

DeepVariant is composed of three programs: `make_examples`, `call_variants`, and
`postprocess_variants`. More details about each program are described in detail
in the [Inputs and outputs](#inputs-and-outputs) section.

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

### make_examples

`make_examples` consumes reads and the reference genome to create TensorFlow
examples for evaluation with our deep learning models. The tf.Example protos are
written out in TFRecord format. To learn more about tf.Example and TFRecord, see
the
[Using TFRecords and tf.Example](https://www.tensorflow.org/tutorials/load_data/tfrecord)
Colab.

`make_examples` is a single-threaded program using 1-2 GB of RAM. Since the
process of generating examples is embarrassingly parallel across the genome,
`make_examples` supports sharding of its input and output via the `--task`
argument with a sharded output specification. For example, if the output is
specified as `--examples examples.tfrecord@10.gz` and `--task 0`, the input to
the program will be 10% of the regions and the output will be written to
`examples.tfrecord-00000-of-00010.gz`.

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
cannot pipe it into DeepVariant. Duplicate marking may be performed, in our
analyses there is almost no difference in accuracy except at lower (<20x)
coverages. Finally, we recommend that you do not perform BQSR. Running BQSR has
a small decrease on accuracy. It is not necessary to do any form of indel
realignment, though there is not a difference in DeepVariant accuracy either
way.

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

`call_variants` consumes TFRecord file(s) of tf.Examples protos created
by `make_examples` and a deep learning model checkpoint and evaluates the model
on each example in the input TFRecord. The output here is a TFRecord of
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
speedup on the call_variants step using a single P100 GPU and 8 CPUs. Note that
currently `call_variants` can only use one GPU at most. So it doesn't improve
the speed if you get a multiple-GPU machine.

### postprocess_variants

`postprocess_variants` reads all of the output TFRecord files from
`call_variants`, sorts them, combines multi-allelic records, and writes out a
VCF file. When [gVCF output](deepvariant-gvcf-support.md) is requested, it also
outputs a gVCF file which merges the VCF with the non-variant sites.

Because `postprocess_variants` combines and sorts the output of `call_variants`,
it needs to see all of the outputs from `call_variants` for a single sample to
merge into a final VCF. `postprocess_variants` is single-threaded and needs a
non-trivial amount of memory to run (20-30 GB), so it is best run on a
single/dual core machine with sufficient memory.

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

You can also see the [Colab example] to see how you can visualize the pileup
images.

## Training data over time

For the models we've released over time, you can find more details about the
training data in
[DeepVariant training data](deepvariant-details-training-data.md).

## CRAM support

As of v0.7, DeepVariant accepts CRAM files as input in addition to BAM files.

As of v0.9.0, we changed the default to use the reference file specified by the
`--ref` flag, instead of the path to the original reference in the CRAM file
(encoded in the file's "UR" tag).

For more information about CRAM, see the
[`samtools` documentation](http://www.htslib.org/doc/samtools.html) in general
but particularly the sections on
[Global Options](http://www.htslib.org/doc/samtools.html#GLOBAL_OPTIONS) and
[reference sequences in CRAM](http://www.htslib.org/doc/samtools.html#REFERENCE_SEQUENCES).

`htslib` also hosts a nice page
[benchmarking CRAM](http://www.htslib.org/benchmarks/CRAM.html) with information
on the effect of different CRAM options on file size and encoding/decoding
performance.

Here are some basic file size and runtime numbers for running a single
`make_examples` job on a 30x whole genome sample in BAM and CRAM format.

Filetype | Size (Gb) | Runtime (min)
-------- | --------- | -------------
BAM      | 66.99     | 79m47.37307s
CRAM     | 37.85     | 96m53.477s
Ratio    | 56.50%    | 121.43%

*   BAM file:
    `gs://deepvariant/performance-testdata/HG002_NIST_150bp_downsampled_30x.bam`
*   CRAM file:
    `gs://deepvariant/performance-testdata/HG002_NIST_150bp_downsampled_30x.cram`

Runtime was measured on
[n1-standard-64](https://cloud.google.com/compute/docs/machine-types#n1_machine_types)
machines.

## Starting from v1.2.0, we include `samtools` and `bcftools`.

Based on user feedback ([GitHub issue #414](https://github.com/google/deepvariant/issues/414)),
we added samtools and bcftools in our Docker image:

```bash
docker run google/deepvariant:"${BIN_VERSION}" samtools
```

and

```bash
docker run google/deepvariant:"${BIN_VERSION}" bcftools
```

You can read more about samtools and bcftools here: http://www.htslib.org/doc/.

## Commands for requesting machines used in case studies

We report runtime in our case studies documentation. In order to make sure the
results we report are reproducible without too much variation, we provide the
commands we used here to show you what kind of machines we ran the case studies
on. This is NOT the fastest or cheapest configuration.

### Command for a CPU-only machine on Google Cloud Platform.

We used a 64-core (vCPU) machine with 240GiB of memory and no GPU, on the Google
Cloud Platform. Specifying the CPU platform also allows us to report the runtime
more consistently.

```shell
gcloud compute instances create "${USER}-cpu"  \
  --scopes "compute-rw,storage-full,cloud-platform" \
  --image-family "ubuntu-2204-lts" \
  --image-project "ubuntu-os-cloud" \
  --machine-type "n1-standard-64" \
  --boot-disk-size "300" \
  --zone "us-west1-b" \
  --min-cpu-platform "Intel Skylake"
```

### Command for a GPU machine on Google Cloud Platform

```shell
gcloud compute instances create "${USER}-gpu" \
  --scopes "compute-rw,storage-full,cloud-platform" \
  --maintenance-policy "TERMINATE" \
  --accelerator=type=nvidia-tesla-p100,count=1 \
  --image-family "ubuntu-2204-lts" \
  --image-project "ubuntu-os-cloud" \
  --machine-type "n1-standard-16" \
  --boot-disk-size "300" \
  --zone "us-west1-b" \
  --min-cpu-platform "Intel Skylake"
```

NOTE: Having an instance up and running could cost you. Remember to delete the
instances you're not using. You can find the instances at:
https://console.cloud.google.com/compute/instances?project=YOUR_PROJECT

[exome case study]: deepvariant-exome-case-study.md
[whole genome case study]: deepvariant-case-study.md
[quick start]: deepvariant-quick-start.md
[Running DeepVariant on Google Cloud Platform]: https://cloud.google.com/life-sciences/docs/tutorials/deepvariant
[TensorFlow]: http://www.tensorflow.org/
[Colab example]: visualizing_examples.ipynb
