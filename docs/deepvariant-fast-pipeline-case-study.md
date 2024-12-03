# DeepVariant Fast Pipeline case study

This document goes over the example of using DeepVariant Fast Pipeline with
PacBio data.

## Background

Fast Pipeline is a DeepVariant feature that allows parallelization of the
make_examples and call_variant stages. It is especially useful for machines with
a GPU. Examples are streamed to call_variants inference, allowing simultaneous
utilization of both the CPU and GPU. Please note that this feature is still
experimental.

This setup requires a machine with a GPU. For this case study, we will use a
`n1-standard-16` compute instance with 1 Nvidia P4 GPU. However, this setup is
not optimal, as 16 cores may not be sufficient to fully utilize the GPU. In a
real-life scenario, allocating 32 cores for `make_examples` would ensure better
GPU utilization and improved runtime.

## Provision the compute instance

Here we create Google Cloud compute instance. You may skip this step if you run
the case study on a local computer with GPU.

```bash
gcloud compute instances create "deepvariant-fast-pipeline" \
  --scopes "compute-rw,storage-full,cloud-platform" \
  --maintenance-policy "TERMINATE" \
  --accelerator=type=nvidia-tesla-p4,count=1 \
  --image-family "ubuntu-2204-lts" \
  --image-project "ubuntu-os-cloud" \
  --machine-type "n1-standard-16" \
  --boot-disk-size "100" \
  --zone "us-central1-a"
```

You can then ssh into the machine by running:

```bash
gcloud compute ssh "deepvariant-fast-pipeline" --zone us-central1-a
```

## Install Nvidia drivers and Nvidia container toolkit (optional)

CUDA drivers and NVIDIA Container toolkit are required to run the case study.
Please refer to the following documentation for more details.
[NVIDIA CUDA Installation Guide for Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/),
[Installing the NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

For this case study we used the
[script](https://github.com/google/deepvariant/blob/r1.8.0/scripts/install_nvidia_docker.sh)
that automates the CUDA and container tools kit installation.

Please note that the script takes about 30 minutes to run.

```bash
wget https://raw.githubusercontent.com/google/deepvariant/refs/heads/r1.8.0/scripts/install_nvidia_docker.sh
chmod +x install_nvidia_docker.sh
./install_nvidia_docker.sh
```

## Get Docker image, models, and test data

### Get DeepVariant Docker image

```bash
BIN_VERSION="1.8.0"
sudo docker pull google/deepvariant:"${BIN_VERSION}-gpu"
```

### Download test data

Before you start running, you need to have the following input files:

1.  A reference genome in [FASTA] format and its corresponding index file
    (.fai).

```bash
mkdir -p reference
gcloud storage cp gs://deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna* reference/
```

1.  An aligned reads file in [BAM] format and its corresponding index file
    (.bai). You get this by aligning the reads from a sequencing instrument,
    using an aligner like [BWA] for example.

```bash
mkdir -p input
gcloud storage cp gs://deepvariant/pacbio-case-study-testdata/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam* input/
```

## Run DeepVariant pipeline on chromosome 20 alignments

`fast_pipeline` binary in DeepVariant docker allows to run `make_examples` and
`call_variant` stages of DeepVariant in stream mode. Here is the command line to
run the `fast_pipeline`

### Prepare files containing command line parameters for all `DeepVariant` binaries

Config files below contain all default command line parameters for PacBio data.
`--examples` and `--gvcf` flags are set with the sharded file names. ==It is
important to ensure that the number of shards matches in all config files and
the `--num_shards` flag in the fast_pipeline binary. In our case it is set to
14==

The machine has 16 virtual cores, but we set the number of shards to 14 to
reserve 2 cores for the input pipeline in call_variants. Insufficient CPU
resources for the inference pipeline can cause an input bottleneck, leading to a
slowdown in the inference stage.

```bash
mkdir -p config
FILE=config/make_examples.ini

cat <<EOM >$FILE
--examples=/tmp/examples.tfrecords@14.gz
--gvcf=/tmp/examples.gvcf.tfrecord@14.gz
--mode=calling
--reads=/input/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam
--ref=/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
--alt_aligned_pileup=diff_channels
--max_reads_per_partition=600
--min_mapping_quality=1
--parse_sam_aux_fields
--partition_size=25000
--phase_reads
--pileup_image_width=147
--norealign_reads
--sort_by_haplotypes
--track_ref_reads
--vsc_min_fraction_indels=0.12
--trim_reads_for_pileup
--call_small_model_examples
--trained_small_model_path=/opt/smallmodels/pacbio
--small_model_snp_gq_threshold=25
--small_model_indel_gq_threshold=30
--small_model_vaf_context_window_size=51
--output_phase_info
--checkpoint=/opt/models/pacbio
--regions=chr20
EOM
```

```bash
FILE=config/call_variants.ini

cat <<EOM >$FILE
--outfile=/output/case_study.cvo.tfrecord.gz
--checkpoint=/opt/models/pacbio
--batch_size=1024
--writer_threads=1
EOM
```

```bash
FILE=config/postprocess_variants.ini

cat <<EOM >$FILE
--ref=/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
--infile=/output/case_study.cvo.tfrecord.gz
--nonvariant_site_tfrecord_path=/tmp/examples.gvcf.tfrecord@14.gz
--outfile=/output/variants.chr20.vcf
--gvcf_outfile=/output/variants.gvcf.chr20.vcf
--small_model_cvo_records=/tmp/examples_call_variant_outputs.tfrecords@14.gz
--cpus=14
EOM
```

```bash
time sudo docker run \
  -v "${PWD}/config":"/config" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v /tmp:/tmp \
  --gpus all \
  -e DV_BIN_PATH=/opt/deepvariant/bin \
  --shm-size=2gb \
  google/deepvariant:"${BIN_VERSION}-gpu" \
    /opt/deepvariant/bin/fast_pipeline \
    --make_example_flags /config/make_examples.ini \
    --call_variants_flags /config/call_variants.ini \
    --postprocess_variants_flags /config/postprocess_variants.ini \
    --shm_prefix dv \
    --num_shards 14 \
    --buffer_size 10485760 \
    2>&1  | tee /tmp/fast_pipeline.docker.log
```

*   `-v` allows to map local directory inside docker container.
*   `-e` we need to set `DV_BIN_PATH` environment variable to point to
    DeepVariant binaries directory inside the container.
*   `--shm-size` sets the size of shared memory available to the container. It
    has to be larger than `--buffer_size` x `--num_shards`. In our case
    buffer_size is 10M and we run 14 shards, so 2gb would be large enough to
    accommodate buffers and all synchronization objects for each shard.

#### `fast_pipeline` command line parameters:

*   `--make_example_flags` - path to the file containing `make_examples` command
    line parameters.
*   `--call_variants_flags` - path to the file containing `call_variants`
    command line parameters.
*   `--postprocess_variants_flags` - path to the file containing
    `postprocess_variants` command line parameters.
*   `--shm_prefix` - prefix for shared memory files. It is an arbitrary name.
*   `--num_shards` - number of parallel processes to run `make_examples`
*   `--buffer_size` - shared memory buffer size for each process.

On a successful completion the `output` directory will contain two VCF files:

```
variants.chr20.vcf
variants.gvcf.chr20.vcf
```

With the same settings the pipeline takes approximately 10 minutes.

```
real    8m15.252s
user    0m0.007s
sys     0m0.035s
```

## Benchmark output

Download benchmark data:

```bash
mkdir -p benchmark
FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

```
HAPPY_VERSION=v0.3.12

time sudo docker run \
  -v ${PWD}/output:/output \
  -v ${PWD}/benchmark:/benchmark \
  -v ${PWD}/reference:/reference \
  jmcdani20/hap.py:${HAPPY_VERSION} \
  /opt/hap.py/bin/hap.py \
    /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    /output/variants.chr20.vcf \
    -f /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    -r /reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -o /output/happy.output \
    --engine=vcfeval \
    --pass-only \
    -l "chr20"
```

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10543        85        22403        74      11375     40     29       0.992002          0.993290        0.507744         0.992646                     NaN                     NaN                   1.748961                   2.138647
INDEL   PASS        10628     10543        85        22403        74      11375     40     29       0.992002          0.993290        0.507744         0.992646                     NaN                     NaN                   1.748961                   2.138647
  SNP    ALL        70166     70101        65       105602        71      35342     12     12       0.999074          0.998989        0.334672         0.999032                2.296566                1.713281                   1.883951                   1.503192
  SNP   PASS        70166     70101        65       105602        71      35342     12     12       0.999074          0.998989        0.334672         0.999032                2.296566                1.713281                   1.883951                   1.503192
```
