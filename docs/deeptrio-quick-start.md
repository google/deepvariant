# DeepTrio quick start

This document explains how to quickly start using
[DeepTrio](deeptrio-details.md) to generate variant calls for trio samples. This
tutorial does not cover all possible settings of DeepTrio. It is intended to be
a starting point for using DeepTrio.

## Background

To get started, we've provided a Docker image, and some test data in a bucket on
Google Cloud Storage. The instructions below show how to download the data
through the corresponding public URLs.

This setup requires a machine with the AVX instruction set. To see if your
machine meets this requirement, you can check the `/proc/cpuinfo` file, which
lists this information under "flags". If you do not have the necessary
instructions, see the next section for more information on how to build your own
Docker image.

### Use Docker to run DeepTrio in one command.

Although DeepTrio can be built from a source, we provide a docker image that
allows to run through all steps in one command to generate VCF/gVCF output files
from input BAM files and the reference.

If you want to compile the binaries for yourself, we also have a [Dockerfile]
that you can use to build your own Docker image. You can read the [docker build]
documentation on how to build.

## Get Docker image, models, and test data

### Get Docker image

```bash
BIN_VERSION="1.6.0"

sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:deeptrio-"${BIN_VERSION}"
```

### Download test data

Before you start, you need to have the following input files:

1.  A reference genome in [FASTA] format and its corresponding index file
    (.fai).

1.  For each sample, one aligned reads file in [BAM] format and its
    corresponding index file (.bai). You get this by aligning the reads from a
    sequencing instrument, using an aligner like [BWA] for example.

We've prepared a small test data bundle for use in this quick start guide that
can be downloaded to your instance from the public URLs.

Download the test bundle:

```bash
INPUT_DIR="${PWD}/quickstart-testdata"
mkdir -p ${INPUT_DIR}

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio

curl ${FTPDIR}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > "${INPUT_DIR}"/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > "${INPUT_DIR}"/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > "${INPUT_DIR}"/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

curl ${FTPDIR}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > "${INPUT_DIR}"/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > "${INPUT_DIR}"/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > "${INPUT_DIR}"/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

curl ${FTPDIR}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > "${INPUT_DIR}"/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > "${INPUT_DIR}"/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > "${INPUT_DIR}"/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

HTTPDIR=https://storage.googleapis.com/deepvariant/quickstart-testdata

wget -P ${INPUT_DIR} "${HTTPDIR}"/HG002.chr20.10_10p1mb.bam
curl "${HTTPDIR}/HG002.chr20.10_10p1mb.bam" > "${INPUT_DIR}/HG002.chr20.10_10p1mb.bam"
curl "${HTTPDIR}/HG002.chr20.10_10p1mb.bam.bai" > "${INPUT_DIR}/HG002.chr20.10_10p1mb.bam.bai"

curl "${HTTPDIR}/HG003.chr20.10_10p1mb.bam" > "${INPUT_DIR}/HG003.chr20.10_10p1mb.bam"
curl "${HTTPDIR}/HG003.chr20.10_10p1mb.bam.bai" > "${INPUT_DIR}/HG003.chr20.10_10p1mb.bam.bai"

curl "${HTTPDIR}/HG004.chr20.10_10p1mb.bam" > "${INPUT_DIR}/HG004.chr20.10_10p1mb.bam"
curl "${HTTPDIR}/HG004.chr20.10_10p1mb.bam.bai" > "${INPUT_DIR}/HG004.chr20.10_10p1mb.bam.bai"

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > "${INPUT_DIR}"/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > "${INPUT_DIR}"/GRCh38_no_alt_analysis_set.fasta.fai
```

This should create a subdirectory in the current directory containing the actual
data files:

```bash
ls -1 ${INPUT_DIR}
```

output:

```
GRCh38_no_alt_analysis_set.fasta
GRCh38_no_alt_analysis_set.fasta.fai
HG002.chr20.10_10p1mb.bam
HG002.chr20.10_10p1mb.bam.bai
HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG003.chr20.10_10p1mb.bam
HG003.chr20.10_10p1mb.bam.bai
HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG004.chr20.10_10p1mb.bam
HG004.chr20.10_10p1mb.bam.bai
HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

## Run DeepTrio with one command

We create one entrypoint that can be directly run as a docker command. If you
want to see the details, you can read through [run_deeptrio.py].

```bash
OUTPUT_DIR="${PWD}/quickstart-output"
mkdir -p "${OUTPUT_DIR}"
```

You can run everything with the following command:

```bash
sudo docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}":"/output" \
  google/deepvariant:deeptrio-"${BIN_VERSION}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=WGS \
  --ref=/input/GRCh38_no_alt_analysis_set.fasta \
  --reads_child=/input/HG002.chr20.10_10p1mb.bam \
  --reads_parent1=/input/HG003.chr20.10_10p1mb.bam \
  --reads_parent2=/input/HG004.chr20.10_10p1mb.bam \
  --output_vcf_child /output/HG002.output.vcf.gz \
  --output_vcf_parent1 /output/HG003.output.vcf.gz \
  --output_vcf_parent2 /output/HG004.output.vcf.gz \
  --sample_name_child 'HG002' \
  --sample_name_parent1 'HG003' \
  --sample_name_parent2 'HG004' \
  --num_shards $(nproc)  \
  --regions "chr20:10,000,000-10,010,000" \
  --intermediate_results_dir /output/intermediate_results_dir \
  --output_gvcf_child /output/HG002.g.vcf.gz \
  --output_gvcf_parent1 /output/HG003.g.vcf.gz \
  --output_gvcf_parent2 /output/HG004.g.vcf.gz
```

NOTE: If you want to look at all the commands being run, you can add
`--dry_run=true` to the command above, which will print out all the commands
but not execute them.

This will generate 15 files and 1 directory in `${OUTPUT_DIR}`:

```bash
ls -1 ${OUTPUT_DIR}
```

output:

```
HG002.g.vcf.gz
HG002.g.vcf.gz.tbi
HG002.output.vcf.gz
HG002.output.vcf.gz.tbi
HG002.output.visual_report.html
HG003.g.vcf.gz
HG003.g.vcf.gz.tbi
HG003.output.vcf.gz
HG003.output.vcf.gz.tbi
HG003.output.visual_report.html
HG004.g.vcf.gz
HG004.g.vcf.gz.tbi
HG004.output.vcf.gz
HG004.output.vcf.gz.tbi
HG004.output.visual_report.html
intermediate_results_dir
```

The directory "intermediate_results_dir" exists because
`--intermediate_results_dir /output/intermediate_results_dir` is specified. This
directory contains the intermediate output of make_examples and call_variants
steps.

For more information about the `HG00*.output.visual_report.html` files, see the
[VCF stats report documentation](deepvariant-vcf-stats-report.md).

## Notes on GPU image

If you are using GPUs, you can pull the GPU version, and make sure you run with
`--gpus 1`. `call_variants` is the only step that uses the GPU, and can only use
one at a time. `make_examples` and `postprocess_variants` do not run on GPU.

For an example to install GPU driver and docker, see [install_nvidia_docker.sh].

```
sudo docker run --gpus 1 \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  google/deepvariant:deeptrio-"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  ...
```

## Notes on Singularity

### CPU version

```
# Pull the image.
singularity pull docker://google/deepvariant:deeptrio-"${BIN_VERSION}"

# Run DeepTrio.
singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
  docker://google/deepvariant:deeptrio-"${BIN_VERSION}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=WGS \
  --ref="${INPUT_DIR}"/GRCh38_no_alt_analysis_set.fasta \
  --reads_child="${INPUT_DIR}"/HG002.chr20.10_10p1mb.bam \
  --reads_parent1="${INPUT_DIR}"/HG003.chr20.10_10p1mb.bam \
  --reads_parent2="${INPUT_DIR}"/HG004.chr20.10_10p1mb.bam \
  --output_vcf_child "${OUTPUT_DIR}"/HG002.output.vcf.gz \
  --output_vcf_parent1 "${OUTPUT_DIR}"/HG003.output.vcf.gz \
  --output_vcf_parent2 "${OUTPUT_DIR}"/HG004.output.vcf.gz \
  --sample_name_child 'HG002' \
  --sample_name_parent1 'HG003' \
  --sample_name_parent2 'HG004' \
  --num_shards $(nproc)  \
  --regions "chr20:10,000,000-10,010,000" \
  --intermediate_results_dir "${OUTPUT_DIR}"/intermediate_results_dir \
  --output_gvcf_child "${OUTPUT_DIR}"/HG002.g.vcf.gz \
  --output_gvcf_parent1 "${OUTPUT_DIR}"/HG003.g.vcf.gz \
  --output_gvcf_parent2 "${OUTPUT_DIR}"/HG004.g.vcf.gz
```

### GPU version

```
# Pull the image.
singularity pull docker://google/deepvariant:deeptrio-"${BIN_VERSION}-gpu"

# Run DeepTrio.
# Using "--nv" and "${BIN_VERSION}-gpu" is important.
singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ \
  docker://google/deepvariant:deeptrio-"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  ...
```

## Evaluating the results

Here we use the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting three VCF files (only covering
10 kb of chr20 for this small example). Here the DeepTrio output VCF for each
sample is compared against the corresponding truth set from GIAB.

```bash
# Pull docker image locally
sudo docker pull jmcdani20/hap.py:v0.3.12

# Evaluate HG002
sudo docker run -it \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /input/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG002.output.vcf.gz \
  -f "/input/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
  -r "/input/GRCh38_no_alt_analysis_set.fasta" \
  -o "/output/HG002.happy" \
  --engine=vcfeval \
  --pass-only \
  -l chr20:10000000-10010000

# Evaluate HG003
sudo docker run -it \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /input/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG003.output.vcf.gz \
  -f "/input/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
  -r "/input/GRCh38_no_alt_analysis_set.fasta" \
  -o "/output/HG003.happy" \
  --engine=vcfeval \
  --pass-only \
  -l chr20:10000000-10010000

# Evaluate HG004
sudo docker run -it \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /input/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG004.output.vcf.gz \
  -f "/input/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
  -r "/input/GRCh38_no_alt_analysis_set.fasta" \
  -o "/output/HG004.happy" \
  --engine=vcfeval \
  --pass-only \
  -l chr20:10000000-10010000
```

You should see output similar to the following.

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL            3         3         0            4         0          0      0      0            1.0               1.0             0.0              1.0                     NaN                     NaN                   1.000000                   3.000000
INDEL   PASS            3         3         0            4         0          0      0      0            1.0               1.0             0.0              1.0                     NaN                     NaN                   1.000000                   3.000000
  SNP    ALL           19        19         0           19         0          0      0      0            1.0               1.0             0.0              1.0                5.333333                5.333333                   2.166667                   2.166667
  SNP   PASS           19        19         0           19         0          0      0      0            1.0               1.0             0.0              1.0                5.333333                5.333333                   2.166667                   2.166667

Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL            3         3         0            3         0          0      0      0            1.0               1.0             0.0              1.0                     NaN                     NaN                   0.500000                   0.500000
INDEL   PASS            3         3         0            3         0          0      0      0            1.0               1.0             0.0              1.0                     NaN                     NaN                   0.500000                   0.500000
  SNP    ALL           21        21         0           21         0          0      0      0            1.0               1.0             0.0              1.0                     6.0                     6.0                   1.333333                   1.333333
  SNP   PASS           21        21         0           21         0          0      0      0            1.0               1.0             0.0              1.0                     6.0                     6.0                   1.333333                   1.333333

Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL            2         2         0            2         0          0      0      0            1.0               1.0             0.0              1.0                     NaN                     NaN                   1.000000                   1.000000
INDEL   PASS            2         2         0            2         0          0      0      0            1.0               1.0             0.0              1.0                     NaN                     NaN                   1.000000                   1.000000
  SNP    ALL           10        10         0           10         0          0      0      0            1.0               1.0             0.0              1.0                     4.0                     4.0                   2.333333                   2.333333
  SNP   PASS           10        10         0           10         0          0      0      0            1.0               1.0             0.0              1.0                     4.0                     4.0                   2.333333                   2.333333
```

[BAM]: http://genome.sph.umich.edu/wiki/BAM
[BWA]: https://academic.oup.com/bioinformatics/article/25/14/1754/225615/Fast-and-accurate-short-read-alignment-with
[docker build]: https://docs.docker.com/engine/reference/commandline/build/
[Dockerfile]: https://github.com/google/deepvariant/blob/r1.6/Dockerfile.deeptrio
[FASTA]: https://en.wikipedia.org/wiki/FASTA_format
[VCF]: https://samtools.github.io/hts-specs/VCFv4.3.pdf
[run_deeptrio.py]: ../scripts/run_deeptrio.py
[install_nvidia_docker.sh]: ../scripts/install_nvidia_docker.sh
