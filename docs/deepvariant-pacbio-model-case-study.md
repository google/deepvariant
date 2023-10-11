# Using DeepVariant for small variant calling from PacBio HiFi reads

#### Author: William Rowell <wrowell@pacificbiosciences.com>

In this case study we describe applying DeepVariant to PacBio HiFi reads to call
variants. We will call small variants from a publicly available whole genome
HiFi dataset from PacBio.

Starting in v1.4.0, PacBio calling uses one-step variant calling. If you're
looking for documentation for the two-step process, please look at v1.3.0.

## Prepare environment

### Tools

[Singularity](https://sylabs.io/docs/) will be used to run DeepVariant and
[hap.py](https://github.com/illumina/hap.py), and we'll use
[miniconda](https://docs.conda.io/en/latest/miniconda.html) and a conda
environment to handle the other dependencies for the case study and samtools.

-   singularity (must be installed by `root` user; outside of the scope of this
    case study)
-   samtools

```bash
# add channels to conda configuration
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the environment and install dependencies
conda create -y -n deepvariant_env
conda activate deepvariant_env
conda install -y samtools==1.10
```

### Download Reference

We will be using GRCh38 for this case study.

```bash
mkdir -p reference

# download and decompress
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta

# index reference
samtools faidx reference/GRCh38_no_alt_analysis_set.fasta
```

### Download Genome in a Bottle Benchmarks

We will benchmark our variant calls against v4.2.1 of the Genome in a Bottle
small variant benchmarks for HG003.

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38

curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

### Download HG003 chr20 HiFi alignments

We'll use HG003 chr20 HiFi reads publicly available from the [PrecisionFDA Truth v2 Challenge](https://precision.fda.gov/challenges/10).

```bash
mkdir -p input
HTTPDIR=https://downloads.pacbcloud.com/public/dataset/HG003/deepvariant-case-study

curl ${HTTPDIR}/HG003.GRCh38.chr20.pFDA_truthv2.bam > input/HG003.GRCh38.chr20.pFDA_truthv2.bam
curl ${HTTPDIR}/HG003.GRCh38.chr20.pFDA_truthv2.bam.bai > input/HG003.GRCh38.chr20.pFDA_truthv2.bam.bai
```

## Run DeepVariant on chromosome 20 alignments

```bash
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
BIN_VERSION="1.6.0"
mkdir -p deepvariant_output

singularity exec --bind /usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads input/HG003.GRCh38.chr20.pFDA_truthv2.bam \
    --output_vcf deepvariant_output/output.vcf.gz \
    --num_shards $(nproc) \
    --regions chr20
```

NOTE: If you want to run each of the steps separately, add `--dry_run=true`
to the command above to figure out what flags you need in each step. Based on
the different model types, different flags are needed in the `make_examples`
step.

## Benchmark output

```bash
mkdir -p happy

singularity exec docker://jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
        --threads $(nproc) \
        -r reference/GRCh38_no_alt_analysis_set.fasta \
        -f benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        -o happy/giab-comparison.v4.2.first_pass \
        --engine=vcfeval \
        --pass-only \
        -l chr20 \
        benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        deepvariant_output/output.vcf.gz
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10551        77        22590        69      11527     39     29       0.992755          0.993763        0.510270         0.993259                     NaN                     NaN                   1.748961                   2.275319
INDEL   PASS        10628     10551        77        22590        69      11527     39     29       0.992755          0.993763        0.510270         0.993259                     NaN                     NaN                   1.748961                   2.275319
  SNP    ALL        70166     70141        25        98780        23      28559      5     11       0.999644          0.999672        0.289117         0.999658                2.296566                1.823452                   1.883951                   1.913585
  SNP   PASS        70166     70141        25        98780        23      28559      5     11       0.999644          0.999672        0.289117         0.999658                2.296566                1.823452                   1.883951                   1.913585
```
