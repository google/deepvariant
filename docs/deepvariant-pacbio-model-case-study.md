# Using DeepVariant for small variant calling from PacBio HiFi reads

#### Author: William Rowell <wrowell@pacificbiosciences.com>

In this case study we describe applying DeepVariant to PacBio HiFi reads to call
variants. We will call small variants from a publicly available whole genome
HiFi dataset from PacBio.

This case study involves a two-step process of variant calling. After the first
round of calling, SNVs are phased and used to haplotag the input BAM. For
the highest accuracy, variants are called again in a second pass. If somewhat
lower Indel accuracy is acceptable in exchange for shorter run-time, the
first-pass calls be used. Accuracy benchmarks for each pass are shown for this
case study.


## Prepare environment

### Tools

[Singularity](https://sylabs.io/docs/) will be used to run DeepVariant and
[hap.py](https://github.com/illumina/hap.py), and we'll use
[miniconda](https://docs.conda.io/en/latest/miniconda.html) and a conda
environment to handle the other dependencies for the case study, samtools and
whatshap.

-   singularity (must be installed by `root` user; outside of the scope of this
    case study)
-   samtools
-   whatshap

```bash
# add channels to conda configuration
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the environment and install dependencies
conda create -y -n deepvariant_whatshap
conda activate deepvariant_whatshap
conda install -y whatshap==1.0 samtools==1.10
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
BIN_VERSION="1.1.0"
mkdir -p deepvariant1

singularity exec --bind /usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads input/HG003.GRCh38.chr20.pFDA_truthv2.bam \
    --output_vcf deepvariant1/output.vcf.gz \
    --num_shards $(nproc) \
    --regions chr20
```

## Phase SNPs on chromosome 20

```bash
mkdir -p whatshap

whatshap phase \
        --output whatshap/deepvariant1.phased.vcf.gz \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        --chromosome chr20 \
        deepvariant1/output.vcf.gz \
        input/HG003.GRCh38.chr20.pFDA_truthv2.bam

tabix -p vcf whatshap/deepvariant1.phased.vcf.gz
```

## Haplotag chromosome 20 alignments

```bash
whatshap haplotag \
        --output whatshap/HG003.GRCh38.chr20.haplotagged.bam \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        whatshap/deepvariant1.phased.vcf.gz \
        input/HG003.GRCh38.chr20.pFDA_truthv2.bam

samtools index whatshap/HG003.GRCh38.chr20.haplotagged.bam
```

## Run DeepVariant on haplotagged chromosome 20 alignments

```bash
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
BIN_VERSION="1.1.0"
mkdir -p deepvariant2

singularity exec --bind /usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads whatshap/HG003.GRCh38.chr20.haplotagged.bam \
    --use_hp_information \
    --output_vcf deepvariant2/output.vcf.gz \
    --num_shards $(nproc) \
    --regions chr20
```

## Clarification of the `--use_hp_information` flag

In order to allow DeepVariant to take advantage of the haplotype (HP) tags in
haplotagged BAMs, three flags must be passed to
`/opt/deepvariant/bin/make_examples`:

```bash
--sort_by_haplotypes --parse_sam_aux_fields --add_hp_channel
```

For the convenient one-step command `/opt/deepvariant/bin/run_deepvariant`, the
`--use_hp_information` flag will provide these three flags to make_examples.

## Benchmark First Pass

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
        deepvariant1/output.vcf.gz
```

First pass output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10493       135        21945       131      10889     77     30       0.987298          0.988151        0.496195         0.987724                     NaN                     NaN                   1.748961                   2.376968
INDEL   PASS        10628     10493       135        21945       131      10889     77     30       0.987298          0.988151        0.496195         0.987724                     NaN                     NaN                   1.748961                   2.376968
  SNP    ALL        70166     70144        22        96093        38      25845      7      1       0.999686          0.999459        0.268958         0.999573                2.296566                1.961428                   1.883951                   2.115808
  SNP   PASS        70166     70144        22        96093        38      25845      7      1       0.999686          0.999459        0.268958         0.999573                2.296566                1.961428                   1.883951                   2.115808
```

## Benchmark Second Pass

```bash
mkdir -p happy

singularity exec docker://jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
        --threads $(nproc) \
        -r reference/GRCh38_no_alt_analysis_set.fasta \
        -f benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        -o happy/giab-comparison.v4.2.second_pass \
        --engine=vcfeval \
        --pass-only \
        -l chr20 \
        benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        deepvariant2/output.vcf.gz
```

Second pass output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10556        72        22270        78      11198     48     23       0.993225          0.992955        0.502829         0.993090                     NaN                     NaN                   1.748961                   2.452495
INDEL   PASS        10628     10556        72        22270        78      11198     48     23       0.993225          0.992955        0.502829         0.993090                     NaN                     NaN                   1.748961                   2.452495
  SNP    ALL        70166     70144        22        96217        38      25971      6      3       0.999686          0.999459        0.269921         0.999573                2.296566                1.952486                   1.883951                   2.068213
  SNP   PASS        70166     70144        22        96217        38      25971      6      3       0.999686          0.999459        0.269921         0.999573                2.296566                1.952486                   1.883951                   2.068213
```
