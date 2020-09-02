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
conda create -n deepvariant_whatshap
conda activate deepvariant_whatshap
conda install whatshap==1.0 samtools==1.10
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

We will benchmark our variant calls against v4.2 of the Genome in a Bottle small
variant benchmarks for HG002.

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020

curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi
```

### Download HG002 chr20 HiFi alignments

We'll use HG002 chr20 HiFi reads publicly available from the [PrecisionFDA Truth v2 Challenge](https://precision.fda.gov/challenges/10).

```bash
mkdir -p input
HTTPDIR=https://downloads.pacbcloud.com/public/dataset/HG002/deepvariant-case-study

curl ${HTTPDIR}/HG002.GRCh38.chr20.pFDA_truthv2.bam > input/HG002.GRCh38.chr20.pFDA_truthv2.bam
curl ${HTTPDIR}/HG002.GRCh38.chr20.pFDA_truthv2.bam.bai > input/HG002.GRCh38.chr20.pFDA_truthv2.bam.bai
```

## Run DeepVariant on chromosome 20 alignments

```bash
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
BIN_VERSION=1.0.0
mkdir -p deepvariant1

singularity exec --bind /usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads input/HG002.GRCh38.chr20.pFDA_truthv2.bam \
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
        input/HG002.GRCh38.chr20.pFDA_truthv2.bam

tabix -p vcf whatshap/deepvariant1.phased.vcf.gz
```

## Haplotag chromosome 20 alignments

```bash
whatshap haplotag \
        --output whatshap/HG002.GRCh38.chr20.haplotagged.bam \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        whatshap/deepvariant1.phased.vcf.gz \
        input/HG002.GRCh38.chr20.pFDA_truthv2.bam

samtools index whatshap/HG002.GRCh38.chr20.haplotagged.bam
```

## Run DeepVariant on haplotagged chromosome 20 alignments

```bash
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
BIN_VERSION=1.0.0
mkdir -p deepvariant2

singularity exec --bind /usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads whatshap/HG002.GRCh38.chr20.haplotagged.bam \
    --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true" \
    --output_vcf deepvariant2/output.vcf.gz \
    --num_shards $(nproc) \
    --regions chr20
```

## Benchmark First Pass

```bash
mkdir -p happy

singularity exec docker://pkrusche/hap.py:latest \
    /opt/hap.py/bin/hap.py \
        --threads $(nproc) \
        -r reference/GRCh38_no_alt_analysis_set.fasta \
        -f benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed \
        -o happy/giab-comparison.v4.2.first_pass \
        --engine=vcfeval \
        -l chr20 \
        benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz \
        deepvariant1/output.vcf.gz
```

First pass output:

```
Benchmarking Summary:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL        11256     11048       208        21583       150       9985     97       0.981521          0.987067        0.462633         0.984286                     NaN                     NaN                   1.561710                   2.063683
 INDEL   PASS        11256     11048       208        21583       150       9985     97       0.981521          0.987067        0.462633         0.984286                     NaN                     NaN                   1.561710                   2.063683
   SNP    ALL        71333     71277        56        95048         7      23684      6       0.999215          0.999902        0.249179         0.999558                2.314904                2.018283                   1.715978                   2.017439
   SNP   PASS        71333     71277        56        95048         7      23684      6       0.999215          0.999902        0.249179         0.999558                2.314904                2.018283                   1.715978                   2.017439
```

## Benchmark Second Pass

```bash
mkdir -p happy

singularity exec docker://pkrusche/hap.py:latest \
    /opt/hap.py/bin/hap.py \
        --threads $(nproc) \
        -r reference/GRCh38_no_alt_analysis_set.fasta \
        -f benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed \
        -o happy/giab-comparison.v4.2.second_pass \
        --engine=vcfeval \
        -l chr20 \
        benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz \
        deepvariant2/output.vcf.gz
```

Second pass output:

```
Benchmarking Summary:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL        11256     11146       110        22043       116      10374     66       0.990227          0.990059        0.470626         0.990143                     NaN                     NaN                   1.561710                   2.225178
 INDEL   PASS        11256     11146       110        22043       116      10374     66       0.990227          0.990059        0.470626         0.990143                     NaN                     NaN                   1.561710                   2.225178
   SNP    ALL        71333     71274        59        94898        10      23533      8       0.999173          0.999860        0.247982         0.999516                2.314904                2.016874                   1.715978                   2.001392
   SNP   PASS        71333     71274        59        94898        10      23533      8       0.999173          0.999860        0.247982         0.999516                2.314904                2.016874                   1.715978                   2.001392
```
