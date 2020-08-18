# DeepVariant + PacBio local Quickstart

## Prepare environment

### Tools

- singularity
- samtools
- whatshap

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n deepvariant_whatshap
conda activate deepvariant_whatshap

conda install singularity==3.5.3 whatshap==1.0 samtools==1.10
```

### Download Reference

```bash
mkdir -p reference

curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
samtools faidx reference/GRCh38_no_alt_analysis_set.fasta
```

### Download Genome in a Bottle Benchmarks

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020

curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi
```

### Download HG002 chr20 HiFi alignments

```bash
mkdir -p input

FTPDIR=http://storage.googleapis.com/deepvariant/pacbio-case-study-testdata

curl ${FTPDIR}/HG002.GRCh38.chr20.bam > input/HG002.GRCh38.chr20.bam
curl ${FTPDIR}/HG002.GRCh38.chr20.bam.bai > input/HG002.GRCh38.chr20.bam.bai
```

## Run DeepVariant on chromosome 20 alignments

```bash
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
BIN_VERSION=rc1.0.0
mkdir -p deepvariant1

singularity exec --bind /usr/lib/locale/ \
  docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads input/HG002.GRCh38.chr20.bam \
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
        input/HG002.GRCh38.chr20.bam

tabix -p vcf whatshap/deepvariant1.phased.vcf.gz
```

## Haplotag chromosome 20 alignments

```bash
whatshap haplotag \
        --output whatshap/HG002.GRCh38.chr20.haplotagged.bam \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        whatshap/deepvariant1.phased.vcf.gz \
        input/HG002.GRCh38.chr20.bam

samtools index whatshap/HG002.GRCh38.chr20.haplotagged.bam
```

## Run DeepVariant on haplotagged chromosome 20 alignments

```bash
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
BIN_VERSION=rc1.0.0
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
 INDEL    ALL        11256     11068       188        21670       139      10059     90       0.983298          0.988029        0.464190         0.985658                     NaN                     NaN                   1.561710                   2.079747
 INDEL   PASS        11256     11068       188        21670       139      10059     90       0.983298          0.988029        0.464190         0.985658                     NaN                     NaN                   1.561710                   2.079747
   SNP    ALL        71333     71279        54        95192         5      23828      5       0.999243          0.999930        0.250315         0.999586                2.314904                2.011701                   1.715978                   2.020670
   SNP   PASS        71333     71279        54        95192         5      23828      5       0.999243          0.999930        0.250315         0.999586                2.314904                2.011701                   1.715978                   2.020670
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
 INDEL    ALL        11256     11161        95        22092        99      10424     62       0.991560          0.991515        0.471845         0.991538                     NaN                     NaN                   1.561710                   2.223964
 INDEL   PASS        11256     11161        95        22092        99      10424     62       0.991560          0.991515        0.471845         0.991538                     NaN                     NaN                   1.561710                   2.223964
   SNP    ALL        71333     71280        53        95017         5      23651      5       0.999257          0.999930        0.248913         0.999593                2.314904                2.010292                   1.715978                   2.005063
   SNP   PASS        71333     71280        53        95017         5      23651      5       0.999257          0.999930        0.248913         0.999593                2.314904                2.010292                   1.715978                   2.005063
```
