# Pangenome-aware DeepVariant Complete Genomics T7+ case study

In this case study, we describe applying Pangenome-aware DeepVariant to a
Complete Genomics T7+ sample. Then we assess the quality of the DeepVariant
variant calls with `hap.py`.

To make it faster to run over this case study, we run only on chromosome 20.

For how to prepare environment, the steps are the same as
[this doc](deepvariant-case-study.md).

## Download Reference

We will be using GRCh38 for this case study.

```bash
mkdir -p reference

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai
```

## Download Complete Genomics T7+ HG002 chr20 BAM

The original FASTQs are from:
https://www.completegenomics.com/demo-data/wgs-dnbseq-t7-plus-pe150-pcr-free/

```bash
mkdir -p input

gcloud storage cp \
  gs://deepvariant/complete-case-study-testdata/complete-t7plus/bam/T7plus_WGS_HG002.chr20.bam* \
  input/
```

## Download Genome in a Bottle v5.0q Benchmarks for HG002

```bash
mkdir -p benchmark

gcloud storage cp \
  gs://deepvariant/GIAB_v5q0/HG002_GRCh38_v5.0q_smvar.vcf.gz* \
  gs://deepvariant/GIAB_v5q0/HG002_GRCh38_v5.0q_smvar.benchmark.bed \
  benchmark/
```

### Download GBZ built for GRCh38

```bash
HTTPDIR=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38

curl ${HTTPDIR}/hprc-v1.1-mc-grch38.gbz > input/hprc-v1.1-mc-grch38.gbz
```

## Download Complete Genomics T7+ model

```bash
mkdir -p input/model

gcloud storage cp \
  gs://deepvariant/complete-case-study-testdata/complete-t7plus/2026/* \
  input/model/
```

## Running DeepVariant with one command

On a CPU-only machine:

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="pangenome_aware_deepvariant-1.10.0"

docker pull google/deepvariant:"${BIN_VERSION}"

docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  --shm-size 12gb \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
  --model_type WGS \
  --ref "/reference/GRCh38_no_alt_analysis_set.fasta" \
  --reads "/input/T7plus_WGS_HG002.chr20.bam" \
  --pangenome "/input/hprc-v1.1-mc-grch38.gbz" \
  --output_vcf "/output/HG002.chr20.output.vcf.gz" \
  --output_gvcf "/output/HG002.chr20.output.g.vcf.gz" \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir \
  --customized_model "/input/model/checkpoint-76800-0.99361-1"
```

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md).

## Benchmark on chr20

```bash
mkdir -p happy

docker pull jmcdani20/hap.py:v0.3.12

docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG002_GRCh38_v5.0q_smvar.vcf.gz \
  /output/HG002.chr20.output.vcf.gz \
  -f /benchmark/HG002_GRCh38_v5.0q_smvar.benchmark.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval \
  --pass-only \
  -l chr20
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        18252     17903       349        22254       231       2601     84    128       0.980879          0.988246        0.116878         0.984549                     NaN                     NaN                   1.699834                   2.102103
INDEL   PASS        18252     17903       349        22254       231       2601     84    128       0.980879          0.988246        0.116878         0.984549                     NaN                     NaN                   1.699834                   2.102103
  SNP    ALL        75944     75666       278        93335       143      17896     52     79       0.996339          0.998104        0.191739         0.997221                2.236314                1.954288                   1.722736                   1.768279
  SNP   PASS        75944     75666       278        93335       143      17896     52     79       0.996339          0.998104        0.191739         0.997221                2.236314                1.954288                   1.722736                   1.768279
```
