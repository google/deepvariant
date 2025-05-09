# DeepVariant Complete Genomics G400 case study

In this case study, we describe applying DeepVariant to a Complete Genomics G400
sample.
Then we assess the quality of the DeepVariant variant calls with `hap.py`.

To make it faster to run over this case study, we run only on chromosome 20.

For how to prepare environment, the steps are the same as
[this doc](deepvariant-case-study.md).


## Download Complete Genomics G400 HG002 chr20 BAM

```bash
mkdir -p input

HTTPDIR=https://storage.googleapis.com/deepvariant/complete-case-study-testdata

curl ${HTTPDIR}/HG002.complete_g400.V350151728.grch38.chr20.bam > input/HG002.complete_g400.V350151728.grch38.chr20.bam

curl ${HTTPDIR}/HG002.complete_g400.V350151728.grch38.chr20.bam.bai > input/HG002.complete_g400.V350151728.grch38.chr20.bam.bai
```

## Download Genome in a Bottle Benchmarks for HG002

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38

curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

## Download Complete Genomics G400 model

```bash
HTTPDIR=https://storage.googleapis.com/deepvariant/complete-case-study-testdata

curl ${HTTPDIR}/complete-g400/weights-60-0.993753.ckpt.data-00000-of-00001 > input/weights-60-0.993753.ckpt.data-00000-of-00001

curl ${HTTPDIR}/complete-g400/weights-60-0.993753.ckpt.index > input/weights-60-0.993753.ckpt.index

curl ${HTTPDIR}/complete-g400/example_info.json > input/example_info.json
```

## Running DeepVariant with one command

On a CPU-only machine:

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="1.9.0"

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG002.complete_g400.V350151728.grch38.chr20.bam \
  --output_vcf /output/HG002.output.vcf.gz \
  --output_gvcf /output/HG002.output.g.vcf.gz \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir \
  --customized_model /input/weights-60-0.993753.ckpt
```

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md).

## Benchmark on chr20

```bash
mkdir -p happy

sudo docker pull jmcdani20/hap.py:v0.3.12

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG002.output.vcf.gz \
  -f /benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
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
INDEL    ALL        11256     11129       127        20905        30       9322     25      4       0.988717          0.997410        0.445922         0.993045                     NaN                     NaN                   1.561710                   2.053139
INDEL   PASS        11256     11129       127        20905        30       9322     25      4       0.988717          0.997410        0.445922         0.993045                     NaN                     NaN                   1.561710                   2.053139
  SNP    ALL        71333     70954       379        85776        52      14722     28      8       0.994687          0.999268        0.171633         0.996972                2.314904                2.098765                   1.715978                   1.753260
  SNP   PASS        71333     70954       379        85776        52      14722     28      8       0.994687          0.999268        0.171633         0.996972                2.314904                2.098765                   1.715978                   1.753260
```

To summarize:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 11129    | 127      | 30       | 0.988717      | 0.997410         | 0.993045        |
| SNP   | 70954    | 379      | 52       | 0.994687      | 0.999268         | 0.996972        |
