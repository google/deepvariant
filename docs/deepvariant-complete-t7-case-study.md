# DeepVariant Complete Genomics T7 case study

In this case study, we describe applying DeepVariant to a Complete Genomics T7
sample.
Then we assess the quality of the DeepVariant variant calls with `hap.py`.

To make it faster to run over this case study, we run only on chromosome 20.

For how to prepare environment, the steps are the same as
[this doc](deepvariant-case-study.md).


## Download Complete Genomics T7 HG001 chr20 BAM

```bash
mkdir -p input

HTTPDIR=https://storage.googleapis.com/deepvariant/complete-case-study-testdata

curl ${HTTPDIR}/HG001.complete_t7.E100030471QC960.grch38.chr20.bam > input/HG001.complete_t7.E100030471QC960.grch38.chr20.bam

curl ${HTTPDIR}/HG001.complete_t7.E100030471QC960.grch38.chr20.bam.bai > input/HG001.complete_t7.E100030471QC960.grch38.chr20.bam.bai
```

## Download Genome in a Bottle Benchmarks for HG001

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38

curl ${FTPDIR}/HG001_GRCh38_1_22_v4.2.1_benchmark.bed > benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
curl ${FTPDIR}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

## Download Complete Genomics T7 model

```bash
HTTPDIR=https://storage.googleapis.com/deepvariant/complete-case-study-testdata

curl ${HTTPDIR}/complete-t7/weights-51-0.995354.ckpt.data-00000-of-00001 > input/weights-51-0.995354.ckpt.data-00000-of-00001

curl ${HTTPDIR}/complete-t7/weights-51-0.995354.ckpt.index > input/weights-51-0.995354.ckpt.index
```

## Running DeepVariant with one command

On a CPU-only machine:

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="1.7.0"

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG001.complete_t7.E100030471QC960.grch38.chr20.bam \
  --output_vcf /output/HG001.output.vcf.gz \
  --output_gvcf /output/HG001.output.g.vcf.gz \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir \
  --customized_model /input/weights-51-0.995354.ckpt
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
  /benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG001.output.vcf.gz \
  -f /benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
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
INDEL    ALL         9974      9947        27        21052         9      10750      3      5       0.997293          0.999126        0.510640         0.998209                     NaN                     NaN                   1.630447                   2.156149
INDEL   PASS         9974      9947        27        21052         9      10750      3      5       0.997293          0.999126        0.510640         0.998209                     NaN                     NaN                   1.630447                   2.156149
  SNP    ALL        69175     68874       301        85030        44      16068      8      2       0.995649          0.999362        0.188969         0.997502                2.288757                2.084645                   1.730097                   1.781789
  SNP   PASS        69175     68874       301        85030        44      16068      8      2       0.995649          0.999362        0.188969         0.997502                2.288757                2.084645                   1.730097                   1.781789
```

To summarize:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 9947     | 27       | 9        | 0.997293      | 0.999126         | 0.998209        |
| SNP   | 68874    | 301      | 44       | 0.995649      | 0.999362         | 0.997502        |
