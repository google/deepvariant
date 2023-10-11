# DeepVariant whole genome sequencing case study

In this case study, we describe applying DeepVariant to a real WGS sample.
Then we assess the quality of the DeepVariant variant calls with `hap.py`.

To make it faster to run over this case study, we run only on chromosome 20.


## Prepare environment

### Tools

[Docker](https://docs.docker.com/get-docker/) will be used to run DeepVariant
and [hap.py](https://github.com/illumina/hap.py),

### Download Reference

We will be using GRCh38 for this case study.

```bash
mkdir -p reference

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai
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

### Download HG003 chr20 BAM

We'll use HG003 Illumina WGS reads publicly available from the
[PrecisionFDA Truth v2 Challenge](https://doi.org/10.1101/2020.11.13.380741).

```bash
mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata

curl ${HTTPDIR}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam > input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
curl ${HTTPDIR}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai > input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai
```


## Running DeepVariant with one command

DeepVariant pipeline consists of 3 steps: `make_examples`, `call_variants`, and
`postprocess_variants`. You can now run DeepVariant with one command using the
`run_deepvariant` script.

### Running on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="1.6.0"

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
  --output_vcf /output/HG003.output.vcf.gz \
  --output_gvcf /output/HG003.output.g.vcf.gz \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir
```

By specifying `--model_type WGS`, you'll be using a model that is best suited
for Illumina Whole Genome Sequencing data.

NOTE: If you want to run each of the steps separately, add `--dry_run=true`
to the command above to figure out what flags you need in each step. Based on
the different model types, different flags are needed in the `make_examples`
step.

`--intermediate_results_dir` flag is optional. By specifying it, the
intermediate outputs of `make_examples` and `call_variants` stages can be found
in the directory. After the command, you can find these files in the directory:

```
call_variants_output.tfrecord.gz
gvcf.tfrecord-?????-of-?????.gz
make_examples.tfrecord-?????-of-?????.gz
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
  /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG003.output.vcf.gz \
  -f /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
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
INDEL    ALL        10628     10588        40        21099        19      10036     15      3       0.996236          0.998283        0.475662         0.997258                     NaN                     NaN                   1.748961                   2.318182
INDEL   PASS        10628     10588        40        21099        19      10036     15      3       0.996236          0.998283        0.475662         0.997258                     NaN                     NaN                   1.748961                   2.318182
  SNP    ALL        70166     69917       249        84796        59      14782     13      3       0.996451          0.999157        0.174324         0.997802                2.296566                2.085786                   1.883951                   1.920577
  SNP   PASS        70166     69917       249        84796        59      14782     13      3       0.996451          0.999157        0.174324         0.997802                2.296566                2.085786                   1.883951                   1.920577
```
