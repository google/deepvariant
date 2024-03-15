# DeepVariant whole exome sequencing (WES) case study

Similar to the [case study on whole genome sequencing data], in this
study we describe applying DeepVariant to a real exome sample using a single
machine.

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

### Download HG003 BAM

```bash
mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/exome-case-study-testdata

curl ${HTTPDIR}/HG003.novaseq.wes_idt.100x.dedup.bam > input/HG003.novaseq.wes_idt.100x.dedup.bam
curl ${HTTPDIR}/HG003.novaseq.wes_idt.100x.dedup.bam.bai > input/HG003.novaseq.wes_idt.100x.dedup.bam.bai
```

### Download capture target BED file

In this case study we'll use `idt_capture_novogene.grch38.bed` as the capture
target BED file. For evaluation, `hap.py` will intersect this BED with the GIAB
confident regions.

```bash
HTTPDIR=https://storage.googleapis.com/deepvariant/exome-case-study-testdata

curl ${HTTPDIR}/idt_capture_novogene.grch38.bed > input/idt_capture_novogene.grch38.bed
```


## Running on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="1.6.1"

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG003.novaseq.wes_idt.100x.dedup.bam \
  --regions /input/idt_capture_novogene.grch38.bed \
  --output_vcf /output/HG003.output.vcf.gz \
  --output_gvcf /output/HG003.output.g.vcf.gz \
  --num_shards $(nproc) \
  --intermediate_results_dir /output/intermediate_results_dir
```

By specifying `--model_type WES`, you'll be using a model that is best suited
for Illumina Whole Exome Sequencing data.

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

## Benchmark on all chromosomes

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
  -T /input/idt_capture_novogene.grch38.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval \
  --pass-only
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL         1051      1022        29         1476        13        418      8      3       0.972407          0.987713        0.283198         0.980000                     NaN                     NaN                   1.747283                   1.859406
INDEL   PASS         1051      1022        29         1476        13        418      8      3       0.972407          0.987713        0.283198         0.980000                     NaN                     NaN                   1.747283                   1.859406
  SNP    ALL        25279     24987       292        27710        59       2662     34      2       0.988449          0.997645        0.096066         0.993025                2.854703                2.749729                   1.623027                   1.636078
  SNP   PASS        25279     24987       292        27710        59       2662     34      2       0.988449          0.997645        0.096066         0.993025                2.854703                2.749729                   1.623027                   1.636078
```

[case study on whole genome sequencing data]: deepvariant-case-study.md
