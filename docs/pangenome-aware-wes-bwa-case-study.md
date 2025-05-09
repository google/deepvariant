# DeepVariant Pangenome-aware WES case study (mapped with BWA)

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

We'll use HG003 Illumina WES reads.

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

### Download GBZ built for GRCh38

```bash
mkdir -p input
HTTPDIR=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38

curl ${HTTPDIR}/hprc-v1.1-mc-grch38.gbz > input/hprc-v1.1-mc-grch38.gbz
```

## Running Pangenome-aware DeepVariant with one command

DeepVariant pipeline consists of 3 steps: `make_examples`, `call_variants`, and
`postprocess_variants`. You can now run DeepVariant with one command using the
`run_pangenome_aware_deepvariant` script.

### Running on a CPU-only machine

In this example, we used a
[n2-standard-96](https://cloud.google.com/compute/docs/general-purpose-machines)
machine.

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="pangenome_aware_deepvariant-1.9.0"

sudo docker pull google/deepvariant:"${BIN_VERSION}"

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  --shm-size 12gb \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
  --model_type WES \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG003.novaseq.wes_idt.100x.dedup.bam  \
  --pangenome /input/hprc-v1.1-mc-grch38.gbz \
  --output_vcf /output/HG003.output.vcf.gz \
  --output_gvcf /output/HG003.output.g.vcf.gz \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir
```


By specifying `--model_type WES`, you'll be using a model that is best
suited for short-read WES data.

NOTE: If you want to run each of the steps separately, add `--dry_run=true`
to the command above to figure out what flags you need in each step.

`--intermediate_results_dir` flag is optional. By specifying it, the
intermediate outputs can be found in the directory. After the command, you can
find these intermediate files in the directory:

```
call_variants_output-?????-of-?????.tfrecord.gz
gvcf.tfrecord-?????-of-?????.gz
make_examples_pangenome_aware_dv.tfrecord-?????-of-?????.gz
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
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
  /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG003.output.vcf.gz \
  -f /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -T /input/idt_capture_novogene.grch38.bed \
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
INDEL    ALL           29        29         0           44         0         14      0      0        1.00000               1.0        0.318182         1.000000                     NaN                     NaN                   3.000000                   3.000000
INDEL   PASS           29        29         0           44         0         14      0      0        1.00000               1.0        0.318182         1.000000                     NaN                     NaN                   3.000000                   3.000000
  SNP    ALL          685       683         2          703         0         20      0      0        0.99708               1.0        0.028450         0.998538                 3.28125                3.260606                   1.795918                   1.834677
  SNP   PASS          685       683         2          703         0         20      0      0        0.99708               1.0        0.028450         0.998538                 3.28125                3.260606                   1.795918                   1.834677
```
