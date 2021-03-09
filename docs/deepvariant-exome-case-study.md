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

We will benchmark our variant calls against v4.2 of the Genome in a Bottle small
variant benchmarks for HG003.

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020

curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG003_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi
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

BIN_VERSION="1.1.0"

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
  /benchmark/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz \
  /output/HG003.output.vcf.gz \
  -f /benchmark/HG003_GRCh38_1_22_v4.2_benchmark.bed \
  -T /input/idt_capture_novogene.grch38.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval
```

Output:

```
Benchmarking Summary:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL         1053      1026        27         1519        23        450     10       0.974359          0.978485        0.296248         0.976417                     NaN                     NaN                   1.752717                   1.906796
 INDEL   PASS         1053      1026        27         1519        23        450     10       0.974359          0.978485        0.296248         0.976417                     NaN                     NaN                   1.752717                   1.906796
   SNP    ALL        25324     25007       317        27947       167       2771     33       0.987482          0.993367        0.099152         0.990416                2.856273                2.759919                   1.625246                   1.665680
   SNP   PASS        25324     25007       317        27947       167       2771     33       0.987482          0.993367        0.099152         0.990416                2.856273                2.759919                   1.625246                   1.665680
```

[case study on whole genome sequencing data]: deepvariant-case-study.md
