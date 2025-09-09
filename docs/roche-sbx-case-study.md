# Roche SBX case study

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

### Download GIAB v4.2.1 truth for benchmarking

We will benchmark our variant calls against GIAB v4.2.1 truth for HG002.

```bash
mkdir -p benchmark

HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata

curl ${HTTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${HTTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl ${HTTPDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

### Download GBZ built for GRCh38

```bash
mkdir -p input
HTTPDIR=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38

curl ${HTTPDIR}/hprc-v1.1-mc-grch38.gbz > input/hprc-v1.1-mc-grch38.gbz
```

### Download HG002 chr20 BAM

We will use Roche SBX HG002 chr20 BAM for this case-study.

```bash
mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/roche-sbx-case-study-testdata

curl ${HTTPDIR}/HG002.roche_sbx.chr20.bam > input/HG002.roche_sbx.chr20.bam
curl ${HTTPDIR}/HG002.roche_sbx.chr20.bam.bai > input/HG002.roche_sbx.chr20.bam.bai
```

### Download the model

In this case study, we're calling variants on HG002 chr20.

```bash
mkdir -p model

HTTPDIR=https://storage.googleapis.com/brain-genomics-public/research/sbx/2025/model/leave-out-HG001

curl ${HTTPDIR}/model.ckpt.data-00000-of-00001 > model/model.ckpt.data-00000-of-00001
curl ${HTTPDIR}/model.ckpt.index > model/model.ckpt.index
curl ${HTTPDIR}/example_info.json > model/example_info.json
```

## Running DeepVariant with one command, on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="pangenome_aware_deepvariant-sbx"

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/model":"/model" \
  --shm-size 12gb \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
  --model_type WGS \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG002.roche_sbx.chr20.bam \
  --pangenome /input/hprc-v1.1-mc-grch38.gbz \
  --output_vcf /output/HG002.chr20.output.vcf.gz \
  --output_gvcf /output/HG002.chr20.output.g.vcf.gz \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir \
  --customized_model /model/model.ckpt \
  --make_examples_extra_args="alt_aligned_pileup=single_row,create_complex_alleles=true,enable_strict_insertion_filter=true,keep_legacy_allele_counter_behavior=true,keep_only_window_spanning_haplotypes=true,keep_supplementary_alignments=true,min_mapping_quality=0,normalize_reads=true,pileup_image_height_pangenome=100,pileup_image_height_reads=100,pileup_image_width=301,sort_by_haplotypes=true,trim_reads_for_pileup=true,vsc_min_fraction_indels=0.08,ws_min_base_quality=25" \
  --postprocess_variants_extra_args="multiallelic_mode=product"
```

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
  /benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG002.chr20.output.vcf.gz \
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
INDEL    ALL        11256     11237        19        22167        22      10474     10     10       0.998312          0.998119        0.472504         0.998215                     NaN                     NaN                   1.561710                   2.089132
INDEL   PASS        11256     11237        19        22167        22      10474     10     10       0.998312          0.998119        0.472504         0.998215                     NaN                     NaN                   1.561710                   2.089132
  SNP    ALL        71333     71286        47        91930        41      20553     12      3       0.999341          0.999426        0.223572         0.999383                2.314904                1.943955                   1.715978                   1.640709
  SNP   PASS        71333     71286        47        91930        41      20553     12      3       0.999341          0.999426        0.223572         0.999383                2.314904                1.943955                   1.715978                   1.640709
```

For all Roche SBX 30x BAMs and VCFs follow this [link](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/sbx/2025/).
