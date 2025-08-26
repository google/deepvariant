# SBX case study for SBX-D and SBX-Fast data

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

### Download T2T v1.1 truth for benchmarking

We will benchmark our variant calls against T2T v1.1 truth for HG002.

```bash
mkdir -p benchmark

HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata

curl ${HTTPDIR}/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz > benchmark/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz
curl ${HTTPDIR}/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz.tbi > benchmark/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz.tbi
curl ${HTTPDIR}/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed > benchmark/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed
```

### Download GBZ built for GRCh38

```bash
mkdir -p input
HTTPDIR=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38

curl ${HTTPDIR}/hprc-v1.1-mc-grch38.gbz > input/hprc-v1.1-mc-grch38.gbz
```

### Download HG002 BAM

Please download a SBX-D (or SBX-Fast) HG002 BAM and put it in your input dir.

```bash
# Download your HG002 BAM
HG002_BAM=/input/your_HG002.bam  # This will be used in the command later.
```

### Download the model

In this case study, we're calling variants on HG002 chr20, and we'll evaluate
with T2T v1.1 truth. We'll use the "leave-out-HG001" model, which we also left
out all chromosome 20 from training or tuning. Refer to the technical white
paper for more details on all experiments.

```bash
mkdir -p model

HTTPDIR=https://storage.googleapis.com/brain-genomics-public/research/sbx/2025/models/leave-out-HG001

curl ${HTTPDIR}/model.ckpt.data-00000-of-00001 > model/model.ckpt.data-00000-of-00001
curl ${HTTPDIR}/model.ckpt.index > model/model.ckpt.index
curl ${HTTPDIR}/example_info.json > model/example_info.json
```

## Running DeepVariant with one command, on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="pangenome_aware_deepvariant-head784362481"

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
  --reads "${HG002_BAM}" \
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
  /benchmark/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz \
  /output/HG002.chr20.output.vcf.gz \
  -f /benchmark/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval \
  --pass-only \
  -l chr20
```
