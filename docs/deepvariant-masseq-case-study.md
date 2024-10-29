# DeepVariant PacBio Iso-Seq/MAS-Seq Case Study

This case study will demonstrate how to run DeepVariant using the
PacBio Iso-Seq/MAS-Seq model, and evaluate the result using `hap.py`.

## Overview

### Data

We will use these data in our analysis. Files will be downloaded in subsequent
steps.

*   HG004 MAS-Seq BAM (chr20)
*   Model Checkpoint Files
*   GRCh38 Reference + Index
*   GIAB benchmark data

## Prepare Data

### Setup directories

Lets first create directories to organize files.

```bash
mkdir -p data benchmark reference output happy
```

### Download the GRCh38 Reference

We will be using GRCh38 for this case study.

```bash
FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai
```

### Download Genome in a Bottle Benchmarks

We will benchmark our variant calls against v4.2.1 of the Genome in a Bottle
small variant benchmarks for HG004.

```bash
FTPDIR=ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38

curl -L ${FTPDIR}/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl -L ${FTPDIR}/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -L ${FTPDIR}/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

### Download HG004 BAM (chr20)

For this case study, we download the chr20 of a HG004 MAS-Seq BAM.

```bash
HTTPDIR=https://storage.googleapis.com/deepvariant/masseq-case-study

curl -L ${HTTPDIR}/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.chr20.bam > data/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.chr20.bam
curl -L ${HTTPDIR}/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.chr20.bam.bai > data/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.chr20.bam.bai
```


### Download a BED file for evaluation

When evaluating, we will use a BED file that restricts to exon regions, and only
include regions where the BAM file has 10x or more coverage.

```bash
HTTPDIR=https://storage.googleapis.com/deepvariant/masseq-case-study

curl -L ${HTTPDIR}/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.depth.10x.exons.bed > data/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.depth.10x.exons.bed
```




### Download the MAS-Seq model

Finally, lets download the MAS-Seq model that we will use to call variants.

```bash
gsutil cp -R gs://deepvariant/models/DeepVariant/1.8.0/savedmodels/deepvariant.masseq.savedmodel .
```

### Running DeepVariant MAS-Seq on a CPU-only machine

The command below will run the DeepVariant MAS-Seq model and produce an output
VCF (`output/out.vcf.gz`).

```bash
BIN_VERSION="head687331500"

sudo docker run \
  -v "$(pwd):$(pwd)" \
  -w $(pwd) \
  google/deepvariant:"${BIN_VERSION}" \
  run_deepvariant \
    --model_type=PACBIO \
    --customized_model=deepvariant.masseq.savedmodel \
    --ref=reference/GRCh38_no_alt_analysis_set.fasta \
    --reads=data/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.chr20.bam \
    --output_vcf=output/HG004.output.vcf.gz \
    --num_shards=$(nproc) \
    --regions=chr20 \
    --make_examples_extra_args="phase_reads=true,sort_by_haplotypes=true,parse_sam_aux_fields=true,realign_reads=false,vsc_min_fraction_indels=0.12,alt_aligned_pileup=diff_channels,trim_reads_for_pileup=true,pileup_image_width=199,min_mapping_quality=1,track_ref_reads=true,partition_size=25000,max_reads_per_partition=0,max_reads_for_dynamic_bases_per_region=1500" \
    --disable_small_model=true \
    --intermediate_results_dir=output/intermediate_results_dir
```

**Flag summary**

*   `--model_type` - Sets the model and options, but we will override the model
    with `--customized model`.
*   `--customized_model` - Points to a model trained using MAS-Seq data.
*   `--ref` - Specifies the reference sequence.
*   `--reads` - Specifies the input bam file.
*   `--output_vcf` - Specifies the output variant file.
*   `--num_shards` - Sets the number of shards to the number of available
    processors (`$(nproc)`). This is used to perform parallelization.
*   `--regions` - Restricts to chr20 to make this case study faster.
*   `--make_examples_extra_args=` - Passes additional arguments to
    make_examples.
*   `--intermediate_results_dir` - Outputs results to an intermediate directory.
    This is optional. If you don't need the intermediate files, no need to
    specify this flag.

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md).

## Benchmark on chr20

```bash
sudo docker run \
  -v $(pwd):$(pwd) \
  -w $(pwd) \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    output/HG004.output.vcf.gz \
    -f benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    -r reference/GRCh38_no_alt_analysis_set.fasta \
    -o happy/happy.output \
    --engine=vcfeval \
    --pass-only \
    -l chr20 \
    --target-regions=data/HG004.giab_na24143.hifi_reads.lima.0--0.lima.IsoSeqX_bc01_5p--IsoSeqX_3p.refined.grch38.mm2.splitN.fc.depth.10x.exons.bed \
    --threads=$(nproc)
```

**Flag summary**

*   `-f` - Sets the benchmark regions (regions of interest that we want to
    benchmark.)
*   `-r` - Sets the reference genome.
*   `-o` - Specifies the output location.
*   `--engine` - Sets the variant comparison engine. See
    [hap.py documentation](https://github.com/Illumina/hap.py) for details.
*   `--pass-only` - Restricts benchmarking to variants that have passed all
    filters.
*   `-l` - Comma-separated list of locations.
*   `--target-regions` - Restricts analysis to given regions only.
*   `--threads` - Level of parallelization to use.

**Output:**

The above command should output the following results:

```
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL          135       106        29          161        13         39      7      2       0.785185          0.893443        0.242236         0.835823                     NaN                     NaN                   2.628571                   2.651163
INDEL   PASS          135       106        29          161        13         39      7      2       0.785185          0.893443        0.242236         0.835823                     NaN                     NaN                   2.628571                   2.651163
  SNP    ALL         1002       935        67          978        12         31      7      1       0.933134          0.987328        0.031697         0.959466                2.795455                2.761538                   2.083077                   2.046729
  SNP   PASS         1002       935        67          978        12         31      7      1       0.933134          0.987328        0.031697         0.959466                2.795455                2.761538                   2.083077                   2.046729
```
