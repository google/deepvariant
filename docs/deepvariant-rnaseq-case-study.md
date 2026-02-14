# DeepVariant RNA-seq Case Study

This case study will demonstrate how to run DeepVariant using the RNA-seq model,
and evaluate the result using `hap.py`.

## Overview

### Tools

We will use the following tools:

*   [Docker](https://docs.docker.com/get-docker/) - Used to run DeepVariant.
*   [mosdepth](https://github.com/brentp/mosdepth) - For calculating coverage.
*   [bedtools](https://bedtools.readthedocs.io) - Used to intersect bedfiles.
*   [hap.py](https://github.com/illumina/hap.py) - Used to evaluate the results.
    We will use Docker to run `hap.py`.

### Data

We will use these data in our analysis. Files will be downloaded in subsequent
steps.

*   HG005 RNA-seq BAM
*   Model Checkpoint Files
*   GRCh38 Reference + Index
*   CDS bedfile (chr20 only)
*   GIAB benchmark data

## Prepare Data

### Setup directories

Lets first create directories to organize files.

```bash
mkdir -p data benchmark reference model output happy
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
small variant benchmarks for HG005. We will also restrict analysis to CDS
regions on chromosome 20 to make this demonstration quicker.

The benchmarks consist of a bedfile containing confident regions, a VCF of
'true' variants, and a VCF index.

```bash
FTPDIR=ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38

curl -L ${FTPDIR}/HG005_GRCh38_1_22_v4.2.1_benchmark.bed > benchmark/HG005_GRCh38_1_22_v4.2.1_benchmark.bed
curl -L ${FTPDIR}/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -L ${FTPDIR}/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

### Download and extract a CDS bedfile

Next, we will download a [gencode](https://www.gencodegenes.org/) gff3
annotation and extract a bed file of chr20 CDS regions.

```bash
curl -L https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gff3.gz > data/gencode.v41.basic.annotation.gff3.gz

# Extract chr20 CDS regions and convert to bed file.
gzip -dc data/gencode.v41.basic.annotation.gff3.gz | \
awk -v OFS='\t' '$1 == "chr20" && $3 == "CDS"  && $4 < $5 { print $1, $4, $5, "CDS" }' | \
awk '!dup[$0]++' > data/chr20_CDS.bed
```

### Download HG005 BAM

We'll use HG005 poly-A selected Illumina RNA-seq reads that are publicly
available.

```bash
HTTPDIR=https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/rna/illumina/mrna

curl -L ${HTTPDIR}/hg005_gm26107.mrna.grch38.bam > data/hg005_gm26107.mrna.grch38.bam
curl -L ${HTTPDIR}/hg005_gm26107.mrna.grch38.bam.bai > data/hg005_gm26107.mrna.grch38.bam.bai
```

### Generate a 3x coverage file

RNA-seq data is only observed in regions that are expressed in a given sample.
Therefore, we will restrict our evaluation to regions of the BAM file that reach
a minimum threshold of 3x in our truth dataset intersected with the confident
GIAB regions. This allows us to better evaluate the accuracy of the model when
it is feasible for a variant to be called from RNA-seq data.

```bash
# Generate a coverage file, and filter for 3x.
sudo docker run \
  -v "$(pwd):$(pwd)" \
  -w $(pwd) \
  -it quay.io/biocontainers/mosdepth:0.3.1--h4dc83fb_1 \
  mosdepth \
    --threads $(nproc) \
    data/hg005_coverage \
    data/hg005_gm26107.mrna.grch38.bam
```

For these next commands, we will run Docker interactively to execute a series of
commands. Run the following command to launch a bedtools container.

```bash
sudo docker run \
  -v "$(pwd):$(pwd)" \
  -w $(pwd) \
  -it quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6 \
  /bin/bash
```

### Extract regions with 3x coverage, and filter out unused contigs.

We will restrict our analysis to regions with a minimum of 3x coverage.

```bash
# (Run within the bedtools container)
min_coverage=3
gzip -dc data/hg005_coverage.per-base.bed.gz | \
egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | \
awk -v OFS="\t" -v min_coverage=${min_coverage} '$4 >= min_coverage { print }' | \
bedtools merge -d 1 -c 4 -o mean -i - > data/hg005_3x.bed
```

### Intersect coverage with CDS regions.

Now we will intersect our 3x bedfile with the CDS bed file:

```bash
# (Run within the bedtools container)
bedtools intersect \
-a data/hg005_3x.bed \
-b data/chr20_CDS.bed > data/chr20_CDS_3x.bed

# We will also intersect this file with confident GIAB regions
bedtools intersect \
-a benchmark/HG005_GRCh38_1_22_v4.2.1_benchmark.bed \
-b data/chr20_CDS_3x.bed > benchmark/chr20_CDS_3x.benchmark_regions.bed
```

Make sure to exit docker environment after running the last command by running:

```
exit
```

We now have a bed file of CDS regions intersected with 3x coverage regions
called `data/chr20_CDS_3x.bed`.

### Directory Structure

After you have run the steps above, your directory structure should look like
this:

```
.
├── benchmark
│   ├── chr20_CDS_3x.benchmark_regions.bed
│   ├── HG005_GRCh38_1_22_v4.2.1_benchmark.bed
│   ├── HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
│   └── HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
├── data
│   ├── chr20_CDS_3x.bed
│   ├── chr20_CDS.bed
│   ├── gencode.v41.basic.annotation.gff3.gz
│   ├── hg005_3x.bed
│   ├── hg005_coverage.mosdepth.global.dist.txt
│   ├── hg005_coverage.mosdepth.summary.txt
│   ├── hg005_coverage.per-base.bed.gz
│   ├── hg005_coverage.per-base.bed.gz.csi
│   ├── hg005_gm26107.mrna.grch38.bam
│   └── hg005_gm26107.mrna.grch38.bam.bai
├── happy
├── model
│   ├── model.ckpt.data-00000-of-00001
│   ├── model.ckpt.index
│   └── model.ckpt.meta
├── output
└── reference
    ├── GRCh38_no_alt_analysis_set.fasta
    └── GRCh38_no_alt_analysis_set.fasta.fai
```

### Running DeepVariant RNA-seq on a CPU-only machine

The command below will run the DeepVariant RNA-seq model and produce an output
VCF (`output/out.vcf.gz`).

```bash
BIN_VERSION="1.10.0"

sudo docker run \
  -v "$(pwd):$(pwd)" \
  -w $(pwd) \
  google/deepvariant:"${BIN_VERSION}" \
  run_deepvariant \
    --model_type=RNASEQ \
    --ref=reference/GRCh38_no_alt_analysis_set.fasta \
    --reads=data/hg005_gm26107.mrna.grch38.bam \
    --output_vcf=output/HG005.output.vcf.gz \
    --disable_small_model \
    --num_shards=$(nproc) \
    --regions=data/chr20_CDS_3x.bed \
    --intermediate_results_dir output/intermediate_results_dir
```

**Flag summary**

*   `--model_type` - Sets the model and options, but we will override the model
    with `--customized model`.
*   `--customized_model` - Points to a model trained using RNA-seq data.
*   `--ref` - Specifies the reference sequence.
*   `--reads` - Specifies the input bam file.
*   `--output_vcf` - Specifies the output variant file.
*   `--num_shards` - Sets the number of shards to the number of available
    processors (`$(nproc)`). This is used to perform parallelization.
*   `--regions` - Restricts analysis to 3x chr20 CDS regions only.
*   `--disable_small_model` - Disables the small model from running.
*   `--make_examples_extra_args=` - Passes additional arguments to
    make_examples.
    *   `split_skip_reads=true` - *Important!* This flag is critical for RNA-seq
        variant calling to work properly. It enables RNA-seq data to be
        processed efficiently.
    *   `channels=''` - Resets the channel list to be appropriate for the
        RNA-seq model. model.
*   `--intermediate_results_dir` - Outputs results to an intermediate directory.

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md).

## Benchmark on chr20

```bash
sudo docker run \
  -v $(pwd):$(pwd) \
  -w $(pwd) \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    benchmark/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    output/HG005.output.vcf.gz \
    -f benchmark/chr20_CDS_3x.benchmark_regions.bed \
    -r reference/GRCh38_no_alt_analysis_set.fasta \
    -o happy/happy.output \
    --engine=vcfeval \
    --pass-only \
    --target-regions=data/chr20_CDS_3x.bed \
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
*   `--target-regions` - Restricts analysis to given regions only.
*   `--threads` - Level of parallelization to use.

**Output:**

The above command should output the following results:

```
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL            9         7         2           13         1          5      0      0       0.777778          0.875000        0.384615         0.823529                     NaN                     NaN                   0.800000                   1.166667
INDEL   PASS            9         7         2           13         1          5      0      0       0.777778          0.875000        0.384615         0.823529                     NaN                     NaN                   0.800000                   1.166667
  SNP    ALL          287       276        11          327         7         44      2      0       0.961672          0.975265        0.134557         0.968421                   4.125                3.954545                   1.141791                   1.194631
  SNP   PASS          287       276        11          327         7         44      2      0       0.961672          0.975265        0.134557         0.968421                   4.125                3.954545                   1.141791                   1.194631
```
