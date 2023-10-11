# Using the DeepVariant hybrid model for calling on PacBio and Illumina data together

In this case study we describe the hybrid DeepVariant model and apply it to the
combination of two datasets:

1.  PacBio HiFi data on HG003 aligned with
    [pbmm2](https://github.com/PacificBiosciences/pbmm2).
2.  Illumina NovaSeq on HG003 aligned with
    [BWA MEM](https://github.com/lh3/bwa).

The FASTQ files come from the
[PrecisionFDA Truth challenge v2](https://precision.fda.gov/challenges/10/view).

They are merged together into a single bam file using `samtools merge`, and then
a new index is created for this hybrid bam using `samtools index`. Note that the
two original bam files must have the same sample name.

Finally, we assess the quality of the DeepVariant variant calls with `hap.py`.

To make it faster to run over this case study, we run only on chromosome 20.

## Background on the hybrid model

This is what the pileup image looks like: The longer PacBio reads are shown at
the top, followed by the shorter Illumina reads at the bottom.

![Example of a hybrid pileup for one variant](images/hybrid_pileup.png)

A DeepVariant hybrid model was first trained for the PrecisionFDA Truth
Challenge V2, and this release model is similar except it has been re-trained
with additional datasets including the HG004 truth set that was held out during
the challenge.

Interestingly, DeepVariant didn't strictly need any code changes to work on
hybrid data -- it worked the first time we tried. But we knew from many previous
experiments that Illumina reads benefit from being realigned to a haplotype
graph, which is too time consuming and unnecessary for the PacBio long reads. We
added a small code change to specifically realign all the short reads to the
haplotype graph, while leaving longer reads with their original alignments. This
created a small but measurable improvement, and was the only code change we made
to enable the hybrid model, aside from training a dedicated hybrid model and
exposing it for easy use through the --model_type parameter in
`run_deepvariant.py`. Much of the work we put into DeepVariant is in
experimenting with different approaches, training on more and better data, and
carefully evaluating the models before releasing them. We did the same with this
hybrid model.

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

We'll use a HG003 BAM file that contains pacbio and illumina data merged
together using `samtools merge`. See the top of this page for more information
on those two datasets.

```bash
mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/hybrid-case-study-testdata

curl ${HTTPDIR}/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.chr20.bam > input/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.chr20.bam
curl ${HTTPDIR}/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.chr20.bam.bai > input/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.chr20.bam.bai
```

## Running DeepVariant

DeepVariant pipeline consists of 3 steps: `make_examples`, `call_variants`, and
`postprocess_variants`. You can run DeepVariant with just one command using the
`run_deepvariant` script.

### Running on a CPU-only machine

Here we specify `--regions chr20` to run on just chromosome 20, saving time so
you can run this case study within about half an hour (tested on 64 CPUs).

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
  --model_type "HYBRID_PACBIO_ILLUMINA" \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads /input/HG003_hybrid_35x_ilmn_35x_pacb.grch38.phased.chr20.bam \
  --output_vcf /output/HG003.output.vcf.gz \
  --output_gvcf /output/HG003.output.g.vcf.gz \
  --num_shards $(nproc) \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir
```

By specifying `--model_type HYBRID_PACBIO_ILLUMINA`, you'll be using a model
that is best suited for (and trained on) the combination of PacBio Hifi long
reads and Illumina short reads.

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

To see the pileup images visually, check out [show_examples](show-examples.md).

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md). Just make sure to use `--model_type
HYBRID_PACBIO_ILLUMINA` when running on combined PacBio and Illumina data.

## Benchmark with hap.py

See [hap.py](https://github.com/illumina/hap.py) documentation for more details
on the parameters and outputs.

```bash
mkdir -p happy

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
INDEL    ALL        10628     10602        26        23385        63      12212     10     51       0.997554          0.994361        0.522215         0.995955                     NaN                     NaN                   1.748961                   2.721448
INDEL   PASS        10628     10602        26        23385        63      12212     10     51       0.997554          0.994361        0.522215         0.995955                     NaN                     NaN                   1.748961                   2.721448
  SNP    ALL        70166     70138        28       105564        43      35354     16     16       0.999601          0.999388        0.334906         0.999494                2.296566                1.812971                   1.883951                   2.187440
  SNP   PASS        70166     70138        28       105564        43      35354     16     16       0.999601          0.999388        0.334906         0.999494                2.296566                1.812971                   1.883951                   2.187440
```

Notice that F1 scores are above 0.999 for SNPs and above 0.995 for indels!
