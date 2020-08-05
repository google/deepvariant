# Using DeepVariant for small variant calling from PacBio HiFi reads

In this case study we describe applying DeepVariant to PacBio HiFi reads to call
variants. We will call small variants from a publicly available whole genome
HiFi dataset from PacBio.

In v0.8, DeepVariant released a model for PacBio HiFi data. Starting from
v0.10.0, sequence from amplified libraries is included in our PacBio HiFi
training set, providing a significant accuracy boost to variant detection from
amplified HiFi data.
In this case study we will apply the PacBio model by specifying `PACBIO` in
the `model_type` parameter in the `run_pacbio_case_study_docker.sh` script.

This case study is run on a standard Google Cloud instance. There are no special
hardware or software requirements for running this case study. For consistency
we use Google Cloud instance with 64 cores and 128 GB of memory. This is NOT the
fastest or cheapest configuration. For more scalable execution of DeepVariant
see the [External Solutions] section.

## Case study overview

Calling small variants using DeepVariant involves multiple steps:

1.  Creating examples. Candidate variants are extracted from an input BAM file
    (previously aligned).
2.  Calling variants. Applying DeepVariant convolutional neural network (CNN)
    model to infer variants.
3.  Exporting results to VCF.

In addition we use Hap.py ([https://github.com/Illumina/hap.py]) to calculate
accuracy metrics.

There are multiple ways to run DeepVariant:

-   Build a binary from the source.
-   Download prebuilt binaries.
-   Download an official DeepVariant Docker image.

This case study is run using the official DeepVariant Docker image.

## Running

For simplicity we provide a script that downloads the input data and runs all
the steps described above using DeepVariant Docker image. **Please note, that if
you create your own script `make_examples` must be called with
`--norealign_reads --vsc_min_fraction_indels 0.12
--alt_aligned_pileup "diff_channels"` flags for PacBio long reads.**

*   Create a Google Cloud virtual instance. This command creates a virtual
    instance with 64 cores and 128 GB of memory.

```shell
gcloud beta compute instances create "${USER}-cpu"  \
  --scopes "compute-rw,storage-full,cloud-platform" \
  --image-family "ubuntu-1604-lts" \
  --image-project "ubuntu-os-cloud" \
  --machine-type "custom-64-131072" \
  --boot-disk-size "300" \
  --zone "us-west1-b" \
  --min-cpu-platform "Intel Skylake"
```

*   Login to a newly created instance:

```shell
gcloud compute ssh "${USER}-cpu" --zone "us-west1-b"
```

*   Download and run the case study script:

```shell
curl https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_pacbio_case_study_docker.sh | bash
```

## Script description

Before running the DeepVariant steps, the following input data is downloaded:

*   BAM file: HG002.pfda_challenge.grch38.phased.bam

    The original FASTQ file comes from the
    [PrecisionFDA Truth Challenge V2](https://precision.fda.gov/challenges/10).
    We mapped with pbmm2.

*   FASTA file: GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

    The original file came from:
    [ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids).
    Because DeepVariant requires bgzip files, we had to unzip and bgzip it, and
    create corresponding index files.

*   Truth VCF and BED

    These come from NIST, as part of the
    [Genome in a Bottle project](http://jimb.stanford.edu/giab/). They are
    downloaded from
    [ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020)

Next the following steps are executed:

*   `make_examples`. This step creates small variant candidates and stores them
    in TensorFlow format.

*   `call_variants`. This step applies DeepVariant DNN to call small variants.

*   `postprocess_variants`. This step converts data from TensorFlow format to
    VCF.

*   `hap.py` ([https://github.com/Illumina/hap.py]) program from Illumina is
    used to evaluate the resulting vcf file. This serves as a check to ensure
    the three DeepVariant commands ran correctly and produced a high-quality
    results.

## Runtime metrics

Step                               | Wall time
---------------------------------- | ---------
`make_examples`                    | ~ 162m
`call_variants`                    | ~ 206m
`postprocess_variants` (with gVCF) | ~ 72m

## Accuracy metrics

The PacBio model was trained using the HG002 genome (the same genome we use for
this case study) with chromosomes 20, 21, 22 excluded. Therefore, we run
evaluation on chr20.

Type  | # TP  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ----- | ---- | ---- | -------- | --------- | ---------
INDEL | 11051 | 205  | 148  | 0.981787 | 0.987240  | 0.984506
SNP   | 71277 | 56   | 6    | 0.999215 | 0.999916  | 0.999565

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[https://github.com/Illumina/hap.py]: https://github.com/Illumina/hap.py
