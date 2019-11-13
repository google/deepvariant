# Using DeepVariant for small variant calling from PacBio CCS reads

In this case study we describe applying DeepVariant to PacBio CCS reads to call
variants. We will call small variants from a publicly available whole genome
from PacBio.

Case study is run on a standard Google Cloud instance. There are no special
hardware or software requirements for running this case study. For consistency
we use Google Cloud instance with 64 cores and 128 GB of memory. This is NOT the
fastest or cheapest configuration. For more scalable execution of DeepVariant
see the [External Solutions] section.

In v0.8 DeepVariant released a special model that works with PacBio data. In
this case study we will apply PacBio model by specifying `PACBIO` in
`model_type` parameter in the `run_pacbio_case_study_docker.sh` script.

## Case study overview

Calling small variants using DeepVariant involves multiple steps:

1.  Creating examples. Variants candidates are extracted from input BAM file
    (previously aligned).
2.  Calling Variants. Applying DeepVariant Convolutional Neural Network (CNN)
    model to infer variants.
3.  Exporting results to VCF.

In addition we use Hap.py ([https://github.com/Illumina/hap.py]) to calculate
accuracy metrics.

There are multiple ways to run DeepVariant:

-   Build a binary from the source.
-   Download prebuilt binaries.
-   Download an official DeepVariant Docker image.

This case study is run using official DeepVariant Docker image.

## Running

For simplicity we provide a script that downloads the input data and runs all
the steps descibed above using DeepVariant Docker image. **Please note, that if
you create your own script make_examples must be called with
`--norealign_reads --vsc_min_fraction_indels 0.12` flag for PacBio long reads.**

1.  Create a Google Cloud virtual instance. This command creates a virtual
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

1.  Login to a newly created instance:

```shell
gcloud compute ssh "${USER}-cpu" --zone "us-west1-b"
```

1.  Download and run the case study script:

```shell
curl https://raw.githubusercontent.com/google/deepvariant/r0.9/scripts/run_pacbio_case_study_docker.sh | bash
```

## Script description

Before running DeepVariant steps following input data is downloaded:

*   BAM file: pacbio.8M.30x.bam. Publicly available PacBio BAM file.
    [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/
    PacBio_SequelII_CCS_11kb/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam)

*   Reference file: hs37d5.fa.gz. The original file came from:
    [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence).
    Because DeepVariant requires bgzip files, we had to unzip and bgzip it, and
    create corresponding index files.

*   Truth VCF and BED. These come from NIST, as part of the
    [Genome in a Bottle project](http://jimb.stanford.edu/giab/). They are
    downloaded from
    [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/)

Next following steps are executed:

*   make_examples. This step creates small variant candidates and stores them in
    TensorFlow format.

*   call_variants. This step applys DeepVariant DNN to call small variants.

*   postprocess_variants. This step converts data from TensorFlow format to VCF.

*   `hap.py` ([https://github.com/Illumina/hap.py]) program from Illumina is
    used to evaluate the resulting vcf file. This serves as a check to ensure
    the three DeepVariant commands ran correctly and produced a high-quality
    results.

## Runtime metrics

Step                               | Wall time
---------------------------------- | ---------
`make_examples`                    | ~ 63m
`call_variants`                    | ~ 3h 19m
`postprocess_variants` (with gVCF) | ~ 60m

## Accuracy metrics

PacBio model was trained using HG002 genome (the same genome we use for this
case study) with chromosomes 20, 21, 22 excluded. Therefore, we run evaluation
on chr20.

Type  | # TP  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ----- | ---- | ---- | -------- | --------- | ---------
INDEL | 9998  | 174  | 160  | 0.982894 | 0.984851  | 0.983872
SNP   | 65200 | 42   | 100  | 0.999356 | 0.998470  | 0.998913

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[https://github.com/Illumina/hap.py]: https://github.com/Illumina/hap.py
