# DeepVariant Nextflow Pipeline

This directory contains the Nextflow pipeline for running DeepVariant. It allows
you to run DeepVariant in a highly scalable and reproducible manner, supporting
both local execution and parallelization.

## Prerequisites

To run this pipeline, you need:

1.  **Nextflow** (version 23.04 or later recommended):
    [Installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2.  **Docker**: [Installation instructions](https://docs.docker.com/get-docker/)
3.  **Java** (version 11 or 17, required by Nextflow).

## Directory Structure

*   `deepvariant.nf`: The main Nextflow pipeline script.
*   `nextflow.config`: Configuration file with default settings (Docker enabled,
    local executor).
*   `modules/`: Pipeline-specific submodules.
*   `case_studies/`: Example sample sheets for running case
    studies with public data.

## Quick Start (Local Run)

1.  **Prepare your input data**: You need a reference genome (FASTA), a BAM/CRAM
    file with aligned reads, and optionally truth VCF/BED files for evaluation.
2.  **Create a sample sheet**: Define your samples in a YAML file. You can use
    the examples in `case_studies/` as a template. Example
    `sample_sheet.yaml`:

    ```yaml
    uid: sample-wgs
    sample: HG003
    bam: /path/to/HG003.bam
    model_type: WGS
    ref: /path/to/GRCh38.fa
    # Optional for evaluation
    truth_vcf: /path/to/truth.vcf.gz
    truth_bed: /path/to/truth.bed
    docker_image: google/deepvariant:1.10.0
    ```

3.  **Run the pipeline**:

    ```bash
    NXF_VER=25.04.2 nextflow run pipelines/deepvariant_nf/deepvariant.nf \
      -config pipelines/nextflow.config \
      --sample_sheet sample_sheet.yaml \
      --output_dir output_dir
    ```

## Parallel Execution

DeepVariant can be computationally intensive. This pipeline supports
parallelizing the `make_examples` step by splitting the genome into chunks and
running them concurrently.

To enable parallel execution, use the `--num_instances` parameter:

```bash
NXF_VER=25.04.2 nextflow run pipelines/deepvariant_nf/deepvariant.nf \
  -config pipelines/nextflow.config \
  --sample_sheet sample_sheet.yaml \
  --num_instances 4 \
  --output_dir output_dir_parallel
```

This will split the work for each uid into 4 parallel instances, significantly
reducing runtime if your machine has sufficient CPU and memory.

## Customizing Configuration

The default `nextflow.config` is tuned for local execution with moderate
resources. You can customize it for your specific environment (e.g., allocating
more memory/CPUs, or configuring it for HPC/Cloud executors).

Refer to the comments in `pipelines/nextflow.config` and the official Nextflow
documentation for details:

*   [Nextflow Configuration Docs](https://www.nextflow.io/docs/latest/config.html)
*   [Nextflow Profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles)

## Running Case Studies

We provide lightweight sample sheets in `case_studies/` that
reference public datasets hosted on Google Cloud Storage. You can use these to
test your installation.

Note: Running these requires internet access and the ability of Docker to access
public GCS buckets.

Example running the WGS case study:

```bash
NXF_VER=25.04.2 nextflow run pipelines/deepvariant_nf/deepvariant.nf \
  -config pipelines/nextflow.config \
  --sample_sheet pipelines/deepvariant_nf/case_studies/dv_release_metrics.yaml \
  --uid_filter wgs \
  --output_dir output_wgs_test
```
