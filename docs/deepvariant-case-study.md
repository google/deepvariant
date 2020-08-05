# DeepVariant whole genome sequencing case study

In this case study, we describe applying DeepVariant to a real WGS sample. We
provide two examples, one running on one CPU-only machine, and one running on a
machine with GPU. And finally we assess the quality of the DeepVariant variant
calls with `hap.py`.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. We report the runtime with
[specific machine type](deepvariant-details.md#commands-for-requesting-machines-used-in-case-studies)
for the sake of consistency in reporting run time. This is NOT the fastest or
cheapest configuration. For more scalable execution of DeepVariant see the
[External Solutions] section.

## Running DeepVariant with one command

DeepVariant pipeline consists of 3 steps: `make_examples`, `call_variants`, and
`postprocess_variants`. You can now run DeepVariant with one command using the
`run_deepvariant` script.

### Running on a CPU-only machine

Here is an example command:

```
sudo docker run \
  -v "${DATA_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/${OUTPUT_VCF} \
  --output_gvcf=/output/${OUTPUT_GVCF} \
  --num_shards=${N_SHARDS}
```

By specifying `--model_type=WGS`, you'll be using a model that is best suited
for Illumina Whole Genome Sequencing data.

The script [run_wgs_case_study_docker.sh] shows a full example of which data to
download, and run DeepVariant with the data.

Before you run the script, you can read through all sections to understand the
details. Here is a quick way to get the script and run it:

```bash
curl https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wgs_case_study_docker.sh | bash
```

### Running on a machine with GPU

Here is an example command:

```bash
sudo docker run --gpus 1 \
  -v "${DATA_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/${OUTPUT_VCF} \
  --output_gvcf=/output/${OUTPUT_GVCF} \
  --num_shards=${N_SHARDS}
```

Note that we add `--gpus 1` to make use of the GPU. `call_variants` is the only
step that uses the GPU, and can only use one at a time.
`make_examples` and `postprocess_variants` do not run on GPU.


The script [run_wgs_case_study_docker_gpu.sh] shows a full example, including
calling a script [install_nvidia_docker.sh] that helps you install.

### Runtime

See
[this page](deepvariant-details.md#commands-for-requesting-machines-used-in-case-studies)
for the commands used to obtain different machine types on Google Cloud.
Currently, the `call_variants` step cannot use more than 1 GPU on the same
machine, and the `postprocess_variants` step is single-process, single-thread.

With the example in [run_wgs_case_study_docker.sh] on a [CPU machine],

Step                               | Hardware | Wall time
---------------------------------- | -------- | ---------
`make_examples`                    | 64 CPUs  | ~ 1.5h
`call_variants`                    | 64 CPUs  | ~ 4.5h
`postprocess_variants` (with gVCF) | 1 CPU    | ~ 1.5h

With the example in [run_wgs_case_study_docker_gpu.sh] on a [GPU machine],

Step                               | Hardware            | Wall time
---------------------------------- | ------------------- | ---------
`make_examples`                    | 16 CPUs             | ~ 5h
`call_variants`                    | 1 P100 GPU, 16 CPUs | ~ 1h 10m
`postprocess_variants` (with gVCF) | 1 CPU               | ~ 1h 10m

Since `make_examples` doesn't utilize GPUs, bringing up one GPU machine for all
steps might not be the most cost-effective solution. For more scalable execution
of DeepVariant see the [External Solutions] section.

## Data description

The original source of data used in this case study:

1.  BAM file: HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam

    The original FASTQ file comes from the
    [PrecisionFDA Truth Challenge V2](https://precision.fda.gov/challenges/10).
    We mapped with BWA-MEM.

1.  FASTA file: GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

    The original file came from:
    [ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids).
    Because DeepVariant requires bgzip files, we had to unzip and bgzip it, and
    create corresponding index files.

1.  Truth VCF and BED

    These come from NIST, as part of the
    [Genome in a Bottle project](http://jimb.stanford.edu/giab/). They are
    downloaded from
    [ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020)

## Variant call quality

We used the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This serves as a check
to ensure the three DeepVariant commands ran correctly and produced high-quality
results.

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1\_Score
----- | ------- | ----- | ---- | -------- | --------- | ---------
INDEL | 522259  | 3207  | 1187 | 0.993897 | 0.997825  | 0.995857
SNP   | 3345988 | 19352 | 3955 | 0.994250 | 0.998820  | 0.996530

[install_nvidia_docker.sh]: ../scripts/install_nvidia_docker.sh
[run_wgs_case_study_docker.sh]: ../scripts/run_wgs_case_study_docker.sh
[run_wgs_case_study_docker_gpu.sh]: ../scripts/run_wgs_case_study_docker_gpu.sh
[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU machine]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
[GPU machine]: deepvariant-details.md#command-for-a-gpu-machine-on-google-cloud-platform
