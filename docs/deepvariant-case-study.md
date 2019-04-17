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
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
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
curl https://raw.githubusercontent.com/google/deepvariant/r0.8/scripts/run_wgs_case_study_docker.sh | bash
```

### Running on a machine with GPU

Here is an example command:

```
sudo nvidia-docker run \
  -v "${DATA_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  gcr.io/deepvariant-docker/deepvariant_gpu:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/${OUTPUT_VCF} \
  --output_gvcf=/output/${OUTPUT_GVCF} \
  --num_shards=${N_SHARDS}
```

Note that instead of using `docker`, we're using `nvidia-docker` to make use of
the GPU. `call_variants` is the only step that uses the GPU.
`make_examples` and `postprocess_variants` do not run on GPU.


The script [run_wgs_case_study_docker_gpu.sh] shows a full example, including
calling a script [install_nvidia_docker.sh] that helps you install
`nvidia-docker`.

### Runtime

See
[this page](deepvariant-details.md#commands-for-requesting-machines-used-in-case-studies)
for the commands used to obtain different machine types on Google Cloud.
Currently, the `call_variants` step cannot use more than 1 GPU on the same
machine, and the `postprocess_variants` step is single-process, single-thread.

With the example in [run_wgs_case_study_docker.sh] on a [CPU machine],

Step                               | Hardware | Wall time
---------------------------------- | -------- | ---------
`make_examples`                    | 64 CPUs  | ~ 1h 10m
`call_variants`                    | 64 CPUs  | ~ 2h 45m
`postprocess_variants` (with gVCF) | 1 CPU    | ~ 35m

With the example in [run_wgs_case_study_docker_gpu.sh] on a [GPU machine],

Step                               | Hardware            | Wall time
---------------------------------- | ------------------- | ---------
`make_examples`                    | 16 CPUs             | ~ 3h 25m
`call_variants`                    | 1 P100 GPU, 16 CPUs | ~    50m
`postprocess_variants` (with gVCF) | 1 CPU               | ~    30m

Since `make_examples` doesn't utilize GPUs, bringing up one GPU machine for all
steps might not be the most cost-effective solution. For more scalable execution
of DeepVariant see the [External Solutions] section.

## Data description

The original source of data used in this case study:

1.  BAM file: HG002_NIST_150bp_50x.bam

    The original FASTQ file comes from the
    [PrecisionFDA Truth Challenge](https://precision.fda.gov/challenges/truth/).
    We ran it through
    [PrecisionFDA's BWA-MEM app](https://precision.fda.gov/apps/app-BpF9YGQ0bxjbjk6Fx1F0KJF0)
    with default setting, and then got the HG002_NIST_150bp_50x.bam file as
    output.

    Since 2019-02-26, the file was further processed using SAMtools 1.9 and
    HTSlib 1.9 to address formatting problems caused by an older version of
    HTSlib:

    ```
    samtools view -bh HG002_NIST_150bp_50x.bam -o HG002_NIST_150bp_50x.bam
    ```

    The FASTQ files are originally from the
    [Genome in a Bottle Consortium](http://jimb.stanford.edu/giab-resources/).
    For more information, see this
    [Scientific Data paper](https://www.nature.com/articles/sdata201625).

1.  FASTA file: hs37d5.fa.gz

    The original file came from:
    [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence).
    Because DeepVariant requires bgzip files, we had to unzip and bgzip it, and
    create corresponding index files.

1.  Truth VCF and BED

    These come from NIST, as part of the
    [Genome in a Bottle project](http://jimb.stanford.edu/giab/). They are
    downloaded from
    [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/)

## Variant call quality

We used the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This serves as a check
to ensure the three DeepVariant commands ran correctly and produced high-quality
results.

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 1488 | 944  | 0.996798 | 0.998048  | 0.997423
SNP   | 1576 | 725  | 0.999483 | 0.999762  | 0.999623

[install_nvidia_docker.sh]: https://github.com/google/deepvariant/blob/r0.8/scripts/install_nvidia_docker.sh
[run_wgs_case_study_docker.sh]: https://github.com/google/deepvariant/blob/r0.8/scripts/run_wgs_case_study_docker.sh
[run_wgs_case_study_docker_gpu.sh]: https://github.com/google/deepvariant/blob/r0.8/scripts/run_wgs_case_study_docker_gpu.sh
[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU machine]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
[GPU machine]: deepvariant-details.md#command-for-a-gpu-machine-on-google-cloud-platform
