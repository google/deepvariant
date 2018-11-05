# DeepVariant whole genome case study

In this case study we describe applying DeepVariant to a real WGS sample.

We provide some guidelines on the computational resources needed for each step.
And finally we assess the quality of the DeepVariant variant calls with
`hap.py`.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. This might not be the fastest or cheapest
configuration for your needs. For more scalable execution of DeepVariant see the
[cost- and speed-optimized, Docker-based
pipelines](https://cloud.google.com/genomics/deepvariant) created for Google
Cloud Platform.

Consult this [document](deepvariant-details.md) for more information about using
GPUs or reading BAM files from Google Cloud Storage (GCS) directly.

## Update since r0.7 : using docker instead of copying binaries.

Starting from the 0.7 release, we use docker to run the binaries instead of
copying binaries to local machines first. You can still read about the previous
approach in
[the Case Study in r0.6](https://github.com/google/deepvariant/blob/r0.6/docs/deepvariant-case-study.md).

We recognize that there might be some overhead of using docker run. But using
docker makes this case study easier to generalize to different versions of Linux
systems. For example, we have verified that you can use docker to run
DeepVariant on other Linux systems such as CentOS 7.

## Update since v0.7.1: Run the case study with one script

Script:
https://github.com/google/deepvariant/blob/r0.7/scripts/run_wgs_case_study_docker.sh

Get the script and run everything. This will install and download everything
needed for this case study. And it will run DeepVariant to generate the output
VCFs, and also run the evaluation as well.

Before you run the script, you can read through all sections to understand the
details.

```bash
wget https://raw.githubusercontent.com/google/deepvariant/r0.7/scripts/run_wgs_case_study_docker.sh -P ${HOME}
chmod +x ${HOME}/run_wgs_case_study_docker.sh
${HOME}/run_wgs_case_study_docker.sh
```

## (Optional) If on Google Cloud Platform, request a machine with this example command

You can run this exercise on any sufficiently capable machine. As a concrete but
simple example of how we performed this study, we used 64-core (vCPU) machine
with 128GiB of memory and no GPU, on the Google Cloud Platform.

We used a command like this to allocate it:

```shell
gcloud beta compute instances create "${USER}-deepvariants-casestudy"  \
--scopes "compute-rw,storage-full,cloud-platform" \
--image-family "ubuntu-1604-lts" \
--image-project "ubuntu-os-cloud" \
--machine-type "custom-64-131072" \
--boot-disk-size "300" \
--boot-disk-type "pd-ssd" \
--zone "us-west1-b"
```

Once the machine is ready, ssh into it:

```
gcloud compute ssh "${USER}-deepvariant-casestudy" --zone "us-west1-b"
```

NOTE: Having an instance up and running could cost you. Remember to delete the
instances you're not using. You can find the instances at:
https://console.cloud.google.com/compute/instances?project=YOUR_PROJECT

## Description for data

The original source of data used in this case study:

1.  BAM file: HG002_NIST_150bp_50x.bam

    The original FASTQ file comes from the [PrecisionFDA Truth
    Challenge](https://precision.fda.gov/challenges/truth/). I ran it through
    [PrecisionFDA's BWA-MEM
    app](https://precision.fda.gov/apps/app-BpF9YGQ0bxjbjk6Fx1F0KJF0) with
    default setting, and then got the HG002_NIST_150bp_50x.bam file as output.
    The FASTQ files are originally from the [Genome in a Bottle
    Consortium](http://jimb.stanford.edu/giab-resources/). For more information,
    see this [Scientific Data
    paper](https://www.nature.com/articles/sdata201625).

1.  FASTA file: hs37d5.fa.gz

    The original file came from:
    [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence).
    Because DeepVariant requires bgzip files, we had to unzip and bgzip it, and
    create corresponding index files.

1.  Truth VCF and BED

    These come from NIST, as part of the [Genome in a Bottle
    project](http://jimb.stanford.edu/giab/). They are downloaded from
    [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/)

## Notes for `call_variants`

We noticed that running multiple `call_variants` on the same machine didn't seem
to save the overall time, because each of the call_variants slowed down when
multiple are running on the same machine.

So if you care about finishing the end-to-end run faster, you could request 64
machines (for example, `n1-standard-8` machines) and run each shard as input
separately, and output to corresponding output shards. For example, the first
machine will run this command:

```shell
time sudo docker run \
  -v "/home/${USER}":"/home/${USER}" \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/call_variants \
  --outfile=${OUTPUT_DIR}/HG002.cvo.tfrecord-00000-of-00064.gz \
  --examples=${OUTPUT_DIR}/HG002.examples.tfrecord-00000-of-00064.gz \
  --checkpoint="${MODEL}"
```

And the rest will process 00001 to 00063. You can also use tools like
[dsub](https://cloud.google.com/genomics/v1alpha2/dsub). In this case study we
only report the runtime on a 1-machine example.

## Notes for `postprocess_variants`

Because this step is single-process, single-thread, if you're orchestrating a
more complicated running pipeline, you might want to request a machine with
fewer cores for this step.

## Resources used by each step

Step                               | wall time
---------------------------------- | -------------------
`make_examples`                    | 113m 19s
`call_variants`                    | 181m 40s
`postprocess_variants` (no gVCF)   | 20m 40s
`postprocess_variants` (with gVCF) | 54m 49s
total time (single machine)        | 315m 39s - 349m 48s

## Variant call quality

We used the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This serves as a check
to ensure the three DeepVariant commands ran correctly and produced high-quality
results.

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 1431 | 916  | 0.996921 | 0.998106  | 0.997513
SNP   | 1329 | 746  | 0.999564 | 0.999755  | 0.999660
