# DeepVariant via Docker

You may run DeepVariant binaries via [Docker](http://www.docker.com). The image
has all DeepVariant binaries and their dependencies preinstalled, which makes it
easier to run DeepVariant locally, on a VM, or at scale using the
[Genomics Pipelines API](https://cloud.google.com/genomics/v1alpha2/pipelines).

## Setup

If you have not already done so, please follow the instructions on the
[quick start](deepvariant-quick-start.md) page.

## Run DeepVariant at scale using the Genomics Pipelines API

Please follow the instructions on
[Google Cloud Genomics](https://cloud.google.com/genomics/deepvariant)
for running DeepVariant at scale using the
[Genomics Pipelines API](https://cloud.google.com/genomics/v1alpha2/pipelines).

## Run the DeepVariant docker image locally or on a dedicated VM

### Prerequisites

You need to install docker to run the image locally or on a dedicated VM.
General installation instructions are
[on the Docker site](https://docs.docker.com/installation/). Ubuntu installation
instructions are
[here](https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/#install-using-the-repository).
To avoid having to run the docker commands using `sudo`, follow the instructions
[here](https://docs.docker.com/engine/installation/linux/linux-postinstall/#manage-docker-as-a-non-root-user).

### Run the pipeline

Copy all required input files to a directory that can be accessed within the
docker image. In this example, we are going to use the quickstart test data:

```bash
mkdir -p input
gsutil -m cp gs://deepvariant/quickstart-testdata/* input/
```

Copy the model:

```bash
mkdir -p models
gsutil -m cp gs://deepvariant/models/DeepVariant/0.4.0/DeepVariant-inception_v3-0.4.0+cl-174375304.data-wgs_standard/* models/
```

Run the docker image mounting the `input` and `models` folders. You may need
to run with `sudo` if you get permission denied errors or follow the
instructions
[here](https://docs.docker.com/engine/installation/linux/linux-postinstall/#manage-docker-as-a-non-root-user)
to manage docker as a non-root user.

```bash
# You may use the 'latest' label to get the latest version.
IMAGE_VERSION=0.4.1
gcloud docker -- pull gcr.io/deepvariant-docker/deepvariant:$IMAGE_VERSION
docker run -it -v $PWD/input:/dv2/input -v $PWD/models:/dv2/models \
    gcr.io/deepvariant-docker/deepvariant:$IMAGE_VERSION
```

Finally, within the docker image, run:

```bash
./opt/deepvariant/bin/make_examples \
    --mode calling \
    --ref /dv2/input/ucsc.hg19.chr20.unittest.fasta.gz \
    --reads /dv2/input/NA12878_S1.chr20.10_10p1mb.bam \
    --examples output.examples.tfrecord \
    --regions "chr20:10,000,000-10,010,000"

./opt/deepvariant/bin/call_variants \
    --outfile call_variants_output.tfrecord \
    --examples output.examples.tfrecord \
    --checkpoint /dv2/models/model.ckpt

./opt/deepvariant/bin/postprocess_variants \
    --ref /dv2/input/ucsc.hg19.chr20.unittest.fasta.gz \
    --infile call_variants_output.tfrecord \
    --outfile output.vcf
```

You may mount an additional `output` folder (using `-v`) when running the docker
image to store the resulting VCF file directly on your local machine.
