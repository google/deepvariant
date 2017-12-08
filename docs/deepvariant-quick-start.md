# DeepVariant quick start

## Set up a Google Cloud account

To get started using DeepVariant on Google Cloud Platform (GCP), you first need
to set up an account and a project to contain your cloud resources.

* If you do not have an account yet, you should create one at
[cloud.google.com](https://cloud.google.com). You should then [enable billing
for your account](https://support.google.com/cloud/answer/6288653?hl=en) but
note that if your account is new, [you receive $300 of free
credit](https://cloud.google.com/free/).  Once your cloud account is set up, you
should be able to log in to the [Cloud
Console](https://console.cloud.google.com) to view or administer your cloud
resources.

* From the Cloud Console, [set up a
project](https://cloud.google.com/resource-manager/docs/creating-managing-projects)
to house all of the cloud resources (storage, compute, services) that you will
associate with your use of DeepVariant. For example, if your organization is
AcmeCorp, you might call your project `acmecorp-deepvariant`.

* Finally, please visit the ["Compute Engine" page on Cloud
Console](https://console.cloud.google.com/compute).  You don't need to create
Compute Engine instances at this time, but simply visiting this page will
initialize your compute engine "service account" so that we can authorize it.

(As you progress in your use of Google Cloud Platform, you will likely find it
useful to create a [Cloud
Organization](https://cloud.google.com/resource-manager/docs/creating-managing-organization)
to house your projects. Here are some [best
practices](https://cloud.google.com/docs/enterprise/best-practices-for-enterprise-organizations)
for organizating cloud projects for an enterprise.)

## Install the Google Cloud SDK

The Google Cloud SDK comes with two very useful command line utilities that you
can use on your local workstation---`gcloud`, which lets you administer your
cloud resources, and `gsutil`, which lets you manage and transfer data to Google
Cloud Storage buckets. We will make use of these tools in the following
instructions. To install the Cloud SDK, [follow the installation instructions
here](https://cloud.google.com/sdk/downloads).

The final step in the installation process (`gcloud init`) will have you
authenticate via your web browser and select a default [zone and
region](https://cloud.google.com/compute/docs/regions-zones/regions-zones) for
your cloud resources, which you can choose based on your location and regional
hardware availability.

NOTE: Not all zones are equipped with GPUs, so if you want to use GPUs for your
project, please take note of the availability listing
[here](https://cloud.google.com/compute/docs/gpus/).

To verify that the installation and authentication succeeded, run

```shell
gcloud auth list
```

and verify that your account email address is printed.

## Starting a Compute Engine instance

A simple way to access compute on GCP is Google Compute Engine. Compute Engine
instances can be sized to meet computational and storage needs for your project.

Before we get started, [ensure you have adequate quota
provisioned](https://cloud.google.com/compute/quotas) so that you can get all
the CPUs/GPUs that you need. To start with, you might want to request quota for
64 CPUs and 2 GPUs in your zone.

DeepVariant can make use of multiple CPU cores and (currently, a single) GPU
device. For this "quick start" guide, let's allocate an 8-core non-preemptible
instance in your default zone with a single GPU, running Ubuntu 16.04, with a
disk of reasonable size for modest work with genomic data. From our local
command line, we do:

```shell
gcloud beta compute instances create "${USER}-deepvariant-quickstart" \
  --scopes "compute-rw,storage-full,cloud-platform"  \
  --image-family ubuntu-1604-lts --image-project ubuntu-os-cloud \
  --machine-type n1-standard-8  \
  --boot-disk-size=200GB \
  --zone us-west1-b \
  --accelerator type=nvidia-tesla-k80,count=1 --maintenance-policy TERMINATE --restart-on-failure
```

NOTE: To create an instance *without GPU*, simply omit the last line from the
command.

Check that the instance has been created and started:

```shell
gcloud compute instances list
```

which should produce output like:

```
NAME                    ZONE        MACHINE_TYPE    PREEMPTIBLE   INTERNAL_IP  EXTERNAL_IP     STATUS
[USER]-deepvariant-quickstart  us-west1-b  n1-standard-8                 10.138.0.4   35.185.203.59   RUNNING
```

Then connect to your instance via SSH:

```shell
gcloud compute ssh --zone us-west1-b "${USER}-deepvariant-quickstart"
```

You should land at a shell prompt in your new instance!

NOTE: All of these steps can also be completed from the Cloud Console, if you
prefer. Consult [this
guide](https://cloud.google.com/compute/docs/quickstart-linux), but be sure to
choose Ubuntu 16.04 as your image, as DeepVariant has not been tested on other
Linux distributions.

For more information about getting started with Compute Engine, see:

*   [Compute Engine instance creation in `gcloud` manual](https://cloud.google.com/sdk/gcloud/reference/compute/instances/create)
*   [Reference to machine sizes/types](https://cloud.google.com/compute/docs/machine-types)

## Download binaries, models, and test data

Before you start running, you need to have the following input files:

1.  A reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
    format and its corresponding index file (.fai). We'll refer to this as
    `${REF}` below.

1.  An aligned reads file in [BAM](http://genome.sph.umich.edu/wiki/BAM) format
    and its corresponding index file (.bai). We'll refer to this as `${BAM}`
    below. You get this by aligning the reads from a sequencing instrument,
    using an aligner like
    [BWA](https://academic.oup.com/bioinformatics/article/25/14/1754/225615/Fast-and-accurate-short-read-alignment-with)
    for example.

1.  A model checkpoint for DeepVariant. We'll refer to this as `${MODEL}` below.

### Download the DeepVariant binaries and install prerequisites

The DeepVariant binaries can be downloaded to your instance with the
[`gsutil`](https://cloud.google.com/storage/docs/gsutil) command.

```bash
BUCKET="gs://deepvariant"
BIN_VERSION="0.4.0"
MODEL_VERSION="0.4.0"
MODEL_CL="174375304"

# Note that we don't specify the CL number for the binary, only the bin version.
BIN_BUCKET="${BUCKET}/binaries/DeepVariant/${BIN_VERSION}/DeepVariant-${BIN_VERSION}+cl-*"
MODEL_NAME="DeepVariant-inception_v3-${MODEL_VERSION}+cl-${MODEL_CL}.data-wgs_standard"
MODEL_BUCKET="${BUCKET}/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}"
DATA_BUCKET="${BUCKET}/quickstart-testdata"

mkdir -p bin
# Download the DeepVariant binaries.
gsutil -m cp -P "${BIN_BUCKET}/*" bin/
```

Which should copy files into the bin subdirectory:

```bash
ls -1 bin/
```

producing

```
call_variants.zip
licenses.zip
make_examples.zip
model_eval.zip
model_train.zip
postprocess_variants.zip
run-prereq.sh
settings.sh
```

Install the prerequisites onto the machine:

```bash
cd bin; bash run-prereq.sh; cd -
```

### Download a model

Models file are stored in the shared Cloud Storage bucket:

`gs://deepvariant/models`

In this bucket models are organized into subdirectories by program name and
version, such as:

```
DeepVariant/0.4.0/DeepVariant-inception_v3-0.4.0+cl-12345.data-wgs_standard/
```

The model files are tagged with the program name and version, model name and
the data used to train the model. The CL number (google's code commit
identifier) may be safely ignored.

IMPORTANT: Models are tied to specific software versions. For example, you can
use model version 0.2.\* with any software version 0.2.\*. We recommend using
the latest software and the latest model.
See [Release Notes](deepvariant-release-notes.md).

Once you've selected an appropriate model directory, you can download it with
the [`gsutil`](https://cloud.google.com/storage/docs/gsutil) command. The path
to these model checkpoint files can then be provided to `call_variants`.

For example, let's download from the repository:

```bash
gsutil cp -R "${MODEL_BUCKET}" .
```

This should create a subdirectory in the current directory containing three
files:

```bash
ls -1 "${MODEL_NAME}/"
```

producing:

```
model.ckpt.data-00000-of-00001
model.ckpt.index
model.ckpt.meta
```

These files are in the standard [TensorFlow checkpoint format](https://www.tensorflow.org/programmers_guide/variables#checkpoint_files).

### Download test data

We've prepared a small test data bundle for use in this quick start guide that
can be downloaded to your instance with the
[`gsutil`](https://cloud.google.com/storage/docs/gsutil) command.

Download the test bundle:

```bash
gsutil cp -R "${DATA_BUCKET}" .
```

This should create a subdirectory in the current directory containing the actual
data files:

```bash
ls -1 quickstart-testdata/
```

outputting:

```
NA12878_S1.chr20.10_10p1mb.bam
NA12878_S1.chr20.10_10p1mb.bam.bai
test_nist.b37_chr20_100kbp_at_10mb.bed
test_nist.b37_chr20_100kbp_at_10mb.vcf.gz
test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi
ucsc.hg19.chr20.unittest.fasta
ucsc.hg19.chr20.unittest.fasta.fai
ucsc.hg19.chr20.unittest.fasta.gz
ucsc.hg19.chr20.unittest.fasta.gz.fai
ucsc.hg19.chr20.unittest.fasta.gz.gzi
```

## Run pipeline

DeepVariant consists of 3 main binaries: `make_examples`, `call_variants`, and
`postprocess_variants`. For this quick start guide we'll store the output in a
new directory on the instance and set up variables to refer to the test data we
downloaded above:

```bash
OUTPUT_DIR=quickstart-output
mkdir -p "${OUTPUT_DIR}"
REF=quickstart-testdata/ucsc.hg19.chr20.unittest.fasta
BAM=quickstart-testdata/NA12878_S1.chr20.10_10p1mb.bam
MODEL="${MODEL_NAME}/model.ckpt"
```

### `make_examples`

`make_examples` is the command used to extract pileup images from your BAM
files, encoding each as `tf.Example` (a kind of protocol buffer that TensorFlow
knows about) in "tfrecord" files. This tool is used as the first step of both
training and inference pipelines.

Here is an example command that would be used in inference (variant "calling")
mode:

```bash
python bin/make_examples.zip \
  --mode calling   \
  --ref "${REF}"   \
  --reads "${BAM}" \
  --regions "chr20:10,000,000-10,010,000" \
  --examples "${OUTPUT_DIR}/examples.tfrecord.gz"
```

If your machine has multiple cores, you can utilize them by running multiple
`make_examples` using the `--task` flag, which helps you split the input and
generates *sharded output*. Here is an example:

First install the GNU `parallel` tool to allow running multiple processes.

```bash
sudo apt-get -y install parallel
```

```bash
LOGDIR=./logs
N_SHARDS=3

mkdir -p "${LOGDIR}"
time seq 0 $((N_SHARDS-1)) | \
  parallel --eta --halt 2 --joblog "${LOGDIR}/log" --res "${LOGDIR}" \
  python bin/make_examples.zip \
    --mode calling \
    --ref "${REF}" \
    --reads "${BAM}" \
    --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
    --regions '"chr20:10,000,000-10,010,000"' \
    --task {}
```

To explain what *sharded output* is, for example, if `N_SHARDS` is 3, you'll get
three output files from the above command:

```bash
${OUTPUT_DIR}/examples.tfrecord-00000-of-00003.gz
${OUTPUT_DIR}/examples.tfrecord-00001-of-00003.gz
${OUTPUT_DIR}/examples.tfrecord-00002-of-00003.gz
```

In this example command above, we also used the
`--regions '"chr20:10,000,000-10,010,000"'` flag, which means we'll only process
10 kilobases of chromosome 20.

### `call_variants`

Once we have constructed the pileup images as "examples", we can then invoke the
variant calling tool to perform inference---identifying and labeling genomic
variants. We need to tell the `call_variants` command where to find the trained
machine learning model, using the `--checkpoint` argument.

Example command:

```bash
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/call_variants_output.tfrecord.gz"

python bin/call_variants.zip \
 --outfile "${CALL_VARIANTS_OUTPUT}" \
 --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
 --checkpoint "${MODEL}"
```

Notice that the output of this command is another "tfrecord.gz" file---this,
again, is serialized protocol buffer data.

### `postprocess_variants`

To convert the tfrecord output of `call_variants` into the
[VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format that is familiar
to bioinformaticists, we need to invoke the `postprocess_variants` tool.

An example command is below. Note that the output file should end with either
`.vcf` or `.vcf.gz`.

```bash
FINAL_OUTPUT_VCF="${OUTPUT_DIR}/output.vcf.gz"

python bin/postprocess_variants.zip \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "${FINAL_OUTPUT_VCF}"
```

## Evaluating the results

Here we use the `hap.py` ([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting 10 kilobase vcf file. This
serves as a quick check to ensure the three DeepVariant commands ran correctly.

```bash
sudo apt-get -y install docker.io

sudo docker pull pkrusche/hap.py
sudo docker run -it -v `pwd`:/data pkrusche/hap.py /opt/hap.py/bin/hap.py \
  /data/quickstart-testdata/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz \
  "/data/${FINAL_OUTPUT_VCF}" \
  -f /data/quickstart-testdata/test_nist.b37_chr20_100kbp_at_10mb.bed \
  -r "/data/${REF}" \
  -o "/data/${OUTPUT_DIR}/happy.output" \
  -l chr20:10000000-10010000
```

You should see output similar to the following.

```
Benchmarking Summary:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL            4         4         0           14         0         10      0              1                 1        0.714286                1                     NaN                     NaN                   0.333333                   1.333333
 INDEL   PASS            4         4         0           14         0         10      0              1                 1        0.714286                1                     NaN                     NaN                   0.333333                   1.333333
   SNP    ALL           44        44         0           60         0         16      0              1                 1        0.266667                1                     1.2                1.307692                   0.333333                   0.395349
   SNP   PASS           44        44         0           60         0         16      0              1                 1        0.266667                1                     1.2                1.307692                   0.333333                   0.395349
```
