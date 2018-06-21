# (BETA) Advanced Case Study: Train a customized SNP and small indel variant caller for BGISEQ-500 data.

redacted
not a standard workflow of DeepVariant itself. To clarify that "DeepVariant" as
a product is the pre-packaged variant caller, without needing the user to run
their own training. Make sure we don't confuse our users.

redacted

NOTE: This case study demonstrates an example of how to train a customized model
on one machine. This is NOT the fastest or cheapest configuration for your
needs. This document is in beta; we don't currently have a suggestion for a
production-grade pipeline for training.

## Request a machine

redacted

```shell
gcloud beta compute instances create "${USER}-training-casestudy"  \
--scopes "compute-rw,storage-full,cloud-platform" \
--image-family "ubuntu-1604-lts" \
--image-project "ubuntu-os-cloud" \
--machine-type "custom-64-131072" \
--boot-disk-size "300" \
--boot-disk-type "pd-ssd" \
--zone "us-west1-b"
```

Set the variables:

```
BASE="${HOME}/training-case-study"
# redacted
DATA_BUCKET=gs://brain-genomics/pichuan/BGISEQ-HG001

INPUT_DIR="${BASE}/input"

BIN_DIR="${INPUT_DIR}/bin"
DATA_DIR="${INPUT_DIR}/data"

REF="${DATA_DIR}/ucsc_hg19.fa"
BAM="${DATA_DIR}/BGISEQ_PE100_NA12878.sorted.bam"
TRUTH_VCF="${DATA_DIR}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer_chrs_FIXED.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_chr.bed"

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"

LOG_DIR="${OUTPUT_DIR}/logs"

# redacted
OUTPUT_BUCKET=gs://brain-genomics/pichuan/customized_training
```

Create directories:

```
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${LOG_DIR}"
```

Copy data

```
gsutil -m cp ${DATA_BUCKET}/BGISEQ_PE100_NA12878.sorted.bam* "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/ucsc_hg19.fa*" "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_*" "${DATA_DIR}"

gunzip "${DATA_BUCKET}/ucsc_hg19.fa.gz"
```

Download extra packages

```
sudo apt-get -y install parallel
```

## Run make_examples in “training” mode to create training, validation, and test sets.

Create examples in "training" mode (which means these `tensorflow.Example`s will
contain a `label` field).

In this tutorial, we create examples on one replicate of HG001 sequenced by
BGISEQ-500 provided on the Genome In a Botton FTP site
[[readme](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/BGISEQ500/standard_library/readme.txt)].

We will create examples in 3 different sets: Training set (everything except for
chr20, 21, and 22), validation set (chr21 and 22), and test sets (chr20).

For the definition of these 3 sets in commonly used machine learning
terminology, please refer to [Machine Learning
Glossary](https://developers.google.com/machine-learning/crash-course/glossary).

### Training set

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    ${BIN_DIR}/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${OUTPUT_DIR}/training_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --exclude_regions "'chr20 chr21 chr22'" \
      --use_fast_pass_aligner \
) >"${LOG_DIR}/training_set.with_label.make_examples.log" 2>&1
```

This took 136m23.478s. We will want to shuffle this on Dataflow later, so I will
copy it to GCS bucket first:

```
gsutil -m cp ${OUTPUT_DIR}/training_set.with_label.tfrecord-?????-of-00064.gz \
  ${OUTPUT_BUCKET}
```

### Validation set

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    ${BIN_DIR}/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${OUTPUT_DIR}/validation_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "'chr21 chr22'" \
      --use_fast_pass_aligner \
) >"${LOG_DIR}/validation_set.with_label.make_examples.log" 2>&1
```

This took 9m15.127s

Validation set is small here. We will just shuffle locally later, so no need to
copy to out GCS bucket.

### Test set

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    ${BIN_DIR}/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${OUTPUT_DIR}/test_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "chr20" \
      --use_fast_pass_aligner \
) >"${LOG_DIR}/test_set.with_label.make_examples.log" 2>&1
```

This took: 6m30.693s

Test set is small here. We will just shuffle locally later, so no need to copy
to out GCS bucket.

## Shuffle each set of examples and generate a data configuration file for each.

Shuffling the `tensorflow.Example`s is an important step for training a model.
In our training logic, we decided to not shuffle our examples within TensorFlow,
but rely a preprocessing step that shuffles the examples globally.

Here we provide an example of running on [Cloud Dataflow
Runner](https://beam.apache.org/documentation/runners/dataflow/) as well locally
using [DirectRunner](https://beam.apache.org/documentation/runners/direct/).
Beam can also use other runners, such as [Spark
Runner](https://beam.apache.org/documentation/runners/spark/).

First, install beam on your machine following the instructions at
https://beam.apache.org/get-started/quickstart-py/

Here is an example. You might or might not need to install everything below:

```
sudo apt-get -y install python-dev python-pip
pip install --upgrade pip
pip install --upgrade virtualenv
```

A virtual environment is a directory tree containing its own Python
distribution. To create a virtual environment, create a directory and run:

```
virtualenv ${HOME}/virtualenv_beam
```

A virtual environment needs to be activated for each shell that is to use it.
Activating it sets some environment variables that point to the virtual
environment’s directories.

To activate a virtual environment in Bash, run:

```
. ${HOME}/virtualenv_beam/bin/activate
```

Once this is activated, install Beam:

```
pip install apache-beam
```

Validation and test sets are small, so we will just shuffle locally using
DirectRunner:

```
time python $HOME/shuffle_tfrecords_beam.py \
  --input_pattern_list="${OUTPUT_DIR}/validation_set.with_label.tfrecord-?????-of-00064.gz" \
  --output_pattern="${OUTPUT_DIR}/validation_set.with_label.shuffled.tfrecord.gz" \
  --output_dataset_name="HG001" \
  --runner=DirectRunner
```

This took 10m57.418s

```
time python $HOME/shuffle_tfrecords_beam.py \
  --input_pattern_list="${OUTPUT_DIR}/test_set.with_label.tfrecord-?????-of-00064.gz" \
  --output_pattern="${OUTPUT_DIR}/test_set.with_label.shuffled.tfrecord.gz" \
  --output_dataset_name="HG001" \
  --runner=DirectRunner
```

This took 8m52.594s

For the training set, it is too large to be running with DirectRunner on this
instance, so we use the DataflowRunner. Before that, please make sure you enable
Dataflow API for your project:
http://console.cloud.google.com/flows/enableapi?apiid=dataflow.

Then, install Dataflow:

```
pip install google-cloud-dataflow
```

Shuffle using Dataflow.

```
time python ${HOME}/shuffle_tfrecords_beam.py \
  --project=YOUR_PROJECT \
  --input_pattern_list="${OUTPUT_BUCKET}"/training_set.with_label.tfrecord-?????-of-00064.gz \
  --output_pattern="${OUTPUT_BUCKET}/training_set.with_label.shuffled.tfrecord.gz" \
  --output_dataset_name="HG001" \
  --output_dataset_config_pbtxt="${OUTPUT_BUCKET}/output.data_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DataflowRunner --project=google.com:brain-genomics \
  --staging_location=$OUTPUT_BUCKET/staging --temp_location=$OUTPUT_BUCKET/tempdir \
  --save_main_session
```

Then, you should be able to see the run on:
https://console.cloud.google.com/dataflow?project=<YOUR_PROJECT_NAME>

Here is an example of my run:

![Dataflow](dataflow.png?raw=true "Dataflow shuffle-examples")

In order to have the best performance, you might need extra resources such as
machines or IPs within a region. That will not be in the scope of this case
study here.

My run took took about 41 min on Dataflow. The output is in
`${OUTPUT_BUCKET}/training_set.with_label.shuffled.tfrecord.gz*`. In this case,
it wrote to 356 shards:
`${OUTPUT_BUCKET}/training_set.with_label.shuffled.tfrecord.gz-?????-of-00356`

redacted

1.  I need to change shuffle_tfrecords_beam.py to actually write on the
    output.data_config.pbtxt file.
1.  Then, proceed with next step of training.
