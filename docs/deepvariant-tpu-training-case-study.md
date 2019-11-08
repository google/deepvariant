# Advanced Case Study: Train a customized SNP and small indel variant caller for BGISEQ-500 data.

DeepVariant is an analysis pipeline that uses a deep neural network to call
genetic variants from next-generation DNA sequencing (NGS) data. While
DeepVariant is highly accurate for [many types of NGS
data](https://doi.org/10.1101/092890), some users may be interested in training
custom deep learning models that have been optimized for very specific data.

This case study describes one way to train such a custom model using a TPU, in
this case for BGISEQ-500 data.

Please note that there is not yet a production-grade training pipeline. This is
just one example of how to train a custom model, and is neither the fastest nor
the cheapest possible configuration. The resulting model also does not represent
the greatest achievable accuracy for BGISEQ-500 data.

## High level summary of result

We demonstrated that by training on 1 replicate of BGISEQ-500 whole genome data
(everything except for chromosome 20-22), we can significantly improve the
accuracy comparing to the WGS model as a baseline: Indel F1 94.39% --> 98.22%;
SNP F1: 99.85% --> 99.90%.

Training for 50k steps took about 1.5 hours on Cloud TPU.

This tutorial is meant as an example for training; all the other processing in
this tutorial were done serially with no pipeline optimization.

## Request a machine

Any sufficiently capable machine will do. For this case study, we used a 64-core
non-preemptible instance with 128GiB and no GPU.

Here is an example on how to get a [CPU machine].

Set the variables:

```
YOUR_PROJECT=REPLACE_WITH_YOUR_PROJECT
OUTPUT_GCS_BUCKET=REPLACE_WITH_YOUR_GCS_BUCKET

BUCKET="gs://deepvariant"
BIN_VERSION="0.9.0"

MODEL_BUCKET="${BUCKET}/models/DeepVariant/${BIN_VERSION}/DeepVariant-inception_v3-${BIN_VERSION}+data-wgs_standard"
GCS_PRETRAINED_WGS_MODEL="${MODEL_BUCKET}/model.ckpt"

OUTPUT_BUCKET="${OUTPUT_GCS_BUCKET}/customized_training"
TRAINING_DIR="${OUTPUT_BUCKET}/training_dir"

BASE="${HOME}/training-case-study"
DATA_BUCKET=gs://deepvariant/training-case-study/BGISEQ-HG001

INPUT_DIR="${BASE}/input"
BIN_DIR="${INPUT_DIR}/bin"
DATA_DIR="${INPUT_DIR}/data"
OUTPUT_DIR="${BASE}/output"
LOG_DIR="${OUTPUT_DIR}/logs"
SHUFFLE_SCRIPT_DIR="${HOME}/deepvariant/tools"

REF="${DATA_DIR}/ucsc_hg19.fa"
BAM="${DATA_DIR}/BGISEQ_PE100_NA12878.sorted.bam"
TRUTH_VCF="${DATA_DIR}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer_chrs_FIXED.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_chr.bed"

N_SHARDS="64"
```

## Download binaries, models, and data

### Create directories:

```
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BIN_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${LOG_DIR}"
```

### Copy data

```
gsutil -m cp ${DATA_BUCKET}/BGISEQ_PE100_NA12878.sorted.bam* "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/ucsc_hg19.fa*" "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_*" "${DATA_DIR}"

gunzip "${DATA_DIR}/ucsc_hg19.fa.gz"
```

### Download extra packages

```
sudo apt -y update
sudo apt -y install parallel
sudo apt -y install docker.io
```

## Run make_examples in “training” mode for training and validation sets.

Create examples in "training" mode (which means these `tensorflow.Example`s will
contain a `label` field).

In this tutorial, we create examples on one replicate of HG001 sequenced by
BGISEQ-500 provided on the
[Genome In a Bottle FTP site](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/BGISEQ500/standard_library/readme.txt).

We will create examples in 2 different sets: Training set (everything except for
chr20, 21, and 22), validation set (chr21 and 22) - These 2 sets will be used in
`model_train` and `model_eval`, so we'll create them in "training" mode so they
have the real labels. We will leave out chr20 for final evaluation for our
trained model at the end.

For the definition of these 3 sets in commonly used machine learning
terminology, please refer to
[Machine Learning Glossary](https://developers.google.com/machine-learning/glossary/).

### Training set

First, to set up,

```
sudo docker pull google/deepvariant:"${BIN_VERSION}"
```

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    sudo docker run \
      -v ${HOME}:${HOME} \
      google/deepvariant:"${BIN_VERSION}" \
      /opt/deepvariant/bin/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${OUTPUT_DIR}/training_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --exclude_regions "'chr20 chr21 chr22'" \
) >"${LOG_DIR}/training_set.with_label.make_examples.log" 2>&1
```

This took ~98min. We will want to shuffle this on Dataflow later, so we copy the
data to GCS bucket first:

```
gsutil -m cp ${OUTPUT_DIR}/training_set.with_label.tfrecord-?????-of-00064.gz \
  ${OUTPUT_BUCKET}
```

### Validation set

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    sudo docker run \
      -v /home/${USER}:/home/${USER} \
      google/deepvariant:"${BIN_VERSION}" \
      /opt/deepvariant/bin/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${OUTPUT_DIR}/validation_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "'chr21 chr22'" \
) >"${LOG_DIR}/validation_set.with_label.make_examples.log" 2>&1
```

This took: ~10min.

Validation set is small here. We will just shuffle locally later, so no need to
copy to out GCS bucket.

## Shuffle each set of examples and generate a data configuration file for each.

Shuffling the `tensorflow.Example`s is an important step for training a model.
In our training logic, we decided to not shuffle our examples within TensorFlow,
but rely a preprocessing step that shuffles the examples globally.

First, if you have run this step before, and want to rerun it, you might want to
consider cleaning up previous data first to avoid confusion:

```
# (Optional) Clean up existing files.
gsutil -m rm -f "${OUTPUT_BUCKET}/training_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
gsutil rm -f "${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt"
```

Here we provide an example of running on [Cloud Dataflow
Runner](https://beam.apache.org/documentation/runners/dataflow/) as well locally
using [DirectRunner](https://beam.apache.org/documentation/runners/direct/).
Beam can also use other runners, such as [Spark
Runner](https://beam.apache.org/documentation/runners/spark/).

First, install beam on your machine following the instructions at
https://beam.apache.org/get-started/quickstart-py/

Here is an example. You might or might not need to install everything below:

```
sudo apt -y install python-dev python-pip
# This will make sure the pip command is bound to python2
python2 -m pip install --user --upgrade --force-reinstall pip

export PATH="$HOME/.local/bin:$PATH"
```

Install Beam:

```
pip install --user apache-beam
```

Get the code that shuffles:

```
mkdir -p ${SHUFFLE_SCRIPT_DIR}
wget https://raw.githubusercontent.com/google/deepvariant/r0.8/tools/shuffle_tfrecords_beam.py -O ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py
```

Validation set is small, so we will just shuffle locally using DirectRunner:

```
# First, clean up existing files.
rm -f "${OUTPUT_DIR}/validation_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
rm -f "${OUTPUT_DIR}/validation_set.dataset_config.pbtxt"


time python ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --input_pattern_list="${OUTPUT_DIR}/validation_set.with_label.tfrecord-?????-of-00064.gz" \
  --output_pattern_prefix="${OUTPUT_DIR}/validation_set.with_label.shuffled" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/validation_set.dataset_config.pbtxt" \
  --output_dataset_name="HG001" \
  --runner=DirectRunner
```

Output is in
`${OUTPUT_DIR}/validation_set.with_label.shuffled-00000-of-00001.tfrecord.gz`

Data config file is in `${OUTPUT_DIR}/validation_set.dataset_config.pbtxt`.

This took ~13min.

For the training set, it is too large to be running with DirectRunner on this
instance, so we use the DataflowRunner. Before that, please make sure you enable
Dataflow API for your project:
http://console.cloud.google.com/flows/enableapi?apiid=dataflow.

Then, install Dataflow:

```
pip install --user google-cloud-dataflow
```

Shuffle using Dataflow.

```
time python ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="${YOUR_PROJECT}" \
  --input_pattern_list="${OUTPUT_BUCKET}"/training_set.with_label.tfrecord-?????-of-00064.gz \
  --output_pattern_prefix="${OUTPUT_BUCKET}/training_set.with_label.shuffled" \
  --output_dataset_name="HG001" \
  --output_dataset_config_pbtxt="${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DataflowRunner \
  --staging_location="${OUTPUT_BUCKET}/staging" \
  --temp_location="${OUTPUT_BUCKET}/tempdir" \
  --save_main_session
```

Then, you should be able to see the run on:
https://console.cloud.google.com/dataflow?project=YOUR_PROJECT

Here is an example of my run:

![Dataflow](images/dataflow.png?raw=true "Dataflow shuffle-examples")

In order to have the best performance, you might need extra resources such as
machines or IPs within a region. That will not be in the scope of this case
study here.

In one test run, this took about 81min on Dataflow.

The output path can be found in the dataset_config file by:

```
gsutil cat "${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt"
```

In the output, the `tfrecord_path` should be valid paths in gs://.

```
# Generated by shuffle_tfrecords_beam.py
#
# --input_pattern_list=YOUR_GCS_BUCKET/training_set.with_label.tfrecord-?????-of-00064.gz
# --output_pattern_prefix=YOUR_GCS_BUCKET/training_set.with_label.shuffled
#

name: "HG001"
tfrecord_path: "YOUR_GCS_BUCKET/training_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
num_examples: 4157931
```

In one test run, it wrote to 365 shards:
`${OUTPUT_BUCKET}/training_set.with_label.shuffled-?????-of-03456.tfrecord.gz`

### Start a Cloud TPU

Before you start, please read about
[Cloud TPU pricing](https://cloud.google.com/tpu/docs/pricing) and consider the
expenses for using TPUs.

In this tutorial, we provide one way to start and use a TPU. This might not be
the easiest way to use TPU. Please refer to the [Cloud TPU documentation] for
more information on Cloud TPUs.

First, check all existing TPUs by running this command:

```
gcloud compute tpus list --zone=us-central1-c
```

In the test run, we didn't have any other existing TPUs.
Then, run the following command to start a TPU:

```
time gcloud compute tpus create ${USER}-demo-tpu \
  --network=default \
  --range=10.240.1.0/29 \
  --version=1.13 \
  --zone=us-central1-c
```

Once you're done with the TPU, you should remember to delete it.

This command took about 3.5min to finish.

After the TPU is created, we can query it by:

```
gcloud compute tpus list --zone=us-central1-c
```

It might look like:

```
NAME              ZONE           ACCELERATOR_TYPE  NETWORK_ENDPOINTS  NETWORK  RANGE          STATUS
USER-demo-tpu  us-central1-c  v2-8              10.240.1.2:8470    default  10.240.1.0/29  READY
```

In this example, set up these variables:

```
export TPU_NAME="${USER}-demo-tpu"
export TPU_IP="10.240.1.2"
```

(One more reminder about
[Cloud TPU pricing](https://cloud.google.com/tpu/docs/pricing) - when you're
done using the TPU, make sure to delete it otherwise you'll continue to be
billed for the TPU.)

### Start `model_train` and `model_eval`

NOTE: all parameters below are used as an example. They were not optimized for
this dataset, and we do not recommend these as the best default either.

```
( time sudo docker run \
  -v /home/${USER}:/home/${USER} \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/model_train \
  --use_tpu \
  --master="grpc://${TPU_IP}:8470" \
  --dataset_config_pbtxt="${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt" \
  --train_dir="${TRAINING_DIR}" \
  --model_name="inception_v3" \
  --number_of_steps=50000 \
  --save_interval_secs=300 \
  --batch_size=512 \
  --learning_rate=0.008 \
  --start_from_checkpoint="${GCS_PRETRAINED_WGS_MODEL}" \
) > "${LOG_DIR}/train.log" 2>&1 &
```

Pointers for common issues or things you can tune:

1.  TPU might not have write access to GCS bucket:

    https://cloud.google.com/tpu/docs/storage-buckets#storage_access

1.  Change `save_interval_secs` to save checkpoints more frequently:

    As you train, you'll notice that not every checkpoint gets saved. If you
    want model checkpoints to be saved more frequently, one flag you can
    consider providing is `save_interval_secs`. With the current default (300)
    when I ran this case study, the checkpoints were saved every 5 minutes. If
    this doesn't provide enough granularity, you can decrease the number of
    seconds, then you'll have more saved checkpoints over time. The tradeoff is
    the space needed on your GCS buckets to keep these checkpoints around.

1.  If you have already run this command before, you might see this in your log:
    "Skipping training since max_steps has already saved."

    This is because the training checkpoints from previous runs are still in
    `${TRAINING_DIR}`. If you're starting a new run, point that variable to a
    new directory, or clean up previous training results.

At the same time, start `model_eval` on CPUs:

```
sudo docker run \
  -v /home/${USER}:/home/${USER} \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/model_eval \
  --dataset_config_pbtxt="${OUTPUT_DIR}/validation_set.dataset_config.pbtxt" \
  --checkpoint_dir="${TRAINING_DIR}" \
  --batch_size=512 > "${LOG_DIR}/eval.log" 2>&1 &
```

`model_eval` will watch the `${TRAINING_DIR}` and start evaluating when there
are newly saved checkpoints. It evaluates the checkpoints on the data specified
in `validation_set.dataset_config.pbtxt`, and saves `*metrics` file to the
directory. These files are used later to pick the best model based on how
accurate they are on the validation set.

When I ran this case study, running `model_eval` on CPUs is fast enough because
`model_train` didn't save checkpoints too frequently.

In my run, `model_train` took about 1.5hr to finish 50k steps (with
batch_size 512). Note that `model_eval` will not stop on its own, so I had to
kill the process after training is no longer producing more checkpoints. And,
**remember to delete the TPU once you're done with training.** You can use this
command to delete the TPU:

```
gcloud compute tpus delete ${TPU_NAME} --zone us-central1-c
```

### Use TensorBoard to visualize progress

You'll want to let model_train and model_eval run for a while before you start a
TensorBoard. (You can start a TensorBoard immediately, but you just won't see
the metrics summary until later.)

We can start a TensorBoard to visualize the progress of training better. We did
this through a Google Cloud Shell from https://console.cloud.google.com , on the
top right:

![Shell](images/ActivateShell.png?raw=true "Activate Google Cloud Shell")

This opens up a terminal at the bottom of the browser page, then run:

```
# Change to your OUTPUT_BUCKET from earlier.
OUTPUT_BUCKET="${OUTPUT_GCS_BUCKET}/customized_training"
TRAINING_DIR="${OUTPUT_BUCKET}/training_dir"
tensorboard --logdir ${TRAINING_DIR} --port=8080
```

This gives some message like:

```
TensorBoard 2.0.1 at http://localhost:8080/ (Press CTRL+C to quit)
```

But that link is not usable directly. I clicked on the “Web Preview” on the top
right of the mini terminal:

![WebPreview](images/WebPreview.png?raw=true "Web Preview")

And clicked on "Preview on port 8080":

![PreviewOnPort](images/PreviewOnPort.png?raw=true "Preview on Port 8080")

Once it starts, you can see many metrics, including accuracy, speed, etc. You
will need to wait for both `model_train` and `model_eval` to run for a while
before the plots will make more sense.

### Clean up the TPU once you're done training

When you are done with training, make sure to clean up the TPU:

```
gcloud compute tpus delete ${TPU_NAME} --zone us-central1-c
```

### Pick a model

You can directly look up the best checkpoint by running:

```
gsutil cat "${TRAINING_DIR}"/best_checkpoint.txt
```

In my run, this showed that the model checkpoint that performs the best on the
validation set was `${TRAINING_DIR}/model.ckpt-36200`.
Based on this result, a few thoughts came into mind:

1.  Training more steps didn't seem to help much. Did the training overfit?
1.  Currently we save checkpoints every 5 minutes (set by `save_interval_secs`).
    We could save checkpoints more frequently so we can observe this curve with
    finer granularity.

But for now, let's use this model to do the final evaluation on the test set and
see how we do. We can use the one-step command to call:

```
sudo docker run \
  -v /home/${USER}:/home/${USER} \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --customized_model "${TRAINING_DIR}/model.ckpt-36200" \
  --ref "${REF}" \
  --reads "${BAM}" \
  --regions "chr20" \
  --output_vcf "${OUTPUT_DIR}/test_set.vcf.gz" \
  --num_shards=${N_SHARDS}
```

Once this is done, we have the final callset in VCF
format here: `${OUTPUT_DIR}/test_set.vcf.gz`. Next step is to run `hap.py` to
complete the evaluation on chromosome 20:

```
sudo docker pull pkrusche/hap.py

time sudo docker run -it \
-v "${DATA_DIR}:${DATA_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${OUTPUT_DIR}/test_set.vcf.gz" \
  -f "${TRUTH_BED}" \
  -r "${REF}" \
  -o "${OUTPUT_DIR}/chr20-calling.happy.output" \
  -l chr20 \
  --engine=vcfeval
```

This takes about 3 minutes. The output of `hap.py` is here:

```
[I] Total VCF records:         3775119
[I] Non-reference VCF records: 3775119
[I] Total VCF records:         130929
[I] Non-reference VCF records: 95922
Benchmarking Summary:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL        10023      9805       218        19078       142       8728     97       0.978250           0.98628        0.457490         0.982249                     NaN                     NaN                   1.547658                   2.058959
 INDEL   PASS        10023      9805       218        19078       142       8728     97       0.978250           0.98628        0.457490         0.982249                     NaN                     NaN                   1.547658                   2.058959
   SNP    ALL        66237     66160        77        77821        57      11568     15       0.998838           0.99914        0.148649         0.998989                2.284397                2.203185                   1.700387                   1.784817
   SNP   PASS        66237     66160        77        77821        57      11568     15       0.998838           0.99914        0.148649         0.998989                2.284397                2.203185                   1.700387                   1.784817
```

To summarize, the accuracy is:

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 218  | 142  | 0.978250 | 0.986280  | 0.982249
SNP   | 77   | 57   | 0.998838 | 0.999140  | 0.998989

The baseline we're comparing to is to directly use the WGS model to make the
calls, using this command:

```bash
sudo docker run \
  -v "${DATA_DIR}:${DATA_DIR}" \
  -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="${REF}" \
  --reads="${BAM}" \
  --output_vcf="${OUTPUT_DIR}/baseline.vcf.gz" \
  --num_shards=${N_SHARDS}
```

Baseline:

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 426  | 745  | 0.957498 | 0.930659  | 0.943888
SNP   | 139  | 50   | 0.997901 | 0.999094  | 0.998497

### Additional things to try

#### Parameters to tune

Starting from the default setting of this tutorial is a good starting point, but
this training case study is by no means the best setting. Training is both a
science and an art. There are many knobs that we could potentially tune. Users
might be able to use different parameters to train a more accurate model even
with the same data, such as `batch_size`, `learning_rate`,
`learning_rate_decay_factor` in modeling.py.

#### Downsampling the BAM file to generate more training examples

When generating the training set, we can make some adjustment to create more
training data. For example, when we train the released WGS model for
DeepVariant, for each BAM file, we created an extra set of training examples
using `--downsample_fraction=0.5`, which downsamples the reads and creates
training examples with lower coverage. We found that this makes the trained
model more robust.

[CPU machine]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
[Cloud TPU documentation]: https://cloud.google.com/tpu/docs/
