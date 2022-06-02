# Advanced Case Study: Train a customized SNP and small indel variant caller for BGISEQ-500 data.

DeepVariant is an analysis pipeline that uses a deep neural network to call
genetic variants from next-generation DNA sequencing (NGS) data. While
DeepVariant is highly accurate for [many types of NGS
data](https://rdcu.be/7Dhl), some users may be interested in training
custom deep learning models that have been optimized for very specific data.

This case study describes one way to train such a custom model using a GPU, in
this case for BGISEQ-500 data.

Please note that there is not yet a production-grade training pipeline. This is
just one example of how to train a custom model, and is neither the fastest nor
the cheapest possible configuration. The resulting model also does not represent
the greatest achievable accuracy for BGISEQ-500 data.

## High level summary of result

We demonstrated that by training on 1 replicate of BGISEQ-500 whole genome data
(everything except for chromosome 20-22), we can significantly improve the
accuracy comparing to the WGS model as a baseline:

* Indel F1 :93.4908% --> 98.1305%;
*   SNP F1: 99.8838% --> 99.9011%.

Training for 50,000 steps took about 1.5 hours on 1 GPU. Currently we cannot
train on multiple GPUs.

This tutorial is meant as an example for training; all the other processing in
this tutorial were done serially with no pipeline optimization.

## Request a machine

For this case study, we use a [GPU machine] with 16 vCPUs.

Set the variables:

```
YOUR_PROJECT=REPLACE_WITH_YOUR_PROJECT
OUTPUT_GCS_BUCKET=REPLACE_WITH_YOUR_GCS_BUCKET

BUCKET="gs://deepvariant"
BIN_VERSION="1.4.0"

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
BAM_CHR1="${DATA_DIR}/BGISEQ_PE100_NA12878.sorted.chr1.bam"
BAM_CHR20="${DATA_DIR}/BGISEQ_PE100_NA12878.sorted.chr20.bam"
BAM_CHR21="${DATA_DIR}/BGISEQ_PE100_NA12878.sorted.chr21.bam"
TRUTH_VCF="${DATA_DIR}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer_chrs_FIXED.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_chr.bed"

N_SHARDS=16
```

## Download binaries and data

### Create directories:

```
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BIN_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${LOG_DIR}"
```

### Copy data

```
gsutil -m cp ${DATA_BUCKET}/BGISEQ_PE100_NA12878.sorted.chr*.bam* "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/ucsc_hg19.fa*" "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_*" "${DATA_DIR}"

gunzip "${DATA_DIR}/ucsc_hg19.fa.gz"
```

### Download extra packages

```
sudo apt -y update
sudo apt -y install parallel
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.4/scripts/install_nvidia_docker.sh
bash -x install_nvidia_docker.sh
```

## Run make_examples in “training” mode for training and validation sets.

Create examples in "training" mode (which means these `tensorflow.Example`s will
contain a `label` field).

In this tutorial, we create examples on one replicate of HG001 sequenced by
BGISEQ-500 provided on the
[Genome In a Bottle FTP site](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/BGISEQ500/standard_library/readme.txt).

In this tutorial, we show how to create examples in 2 different sets:
Training set (chr1), validation set (chr21) - These 2 sets are used in
`model_train` and `model_eval`, so we create them in "training" mode so they
have the real labels. We use chr20 for final evaluation for our trained model at
the end.

For the definition of these 3 sets in commonly used machine learning
terminology, please refer to
[Machine Learning Glossary](https://developers.google.com/machine-learning/glossary/).

### Training set

First, to set up,

```
sudo docker pull google/deepvariant:"${BIN_VERSION}"
sudo docker pull google/deepvariant:"${BIN_VERSION}-gpu"
```

`make_examples` step doesn't use GPU:

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --line-buffer \
    sudo docker run \
      -v ${HOME}:${HOME} \
      google/deepvariant:"${BIN_VERSION}" \
      /opt/deepvariant/bin/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM_CHR1}" \
      --examples "${OUTPUT_DIR}/training_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "'chr1'" \
      --channels "insert_size" \
) 2>&1 | tee "${LOG_DIR}/training_set.with_label.make_examples.log"
```

This took about 22min.

Starting in v1.4.0, we added an extra channel in our WGS setting using the
`--channels "insert_size"` flag. And, the make_examples step creates
`*.example_info.json` files. For example, you can see it here:

```
$ cat "${OUTPUT_DIR}/training_set.with_label.tfrecord-00000-of-00016.gz.example_info.json"
{"version": "1.4.0", "shape": [100, 221, 7], "channels": [1, 2, 3, 4, 5, 6, 19]}
```

Depending on your data type, you might want to tweak the flags for the
`make_examples` step, which can result in different shape of the output
examples.

We will want to shuffle this on Dataflow later, so we copy the data to GCS
bucket first:

```
gsutil -m cp ${OUTPUT_DIR}/training_set.with_label.tfrecord-?????-of-00016.gz* \
  ${OUTPUT_BUCKET}
```

NOTE: If you prefer shuffling locally, please take a look at this user-provided
shuffler option:
https://github.com/google/deepvariant/issues/360#issuecomment-1019990366

### Validation set

```
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --line-buffer \
    sudo docker run \
      -v /home/${USER}:/home/${USER} \
      google/deepvariant:"${BIN_VERSION}" \
      /opt/deepvariant/bin/make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM_CHR21}" \
      --examples "${OUTPUT_DIR}/validation_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "'chr21'" \
      --channels "insert_size" \
) 2>&1 | tee "${LOG_DIR}/validation_set.with_label.make_examples.log"
```

This took: ~5min.

Copy to GCS bucket:

```
gsutil -m cp ${OUTPUT_DIR}/validation_set.with_label.tfrecord-?????-of-00016.gz* \
  ${OUTPUT_BUCKET}
```

## Shuffle each set of examples and generate a data configuration file for each.

Shuffling the `tensorflow.Example`s is an important step for training a model.
In our training logic, we shuffle examples globally using a preprocessing step.

First, if you have run this step before, and want to rerun it, you might want to
consider cleaning up previous data first to avoid confusion:

```
# (Optional) Clean up existing files.
gsutil -m rm -f "${OUTPUT_BUCKET}/training_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
gsutil rm -f "${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt"
gsutil -m rm -f "${OUTPUT_BUCKET}/validation_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
gsutil rm -f "${OUTPUT_BUCKET}/validation_set.dataset_config.pbtxt"
gsutil rm -f "${OUTPUT_BUCKET}/example_info.json"
```

Here we provide examples for running on [Cloud Dataflow
Runner](https://beam.apache.org/documentation/runners/dataflow/) and also
[DirectRunner](https://beam.apache.org/documentation/runners/direct/).
Beam can also use other runners, such as [Spark
Runner](https://beam.apache.org/documentation/runners/spark/).

First, activate a virtual environment to install beam on your machine following
the instructions at https://beam.apache.org/get-started/quickstart-py/.

Then, get the code that shuffles:

```
mkdir -p ${SHUFFLE_SCRIPT_DIR}
wget https://raw.githubusercontent.com/google/deepvariant/r1.4/tools/shuffle_tfrecords_beam.py -O ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py
```

Next, we shuffle the data using DataflowRunner. Before that, please make sure
you enable Dataflow API for your project:
http://console.cloud.google.com/flows/enableapi?apiid=dataflow.

To access `gs://` path, make sure you run this in your virtual environment:

```
sudo apt -y update && sudo apt -y install python3-pip
pip3 install --upgrade pip
pip3 install setuptools --upgrade
pip3 install apache_beam[gcp]
pip3 install tensorflow  # For parsing tf.Example in shuffle_tfrecords_beam.py.
```

Shuffle using Dataflow.

```
time python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="${YOUR_PROJECT}" \
  --input_pattern_list="${OUTPUT_BUCKET}"/training_set.with_label.tfrecord-?????-of-00016.gz \
  --output_pattern_prefix="${OUTPUT_BUCKET}/training_set.with_label.shuffled" \
  --output_dataset_name="HG001" \
  --output_dataset_config_pbtxt="${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DataflowRunner \
  --staging_location="${OUTPUT_BUCKET}/staging" \
  --temp_location="${OUTPUT_BUCKET}/tempdir" \
  --save_main_session \
  --region us-east1
```

Then, you should be able to see the run on:
https://console.cloud.google.com/dataflow?project=YOUR_PROJECT

In order to have the best performance, you might need extra resources such as
machines or IPs within a region. That will not be in the scope of this case
study here.

The output path can be found in the dataset_config file by:

```
gsutil cat "${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt"
```

In the output, the `tfrecord_path` should be valid paths in gs://.

```
# Generated by shuffle_tfrecords_beam.py
# class2: 124564
# class1: 173668
# class0: 44526
#
# --input_pattern_list=OUTPUT_BUCKET/training_set.with_label.tfrecord-?????-of-00016.gz
# --output_pattern_prefix=OUTPUT_BUCKET/training_set.with_label.shuffled
#

name: "HG001"
tfrecord_path: "OUTPUT_GCS_BUCKET/training_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
num_examples: 342758
```

We can shuffle the validation set locally using
[DirectRunner](https://beam.apache.org/documentation/runners/direct/). Adding
`--direct_num_workers=0` sets the number of threads/subprocess to the number
of cores of the machine where the pipeline is running.

```
time python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="${YOUR_PROJECT}" \
  --input_pattern_list="${OUTPUT_DIR}"/validation_set.with_label.tfrecord-?????-of-00016.gz \
  --output_pattern_prefix="${OUTPUT_DIR}/validation_set.with_label.shuffled" \
  --output_dataset_name="HG001" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/validation_set.dataset_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=0
```

Here is the validation_set:

```
cat "${OUTPUT_DIR}/validation_set.dataset_config.pbtxt"
```

```
# Generated by shuffle_tfrecords_beam.py
# class0: 5595
# class1: 31852
# class2: 21954
#
# --input_pattern_list=OUTPUT_DIR/validation_set.with_label.tfrecord-?????-of-00016.gz
# --output_pattern_prefix=OUTPUT_DIR/validation_set.with_label.shuffled
#

name: "HG001"
tfrecord_path: "OUTPUT_DIR/validation_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
num_examples: 59401
```

### Start `model_train` and `model_eval`

NOTE: all parameters below are used as an example. They are not optimized for
this dataset, and are not recommended as the best default either.

```
( time sudo docker run --gpus 1 \
  -v /home/${USER}:/home/${USER} \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/model_train \
  --dataset_config_pbtxt="${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt" \
  --train_dir="${TRAINING_DIR}" \
  --model_name="inception_v3" \
  --number_of_steps=50000 \
  --save_interval_secs=300 \
  --batch_size=32 \
  --learning_rate=0.0005 \
  --start_from_checkpoint="${GCS_PRETRAINED_WGS_MODEL}" \
) > "${LOG_DIR}/train.log" 2>&1 &
```

At the same time, we start `model_eval` on the same machine. Given we only have
1 GPU in this example and is being used in `model_train`, we run `model_eval`
on CPUs instead (without `--gpus 1`).

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
batch_size 32). Note that `model_eval` will not stop on its own, so I had to
kill the process after training is no longer producing more checkpoints.

### (Optional) Use TensorBoard to visualize progress

We can start a TensorBoard to visualize the progress of training better. This
step is optional.

You'll want to let model_train and model_eval run for a while before you start a
TensorBoard. (You can start a TensorBoard immediately, but you just won't see
the metrics summary until later.)

We did this through a Google Cloud Shell from https://console.cloud.google.com,
on the top right:

![Shell](images/ActivateShell.png?raw=true "Activate Google Cloud Shell")

This opens up a terminal at the bottom of the browser page, then run:

```
# Change to your OUTPUT_BUCKET from earlier.
OUTPUT_BUCKET="${OUTPUT_GCS_BUCKET}/customized_training"
TRAINING_DIR="${OUTPUT_BUCKET}/training_dir"
tensorboard --logdir ${TRAINING_DIR} --port=8080
```

After it started, I clicked on the “Web Preview” on the top right of the mini
terminal:

![WebPreview](images/WebPreview.png?raw=true "Web Preview")

And clicked on "Preview on port 8080":

![PreviewOnPort](images/PreviewOnPort.png?raw=true "Preview on Port 8080")

Once it starts, you can see many metrics, including accuracy, speed, etc. You
will need to wait for both `model_train` and `model_eval` to run for a while
before the plots will make more sense.

### Pick a model

You can directly look up the best checkpoint by running:

```
gsutil cat "${TRAINING_DIR}"/best_checkpoint.txt
```

In my run, this showed that the model checkpoint that performs the best on the
validation set was `${TRAINING_DIR}/model.ckpt-33739`.

It's possible that training more steps can result in better accuracy. For now
let's use this model to do the final evaluation on the test set and see how we
do. We can use the one-step command to call:

```
sudo docker run --gpus 1 \
  -v /home/${USER}:/home/${USER} \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --customized_model "${TRAINING_DIR}/model.ckpt-33739" \
  --ref "${REF}" \
  --reads "${BAM_CHR20}" \
  --regions "chr20" \
  --output_vcf "${OUTPUT_DIR}/test_set.vcf.gz" \
  --num_shards=${N_SHARDS}
```

In v1.4.0, by using `--model_type WGS`, `run_deepvariant` will
automatically add `insert_size` as an extra channel in the `make_examples` step.
So we don't need to add it in `--make_examples_extra_args`.

Once this is done, we have the final callset in VCF
format here: `${OUTPUT_DIR}/test_set.vcf.gz`. Next step is to run `hap.py` to
complete the evaluation on chromosome 20:

```
sudo docker pull jmcdani20/hap.py:v0.3.12

time sudo docker run -it \
-v "${DATA_DIR}:${DATA_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${OUTPUT_DIR}/test_set.vcf.gz" \
  -f "${TRUTH_BED}" \
  -r "${REF}" \
  -o "${OUTPUT_DIR}/chr20-calling.happy.output" \
  -l chr20 \
  --engine=vcfeval \
  --pass-only
```

The output of `hap.py` is here:

```
[I] Total VCF records:         3775119
[I] Non-reference VCF records: 3775119
[W] overlapping records at chr20:35754687 for sample 0
[W] Variants that overlap on the reference allele: 3
[I] Total VCF records:         132914
[I] Non-reference VCF records: 96580
2022-06-02 00:36:08,582 WARNING  Creating template for vcfeval. You can speed this up by supplying a SDF template that corresponds to /home/pichuan_google_com/training-case-study/input/data/ucsc_hg19.fa
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10023      9806       217        19266       163       8898    107     33       0.978350          0.984279        0.461850         0.981305                     NaN                     NaN                   1.547658                   2.046311
INDEL   PASS        10023      9806       217        19266       163       8898    107     33       0.978350          0.984279        0.461850         0.981305                     NaN                     NaN                   1.547658                   2.046311
  SNP    ALL        66237     66160        77        78315        54      12065     15      4       0.998838          0.999185        0.154057         0.999011                2.284397                2.200204                   1.700387                   1.798656
  SNP   PASS        66237     66160        77        78315        54      12065     15      4       0.998838          0.999185        0.154057         0.999011                2.284397                2.200204                   1.700387                   1.798656
```

To summarize, the accuracy is:

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 217  | 163  | 0.978350 | 0.984279  | 0.981305
SNP   | 77   | 54   | 0.998838 | 0.999185  | 0.999011

The baseline we're comparing to is to directly use the WGS model to make the
calls, using this command:

```bash
sudo docker run --gpus 1 \
  -v /home/${USER}:/home/${USER} \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --ref "${REF}" \
  --reads "${BAM_CHR20}" \
  --regions "chr20" \
  --output_vcf "${OUTPUT_DIR}/baseline.vcf.gz" \
  --num_shards=${N_SHARDS}
```

Baseline:

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 457  | 912  | 0.954405 | 0.916192  | 0.934908
SNP   |  69  | 85   | 0.998958 | 0.998718  | 0.998838

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

[GPU machine]: deepvariant-details.md#command-for-a-gpu-machine-on-google-cloud-platform
