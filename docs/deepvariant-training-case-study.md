# Advanced Case Study: Train a customized SNP and small indel variant caller for BGISEQ-500 data.

DeepVariant is an analysis pipeline that uses a deep neural network to call
genetic variants from next-generation DNA sequencing (NGS) data. While
DeepVariant is highly accurate for
[many types of NGS data](https://rdcu.be/7Dhl), some users may be interested in
training custom deep learning models that have been optimized for very specific
data.

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

*   Indel F1 `94.1957%` --> `98.1289%`
*   SNP F1: `99.8725%` --> `99.9162%`

This tutorial is meant as an example for training; all the other processing in
this tutorial were done serially with no pipeline optimization.

## Request a machine

For this case study, we use a [GPU machine] with 16 vCPUs. You can request this
machine on Google Cloud using the following command:

```bash
host="${USER}-deepvariant-vm"
zone="us-west1-b"

gcloud compute instances create ${host} \
    --scopes "compute-rw,storage-full,cloud-platform" \
    --maintenance-policy "TERMINATE" \
    --accelerator=type=nvidia-tesla-p100,count=1 \
    --image-family "ubuntu-2204-lts" \
    --image-project "ubuntu-os-cloud" \
    --machine-type "n1-standard-16" \
    --boot-disk-size "300" \
    --zone "${zone}" \
    --min-cpu-platform "Intel Skylake"
```

After a minute or two, your VM should be ready and you can ssh into it using the
following command:

```bash
gcloud compute ssh ${host} --zone ${zone}
```

Once you have logged in, set the variables:

```bash
YOUR_PROJECT=REPLACE_WITH_YOUR_PROJECT
OUTPUT_GCS_BUCKET=REPLACE_WITH_YOUR_GCS_BUCKET

BUCKET="gs://deepvariant"
VERSION="1.9.0"
DOCKER_IMAGE="google/deepvariant:${VERSION}"

GCS_PRETRAINED_WGS_MODEL="${BUCKET}/models/DeepVariant/1.9.0/checkpoints/wgs/deepvariant.wgs.ckpt"

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

N_SHARDS=$(nproc)
```

## Download binaries and data

### Create directories:

```bash
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BIN_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${LOG_DIR}"
```

### Copy data

```bash
gsutil -m cp ${DATA_BUCKET}/BGISEQ_PE100_NA12878.sorted.chr*.bam* "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/ucsc_hg19.fa*" "${DATA_DIR}"
gsutil -m cp -r "${DATA_BUCKET}/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_*" "${DATA_DIR}"
```

### Download extra packages

```bash
sudo apt -y update
sudo apt -y install parallel
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.9/scripts/install_nvidia_docker.sh
bash -x install_nvidia_docker.sh
```

## Run make_examples in “training” mode for training and validation sets.

Create examples in "training" mode (which means these `tensorflow.Example`s will
contain a `label` field).

In this tutorial, we create examples on one replicate of HG001 sequenced by
BGISEQ-500 provided on the
[Genome In a Bottle FTP site](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/BGISEQ500/standard_library/readme.txt).

In this tutorial, we will split the genome up into the following datasets:

| chrom | Name                  | Description                                  |
| ----- | --------------------- | -------------------------------------------- |
| chr1  | Training Set          | Examples used to train our model.            |
| chr21 | Validation / Tune Set | Examples used to evaluate the performance of our model during training |
| chr20 | Test Set              | Examples reserved for testing performance of our trained model |

Note that normally, the training dataset will be much larger (e.g. chr1-19),
rather than just a single chromosome. We use just chr1 here to demonstrate how
customized training works.

For the definition of these 3 sets in commonly used machine learning
terminology, please refer to
[Machine Learning Glossary](https://developers.google.com/machine-learning/glossary/).

### Training set

First, to set up, lets pull the docker images.

```bash
sudo docker pull ${DOCKER_IMAGE}     # Standard CPU Docker Image.
sudo docker pull ${DOCKER_IMAGE}-gpu # GPU-enabled Docker image.
```

The `make_examples` step doesn't use GPU, so we will not require the GPU-enabled
image.

```bash
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --line-buffer \
    sudo docker run \
      -v ${HOME}:${HOME} \
      ${DOCKER_IMAGE} \
      make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM_CHR1}" \
      --examples "${OUTPUT_DIR}/training_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "'chr1'" \
      --channel_list "BASE_CHANNELS,insert_size" \
) 2>&1 | tee "${LOG_DIR}/training_set.with_label.make_examples.log"
```

This took 18m12.709s.

Starting in v1.4.0, we added an extra channel in our WGS setting. The
`--channel_list` sets the channels that our model will use. You can use
`BASE_CHANNELS` as a shorthand for specifying six commonly used channels, and
append additional channels to this list (e.g. `insert_size`). The make_examples
step creates `*.example_info.json` files. For example, you can see it here:

```
cat "${OUTPUT_DIR}/training_set.with_label.tfrecord-00000-of-000${N_SHARDS}.gz.example_info.json"
```

```json
{"version": "1.9.0", "shape": [100, 221, 7], "channels": [1, 2, 3, 4, 5, 6, 19]}
```

Depending on your data type, you might want to tweak the flags for the
`make_examples` step, which can result in different shape of the output
examples.

We will want to shuffle this on Dataflow later, so we copy the data to GCS
bucket first:

```
gsutil -m cp ${OUTPUT_DIR}/training_set.with_label.tfrecord-?????-of-000${N_SHARDS}.gz* \
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
      ${DOCKER_IMAGE} \
      make_examples \
      --mode training \
      --ref "${REF}" \
      --reads "${BAM_CHR21}" \
      --examples "${OUTPUT_DIR}/validation_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions "'chr21'" \
      --channel_list "BASE_CHANNELS,insert_size" \
) 2>&1 | tee "${LOG_DIR}/validation_set.with_label.make_examples.log"
```

This took: 4m33.081s.

Copy to GCS bucket:

```bash
gsutil -m cp ${OUTPUT_DIR}/validation_set.with_label.tfrecord-?????-of-000${N_SHARDS}.gz* \
  ${OUTPUT_BUCKET}
```

## Shuffle each set of examples and generate a data configuration file for each.

Shuffling the `tensorflow.Example`s is an important step for training a model.
In our training logic, we shuffle examples globally using a preprocessing step.

First, if you have run this step before, and want to rerun it, you might want to
consider cleaning up previous data first to avoid confusion:

```bash
# (Optional) Clean up existing files.
gsutil -m rm -f "${OUTPUT_BUCKET}/training_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
gsutil rm -f "${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt"
gsutil -m rm -f "${OUTPUT_BUCKET}/validation_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
gsutil rm -f "${OUTPUT_BUCKET}/validation_set.dataset_config.pbtxt"
gsutil rm -f "${OUTPUT_BUCKET}/example_info.json"
```

Here we provide examples for running on
[Cloud Dataflow Runner](https://beam.apache.org/documentation/runners/dataflow/)
and also [DirectRunner](https://beam.apache.org/documentation/runners/direct/).
Beam can also use other runners, such as
[Spark Runner](https://beam.apache.org/documentation/runners/spark/).

First, create a virtual environment to install beam on your machine.

```bash
sudo apt install -y python3.10-venv
# Create a virtualenv
python3 -m venv beam

# Activate the virtualenv
. beam/bin/activate
```

Consult the instructions at https://beam.apache.org/get-started/quickstart-py/
if you run into any issues.

Then, get the script that performs shuffling:

```bash
mkdir -p ${SHUFFLE_SCRIPT_DIR}
wget https://raw.githubusercontent.com/google/deepvariant/r1.9/tools/shuffle_tfrecords_beam.py -O ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py
```

Next, we shuffle the data using DataflowRunner. Before that, please make sure
you enable Dataflow API for your project:
http://console.cloud.google.com/flows/enableapi?apiid=dataflow.

To access `gs://` path, make sure you run this in your virtual environment:

```bash
sudo apt -y update && sudo apt -y install python3-pip
pip3 install --upgrade pip
pip3 install setuptools --upgrade
pip3 install tensorflow  # For parsing tf.Example in shuffle_tfrecords_beam.py.
pip3 install apache_beam[gcp]==2.50.0  # 2.51.0 didn't work in my run.
```

Shuffle using Dataflow.

```bash
time python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="${YOUR_PROJECT}" \
  --input_pattern_list="${OUTPUT_BUCKET}"/training_set.with_label.tfrecord-?????-of-000${N_SHARDS}.gz \
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

```bash
gsutil cat "${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt"
```

In the output, the `tfrecord_path` should be valid paths in gs://.

```
# Generated by shuffle_tfrecords_beam.py
# class2: 124571
# class1: 173664
# class0: 43803

name: "HG001"
tfrecord_path: "OUTPUT_BUCKET/training_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
num_examples: 342038
#
# --input_pattern_list=OUTPUT_BUCKET/training_set.with_label.tfrecord-?????-of-00016.gz
# --output_pattern_prefix=OUTPUT_BUCKET/training_set.with_label.shuffled
#
```


We can shuffle the validation set locally using
[DirectRunner](https://beam.apache.org/documentation/runners/direct/). Adding
`--direct_num_workers=0` sets the number of threads/subprocess to the number of
cores of the machine where the pipeline is running.

```bash
time python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="${YOUR_PROJECT}" \
  --input_pattern_list="${OUTPUT_DIR}"/validation_set.with_label.tfrecord-?????-of-000${N_SHARDS}.gz \
  --output_pattern_prefix="${OUTPUT_DIR}/validation_set.with_label.shuffled" \
  --output_dataset_name="HG001" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/validation_set.dataset_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=0
```

Here is the validation_set:

```bash
cat "${OUTPUT_DIR}/validation_set.dataset_config.pbtxt"
```

```
# Generated by shuffle_tfrecords_beam.py
# class0: 5491
# class1: 31856
# class2: 21955

name: "HG001"
tfrecord_path: "OUTPUT_DIR/validation_set.with_label.shuffled-?????-of-?????.tfrecord.gz"
num_examples: 59302
#
# --input_pattern_list=OUTPUT_DIR/validation_set.with_label.tfrecord-?????-of-00016.gz
# --output_pattern_prefix=OUTPUT_DIR/validation_set.with_label.shuffled
#
```

### Fetch a config file

Before we can begin training, we will need a configuration file containing
training parameters. Parameters within this training file can be overridden when
we run `train` by passing `--config.<param>=<value>`.

```bash
curl https://raw.githubusercontent.com/google/deepvariant/r1.9/deepvariant/dv_config.py > dv_config.py
```

### Start `train`

NOTE: all parameters below are used as an example. They are not optimized for
this dataset, and are not recommended as the best default either.

```bash
( time sudo docker run --gpus 1 \
    -v /home/${USER}:/home/${USER} \
    -w /home/${USER} \
    ${DOCKER_IMAGE}-gpu \
    train \
    --config=dv_config.py:base \
    --config.train_dataset_pbtxt="${OUTPUT_BUCKET}/training_set.dataset_config.pbtxt" \
    --config.tune_dataset_pbtxt="${OUTPUT_DIR}/validation_set.dataset_config.pbtxt" \
    --config.init_checkpoint=${GCS_PRETRAINED_WGS_MODEL} \
    --config.num_epochs=10 \
    --config.learning_rate=0.0001 \
    --config.num_validation_examples=0 \
    --experiment_dir=${TRAINING_DIR} \
    --strategy=mirrored \
    --config.batch_size=384 \
) > "${LOG_DIR}/train.log" 2>&1 &
```

Once training starts, you should see a summary of your datasets:

```
Training Examples: 342038
Tune Examples: 59302
Batch Size: 384
Epochs: 10
Steps per epoch: 890
Steps per tune: 154
Num train steps: 8900
Steps per iter: 128
```

As training runs, the validation/tune dataset will be evaluated at the end of
each epoch, and every n training steps specified by `--config.tune_every_steps`.
You can lower `--config.tune_every_steps` to perform evaluation more frequently.

Checkpoints are stored whenever the `tune/f1_weighted` metric improves when
evaluating the tune dataset. In this way, the last checkpoint stored will always
be the best performing checkpoint. The best performing checkpoint metric can be
configured using `--config.best_checkpoint_metric`.

We tested with 1 GPU - it took 92m6.049s.

You can run with more GPUs to speed up.

Once training is complete, the following command can be used list checkpoints:

```bash
gsutil ls ${TRAINING_DIR}/checkpoints/ema/
```

The best checkpoint can be retrieved using the following command:

```bash
# Get the list of files from gsutil
files=$(gsutil ls ${TRAINING_DIR}/checkpoints/ema)

# Initialize variables to store the best filename and score
BEST_CHECKPOINT=""
best_score=0
best_checkpoint_number=0

# Iterate over the files
for file in $files; do
  # Extract the checkpoint number and score using regular expression
  if [[ $file =~ checkpoint-([0-9]+)-([0-9.]+)- ]]; then
    checkpoint=${BASH_REMATCH[1]}
    score=${BASH_REMATCH[2]}

    # Compare the score and checkpoint number with the current best
    if (( $(echo "$score > $best_score" | bc -l) || \
         ( $(echo "$score == $best_score" | bc -l) && ((checkpoint > best_checkpoint_number)) ) )); then
      # Construct the desired filename format
      BEST_CHECKPOINT="${TRAINING_DIR}/checkpoints/ema/checkpoint-${checkpoint}-${score}-1"
      best_score=$score
      best_checkpoint_number=$checkpoint
    fi
  fi
done

echo $BEST_CHECKPOINT
```

### (Optional) Use TensorBoard to visualize progress

We can start a TensorBoard to visualize the progress of training better. This
step is optional.

You'll want to let `train` run for a while before you start a TensorBoard. (You
can start a TensorBoard immediately, but you just won't see the metrics summary
until later.) We did this through a Google Cloud Shell from
https://console.cloud.google.com, on the top right:

![Shell](images/ActivateShell.png?raw=true "Activate Google Cloud Shell")

This opens up a terminal at the bottom of the browser page, then run:

```bash
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
will need to wait for `train` to run for a while before the plots will appear.

### Test the model

Now that we have performed training, we can test the performance of the new
model using our holdout dataset (chr20).

The following one-step command can be used to call DeepVariant and run our newly
trained model:

```bash
sudo docker run --gpus 1 \
  -v /home/${USER}:/home/${USER} \
  "${DOCKER_IMAGE}-gpu" \
  run_deepvariant \
  --model_type WGS \
  --customized_model "${BEST_CHECKPOINT}" \
  --ref "${REF}" \
  --reads "${BAM_CHR20}" \
  --regions "chr20" \
  --output_vcf "${OUTPUT_DIR}/test_set.vcf.gz" \
  --num_shards=${N_SHARDS} \
  --disable_small_model
```

We use a small model to classify some
candidates. In this example, we set `--disable_small_model` so
that small model is disabled. This allows us to run all examples
through the model we just trained.

Starting in v1.4.0, by using `--model_type WGS`, `run_deepvariant` will
automatically add `insert_size` as an extra channel in the `make_examples` step.
So we don't need to add it in `--make_examples_extra_args`.

When the `make_examples` step is run, you might see messages like:

```
E tensorflow/compiler/xla/stream_executor/cuda/cuda_driver.cc:268] failed call to cuInit: CUDA_ERROR_NO_DEVICE: no CUDA-capable device is detected
```

which should not be a problem.

You can use `nvidia-smi` to confirm whether the GPUs are used. If so, you can
ignore the message.

Once this is done, we have the final callset in VCF format here:
`${OUTPUT_DIR}/test_set.vcf.gz`. Next step is to run `hap.py` to complete the
evaluation on chromosome 20:

```bash
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
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10023      9801       222        19428       158       9067    107     31       0.977851          0.984751        0.466698         0.981289                     NaN                     NaN                   1.547658                   2.155273
INDEL   PASS        10023      9801       222        19428       158       9067    107     31       0.977851          0.984751        0.466698         0.981289                     NaN                     NaN                   1.547658                   2.155273
  SNP    ALL        66237     66178        59        78167        52      11899     11      1       0.999109          0.999215        0.152225         0.999162                2.284397                  2.1969                   1.700387                   1.794365
  SNP   PASS        66237     66178        59        78167        52      11899     11      1       0.999109          0.999215        0.152225         0.999162                2.284397                  2.1969                   1.700387                   1.794365
```

To summarize, the accuracy is:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 9801     | 222      | 158      | 0.977851      | 0.984751         | 0.981289        |
| SNP   | 66178    | 59       | 52       | 0.999109      | 0.999215         | 0.999162        |

The baseline we're comparing to is to directly use the WGS model to make the
calls, using this command:

```bash
sudo docker run --gpus 1 \
  -v /home/${USER}:/home/${USER} \
  ${DOCKER_IMAGE}-gpu \
  run_deepvariant \
  --model_type WGS \
  --ref "${REF}" \
  --reads "${BAM_CHR20}" \
  --regions "chr20" \
  --output_vcf "${OUTPUT_DIR}/baseline.vcf.gz" \
  --num_shards=${N_SHARDS} \
  --disable_small_model
```

Baseline:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 9648     | 375      | 685      | 0.962586      | 0.936172         | 0.949195        |
| SNP   | 66155    | 82       | 87       | 0.998762      | 0.998687         | 0.998725        |

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
