# Training DeepVariant models

It is possible to train your own DeepVariant models with the released
DeepVariant code. By training your own models you can potentially teach
DeepVariant to call variants better on your *specific* sequencing data than the
stock released model. It is often possible to greatly increase the accuracy of
DeepVariant on new data types, particularly if these look significantly
different from the deep whole-genome [Illumina] sequencing data used to train
the stock model, such as exomes or targeted capture or other sequencing
technologies.

Unfortunately, the Genomics team does not have the capacity to support setting
up, running, or evaluating training new models for DeepVariant. That said, the
following high-level overview should help bootstrap training your own custom
models or fine-tune the released model on new datasets:

*   Run `make_examples` in training mode, providing the `--confident_regions`
    and `--truth_variants` arguments so that `make_examples` produces the
    labeled TensorFlow examples uses during training. Confident regions and
    truth variants should be a BED file of regions where variants (or not) have
    been confidently determined and a VCF file containing all of the variants
    within those regions. [Genome in a
    Bottle](https://github.com/genome-in-a-bottle) provides examples of such
    files.
*   Prepare a config text protobuf pointing to the output of `make_examples`,
    which looks something like:

```
name: "my-training-dataset"
tfrecord_path: "/path/to/my/training_examples.tfrecord"
num_examples: 42
```

*   Setup a machine or machines to use with TensorFlow for training. We
    recommend either a distributed configuration or a single machine with a
    large training capacity (multiple GPUs or TPUs) as many millions of training
    steps are required to produce high-quality models. Less are required for
    fine-tuning an existing model to a new dataset, though.
*   Training the model is accomplished with the `model_train.py` script; check
    the source code of this file for required arguments and options. A key
    argument is `--dataset_config_pbtxt`, which should point to the dataset
    config proto create above, and is used by `model_train.py` to feed labeled
    data into TensorFlow for model training. Another key argument is
    `--start_from_checkpoint`, which tells the training system to start not from
    a random parameter set but to initialize the model with the parameters from
    the provided model.

    We often use two pre-trained models to start DeepVariant training. (a) The
    imagenet `inception_v3` checkpoint from
    [Slim](https://github.com/tensorflow/models/tree/master/research/slim#Data),
    which we use for training models "from scratch". (b) The released
    DeepVariant model itself, in order to "fine-tune" the DeepVariant model for
    your own dataset.

It should be possible, to verify everything is setup correctly, to (1) make some
labeled examples using the [WGS case study] data on chr20 (2) start up
`model_train.py`, it's companion `model_eval.py`, and [TensorBoard] on a single
machine with a GPU and (3) interact with [TensorBoard] to get a sense of how
training DeepVariant models works.

Once this is working, it's a matter of scaling up the training to realistic
datasets and steps. The Genomics team trains DeepVariant using a distributed
TensorFlow training cluster to roughly ~10 million steps (with a batch size of
64, meaning we train on 640M pileup tensors in a standard training run). We then
select a few high-quality models within those steps using [TensorBoard] and QC
them for their actual calling performance.

We recommend you read through the general TensorFlow documentation on training
models to better understand how this all works and fits together:

*   [TensorFlow Mechanics
    101](https://www.tensorflow.org/get_started/mnist/mechanics)
*   [Training inception](https://www.tensorflow.org/tutorials/image_retraining)
*   [Distributed TensorFlow](https://www.tensorflow.org/deploy/distributed)

[WGS case study]: deepvariant-case-study.md
[Illumina]: http://www.illumina.com/
[TensorBoard]: https://www.tensorflow.org/get_started/summaries_and_tensorboard
