# Copyright 2023 Google LLC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""Config for use with the custom training DeepVariant Loop."""

import ml_collections

# =============#
# Test Config #
# =============#


def get_test_config(config: ml_collections.ConfigDict):
  """Config parameters for test training."""

  # Config values to reduce the size and time for training
  config.batch_size = 4
  config.num_epochs = 2
  config.num_validation_examples = 1
  config.warmup_steps = 0
  config.limit = 50
  config.steps_per_iter = 4
  config.shuffle_buffer_elements = 50
  config.init_checkpoint = ''


# =====================#
# DeepVariant Configs #
# =====================#


def get_wgs_config(config: ml_collections.ConfigDict):
  """Config parameters for wgs training."""

  # Exome Dataset
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = ''
  # If set to 0, use full validation dataset.
  config.num_validation_examples = 150_000

  config.best_checkpoint_metric = 'tune/f1_weighted'
  config.batch_size = 16384
  config.num_epochs = 10
  config.optimizer = 'adam'
  config.beta_1 = 0.9651804083266324
  config.beta_2 = 0.9665259112630292
  config.weight_decay = 0.00004
  config.adaptive_epsilon = True
  config.optimizer_weight_decay = 0.0

  config.early_stopping_patience = 100
  config.learning_rate = 0.0000796142074327502
  config.learning_rate_num_epochs_per_decay = 2.25
  config.learning_rate_decay_rate = 0.9999
  config.warmup_steps = 0

  config.backbone_dropout_rate = 0.2

  # Exponential Moving Average
  config.use_ema = True
  config.ema_momentum = 0.991463134331829


def get_exome_config(config: ml_collections.ConfigDict):
  """Config parameters for exome training."""

  # Exome Dataset
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = '/path/to/warmstart/checkpoint'
  # If set to 0, use full validation dataset.
  config.num_validation_examples = 0

  config.best_checkpoint_metric = 'tune/f1_weighted'
  config.batch_size = 16384
  config.num_epochs = 4
  config.optimizer = 'adam'
  config.beta_1 = 0.9575875572181167
  config.beta_2 = 0.9475272158875401
  config.adaptive_epsilon = True
  config.weight_decay = 0.00004
  config.optimizer_weight_decay = 0.0

  config.early_stopping_patience = 250
  config.learning_rate = 0.00008663001151624387
  config.learning_rate_num_epochs_per_decay = 2.25
  config.learning_rate_decay_rate = 0.5
  config.warmup_steps = 27249

  config.backbone_dropout_rate = 0.0

  # Exponential Moving Average
  config.use_ema = True
  config.ema_momentum = 0.9725254942104883


def get_pacbio_config(config: ml_collections.ConfigDict):
  """Training parameters."""
  config.num_epochs = 8
  config.num_validation_examples = 150_000
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = '/path/to/warmstart/checkpoint'

  config.best_checkpoint_metric = 'tune/categorical_accuracy'
  config.batch_size = 16384
  config.num_epochs = 5
  config.optimizer = 'adam'
  config.beta_1 = 0.9651804083266324
  config.beta_2 = 0.9665259112630292
  config.adaptive_epsilon = True
  config.weight_decay = 0.00004
  config.optimizer_weight_decay = 0.0

  config.early_stopping_patience = 100
  config.learning_rate = 0.00008663001151624387
  config.learning_rate_num_epochs_per_decay = 2.66
  config.learning_rate_decay_rate = 0.8514735277962562
  config.warmup_steps = 0

  config.backbone_dropout_rate = 0.0

  # Exponential Moving Average
  config.use_ema = True
  config.ema_momentum = 0.991463134331829


# =====================#
# DeepSomatic Configs #
# =====================#


def get_deepsomatic_wgs_config(config: ml_collections.ConfigDict):
  get_wgs_config(config)
  config.learning_rate = 0.00009483389877395854
  config.learning_rate_decay_rate = 0.5
  config.warmup_steps = 10000
  config.num_epochs = 10
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = '/path/to/warmstart/checkpoint'


def get_deepsomatic_wes_config(
    config: ml_collections.ConfigDict,
):
  """Config parameters for wgs training."""
  get_wgs_config(config)
  # Exome Dataset
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = ''
  config.num_epochs = 80
  config.early_stopping_patience = 10


def get_deepsomatic_pacbio_tumor_normal_config(
    config: ml_collections.ConfigDict,
):
  """Training parameters."""
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = ''

  config.batch_size = 8192
  config.best_checkpoint_metric = 'tune/f1_homalt'
  config.class_weights = '1,1,1'
  config.num_epochs = 4
  config.num_validation_examples = 15_000_000
  config.optimizer = 'adam'
  config.beta_1 = 0.9651804083266324
  config.beta_2 = 0.9665259112630292
  config.adaptive_epsilon = True
  config.weight_decay = 0.00004
  config.optimizer_weight_decay = 0.0

  config.early_stopping_patience = 100
  config.learning_rate = 0.00004850252247279374
  config.learning_rate_num_epochs_per_decay = 2.66
  config.learning_rate_decay_rate = 0.8011670994484517
  config.warmup_steps = 8607

  config.backbone_dropout_rate = 0.0

  # Exponential Moving Average
  config.use_ema = True
  config.ema_momentum = 0.991463134331829


def get_deepsomatic_ont_tumor_normal_config(config: ml_collections.ConfigDict):
  get_deepsomatic_pacbio_tumor_normal_config(config)
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = ''


def get_deepsomatic_pacbio_tumor_only_config(config: ml_collections.ConfigDict):
  """Training parameters."""
  get_pacbio_config(config)
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = ''


def get_deepsomatic_wgs_tumor_only_config(config: ml_collections.ConfigDict):
  get_wgs_config(config)
  config.class_weights = '1,1,10'
  config.train_dataset_pbtxt = '/placer/prod/home/brain-genomics/danielecook/deepsomatic/ds_wgs_tumor_only/ds_wgs_tumor_only_train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/placer/prod/home/brain-genomics/danielecook/deepsomatic/ds_wgs_tumor_only/ds_wgs_tumor_only_tune.dataset_config.pbtxt'


# FFPE Configs
def get_deepsomatic_wgs_ffpe_config(
    config: ml_collections.ConfigDict,
):
  get_wgs_config(config)
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = '/path/to/warmstart/checkpoint'


def get_deepsomatic_wes_ffpe_config(
    config: ml_collections.ConfigDict,
):
  """Config parameters for wgs training."""
  get_wgs_config(config)
  # Exome Dataset
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  config.init_checkpoint = ''
  config.num_epochs = 80
  config.early_stopping_patience = 10


# =============#
# Base Config #
# =============#


def get_config(config_name: str) -> ml_collections.ConfigDict:
  """Training parameters."""
  config = ml_collections.ConfigDict()

  config.model_type = 'inception_v3'
  config.trial = 0  # Used to allow for replicates during training.

  # Default Dataset
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'

  config.best_checkpoint_metric = 'tune/f1_weighted'
  config.batch_size = 16384
  config.num_epochs = 10
  config.num_validation_examples = 1500000
  config.optimizer = 'rmsprop'

  config.use_ema = True

  config.denovo_enabled = False
  # class_weights can be specified as a comma-delimited string of weights
  # for each class: e.g. 1,1,10 or 1,1,1
  config.class_weights = ''

  # Training hyperparameters
  config.learning_rate = 0.001
  config.learning_rate_num_epochs_per_decay = 2.0
  config.learning_rate_decay_rate = 0.947
  config.average_decay = 0.999
  config.label_smoothing = 1e-6
  config.rho = 0.9
  config.momentum = 0.9
  config.epsilon = 1.0
  config.warmup_steps = 10_000
  config.init_checkpoint = ''
  config.init_backbone_with_imagenet = False
  config.best_metrics = 'tune/f1_weighted'
  config.weight_decay = 0.00004
  config.backbone_dropout_rate = 0.2
  # Stop training when this many consecutive evaluations yield no improvement.
  config.early_stopping_patience = 250

  # Weight decay of optimizer
  config.optimizer_weight_decay = 0.0

  # An 'iter' refers to a group of train/tune steps run in succession.
  config.steps_per_iter = 128

  # TensorBoard Options
  config.log_every_steps = 1_280
  # Tuning happens at every epoch. The frequency can be increased here.
  config.tune_every_steps = 12_800

  # Data Pipeline Options
  config.prefetch_buffer_bytes = 16 * 1000 * 1000
  config.shuffle_buffer_elements = 50_000
  config.input_read_threads = 8

  # Placeholder value for limiting training examples. 0=No limit.
  config.limit = 0

  if config_name and '+' in config_name:
    config_name, alt_mode = config_name.split('+')
  else:
    alt_mode = None

  config.alt_mode = alt_mode

  if config_name == 'wgs':
    get_wgs_config(config)
  elif config_name == 'exome':
    get_exome_config(config)
  elif config_name == 'pacbio':
    get_pacbio_config(config)
  elif config_name == 'deepsomatic_wgs':
    get_deepsomatic_wgs_config(config)
  elif config_name == 'deepsomatic_wgs_tumor_only':
    get_deepsomatic_wgs_tumor_only_config(config)
  elif config_name == 'deepsomatic_wes':
    get_deepsomatic_wes_config(config)
  elif config_name == 'deepsomatic_wgs_ffpe':
    get_deepsomatic_wgs_ffpe_config(config)
  elif config_name == 'deepsomatic_wes_ffpe':
    get_deepsomatic_wes_ffpe_config(config)
  elif config_name == 'deepsomatic_pacbio_tumor_normal':
    get_deepsomatic_pacbio_tumor_normal_config(config)
  elif config_name == 'deepsomatic_ont_tumor_normal':
    get_deepsomatic_ont_tumor_normal_config(config)
  elif config_name == 'deepsomatic_pacbio_tumor_only':
    get_deepsomatic_pacbio_tumor_only_config(config)
  elif config_name == 'base':
    # Use the base config.
    pass
  else:
    raise ValueError(f'Unknown config_name: {config_name}')

  if config.alt_mode == 'test':
    config.alt_mode = 'test'
    get_test_config(config)

  return config
