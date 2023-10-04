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


def get_exome_config(
    config: ml_collections.ConfigDict,
) -> ml_collections.ConfigDict:
  """Config parameters for exome training."""

  # Exome Dataset
  config.train_dataset_pbtxt = '/path/to/your/train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/path/to/your/tune.dataset_config.pbtxt'
  # If set to 0, use full validation dataset.
  config.num_validation_examples = 0

  config.num_epochs = 80
  config.learning_rate = 0.01
  config.learning_rate_num_epochs_per_decay = 2.0
  config.learning_rate_decay_rate = 0.9999
  config.rho = 0.9763046740422171
  config.momentum = 0.9848544529312561
  config.epsilon = 0.8696723762650027
  config.warmup_steps = 718
  config.weight_decay = 0.1
  config.backbone_dropout_rate = 0.22517227651098964

  config.init_checkpoint = ''

  return config


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
  config.early_stopping_patience = 10

  # TensorBoard Options
  config.log_every_steps = 100
  # Tuning happens at every epoch. The frequency can be increased here.
  config.tune_every_steps = 100_000

  # Data Pipeline Options
  config.prefetch_buffer_bytes = 16 * 1000 * 1000
  config.shuffle_buffer_elements = 10_000
  config.input_read_threads = 32

  # Placeholder value for limiting training examples. 0=No limit.
  config.limit = 0

  if config_name == 'exome':
    config = get_exome_config(config)
  elif config_name == 'base':
    # Use the base config.
    pass
  else:
    raise ValueError(f'Unknown config_name: {config_name}')

  return config
