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
  config.train_dataset_pbtxt = '/placer/prod/home/brain-genomics/pichuan/deepvariant_exome/dv-wes-keras-cl522895908/dv-wes-keras-cl522895908_train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/placer/prod/home/brain-genomics/pichuan/deepvariant_exome/dv-wes-keras-cl522895908/dv-wes-keras-cl522895908_tune.dataset_config.pbtxt'
  config.num_validation_examples = None

  config.num_epochs = 19
  config.learning_rate = 0.006711579692248813
  config.learning_rate_num_epochs_per_decay = 2.0
  config.learning_rate_decay_rate = 0.9913758960978765
  config.rho = 0.9763046740422171
  config.momentum = 0.9848544529312561
  config.epsilon = 0.8696723762650027
  config.warmup_steps = 718
  config.weight_decay = 0.1
  config.backbone_dropout_rate = 0.22517227651098964

  config.init_weights_file = '/readahead/512M/cns/oz-d/home/brain-genomics/kishwar/deepvariant_1.6_release_models/wgs_rc0/wgs.rc0.ckpt'

  return config


def get_config(config_name: str) -> ml_collections.ConfigDict:
  """Training parameters."""
  config = ml_collections.ConfigDict()
  config.early_stopping = ml_collections.ConfigDict()

  config.model_type = 'inception_v3'
  config.trial = 0  # Used to allow for replicates during training.

  # Default Dataset
  config.train_dataset_pbtxt = '/placer/prod/home/brain-genomics/pichuan/deepvariant_wgs/cl511097225-element-vg/cl511097225-element-vg_train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/placer/prod/home/brain-genomics/pichuan/deepvariant_wgs/cl511097225-element-vg/cl511097225-element-vg_tune.dataset_config.pbtxt'

  config.best_checkpoint_metric = 'tune/f1_weighted'
  config.batch_size = 16384
  config.num_epochs = 10
  config.num_validation_examples = 1500000

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
  config.init_weights_file = None
  config.init_backbone_with_imagenet = True
  config.best_metrics = 'val_categorical_accuracy'
  config.weight_decay = 0.00004
  config.backbone_dropout_rate = 0.2

  # TensorBoard Options
  config.log_every_steps = 100
  config.tune_every_steps = 100_000

  # Data Pipeline Options
  config.prefetch_buffer_bytes = 16 * 1000 * 1000
  config.shuffle_buffer_elements = 1024
  config.input_read_threads = 32

  if config_name == 'exome':
    config = get_exome_config(config)
  elif config_name == 'base':
    raise ValueError('No model specified for config type.')
  else:
    raise ValueError(f'Unknown config_name: {config_name}')

  return config
