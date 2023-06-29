# Copyright 2017 Google LLC.
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


def get_config() -> ml_collections.ConfigDict:
  """Training parameters."""
  config = ml_collections.ConfigDict()
  config.early_stopping = ml_collections.ConfigDict()
  # Reduce steps per epoch so that we log more frequently, to workaround
  # the infrequent logging due to the TensorBoard update_freq arg being
  # deprecated (internal): https://datascience.stackexchange.com/a/47416
  config.num_mini_epochs_per_epoch = 25
  config.batch_size = 512
  config.num_epochs = 10
  config.num_validation_examples = 1500000
  config.train_dataset_pbtxt = '/placer/prod/home/brain-genomics/pichuan/deepvariant_wgs/cl511097225-element-vg/cl511097225-element-vg_train.dataset_config.pbtxt'
  config.tune_dataset_pbtxt = '/placer/prod/home/brain-genomics/pichuan/deepvariant_wgs/cl511097225-element-vg/cl511097225-element-vg_tune.dataset_config.pbtxt'
  config.learning_rate = 0.256
  config.learning_rate_num_epochs_per_decay = 2.0
  config.learning_rate_decay_rate = 0.947
  config.average_decay = 0.999
  config.label_smoothing = 1e-6
  config.rho = 0.9
  config.momentum = 0.9
  config.epsilon = 1.0
  config.warmup_steps = 10_000
  config.num_mini_epochs_per_epoch = 25
  config.init_weights_file = None
  config.init_backbone_with_imagenet = True
  # Currently the code in train_inceptionv3 assumes the best_metrics uses 'max'.
  # If you use a best_metrics that uses 'min', you need to change
  # train_inceptionv3.py.
  config.best_metrics = 'val_categorical_accuracy'
  # Number of epochs for EarlyStopping. In our code, this will be multiplied
  # num_mini_epochs_per_epoch to convert into number of mini epochs.
  config.early_stopping.patience = 5

  config.log_every_steps = 100
  config.eval_every_steps = 100_000

  config.prefetch_buffer_bytes = 16 * 1000 * 1000
  config.shuffle_buffer_elements = 1024
  config.input_read_threads = 32

  return config
