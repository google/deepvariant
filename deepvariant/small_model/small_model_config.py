# Copyright 2024 Google LLC.
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
"""Model configs for small model training."""

import enum
import ml_collections


class PresetConfig(enum.Enum):
  """Enum of preset configs."""

  WGS = "wgs"
  PACBIO = "pacbio"
  ONT = "ont"


def set_wgs_config(config: ml_collections.ConfigDict) -> None:
  config.train_tfrecord_directory = "/placer/prod/home/brain-genomics/lucasbrambrink/deepvariant_wgs/b393196659_wgs/b393196659_wgs_train_small_model.examples*.tfrecord.gz"
  config.tune_tfrecord_directory = "/placer/prod/home/brain-genomics/lucasbrambrink/deepvariant_wgs/b393196659_wgs/b393196659_wgs_tune_small_model.examples*.tfrecord.gz"
  config.num_train_samples = 941_401_617
  config.num_tune_samples = 27_693_525
  config.model_params.expand_by_haplotype = False


def set_pacbio_config(config: ml_collections.ConfigDict) -> None:
  config.train_tfrecord_directory = "/placer/prod/home/brain-genomics/lucasbrambrink/deepvariant_pacbio/b393196659_pacbio/b393196659_pacbio_train_small_model.examples*.tfrecord.gz"
  config.tune_tfrecord_directory = "/placer/prod/home/brain-genomics/lucasbrambrink/deepvariant_pacbio/b393196659_pacbio/b393196659_pacbio_tune_small_model.examples*.tfrecord.gz"
  config.num_train_samples = 592_876_080
  config.num_tune_samples = 18_221_432
  config.model_params.expand_by_haplotype = True


def set_ont_config(config: ml_collections.ConfigDict) -> None:
  config.train_tfrecord_directory = "/placer/prod/home/brain-genomics/lucasbrambrink/deepvariant_ont/b393196659_ont/b393196659_ont_train_small_model.examples*.tfrecord.gz"
  config.tune_tfrecord_directory = "/placer/prod/home/brain-genomics/lucasbrambrink/deepvariant_ont/b393196659_ont/b393196659_ont_tune_small_model.examples*.tfrecord.gz"
  config.num_train_samples = 1_667_298_556
  config.num_tune_samples = 47_022_920
  config.model_params.expand_by_haplotype = True


def set_preset_config(
    preset_config: PresetConfig, config: ml_collections.ConfigDict
) -> None:
  """Sets the config for a given preset."""
  if preset_config == PresetConfig.WGS:
    set_wgs_config(config)
  elif preset_config == PresetConfig.PACBIO:
    set_pacbio_config(config)
  elif preset_config == PresetConfig.ONT:
    set_ont_config(config)


def get_config(config_name: str) -> ml_collections.ConfigDict:
  """Returns the default configuration as instance of ConfigDict."""
  # Model hyperparameters
  model_params = ml_collections.ConfigDict()
  model_params.activation = "relu"
  model_params.hidden_layer_sizes = (750, 750)
  model_params.optimizer = "adam"
  model_params.learning_rate = 1e-04
  model_params.learning_rate_num_epochs_per_decay = 1
  model_params.learning_rate_decay_rate = 0.99
  model_params.weight_decay = 0.0000001
  model_params.steps_per_execution = 128
  model_params.features = ()
  model_params.vaf_context_window_size = 51
  model_params.expand_by_haplotype = False

  # Training parameters
  config = ml_collections.ConfigDict()
  config.config_name = config_name
  config.epochs = 50
  config.batch_size = 512
  config.logging_frequency = 16384
  config.interleave_cycle_length = 4
  config.interleave_parallel_calls = 64
  config.shuffle_buffer_elements = 100_000
  config.tfrecord_buffer_size = 64 * 1000 * 1000
  config.tfrecord_num_parallel_calls = 12
  config.prefetch_buffer_size = 12
  config.read_ahead_buffer = "128M"
  config.map_parallel_calls = 100

  # Training data
  config.train_tfrecord_directory = ""
  config.num_train_samples = 0
  config.tune_tfrecord_directory = ""
  config.num_tune_samples = 0

  # Training environment
  config.trial_id = 1
  config.tensorboard_directory = None
  config.is_xmanager_run = False

  config.model_params = model_params
  try:
    preset_config = PresetConfig(config_name)
    set_preset_config(preset_config, config)
  except ValueError as e:
    raise ValueError(f"Unknown preset config name: {config_name}") from e

  return config
