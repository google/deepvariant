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
  WES = "wes"
  PACBIO = "pacbio"
  ONT = "ont"
  HYBRID = "hybrid"


def set_wgs_config(config: ml_collections.ConfigDict) -> None:
  config.train_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_wgs/b345306347_wgs/ttl=14d/b345306347_wgs_train_ME_*.tfrecord.gz.small_model.tsv"
  config.tune_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_wgs/b345306347_wgs/ttl=14d/b345306347_wgs_tune_ME_*.tfrecord.gz.small_model.tsv"
  config.num_train_samples = 784_569_482
  config.num_tune_samples = 22_245_399


def set_wes_config(config: ml_collections.ConfigDict) -> None:
  config.train_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_exome/b345306347_exome/ttl=14d/b345306347_exome_train_ME_*.tfrecord.gz.small_model.tsv"
  config.tune_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_exome/b345306347_exome/ttl=14d/b345306347_exome_tune_ME_*.tfrecord.gz.small_model.tsv"
  config.num_train_samples = 19_217_144
  config.num_tune_samples = 726_013


def set_pacbio_config(config: ml_collections.ConfigDict) -> None:
  config.train_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_pacbio/b345306347_pacbio/ttl=14d/b345306347_pacbio_train_ME_*.tfrecord.gz.small_model.tsv"
  config.tune_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_pacbio/b345306347_pacbio/ttl=14d/b345306347_pacbio_tune_ME_*.tfrecord.gz.small_model.tsv"
  config.num_train_samples = 2_015_187_426
  config.num_tune_samples = 53_932_522


def set_ont_config(config: ml_collections.ConfigDict) -> None:
  config.train_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_ont/b352354037-base/ttl=14d/b352354037-base_train_ME_*.tfrecord.gz.small_model.tsv"
  config.tune_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_ont/b352354037-base/ttl=14d/b352354037-base_tune_ME_*.tfrecord.gz.small_model.tsv"
  config.num_train_samples = 447_500_870
  config.num_tune_samples = 12_288_308


def set_hybrid_config(config: ml_collections.ConfigDict) -> None:
  config.train_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_hybrid/b345306347_hybrid_pacbio_illumina/ttl=14d/b345306347_hybrid_pacbio_illumina_train_ME_*.tfrecord.gz.small_model.tsv"
  config.tune_tsv_directory = "/cns/oz-d/home/brain-genomics/lucasbrambrink/deepvariant_hybrid/b345306347_hybrid_pacbio_illumina/ttl=14d/b345306347_hybrid_pacbio_illumina_tune_ME_*.tfrecord.gz.small_model.tsv"
  config.num_train_samples = 142_447_287
  config.num_tune_samples = 3_916_034


def set_preset_config(
    preset_config: PresetConfig, config: ml_collections.ConfigDict
) -> None:
  """Sets the config for a given preset."""
  if preset_config == PresetConfig.WGS:
    set_wgs_config(config)
  elif preset_config == PresetConfig.WES:
    set_wes_config(config)
  elif preset_config == PresetConfig.PACBIO:
    set_pacbio_config(config)
  elif preset_config == PresetConfig.ONT:
    set_ont_config(config)
  elif preset_config == PresetConfig.HYBRID:
    set_hybrid_config(config)


def get_config(config_name: str) -> ml_collections.ConfigDict:
  """Returns the default configuration as instance of ConfigDict."""
  # Model hyperparameters
  model_params = ml_collections.ConfigDict()
  model_params.activation = "relu"
  model_params.hidden_layer_sizes = (100, 100)
  model_params.optimizer = "adam"

  # Training parameters
  config = ml_collections.ConfigDict()
  config.epochs = 20
  config.batch_size = 64
  config.logging_frequency = 8192
  config.train_tsv_directory = ""
  config.num_train_samples = 0
  config.tune_tsv_directory = ""
  config.num_tune_samples = 0

  # Local training only
  config.test_fraction = 0.2

  # Training environment
  config.trial_id = 1
  config.tensorboard_directory = None
  config.is_xmanager_run = False

  try:
    preset_config = PresetConfig(config_name)
    set_preset_config(preset_config, config)
  except ValueError as e:
    raise ValueError(f"Unknown preset config name: {config_name}") from e

  config.model_params = model_params

  return config
