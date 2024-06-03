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

import ml_collections


def set_wgs_config(config: ml_collections.ConfigDict) -> None:
  # TODO: update training data
  config.train_tsv_directory = "/cns/oz-d/home/brain-genomics/shafin/deepvariant_wgs/b338442845_wgs_small_model/ttl=14d/"
  config.num_train_samples = 10_000_000


def set_ont_config(config: ml_collections.ConfigDict) -> None:
  # TODO: update training data
  config.train_tsv_directory = (
      "/cns/oz-d/home/brain-genomics/lucasbrambrink/smallmodel/examples/ont/"
  )
  config.num_train_samples = 1_000_000


def set_pacbio_config(config: ml_collections.ConfigDict) -> None:
  # TODO: update training data
  config.train_tsv_directory = (
      "/cns/oz-d/home/brain-genomics/lucasbrambrink/smallmodel/examples/pacbio/"
  )
  config.num_train_samples = 1_000_000


def get_config(config_name: str) -> ml_collections.ConfigDict:
  """Returns the default configuration as instance of ConfigDict."""
  # Model hyperparameters
  model_params = ml_collections.ConfigDict()
  model_params.activation = "relu"
  model_params.hidden_layer_sizes = (100, 100)
  model_params.optimizer = "adam"

  # Training parameters
  config = ml_collections.ConfigDict()
  config.epochs = 5
  config.batch_size = 32
  config.logging_frequency = 8192
  config.train_tsv_directory = ""
  config.num_train_samples = 1_000_000
  config.tune_tsv_directory = ""
  config.num_tune_samples = 10_000
  config.test_fraction = 0.2

  # Training environment
  config.tensorboard_directory = None
  config.is_xmanager_run = False

  if config_name == "wgs":
    set_wgs_config(config)
  elif config_name == "ont":
    set_ont_config(config)
  elif config_name == "pacbio":
    set_pacbio_config(config)
  else:
    raise ValueError(f"Unknown config name: {config_name}")

  config.model_params = model_params

  return config
