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
"""Module for training a smaller DeepVariant model for SNPs."""

from absl import app
from absl import flags
from ml_collections.config_flags import config_flags

from deepvariant import logging_level
from deepvariant.small_model import train_utils


_CONFIG = config_flags.DEFINE_config_file(
    "config", None, "Training configuration."
)
_MODEL_PATH = flags.DEFINE_string(
    "model_path",
    None,
    "Path to where keras model will be saved.",
    required=True,
)
_TRUTH_VCF = flags.DEFINE_string("truth_vcf", None, "Path to truth VCF file.")
_CALCULATE_FEATURE_IMPORTANCE = flags.DEFINE_boolean(
    "calculate_feature_importance",
    False,
    "Whether to calculate feature importance.",
)
_SHOW_MISTAKES_ABOVE_GQ_THRESHOLD = flags.DEFINE_integer(
    "show_mistakes_above_threshold",
    -1,
    "If set, log mistakes above this GQ threshold.",
)


def main(unused_argv):
  logging_level.set_from_flag()
  config = _CONFIG.value

  # Load training data
  candidates = train_utils.load_training_data(
      config.tsv_directory, config.max_training_samples
  )
  # Remove mislabels if truth VCF is provided.
  if _TRUTH_VCF.value:
    candidates = train_utils.sanitize_mislabels(_TRUTH_VCF.value, candidates)

  # Split into training and test data
  df_train, df_test = train_utils.split_training_data(
      candidates, config.test_fraction
  )

  # Init keras model
  model = train_utils.keras_mlp_model(config.model_params)

  # Run training
  model_runner = train_utils.ModelRunner(
      model,
      df_train,
      df_test,
      epochs=config.epochs,
      batch_size=config.batch_size,
  )
  model_runner.run()

  # Save model
  model.save(f"{_MODEL_PATH.value}.keras", save_format="tf")

  # Optional: run post-training analysis.
  if _CALCULATE_FEATURE_IMPORTANCE.value:
    model_runner.calculate_feature_importance()
  if _SHOW_MISTAKES_ABOVE_GQ_THRESHOLD.value > 0:
    model_runner.show_mistakes_above_threshold(
        _SHOW_MISTAKES_ABOVE_GQ_THRESHOLD.value
    )


if __name__ == "__main__":
  app.run(main)
