# Copyright 2026 Google LLC.
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
"""Functions for writing small model JSON files."""

import json
import os
from typing import Any, Sequence

from etils import epath
import ml_collections

from deepvariant import dv_vcf_constants

MODEL_INFO_FILENAME = "small_model_info.json"


def write_model_info_json(
    small_model_example_info_filename: str,
    model_features: Sequence[str],
    config: ml_collections.ConfigDict | None = None,
) -> None:
  """Writes the model config in the working directory."""
  working_dir = os.path.dirname(small_model_example_info_filename)
  epath.Path(working_dir).mkdir(parents=True, exist_ok=True)
  serialized_config = {}
  if config:
    serialized_config = config.to_dict()
  with epath.Path(small_model_example_info_filename).open("w") as fout:
    json.dump(
        {
            "version": dv_vcf_constants.DEEP_VARIANT_VERSION,
            "shape": len(model_features),
            "model_features": model_features,
            "config": serialized_config,
        },
        fout,
        indent=2,
    )


def write_model_info_from_config(
    checkpoint_directory: str,
    config: ml_collections.ConfigDict,
    model_features: Sequence[str],
):
  small_model_checkpoint_info_json_filename = os.path.join(
      checkpoint_directory, MODEL_INFO_FILENAME
  )
  write_model_info_json(
      small_model_checkpoint_info_json_filename,
      model_features=model_features,
      config=config,
  )


def write_model_info_from_model_features(
    small_model_examples_filename: str, model_features: Sequence[str]
):
  small_model_example_info_filename = (
      f"{small_model_examples_filename}.{MODEL_INFO_FILENAME}"
  )
  write_model_info_json(
      small_model_example_info_filename, model_features=model_features
  )


def read_model_config(
    checkpoint_path: str, optional: bool = True
) -> dict[str, Any]:
  """Reads the model config from the working directory."""
  small_model_example_info_filename = os.path.join(
      os.path.dirname(checkpoint_path), MODEL_INFO_FILENAME
  )
  if not epath.Path(small_model_example_info_filename).exists():
    if not optional:
      raise FileNotFoundError(
          f"Model config file {small_model_example_info_filename} does not"
          " exist."
      )
    return {}
  with epath.Path(small_model_example_info_filename).open("r") as fin:
    return json.load(fin)


def get_model_features_from_model_config(
    checkpoint_path: str, optional: bool = True
) -> Sequence[str]:
  """Returns the model features from the checkpoint."""
  model_config = read_model_config(checkpoint_path, optional)
  return model_config.get("model_features", [])
