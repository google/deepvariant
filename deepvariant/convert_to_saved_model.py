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
r"""Convert a Checkpoint to a SavedModel.

Example command:
  convert_to_saved_model --checkpoint=/path/to/checkpoint --output=/tmp/output
"""

import os
from typing import Optional

from absl import app
from absl import flags
from absl import logging
import tensorflow as tf

from deepvariant import dv_utils
from deepvariant import keras_modeling as modeling
from tensorflow.python.platform import gfile

# Outputs:
_OUTPUT = flags.DEFINE_string('output', None, 'Output SavedModel name.')

# Model checkpoint:
_CHECKPOINT = flags.DEFINE_string(
    'checkpoint',
    None,
    (
        'Path to checkpoint directory + prefix. '
        'For example: <path/to/model>/checkpoint-50.'
    ),
)

_EXAMPLE_INFO_JSON = flags.DEFINE_string(
    'example_info_json',
    None,
    (
        'Path to json file containing example information. '
        'For example: <path/to/model>/example.info.json.'
    ),
)


def register_required_flags():
  flags.mark_flags_as_required([
      'checkpoint',
      'example_info_json',
      'output',
  ])


def initialize_model(
    example_info_json: str, checkpoint_path: str
) -> Optional[tf.keras.Model]:
  """Initializes the model and gathers parameters.

  Args:
    example_info_json: Path to json file containing example shape.
    checkpoint_path: Path to model checkpoint.

  Returns:
    An initialized model.
  """
  logging.info('Reading example shape from %s', example_info_json)
  example_shape = dv_utils.get_shape_and_channels_from_json(example_info_json)[
      0
  ]
  logging.info('Loading %s', checkpoint_path)
  logging.info('Example shape %s', example_shape)
  # Load model
  model = modeling.inceptionv3(example_shape, init_backbone_with_imagenet=False)
  # model.load_weights(checkpoint_path).expect_partial()
  # checkpoint = tf.train.Checkpoint(model=model)
  # Note that the `print_model_summary` is necessary because we need to run a
  # forward pass with the model in order for assert_existing_objects_matched to
  # work as expected. If you don't do this, then assert_existing_objects_matched
  # will not raise an error even if the wrong checkpoint is used.
  # Some context here: internal.
  input_shape = (1, example_shape[0], example_shape[1], example_shape[2])
  modeling.print_model_summary(model, input_shape)
  model.load_weights(
      checkpoint_path
  ).expect_partial().assert_existing_objects_matched()

  # checkpoint.restore(
  #     checkpoint_path
  # ).expect_partial().assert_existing_objects_matched()

  logging.info('Finished initialize_model.')
  return model


def main(_):
  """Main entry point."""
  loaded_model = initialize_model(
      example_info_json=_EXAMPLE_INFO_JSON.value,
      checkpoint_path=_CHECKPOINT.value,
  )
  tf.saved_model.save(loaded_model, _OUTPUT.value)
  # Copy over the example_info.json.
  gfile.Copy(
      _EXAMPLE_INFO_JSON.value,
      os.path.join(_OUTPUT.value, 'example_info.json'),
      overwrite=True,
  )


if __name__ == '__main__':
  register_required_flags()
  app.run(main)
