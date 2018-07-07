# Copyright 2017 Google Inc.
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
"""Utility functions for creating TensforFlow mock models.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

from absl import flags
import tensorflow as tf
from deepvariant import dv_constants
from deepvariant import modeling

FLAGS = flags.FLAGS
# Constants dictating moving average.
_MOVING_AVERAGE_DECAY = 0.995


def write_fake_checkpoint(model_name,
                          session,
                          checkpoint_dir,
                          moving_average_decay=_MOVING_AVERAGE_DECAY,
                          name='model'):
  """Writes a fake TensorFlow checkpoint to checkpoint_dir."""
  path = os.path.join(checkpoint_dir, name)
  with session as sess:
    model = modeling.get_model(model_name)
    # Needed to protect ourselves for models without an input image shape.
    h, w = getattr(
        model, 'input_image_shape',
        (dv_constants.PILEUP_DEFAULT_HEIGHT, dv_constants.PILEUP_DEFAULT_WIDTH))
    images = tf.placeholder(
        tf.float32, shape=(4, h, w, dv_constants.PILEUP_NUM_CHANNELS))
    model.create(images, num_classes=3, is_training=True)
    # This is gross, but necessary as model_eval assumes the model was trained
    # with model_train which uses exp moving averages. Unfortunately we cannot
    # just call into model_train as it uses FLAGS which conflict with the
    # flags in use by model_eval. So we inline the creation of the EMA here.
    variable_averages = tf.train.ExponentialMovingAverage(
        moving_average_decay, tf.train.get_or_create_global_step())
    tf.add_to_collection(
        tf.GraphKeys.UPDATE_OPS,
        variable_averages.apply(tf.contrib.framework.get_model_variables()))
    sess.run(tf.global_variables_initializer())
    save = tf.train.Saver(tf.contrib.framework.get_variables())
    save.save(sess, path)
  return path


def test_tmpdir(name):
  """Returns a path to a temp directory name in the test tmpdir.

  Args:
    name: str; the name of the file, should not contain any slashes.

  Returns:
    str path to a tmpfile with filename name in our test tmpfile directory.
  """
  dirname = os.path.join(tf.test.get_temp_dir(), name)
  tf.gfile.MakeDirs(dirname)
  return dirname


# redacted
# taking a required_variables_regexps and have a function that looks like:
# return all(any(re.match(var, pat) for var in var_to_shape_map.keys()
#   for pat in required_variable_regexps))
def check_equals_checkpoint_top_scopes(checkpoint, list_of_scopes):
  """Returns true if list_of_scopes equals the top scopes in checkpoint."""
  reader = tf.train.NewCheckpointReader(checkpoint)
  var_to_shape_map = reader.get_variable_to_shape_map()
  top_scopes_in_ckpt = set([x.split('/')[0] for x in var_to_shape_map.keys()])
  return top_scopes_in_ckpt == set(list_of_scopes)
