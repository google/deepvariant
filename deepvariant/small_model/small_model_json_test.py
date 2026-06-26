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
import json
import os

from absl.testing import absltest
import ml_collections

from deepvariant import dv_vcf_constants
from deepvariant.small_model import small_model_json


class SmallModelJsonTest(absltest.TestCase):

  def test_write_and_read_model_info(self):
    temp_dir = self.create_tempdir().full_path
    json_path = os.path.join(temp_dir, 'test_info.json')
    features = ['feat1', 'feat2']
    config = ml_collections.ConfigDict({'param1': 'value1'})

    small_model_json.write_model_info_json(json_path, features, config)

    # Verify file content directly
    with open(json_path, 'r') as f:
      data = json.load(f)
    self.assertEqual(data['version'], dv_vcf_constants.DEEP_VARIANT_VERSION)
    self.assertEqual(data['shape'], 2)
    self.assertEqual(data['model_features'], features)
    self.assertEqual(data['config'], {'param1': 'value1'})

    # Test reading back
    checkpoint_dir = self.create_tempdir().full_path
    checkpoint_path = os.path.join(checkpoint_dir, 'model.ckpt')
    info_path = os.path.join(
        checkpoint_dir, small_model_json.MODEL_INFO_FILENAME
    )

    small_model_json.write_model_info_json(info_path, features, config)

    read_config = small_model_json.read_model_config(checkpoint_path)
    self.assertEqual(read_config['model_features'], features)
    self.assertEqual(read_config['config'], {'param1': 'value1'})

    read_features = small_model_json.get_model_features_from_model_config(
        checkpoint_path
    )
    self.assertEqual(read_features, features)

  def test_write_model_info_from_config(self):
    checkpoint_dir = self.create_tempdir().full_path
    config = ml_collections.ConfigDict({'param1': 'value1'})
    features = ['feat1', 'feat2']

    small_model_json.write_model_info_from_config(
        checkpoint_dir, config, features
    )

    info_path = os.path.join(
        checkpoint_dir, small_model_json.MODEL_INFO_FILENAME
    )
    self.assertTrue(os.path.exists(info_path))
    with open(info_path, 'r') as f:
      data = json.load(f)
    self.assertEqual(data['model_features'], features)
    self.assertEqual(data['config'], {'param1': 'value1'})

  def test_write_model_info_from_model_features(self):
    temp_dir = self.create_tempdir().full_path
    examples_path = os.path.join(temp_dir, 'examples.tfrecord')
    features = ['feat1', 'feat2']

    small_model_json.write_model_info_from_model_features(
        examples_path, features
    )

    info_path = f'{examples_path}.{small_model_json.MODEL_INFO_FILENAME}'
    self.assertTrue(os.path.exists(info_path))
    with open(info_path, 'r') as f:
      data = json.load(f)
    self.assertEqual(data['model_features'], features)
    self.assertEqual(data['config'], {})

  def test_read_model_config_missing_optional(self):
    checkpoint_dir = self.create_tempdir().full_path
    checkpoint_path = os.path.join(checkpoint_dir, 'model.ckpt')
    self.assertEqual(
        small_model_json.read_model_config(checkpoint_path, optional=True), {}
    )

  def test_read_model_config_missing_required(self):
    checkpoint_dir = self.create_tempdir().full_path
    checkpoint_path = os.path.join(checkpoint_dir, 'model.ckpt')
    with self.assertRaises(FileNotFoundError):
      small_model_json.read_model_config(checkpoint_path, optional=False)


if __name__ == '__main__':
  absltest.main()
