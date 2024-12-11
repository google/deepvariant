# Copyright 2022 Google LLC.
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
"""Setting up metrics for use in the DeepVariant model."""

import enum
import ml_collections
import tensorflow as tf
import tensorflow_addons as tfa


class MetricSubset(enum.Enum):
  ALL = ''
  SNP = 'snp'
  INDEL = 'indel'


class F1ScorePerClass(tfa.metrics.F1Score):
  """Reports F1 Score for a target class."""

  def __init__(self, num_classes: int, target_class: int, name: str):
    self.target_class = target_class
    super().__init__(num_classes=num_classes, name=name)

  def result(self) -> tf.Tensor:
    return super().result()[self.target_class]


def create_metrics(
    config: ml_collections.ConfigDict,
) -> list[tf.keras.metrics.Metric]:
  metrics = create_full_metrics()
  if config.include_snp_indel_metrics:
    metrics.extend(create_limited_metrics(MetricSubset.SNP))
    metrics.extend(create_limited_metrics(MetricSubset.INDEL))
  # Leave mean loss as the last metric as it is updated differently.
  metrics.append(tf.keras.metrics.Mean(name='loss'))
  return metrics


def create_full_metrics() -> list[tf.keras.metrics.Metric]:
  return [
      tf.keras.metrics.CategoricalAccuracy(name='categorical_accuracy'),
      tf.keras.metrics.CategoricalCrossentropy(
          name='categorical_cross_entropy'
      ),
      tf.keras.metrics.TruePositives(name='true_positives'),
      tf.keras.metrics.TrueNegatives(name='true_negatives'),
      tf.keras.metrics.FalsePositives(name='false_positives'),
      tf.keras.metrics.FalseNegatives(name='false_negatives'),
      tf.keras.metrics.Precision(name='precision'),
      tf.keras.metrics.Precision(name='precision_homref', class_id=0),
      tf.keras.metrics.Precision(name='precision_het', class_id=1),
      tf.keras.metrics.Precision(name='precision_homalt', class_id=2),
      tf.keras.metrics.Recall(name='recall'),
      tf.keras.metrics.Recall(name='recall_homref', class_id=0),
      tf.keras.metrics.Recall(name='recall_het', class_id=1),
      tf.keras.metrics.Recall(name='recall_homalt', class_id=2),
      tfa.metrics.F1Score(
          num_classes=3,
          average='weighted',
          name='f1_weighted',
      ),
      tfa.metrics.F1Score(num_classes=3, average='micro', name='f1_micro'),
      tfa.metrics.F1Score(num_classes=3, average='macro', name='f1_macro'),
      F1ScorePerClass(num_classes=3, target_class=0, name='f1_homref'),
      F1ScorePerClass(num_classes=3, target_class=1, name='f1_het'),
      F1ScorePerClass(num_classes=3, target_class=2, name='f1_homalt'),
  ]


def create_limited_metrics(
    metric_subset: MetricSubset,
) -> list[tf.keras.metrics.Metric]:
  metric_name_suffix = (
      '' if metric_subset == MetricSubset.ALL else f'_{metric_subset.value}'
  )
  return [
      tf.keras.metrics.CategoricalAccuracy(
          name=f'categorical_accuracy{metric_name_suffix}'
      ),
      tfa.metrics.F1Score(
          num_classes=3,
          average='weighted',
          name=f'f1_weighted{metric_name_suffix}',
      ),
  ]
