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
"""Module for calling variants on examples using a trained keras model."""

from typing import List, Sequence, Tuple

import numpy as np
import tensorflow as tf

from deepvariant.protos import deepvariant_pb2
from deepvariant.small_model import keras_config
from deepvariant.small_model import make_small_model_examples
from third_party.nucleus.util import genomics_math
from third_party.nucleus.util import variantcall_utils


SMALL_MODEL_ID = "small_model"

GENOTYPE_BY_ENCODING = {
    make_small_model_examples.GenotypeEncoding.REF.value: (0, 0),
    make_small_model_examples.GenotypeEncoding.HET.value: (0, 1),
    make_small_model_examples.GenotypeEncoding.HOM_ALT.value: (1, 1),
}


class SmallModelVariantCaller:
  """Class responsible for calling variants on candidate examples.

  It is a stateful object that loads the trained model from `model_path`
  and holds a set of quality predicates that are applied to each call result to
  determine if we can accept to write it as CVO directly or if it should be
  passed to the regular DeepVariant model as a pileup image.
  """

  def __init__(
      self,
      classifier: tf.keras.Model,
      gq_threshold: float,
      batch_size: int,
  ):
    self.classifier = classifier
    self.gq_threshold = gq_threshold
    self.batch_size = batch_size

  @classmethod
  def from_model_path(
      cls, model_path: str, gq_threshold: float, batch_size: int
  ) -> "SmallModelVariantCaller":
    """Init class with a path to a pickled model."""
    return cls(
        keras_config.load_keras_model(model_path), gq_threshold, batch_size
    )

  def _accept_call_result(self, probability) -> bool:
    """Determine if the given probability is above the GQ threshold."""
    return (
        genomics_math.ptrue_to_bounded_phred(max(probability))
        >= self.gq_threshold
    )

  def _get_call_variant_outputs(
      self, candidate, prediction, probability
  ) -> deepvariant_pb2.CallVariantsOutput:
    """Returns a CallVariantsOutput for the given candidate and prediction."""
    variant = candidate.variant
    variant_call = candidate.variant.calls[0]
    del variant_call.genotype[:]
    variant_call.genotype.extend(GENOTYPE_BY_ENCODING[prediction])
    variantcall_utils.set_model_id(variant_call, SMALL_MODEL_ID)
    del variant.calls[:]
    variant.calls.append(variant_call)
    variant.quality = genomics_math.ptrue_to_bounded_phred(max(probability))
    return deepvariant_pb2.CallVariantsOutput(
        variant=variant,
        genotype_probabilities=probability,
        debug_info=None,
        alt_allele_indices=deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
            indices=[0]
        ),
    )

  def classify(self, examples: np.ndarray) -> List[Tuple[float, ...]]:
    """Classifies the given example."""
    predictions = []
    for i in range(0, len(examples), self.batch_size):
      predictions.extend(
          self.classifier.predict_on_batch(examples[i : i + self.batch_size])
      )
    return predictions

  def call_variants(
      self,
      candidates: Sequence[deepvariant_pb2.DeepVariantCall],
      examples: Sequence[Sequence[int]],
  ) -> Tuple[
      List[deepvariant_pb2.CallVariantsOutput],
      List[deepvariant_pb2.DeepVariantCall],
  ]:
    """Calls variants on the given examples."""
    probabilities = self.classify(np.array(examples))
    call_variant_outputs = []
    filtered_candidates = []
    for candidate, probability in zip(candidates, probabilities):
      if self._accept_call_result(probability):
        prediction = np.argmax(probability)
        call_variant_outputs.append(
            self._get_call_variant_outputs(candidate, prediction, probability)
        )
      else:
        filtered_candidates.append(candidate)

    return call_variant_outputs, filtered_candidates
