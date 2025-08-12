# Copyright 2025 Google LLC.
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

"""Compares labels in two sets of tf.Example protos.

This utility takes two TFRecord files of labeled examples, presumably one
generated with HaplotypeLabeler and the other with CombinedLabeler. It compares
the labels for each variant and writes any differences to an output file.

This is useful for debugging and understanding the differences between the two
labeling algorithms.
"""

from absl import flags

from absl import app
import multiprocessing
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import variant_utils
from tensorflow.core.example import example_pb2

FLAGS = flags.FLAGS

_EXAMPLES_LABELER_A = flags.DEFINE_string(
    'labeled_examples_a',
    None,
    'TFRecord file of examples labeld with labeler A.',
)
_EXAMPLES_LABELER_B = flags.DEFINE_string(
    'labeled_examples_b',
    None,
    'TFRecord file of examples labeld with labeler B.',
)
_OUTPUT = flags.DEFINE_string(
    'output', None, 'Output file to write differently labeled variants.'
)


def get_variant_and_genotype_from_example(example_proto):
  """Extracts the variant and genotype from a tf.Example proto."""
  variant = variants_pb2.Variant()
  variant.ParseFromString(
      example_proto.features.feature['variant/encoded'].bytes_list.value[0]
  )
  genotype = tuple(variant.calls[0].genotype)
  return variant, genotype


def read_examples(tfrecord_path):
  """Reads examples from a TFRecord file and returns a label map."""
  labels = {}
  for example in tfrecord.read_tfrecords(
      tfrecord_path, proto=example_pb2.Example
  ):
    variant, genotype = get_variant_and_genotype_from_example(example)
    key = variant_utils.variant_key(variant)
    labels[key] = genotype
  return labels


def main(argv):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Too many command-line arguments.', errors.CommandLineError
      )

    print('Reading haplotype examples...')
    labeler_a_labels = read_examples(_EXAMPLES_LABELER_A.value)
    print(f'Found {len(labeler_a_labels)} haplotype examples.')

    print('Reading combined examples...')
    labeler_b_labels = read_examples(_EXAMPLES_LABELER_B.value)
    print(f'Found {len(labeler_b_labels)} combined examples.')

    with open(_OUTPUT.value, 'w') as output_file:
      all_keys = sorted(
          set(labeler_a_labels.keys()) | set(labeler_b_labels.keys())
      )
      diff_count = 0
      for key in all_keys:
        labels_a_gt = labeler_a_labels.get(key)
        labels_b_gt = labeler_b_labels.get(key)

        haplotype_gt_for_sort = (
            sorted(labels_a_gt) if labels_a_gt is not None else None
        )
        combined_gt_for_sort = (
            sorted(labels_b_gt) if labels_b_gt is not None else None
        )

        if haplotype_gt_for_sort != combined_gt_for_sort:
          diff_count += 1
          output_file.write(f'Variant: {key}\n')
          output_file.write(f'  Labeler A Genotype: {labels_a_gt}\n')
          output_file.write(f'  Labeler B Genotype:  {labels_b_gt}\n')

    print(f'Found {diff_count} differences.')
    print(f'Wrote differences to {_OUTPUT.value}')


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'labeled_examples_a',
      'labeled_examples_b',
      'output',
  ])
  app.run(main)
