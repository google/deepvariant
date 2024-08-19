# Copyright 2021 Google LLC.
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
r"""Experimental multi-sample make_examples for DeepVariant.

This is a prototype for experimentation with multiple samples in DeepVariant, a
proof of concept enabled by a refactoring to join together DeepVariant and
DeepTrio, generalizing the functionality of make_examples to work with multiple
samples.

The output of this script is not compatible with any of DeepVariant's public
models, and the DeepVariant team does not intend to provide support for users
of this script.

Example usage:

multisample_make_examples \
  --mode calling \
  --ref "reference.fa" \
  --reads "sample1.bam;sample2.bam;sample3.bam;sample4.bam;sample5.bam" \
  --sample_names "sample1;sample2;sample3;sample4;sample5" \
  --examples "examples.tfrecord.gz" \
  --pileup_image_heights "20;20;20;20;20" \ # optional
  --downsample_fractions "0.5;0.5;0.5;0.5;0.5" # optional
"""

import os

from absl import app
from absl import flags

from deepvariant import logging_level
from deepvariant import make_examples_core
from deepvariant import make_examples_options
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.util import errors
from third_party.nucleus.util import proto_utils

MAIN_SAMPLE_INDEX = 0  # This is the first sample listed in --reads.

FLAGS = flags.FLAGS

# Adopt more general flags from make_examples_options.
flags.adopt_module_key_flags(make_examples_options)

# Define flags specific to multi-sample make_examples.
flags.DEFINE_string(
    'reads',
    None,
    (
        'Required. A list of BAM/CRAM files, with different samples separated'
        ' by semi-colons. At least one aligned, sorted, indexed BAM/CRAM file'
        ' is required for each sample. All must be aligned to the same'
        ' reference genome compatible with --ref. Can provide multiple BAMs'
        ' (comma-separated) for each sample. Format is, for example:'
        ' sample1;sample2_BAM1,sample2_BAM2;sample3 '
    ),
)
flags.DEFINE_string(
    'sample_names',
    'DEFAULT',
    (
        'Sample names corresponding to the samples (must match order and length'
        ' of samples in --reads). Separate names for each sample with'
        ' semi-colons, e.g. "sample1;sample2;sample3". If not specified, (i.e.'
        ' "sample1;;sample3" or even ";;") any sample without a sample name'
        ' from this flag the will be inferred from the header information from'
        ' --reads.'
    ),
)
flags.DEFINE_string(
    'downsample_fractions',
    'DEFAULT',
    (
        'If not empty string ("") must be a value between 0.0 and 1.0. Reads'
        ' will be kept (randomly) with a probability of downsample_fraction'
        ' from the input sample BAMs. This argument makes it easy to create'
        ' examples as though the input BAM had less coverage. Similar to'
        ' --reads and --sample_name, supply different values for each sample by'
        ' separating them with semi-colons, where the order of samples is the'
        ' same as in --reads.'
    ),
)
flags.DEFINE_string(
    'pileup_image_heights',
    'DEFAULT',
    (
        'Height for the part of the pileup image showing reads from each'
        ' sample. By default, use a height of 100 for all samples. Similar to'
        ' --reads and --sample_name, supply different values for each sample by'
        ' separating them with semi-colons, where the order of samples is the'
        ' same as in --reads.'
    ),
)


def n_samples_from_flags(add_flags=True, flags_obj=None):
  """Collects sample-related options into a list of samples."""

  n_reads = flags_obj.reads.split(';')

  num_samples = len(n_reads)
  flags_organized = {}
  for flag_name in [
      'reads',
      'sample_names',
      'downsample_fractions',
      'pileup_image_heights',
  ]:
    if flags_obj[flag_name].value != 'DEFAULT':
      flags_organized[flag_name] = flags_obj[flag_name].value.split(';')
      if len(flags_organized[flag_name]) != num_samples:
        raise ValueError(
            f'--{flag_name} has {len(flags_organized[flag_name])} '
            'samples, but it should be matching the number of '
            f'samples in --reads, which was {num_samples}.'
        )
    else:
      flags_organized[flag_name] = [''] * num_samples

  n_sample_options = []
  for i in range(num_samples):
    sample_name = make_examples_core.assign_sample_name(
        sample_name_flag=flags_organized['sample_names'][i],
        reads_filenames=flags_organized['reads'][i],
    )

    n_sample_options.append(
        deepvariant_pb2.SampleOptions(
            role=str(i),
            name=sample_name,
            variant_caller_options=make_examples_core.make_vc_options(
                sample_name=sample_name, flags_obj=flags_obj
            ),
            order=range(num_samples),
            pileup_height=100,
        )
    )

  if add_flags:
    for i in range(num_samples):
      n_sample_options[i].reads_filenames.extend(
          flags_organized['reads'][i].split(',')
      )
      if flags_organized['downsample_fractions'][i]:
        n_sample_options[i].downsample_fraction = float(
            flags_organized['downsample_fractions'][i]
        )
      if flags_organized['pileup_image_heights'][i]:
        n_sample_options[i].pileup_height = int(
            flags_organized['pileup_image_heights'][i]
        )

  # Ordering here determines the default order of samples, and when a sample
  # above has a custom .order, then this is the list those indices refer to.
  samples_in_order = n_sample_options
  sample_role_to_train = '0'
  return samples_in_order, sample_role_to_train


def default_options(add_flags=True, flags_obj=None):
  """Creates a MakeExamplesOptions proto populated with reasonable defaults.

  Args:
    add_flags: bool. defaults to True. If True, we will push the value of
      certain FLAGS into our options. If False, those option fields are left
      uninitialized.
    flags_obj: object.  If not None, use as the source of flags, else use global
      FLAGS.

  Returns:
    deepvariant_pb2.MakeExamplesOptions protobuf.

  Raises:
    ValueError: If we observe invalid flag values.
  """
  if not flags_obj:
    flags_obj = FLAGS

  samples_in_order, sample_role_to_train = n_samples_from_flags(
      add_flags=add_flags, flags_obj=flags_obj
  )

  options = make_examples_options.shared_flags_to_options(
      add_flags=add_flags,
      flags_obj=flags_obj,
      samples_in_order=samples_in_order,
      sample_role_to_train=sample_role_to_train,
      main_sample_index=MAIN_SAMPLE_INDEX,
  )

  if add_flags:
    options.bam_fname = '|'.join(
        [os.path.basename(x) for x in flags_obj.reads.split(';')]
    )

  return options


def check_options_are_valid(options):
  """Checks that all the options chosen make sense together."""

  # Check for general flags (shared for DeepVariant and DeepTrio).
  make_examples_options.check_options_are_valid(
      options, main_sample_index=MAIN_SAMPLE_INDEX
  )

  sample_names = [s.name for s in options.sample_options]
  if len(sample_names) != len(set(sample_names)):
    raise ValueError('--sample_names cannot contain duplicate names.')


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: make_examples does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)),
          errors.CommandLineError,
      )
    del argv  # Unused.

    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()
    hts_verbose.set(hts_verbose.htsLogLevel.HTS_LOG_WARNING)

    # Set up options; may do I/O.
    options = default_options(add_flags=True, flags_obj=FLAGS)
    check_options_are_valid(options)

    # Run!
    make_examples_core.make_examples_runner(options)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'mode',
      'reads',
      'ref',
  ])
  app.run(main)
