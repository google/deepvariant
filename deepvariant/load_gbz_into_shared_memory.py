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
"""Load the sequneces of a GBZ file into shared memory, usable by multiple processes."""

from absl import app
from absl import flags

from deepvariant import logging_level
from third_party.nucleus.io.python import gbz_reader
from third_party.nucleus.io.python import hts_verbose
from third_party.nucleus.util import errors


# Flags related to loading GBZ into shared memory.
_PANGENOME_GBZ = flags.DEFINE_string(
    'pangenome_gbz',
    None,
    (
        'Required. Pangenome GBZ file to load into shared memory'
        '(Only the sequences are loaded into shared memory.)'
    ),
)

_NUM_SHARDS = flags.DEFINE_integer(
    'num_shards',
    None,
    (
        'Required. Number of shards that will use the shared memory. It is'
        'important to set this number correctly to make sure that the'
        'shared memory is not deleted before all processes are done using it.'
        'If this value is greater than required then the shared memory will'
        'exist even after all processes are done using it,'
        'which is not desired.'
    ),
)

_SHARED_MEMORY_NAME = flags.DEFINE_string(
    'shared_memory_name',
    'GBZ_SHARED_MEMORY',
    (
        'Name of the shared memory segment.'
    ),
)

_SHARED_MEMORY_SIZE_GB = flags.DEFINE_integer(
    'shared_memory_size_gb',
    12,
    'Size of the shared memory in GB.',
)


def load_gbz_into_shared_memory(
    pangenome_gbz: str,
    shared_memory_name: str,
    shared_memory_size_gb: int,
    num_shards: int,
):
  """Loads a GBZ file into a shared memory segment."""
  sample_name = 'GRCh38'
  context = 1000
  chrom_prefix = ''
  create_shared_memory = True
  use_loaded_shared_memory = False
  num_processes = num_shards

  ## The only parameters that are important are:
  ## - shared_memory_name,
  ## - create_shared_memory,
  ## - use_loaded_shared_memory,
  ## - shared_memory_size_gb,
  ## - num_processes
  ## The rest are set to default values. This is because the shared memory will
  ## keep only the sequences of the pangenome and not the rest of the
  ## information.
  ## Just instantating a gbz reader with create_shared_memory=True
  ## will load the sequences into shared memory.
  ## If num_processes is greater than 0 then the shared memory
  ## will NOT be deleted after this call (to be more precise
  ## after calling ~GbzReader() here).
  gbz_reader.GbzReader(
      pangenome_gbz,
      sample_name,
      context,
      chrom_prefix,
      shared_memory_name,
      create_shared_memory,
      use_loaded_shared_memory,
      shared_memory_size_gb,
      num_processes,
  )


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: load_gbz_into_shared_memory does not'
          'accept positional arguments but some are present on'
          'the command line:'
          '"{}".'.format(str(argv)),
          errors.CommandLineError,
      )
    del argv  # Unused.

    logging_level.set_from_flag()
    hts_verbose.set(hts_verbose.htsLogLevel.HTS_LOG_WARNING)

    load_gbz_into_shared_memory(
        _PANGENOME_GBZ.value,
        _SHARED_MEMORY_NAME.value,
        _SHARED_MEMORY_SIZE_GB.value,
        _NUM_SHARDS.value,
    )


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'pangenome_gbz',
      'num_shards',
  ])
  app.run(main)
