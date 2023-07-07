# Copyright 2018 Google LLC.
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
"""Utility functions for working with sharded files.

A sharded file is a single conceptual file that is broken into a collection
of files to make parallelization easier.  A sharded file spec is like a
filename for a sharded file; the file spec "/some/path/prefix@200.txt"
says that the sharded file consists of 200 actual files that have names like
"/some/path/prefix-00000-of-00200.txt", "/some/path/prefix-00001-of-00200.txt",
etc.  This module contains functions for parsing, generating, detecting and
resolving sharded file specs.
"""

# Important: Please keep this module free of TensorFlow c++ extensions.
# This makes it easy to build pure python packages for training that work with
# CMLE.

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math
import os
import re
from typing import List, Tuple

from etils import epath


import six

SHARD_SPEC_PATTERN = re.compile(R'((.*)\@(\d*[1-9]\d*)(?:\.(.+))?)')
SHARD_FILE_PATTERN = re.compile(R'(.*)-(\d+)-of-(\d*[1-9]\d*)([^/]+)?$')


class ShardError(Exception):
  """An I/O error."""


def parse_sharded_file_spec(spec):
  """Parse a sharded file specification.

  Args:
    spec: str. The sharded file specification. A sharded file spec is one like
      'gs://some/file@200.txt'. Here, '@200' specifies the number of shards.

  Returns:
    basename: str. The basename for the files.
    num_shards: int >= 0. The number of shards.
    suffix: str. The suffix if there is one, or '' if not.
  Raises:
    ShardError: If the spec is not a valid sharded specification.
  """
  m = SHARD_SPEC_PATTERN.match(spec)
  if not m:
    raise ShardError(('The file specification {0} is not a sharded file '
                      'specification because it did not match the regex '
                      '{1}').format(spec, SHARD_SPEC_PATTERN.pattern))

  # If there's a non-empty suffix, we need to prepend '.' so we get files like
  # foo@20.ext instead of foo@ext. The original C++ parser version has:
  # string ext = StrCat(suff.empty() ? "" : ".", suff);
  suffix = '.' + m.group(4) if m.group(4) else ''

  return m.group(2), int(m.group(3)), suffix


def _shard_width(num_shards):
  """Return the width of the shard matcher based on the number of shards."""
  return max(5, int(math.floor(math.log10(num_shards)) + 1))


def generate_sharded_filenames(spec: str) -> List[str]:
  """Generate the list of filenames corresponding to the sharding path.

  Args:
    spec: str. Represents a filename with a sharding specification.
      e.g., 'gs://some/file@200.txt' represents a file sharded 200 ways.

  Returns:
    List of filenames.

  Raises:
    ShardError: If spec is not a valid sharded file specification.
  """
  basename, num_shards, suffix = parse_sharded_file_spec(spec)
  files = []
  width = _shard_width(num_shards)
  format_str = '{{0}}-{{1:0{0}}}-of-{{2:0{0}}}{{3}}'.format(width)
  for i in range(num_shards):
    files.append(format_str.format(basename, i, num_shards, suffix))

  return files


def glob_list_sharded_file_patterns(
    comma_separated_patterns: str, sep: str = ','
) -> List[str]:
  """Generate list of filenames corresponding to `comma_separated_patterns`.

  Args:
    comma_separated_patterns: str. A pattern or a comma-separated list of
      patterns that represent file names.
    sep: char. Separator character.

  Returns:
    List of filenames, sorted and dedupped.
  """
  return sorted(set([
      os.fspath(f)
      for pattern in comma_separated_patterns.split(sep)
      # Note: epath.Path().glob(normalize_to_sharded_file_pattern(pattern))
      # does not work with gs://.
      for f in epath.Path(
          os.path.dirname(normalize_to_sharded_file_pattern(pattern))).glob(
              os.path.basename(normalize_to_sharded_file_pattern(pattern)))
  ]))


def generate_sharded_file_pattern(basename, num_shards, suffix):
  """Generate a sharded file pattern.

  Args:
    basename: str. The basename for the files.
    num_shards: int. The number of shards.
    suffix: str. The suffix if there is one or ''.
  Returns:
    pattern:
  """
  width = _shard_width(num_shards)
  specifier = '?' * width
  format_str = '{{0}}-{{1}}-of-{{2:0{0}}}{{3}}'.format(width)
  return format_str.format(basename, specifier, num_shards, suffix)


def normalize_to_sharded_file_pattern(spec_or_pattern):
  """Take a sharding spec or sharding file pattern and return a sharded pattern.

  The input can be a sharding spec(e.g '/some/file@10') or a sharded file
  pattern (e.g. '/some/file-?????-of-00010)

  Args:
    spec_or_pattern: str. A sharded file specification or sharded file pattern.

  Returns:
    A sharded file pattern.
  """
  try:
    basename, num_shards, suffix = parse_sharded_file_spec(spec_or_pattern)
  except ShardError:
    return spec_or_pattern
  return generate_sharded_file_pattern(basename, num_shards, suffix)


def is_sharded_filename(filename: str) -> bool:
  """Returns True if filename is a sharded filename."""
  m = SHARD_FILE_PATTERN.match(filename)
  return m is not None


def is_sharded_file_spec(spec: str) -> bool:
  """Returns True if spec is a sharded file specification."""
  m = SHARD_SPEC_PATTERN.match(spec)
  return m is not None


# TODO: retire when GenerateShardedFilename is added to library.
def sharded_filename(spec: str, i: int) -> str:
  """Gets a path appropriate for writing the ith file of a sharded spec."""
  return generate_sharded_filenames(spec)[i]


# TODO: Improve the return value (instead of using tuple). It hurts
# readability when there are multiple input filespecs.
def resolve_filespecs(
    shard: int, *filespecs: str
):
  """Transforms potentially sharded filespecs into their paths for single shard.

  This function takes a shard number and a varargs of potentially-sharded
  filespecs, and returns a list where the filespecs have been resolved into
  concrete file paths for a single shard.

  This function has a concept of a master filespec, which is used to constrain
  and check the validity of other filespecs. The first filespec is considered
  the master, and it cannot be None. For example, if master is not sharded, none
  of the other specs can be sharded, and vice versa. They must all also have a
  consistent sharding (e.g., master is @10, then all others must be @10).

  Note that filespecs (except the master) may be None or any other False value,
  which are returned as-is in the output list.

  Args:
    shard: int >= 0. Our shard number.
    *filespecs: list[str]. Contains all of the filespecs we want to resolve into
      shard-specific file paths.

  Returns:
    A list. The first element is the number of shards, which is an int >= 1 when
    filespecs contains sharded paths and 0 if none do. All subsequent
    returned values follow the shard-specific paths for each filespec, in order.

  Raises:
    ValueError: if any filespecs are inconsistent.
  """
  if not filespecs:
    raise ValueError('filespecs must have at least one element.')

  master = filespecs[0]
  master_is_sharded = is_sharded_file_spec(master)

  master_num_shards = 0
  if master_is_sharded:
    _, master_num_shards, _ = parse_sharded_file_spec(master)
    if shard >= master_num_shards or shard < 0:
      raise ValueError('Invalid shard={} value with master={} sharding'.format(
          shard, master))
  elif shard > 0:
    raise ValueError('Output is not sharded but shard > 0: {}'.format(shard))

  def resolve_one(filespec):
    """Resolves a single filespec into a concrete filepath."""
    if not filespec:
      return filespec

    is_sharded = is_sharded_file_spec(filespec)
    if master_is_sharded != is_sharded:
      raise ValueError('Master={} and {} have inconsistent sharding'.format(
          master, filespec))

    if not is_sharded:  # Not sharded => filespec is the actual filename.
      return filespec

    _, filespec_num_shards, _ = parse_sharded_file_spec(filespec)
    if filespec_num_shards != master_num_shards:
      raise ValueError('Master={} and {} have inconsistent sharding'.format(
          master, filespec))
    return sharded_filename(filespec, shard)

  return [master_num_shards] + [resolve_one(spec) for spec in filespecs]


def maybe_generate_sharded_filenames(filespec):
  """Potentially expands sharded filespec into a list of paths.

  This function takes in a potentially sharded filespec and expands it into a
  list containing the full set of corresponding concrete sharded file paths. If
  the input filespec is not sharded then a list containing just that file path
  is returned. This function is useful, for example, when the input to a binary
  can either be sharded or not.

  Args:
    filespec: String. A potentially sharded filespec to expand.

  Returns:
    A list of file paths.

  Raises:
    TypeError: if filespec is not in valid string_types.
  """
  if not isinstance(filespec, six.string_types):
    raise TypeError('Invalid filespec: %s' % filespec)
  if is_sharded_file_spec(filespec):
    return generate_sharded_filenames(filespec)
  else:
    return [filespec]


def parse_sharded_filename(filename):
  """Parse a sharded name specification.

  Args:
    filename: str. The sharded filename. A sharded filename is one like
      'gs://some/file-00000-of-00200.txt'.

  Returns:
    basename: str. The first part before the -ddddd-of-ddddd part.
    shard_index: str: The first part of -ddddd-of-ddddd
    num_shards: str: The number of shards. (Second part after "of")
    suffix: str. Any string after the -ddddd-of-ddddd part.
  Raises:
    ShardError: If the spec is not a valid sharded specification.
  """
  m = SHARD_FILE_PATTERN.match(filename)
  if not m:
    raise ShardError(('The file specification {0} is not a sharded file '
                      'specification because it did not match the regex '
                      '{1}').format(filename, SHARD_FILE_PATTERN.pattern))
  suffix = '' if m.group(4) is None else m.group(4)
  return m.group(1), m.group(2), m.group(3), suffix
