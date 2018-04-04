# Copyright 2018 Google Inc.
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

"""Utilities to help with testing code."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os


from tensorflow import flags
from absl.testing import absltest

from third_party.nucleus.protos import position_pb2
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import cigar as _cigar
from tensorflow.python.platform import gfile

FLAGS = flags.FLAGS

# In the OSS version these will be ''.
DATADIR = ''
DEFAULT_WORKSPACE = ''

# In the OSS version this becomes 'deepvariant/testdata'
RELATIVE_TESTDATA_PATH = 'third_party/nucleus/testdata'


def genomics_testdata(path, datadir=DATADIR):
  """Gets the path to a testdata file in genomics at relative path.

  Args:
    path: A path to a testdata file *relative* to the genomics root
      directory. For example, if you have a test file in
      "datadir/deepvariant/core/testdata/foo.txt", path should be
      "deepvariant/core/testdata/foo.txt" to get a path to it.
    datadir: The path of the genomics root directory *relative* to
      the testing source directory.

  Returns:
    The absolute path to a testdata file.
  """
  if hasattr(FLAGS, 'test_srcdir'):
    # Google code uses FLAG.test_srcdir
    # TensorFlow uses a routine googletest.test_src_dir_path.
    test_workspace = os.environ.get('TEST_WORKSPACE', DEFAULT_WORKSPACE)
    test_srcdir = os.path.join(FLAGS.test_srcdir, test_workspace)
  else:
    # In bazel TEST_SRCDIR points at the runfiles directory, and
    # TEST_WORKSPACE names the workspace.  We need to append to the
    # path the name of the workspace in order to get to the root of our
    # source tree.
    test_workspace = os.environ['TEST_WORKSPACE']
    test_srcdir = os.path.join(os.environ['TEST_SRCDIR'], test_workspace)
  return os.path.join(test_srcdir, datadir, path)


# redacted
def genomics_core_testdata(filename):
  """Gets the path to a testdata named filename in util/testdata.

  Args:
    filename: The name of a testdata file in the core genomics testdata
      directory. For example, if you have a test file in
      "third_party/nucleus/util/testdata/foo.txt", filename should be
      "foo.txt" to get a path to it.

  Returns:
    The absolute path to a testdata file.
  """
  return genomics_testdata(os.path.join(RELATIVE_TESTDATA_PATH, filename))


def test_tmpfile(name, contents=None):
  """Returns a path to a tempfile named name in the test_tmpdir.

  Args:
    name: str; the name of the file, should not contain any slashes.
    contents: bytes, or None. If not None, tmpfile's contents will be set to
      contents before returning the path.

  Returns:
    str path to a tmpfile with filename name in our test tmpfile directory.
  """
  path = os.path.join(absltest.get_default_test_tmpdir(), name)
  if contents is not None:
    with gfile.FastGFile(path, 'wb') as fout:
      fout.write(contents)
  return path


def set_list_values(list_value, values):
  """Sets a ListValue to have the values in values."""

  def format_one(value):
    if isinstance(value, str):
      return struct_pb2.Value(string_value=value)
    elif isinstance(value, float):
      return struct_pb2.Value(number_value=value)
    elif isinstance(value, (int, long)):
      return struct_pb2.Value(int_value=value)
    else:
      raise ValueError('Unsupported type ', value)

  del list_value.values[:]
  list_value.values.extend([format_one(value) for value in values])
  # list_value.values.extend(vals)


def make_variant(chrom='chr1',
                 start=10,
                 alleles=None,
                 end=None,
                 filters=None,
                 qual=None,
                 gt=None,
                 gq=None,
                 sample_name=None,
                 gls=None):
  """Creates a new Variant proto from args.

  Args:
    chrom: str. The reference_name for this variant. Defaults to 'chr1'.
    start: int. The starting position of this variant. Defaults to 10.
    alleles: list of str with at least one element. alleles[0] is the reference
      bases and alleles[1:] will be set to alternate_bases of variant. If None,
      defaults to ['A', 'C'].
    end: int or None. If not None, the variant's end will be set to this value.
      If None, will be set to the start + len(reference_bases).
    filters: str, list of str, or None. Sets the filters field of the variant to
      this value if not None. If filters is a string `value`, this is equivalent
      to an argument [`value`]. If None, no value will be assigned to the
      filters field.
    qual: int or None. The quality score for this variant. If None, no quality
      score will be written in the Variant.
    gt: A list of ints, or None. If present, creates a VariantCall in Variant
      with genotype field set to this value. The special 'DEFAULT' value, if
      provided, will set the genotype to [0, 1]. This is the default behavior.
    gq: int or None. If not None and gt is not None, we will add an this GQ
      value to our VariantCall.
    sample_name: str or None. If not None and gt is not None, sets the
      call_set_name of our VariantCall to this value.
    gls: array-list of float, or None. If not None and gt is not None, sets the
      genotype_likelihoods of our VariantCall to this value.

  Returns:
    third_party.nucleus.protos.Variant proto.
  """
  if alleles is None:
    alleles = ['A', 'C']

  if not end:
    end = start + len(alleles[0])

  variant = variants_pb2.Variant(
      reference_name=chrom,
      start=start,
      end=end,
      reference_bases=alleles[0],
      alternate_bases=alleles[1:],
      quality=qual,
  )

  if filters is not None:
    if not isinstance(filters, (list, tuple)):
      filters = [filters]
    variant.filter[:] = filters

  if gt:
    call = variant.calls.add(genotype=gt)

    if sample_name:
      call.call_set_name = sample_name

    if gq:
      set_list_values(call.info['GQ'], [gq])

    if gls:
      call.genotype_likelihood.extend(gls)

  return variant


def make_read(bases,
              quals=None,
              cigar=None,
              mapq=50,
              chrom='chr1',
              start=1,
              name='read'):
  """Makes a third_party.nucleus.protos.Read for testing."""
  if quals and len(bases) != len(quals):
    raise ValueError('Incompatable bases and quals', bases, quals)
  read = reads_pb2.Read(
      fragment_name=name,
      proper_placement=True,
      read_number=1,
      number_reads=2,
      aligned_sequence=bases,
      aligned_quality=quals,
      alignment=reads_pb2.LinearAlignment(
          position=position_pb2.Position(reference_name=chrom, position=start),
          mapping_quality=mapq,
          cigar=_cigar.to_cigar_units(cigar) if cigar else []))
  return read


def cc_iterable_len(cc_iterable):
  """Count the number of elements in an Iterable object.

  Args:
    cc_iterable: a CLIF-wrap of a subclass of the C++ Iterable class.

  Returns:
    integer count
  """
  count = 0
  while True:
    not_done, _ = cc_iterable.Next()
    if not not_done:
      break
    count += 1
  return count


def iterable_len(iterable):
  """Returns the length of a Python iterable, by advancing it."""
  return sum(1 for _ in iterable)


# redacted
def assert_not_called_workaround(mock):
  """Asserts that a mock has not been called.

  There's a bug in mock where some of the assert functions on a mock are being
  dropped when that mock is created with an autospec:

    https://bugs.python.org/issue28380

  The mock 2.0.0 backport doesn't have the fix yet. The required patch is:

    https://bugs.python.org/file44991/fix_autospecced_mock_functions.patch

  but the current mock (checked 07/22/17) backport code is missing the fix:

    https://github.com/testing-cabal/mock/blob/master/mock/mock.py#L315

  This is an open issue on the mock github repo:

    https://github.com/testing-cabal/mock/issues/398

  And they claim that it'll be a few months (as of April 2017) before it is
  incorporated into the backport.

  Args:
    mock: The mock to assert hasn't been called.

  Raises:
    AssertionError if mock has been called.
  """
  if mock.call_count != 0:
    raise AssertionError("Expected no calls to '{}' but was called {} times"
                         .format(mock.name, mock.call_count))


# redacted
def assert_called_once_workaround(mock):
  """Asserts that a mock has been called exactly once.

  See assert_not_called_workaround for the backstory on why this function
  exists.

  Args:
    mock: The mock that should have been called exactly once.

  Raises:
    AssertionError if mock wasn't called exactly once.
  """
  if mock.call_count != 1:
    raise AssertionError(
        "Expected exactly one call to '{}' but was called {} times".format(
            mock.name, mock.call_count))
