# Copyright 2020 Google LLC.
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
"""Beam sources for genomics file formats."""

from third_party.nucleus.io import bed
from third_party.nucleus.io import sam

# pylint:disable=g-bad-import-order
# Beam needs to be imported after nucleus for this to work in open source.
from apache_beam.io import filebasedsource
from apache_beam.io.iobase import Read
from apache_beam.transforms import PTransform
# pylint:enable=g-bad-import-order


class _GenomicsSource(filebasedsource.FileBasedSource):
  """A base source for reading genomics files.

  Do not use this class directly. Instead, use the subclass for the specific
  file type.
  """

  def __init__(self, file_pattern, validate, **nucleus_kwargs):
    """Initializes a _GenomicsSource for use with readers for genomics files."""

    super(_GenomicsSource, self).__init__(
        file_pattern=file_pattern, splittable=False, validate=validate)
    self.nucleus_kwargs = nucleus_kwargs

  def read_records(self, input_path, offset_range_tracker):
    """Yields records returned by nucleus_reader."""
    if offset_range_tracker.start_position():
      raise ValueError('Start position not 0: %d' %
                       offset_range_tracker.start_position())
    current_offset = offset_range_tracker.start_position()
    reader = self.nucleus_reader(input_path, **self.nucleus_kwargs)

    with reader:
      for record in reader:
        if not offset_range_tracker.try_claim(current_offset):
          raise RuntimeError('Unable to claim position: %d' % current_offset)
        yield record
        current_offset += 1

  @property
  def nucleus_reader(self):
    raise NotImplementedError


class _BedSource(_GenomicsSource):
  """A source for reading BED files."""

  nucleus_reader = bed.BedReader


class _SamSource(_GenomicsSource):
  """A source for reading SAM/BAM files."""

  nucleus_reader = sam.SamReader


class ReadGenomicsFile(PTransform):
  """For reading one or more genomics files.

  Do not use this class directly. Instead, use the subclass for the specific
  file type.
  """

  def __init__(self, file_pattern, validate=True, **nucleus_kwargs):
    """Initializes the ReadSam transform."""

    super(ReadGenomicsFile, self).__init__()
    self._source = self._source_class(
        file_pattern, validate=validate, **nucleus_kwargs)

  def expand(self, pvalue):
    return pvalue.pipeline | Read(self._source)

  @property
  def _source_class(self):
    raise NotImplementedError


class ReadBed(ReadGenomicsFile):
  """For reading records from one or more BED files."""

  _source_class = _BedSource


class ReadSam(ReadGenomicsFile):
  """For reading records from one or more SAM/BAM files."""

  _source_class = _SamSource
