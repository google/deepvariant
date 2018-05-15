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
"""Classes for reading FASTA files.

The FASTA format is described at
https://en.wikipedia.org/wiki/FASTA_format

API for reading:

```python
from third_party.nucleus.io import fasta
from third_party.nucleus.protos import range_pb2

with fasta.RefFastaReader(input_path) as reader:
  region = range_pb2.Range(reference_name='chrM', start=1, end=6)
  basepair_string = reader.query(region)
  print(basepair_string)
```

If `input_path` ends with '.gz', it is assumed to be compressed.  All FASTA
files are assumed to be indexed with the index file located at
`input_path + '.fai'`.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections

import six

from third_party.nucleus.io import genomics_reader
from third_party.nucleus.io.python import fasta_reader
from third_party.nucleus.io.python import reference_fai
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.util import ranges

# redacted
RefFastaHeader = collections.namedtuple(
    'RefFastaHeader', ['contigs'])


class RefFastaReader(genomics_reader.GenomicsReader):
  """Class for reading from FASTA files containing a reference genome."""

  def __init__(self, input_path, cache_size=None):
    """Initializes a RefFastaReader.

    Args:
      input_path: string. A path to a resource containing FASTA records.
      cache_size: integer. Number of bases to cache from previous queries.
        Defaults to 64K.  The cache can be disabled using cache_size=0.
    """
    super(RefFastaReader, self).__init__()

    fasta_path = input_path
    fai_path = fasta_path + '.fai'
    if cache_size is None:
      # Use the C++-defined default cache size.
      self._reader = reference_fai.GenomeReferenceFai.from_file(
          fasta_path, fai_path)
    else:
      self._reader = reference_fai.GenomeReferenceFai.from_file(
          fasta_path, fai_path, cache_size)

    # redacted
    self.header = RefFastaHeader(contigs=self._reader.contigs)

  def iterate(self):
    """Returns an iterator for going through all bases in the file.

    NOTE: This function is not implemented for this data type. Retrieving all
    base pairs in the genome can be performed by querying for the full range of
    each contig defined in the header.
    """
    raise NotImplementedError('Can not iterate through a FASTA file')

  def query(self, region):
    """Returns the base pairs (as a string) in the given region."""
    return self._reader.bases(region)

  def is_valid(self, region):
    """Returns whether the region is contained in this FASTA file."""
    return self._reader.is_valid_interval(region)

  def contig(self, contig_name):
    """Returns a ContigInfo proto for contig_name."""
    return self._reader.contig(contig_name)

  @property
  def c_reader(self):
    """Returns the underlying C++ reader."""
    return self._reader

  def __exit__(self, exit_type, exit_value, exit_traceback):
    self._reader.__exit__(exit_type, exit_value, exit_traceback)


class InMemoryRefReader(genomics_reader.GenomicsReader):
  """A `RefFastaReader` getting its bases from an in-memory data structure.

  An `InMemoryRefReader` provides the same API as `RefFastaReader` but doesn't
  fetch its data from an on-disk FASTA file but rather fetches the bases from an
  in-memory cache containing (chromosome, start, bases) tuples.

  In particular, the `query(Range(chrom, start, end))` operation fetches bases
  from the tuple where `chrom` == chromosome, and then from the bases where the
  first base of bases starts at start. If start > 0, then the bases string is
  assumed to contain bases starting from that position in the region. For
  example, the record ('1', 10, 'ACGT') implies that
  `query(ranges.make_range('1', 11, 12))` will return the base 'C', as the 'A'
  base is at position 10. This makes it straightforward to cache a small region
  of a full chromosome without having to store the entire chromosome sequence in
  memory (potentially big!).
  """

  def __init__(self, chromosomes):
    """Initializes an InMemoryRefReader using data from chromosomes.

    Args:
      chromosomes: list[tuple]. The chromosomes we are caching in memory as a
        list of tuples. Each tuple must be exactly three elements in length,
        containing (chromosome name [str], start [int], bases [str]).

    Raises:
      ValueError: If any of the chromosomes tuples are invalid.
    """
    super(InMemoryRefReader, self).__init__()

    ref_seqs = []
    contigs = []
    for i, (contig_name, start, bases) in enumerate(chromosomes):
      if start < 0:
        raise ValueError('start={} must be >= for chromosome={}'.format(
            start, contig_name))
      if not bases:
        raise ValueError(
            'Bases must contain at least one base, but got "{}"'.format(bases))

      end = start + len(bases)
      ref_seqs.append(reference_pb2.ReferenceSequence(
          region=ranges.make_range(contig_name, start, end), bases=bases))
      contigs.append(
          reference_pb2.ContigInfo(
              name=contig_name, n_bases=end, pos_in_fasta=i))

    self._reader = fasta_reader.InMemoryGenomeReference.create(
        contigs, ref_seqs)
    self.header = RefFastaHeader(contigs=self._reader.contigs)

  def iterate(self):
    """Returns an iterator for going through all bases in the file.

    NOTE: This function is not implemented for this data type. Retrieving all
    base pairs in the genome can be performed by querying for the full range of
    each contig defined in the header.
    """
    raise NotImplementedError('Can not iterate through a FASTA file')

  def query(self, region):
    """Returns the base pairs (as a string) in the given region."""
    return self._reader.bases(region)

  def is_valid(self, region):
    """Returns whether the region is contained in this FASTA file."""
    return self._reader.is_valid_interval(region)

  def contig(self, contig_name):
    """Returns a ContigInfo proto for contig_name."""
    return self._reader.contig(contig_name)

  @property
  def c_reader(self):
    """Returns the underlying C++ reader."""
    return self._reader

  def __str__(self):

    def _format_refseq(refseq):
      bases = refseq.bases
      if len(bases) >= 50:
        bases = bases[0:50] + '...'
      return 'Contig(chrom={} start={}, end={}, bases={})'.format(
          refseq.region.reference_name, refseq.region.start, refseq.region.end,
          bases)

    contigs_strs = [
        _format_refseq(refseq)
        for refseq in six.itervalues(self._reader.reference_sequences)
    ]
    return 'InMemoryRefReader(contigs={})'.format(''.join(contigs_strs))

  __repr__ = __str__
