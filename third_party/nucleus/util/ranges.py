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

"""Utilities for Range overlap detection."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import gzip
import re



import intervaltree

from third_party.nucleus.protos import position_pb2
from third_party.nucleus.protos import range_pb2
from tensorflow.python.platform import gfile


# Regular expressions for matching literal chr:start-stop strings.
_REGION_LITERAL_REGEXP = re.compile(r'^(\S+):([0-9,]+)-([0-9,]+)$')
# Regular expressions for matching literal chr:start strings.
_POSITION_LITERAL_REGEXP = re.compile(r'^(\S+):([0-9,]+)$')


class RangeSet(object):
  """Fast overlap detection of a genomic position against a db of Ranges.

  Enables very fast O(log n) computation of whether a point chr:pos falls within
  one of a large number of genomic ranges.

  This class does not supports overlapping or adjacent intervals. Any such
  intervals will be automatically merged together in the constructor.

  This class is immutable. No methods should be added that directly modify the
  ranges held by the class.
  """

  # redacted
  def __init__(self, ranges=None):
    """Creates an RangeSet backed by ranges.

    Note that the Range objects in ranges are *not* stored directly here, so
    they can safely be modified after they are passed to this RangeSet.

    Args:
      ranges: A list of genomics.Range objects (or anything with
        reference_name, start, and end properties following the
        genomics.Range convention). If None, no ranges will be used,
        and overlaps() will always return False.
    """
    if ranges is None:
      ranges = []

    # Add each range to our contig-specific intervaltrees.
    self._by_chr = collections.defaultdict(intervaltree.IntervalTree)
    for range_ in ranges:
      self._by_chr[range_.reference_name].addi(range_.start, range_.end, None)

    # Merge overlapping / adjacent intervals in each tree.
    for tree in self._by_chr.itervalues():
      tree.merge_overlaps()

  def __iter__(self):
    """Iterate over the ranges in this RangeSet.

    Yields:
      Each range of this RangeSet, in an arbitrary order. These objects are new
      range protos so can be freely modified.
    """
    for refname, chr_ranges in self._by_chr.iteritems():
      for start, end, _ in chr_ranges:
        yield make_range(refname, start, end)

  @classmethod
  def from_regions(cls, regions, contig_map=None):
    """Parses a command-line style literal regions flag into a RangeSet.

    Args:
      regions: An iterable or None. If not None, regions will be parsed with
        ranges.from_regions.
      contig_map: An optional dictionary mapping from contig names to ContigInfo
        protobufs. If provided, allows literals of the format "contig_name",
        which will be parsed into a Range with reference_name=contig_name,
        start=0, end=n_bases where n_bases comes from the ContigInfo.

    Returns:
      A RangeSet object.
    """
    if regions is None:
      return cls(ranges=[])
    else:
      return cls(ranges=from_regions(regions, contig_map=contig_map))

  @classmethod
  def from_contigs(cls, contigs):
    """Creates a RangeSet with an interval covering each base of each contig."""
    return cls(make_range(contig.name, 0, contig.n_bases) for contig in contigs)

  @classmethod
  def from_bed(cls, source):
    """Creates a RangeSet containing the intervals from source.

    Args:
      source: A path to a BED (or equivalent) file of intervals.

    Returns:
      A RangeSet.
    """
    with gfile.GFile(source) as fin:
      if source.endswith('.gz'):
        fin = gzip.GzipFile(fileobj=fin)
      return cls(bed_parser(fin))

  def intersection(self, *others):
    """Computes the intersection among this RangeSet and *others RangeSets.

    This function computes the intersection of all of the intervals in self and
    *others, returning a RangeSet containing only intervals common to all. The
    intersection here is an ranged intersection, not an identity intersection,
    so the resulting set of intervals may not contain any of the original
    intervals in any of the sets.

    To be concrete, suppose we have three sets to intersect, each having two
    intervals:

      self   : chr1:1-10, chr2:20-30
      other1 : chr1:5-8, chr3:10-40
      other2 : chr1:3-7, chr3:10-30

    self.intersection(other1, other2) produces a RangeSet with one interval
    chr1:3-7, the common bases on chr1 in self, other1, and other2. No intervals
    on chr2 or chr3 are included since the chr2 only occurs in self and the two
    intervals on chr3, despite having some shared bases, don't have an
    overlapping interval in self.

    Args:
      *others: A list of RangeSet objects to intersect with the intervals in
        this RangeSet.

    Returns:
      A RangeSet. If *others is empty, this function returns self rather than
      making an unnecessary copy. In all other cases, the returned value will be
      a freshly allocated RangeSet.
    """

    def _intersect2(refname, tree1, tree2):
      """Intersects the intervals of two IntervalTrees."""
      # Yields all of the overlapping intervals from each interval of tree1
      # found in tree2. Since each tree has only non-adjacent, non-overlapping,
      # intervals this calculation is straightforward and safe and produces only
      # non-adjacent, non-overlapping intervals.
      return (make_range(refname, max(interval1.begin, overlapping.begin),
                         min(interval1.end, overlapping.end))
              for interval1 in tree1
              for overlapping in tree2[interval1])

    # Iteratively intersect each of our *other RangeSets with this RangeSet.
    # Sort by size so we do the smallest number of element merge first.
    # redacted
    # common contigs upfront across all others and only looping over those.
    intersected = self
    for other in sorted(others, key=len):
      intersected_intervals = []
      # pylint: disable=protected-access
      # So we can intersect intervals within each contig separately.
      for refname, intervals in intersected._by_chr.iteritems():
        # If refname is present in other, intersect those two IntervalTrees
        # directly and add those contigs to our growing list of intersected
        # intervals. If refname isn't present, all of the intervals on refname
        # should be dropped as there are no intervals to overlap.
        other_chr = other._by_chr.get(refname, None)
        if other_chr:
          intersected_intervals.extend(
              _intersect2(refname, intervals, other_chr))

      # Update our intersected RangeSet with the new intervals.
      intersected = RangeSet(intersected_intervals)

    return intersected

  def exclude_regions(self, other):
    """Chops out all of the intervals in other from this this RangeSet.

    This is a *MUTATING* operation for performance reasons. Make a copy of self
    if you want to avoid modifying the RangeSet.

    Args:
      other: A RangeSet object whose intervals will be removed from this
        RangeSet.
    """
    # pylint: disable=protected-access
    for chrname, chr_intervals in other._by_chr.iteritems():
      # If refname is present in self, difference those two IntervalTrees.
      self_intervals = self._by_chr.get(chrname, None)
      if self_intervals:
        for begin, end, _ in chr_intervals:
          self_intervals.chop(begin, end)
        if self_intervals.is_empty():
          # Cleanup after ourselves by removing empty trees from out map.
          del self._by_chr[chrname]

  def __len__(self):
    """Gets the number of ranges used by this RangeSet."""
    return sum(len(for_chr) for for_chr in self._by_chr.itervalues())

  def __nonzero__(self):
    """Returns True if this RangeSet is not empty."""
    return bool(self._by_chr)

  __bool__ = __nonzero__  # Python 3 compatibility.

  def variant_overlaps(self, variant, empty_set_return_value=True):
    """Returns True if the variant's range overlaps with any in this set."""
    if not self:
      return empty_set_return_value
    else:
      return self.overlaps(variant.reference_name, variant.start)

  def overlaps(self, chrom, pos):
    """Returns True if chr:pos overlaps with any range in this RangeSet.

    Uses a fast bisection algorithm to determine the overlap in O(log n) time.

    Args:
      chrom: The chromosome name as a string.
      pos: The position (0-based) as an integer.

    Returns:
      True if chr:pos overlaps with a range.
    """
    chr_ranges = self._by_chr.get(chrom, None)
    if chr_ranges is None:
      return False
    return chr_ranges.overlaps(pos)

  def partition(self, max_size):
    """Splits our intervals so that none are larger than max_size.

    Slices up the intervals in this RangeSet into a equivalent set of interval (
    i.e., spanning the same set of bases), each of which is at most max_size in
    length.

    This function does not modify this RangeSet.

    Because RangeSet merges adjacent intervals, this function cannot return use
    a RangeSet to represent the partitioned intervals and so instead generates
    these intervals via a yield statement.

    Args:
      max_size: A positive integer (> 0) indicating the maximum size of any
        interval.

    Yields:
      third_party.nucleus.protos.Range protos
      in an arbitrary (but not necessarily random) order.

    Raises:
      ValueError: if max_size <= 0.
    """
    if max_size <= 0:
      raise ValueError('max_size must be > 0', max_size)

    for interval in self:
      refname = interval.reference_name
      for pos in range(interval.start, interval.end, max_size):
        yield make_range(refname, pos, min(interval.end, pos + max_size))


def make_position(chrom, position, reverse_strand=False):
  """Makes a third_party.nucleus.protos.Position."""
  return position_pb2.Position(
      reference_name=chrom, position=position, reverse_strand=reverse_strand)


def make_range(chrom, start, end):
  """Creates a genomics.Range object chr:start-end.

  Args:
    chrom: The chromosome name as a string.
    start: The start position (0-based, inclusive, integer) of this range.
    end: The end position (0-based, exclusive, integer) of this range.

  Returns:
    A third_party.nucleus.protos.Range.
  """
  return range_pb2.Range(reference_name=chrom, start=start, end=end)


def position_overlaps(chrom, pos, interval):
  """Does interval overlap the position chr:pos?

  Args:
    chrom: The chromosome name as a string.
    pos: The position (0-based, integer).
    interval: third_party.nucleus.protos.Range object.

  Returns:
    True if interval overlaps chr:pos.
  """
  return (chrom == interval.reference_name and
          interval.start <= pos < interval.end)


def ranges_overlap(i1, i2):
  """Checks whether ranges i1 and i2 overlap.

  Args:
    i1: third_party.nucleus.protos.Range object.
    i2: third_party.nucleus.protos.Range object.

  Returns:
    True if i1 and i2 overlap.
  """
  return (i1.reference_name == i2.reference_name and i1.end > i2.start and
          i1.start < i2.end)


def bedpe_parser(fd):
  """Parses Range objects from a BEDPE-formatted file object.

  See http://bedtools.readthedocs.org/en/latest/content/general-usage.html
  for more information on the BEDPE format.

  Skips events that span across chromosomes. For example, if the starting
  location is on chr1 and the ending location is on chr2, that record will
  not appear in the output.

  Args:
    fd: An iterable producing string, one per line in BEDPE format.

  Yields:
    third_party.nucleus.protos.Range protobuf objects.
  """
  for line in fd:
    parts = line.split('\t')
    if parts[0] == parts[3]:
      # only keep events on the same chromosome
      yield make_range(parts[0], int(parts[1]), int(parts[5]))


def bed_parser(fd):
  """Parses Range objects from a BED-formatted file object.

  See http://bedtools.readthedocs.org/en/latest/content/general-usage.html
  for more information on the BED format.

  Args:
    fd: An iterable producing string, one per line in BED format.

  Yields:
    third_party.nucleus.protos.Range protobuf objects.
  """
  for line in fd:
    parts = line.split('\t')
    yield make_range(parts[0], int(parts[1]), int(parts[2]))


def from_regions(regions, contig_map=None):
  """Parses each region of `regions` into a Range proto.

  This function provides a super high-level interface for
  reading/parsing/converting objects into Range protos. Each `region` of
  `regions` is processed in turn, yielding one or more Range protos. This
  function inspects the contents of `region` to determine how to convert it to
  Range(s) protos. The following types of `region` strings are supported:

    * If region ends with an extension known in _get_parser_for_file, we treat
      region as a file and read the Range protos from it with the corresponding
      reader from _get_parser_for_file, yielding each Range from the file in
      order.
    * Otherwise we parse region as a region literal (`chr20:1-10`) and return
      the Range proto.

  Args:
    regions: iterable[str]. Converts each element of this iterable into
      region(s).
    contig_map: An optional dictionary mapping from contig names to ContigInfo
      protobufs. If provided, allows literals of the format "contig_name",
      which will be parsed into a Range with reference_name=contig_name,
      start=0, end=n_bases where n_bases comes from the ContigInfo.

  Yields:
    A Range proto.
  """
  for region in regions:
    reader = _get_parser_for_file(region)
    if reader:
      with gfile.GFile(region) as fin:
        for elt in reader(fin):
          yield elt
    else:
      yield parse_literal(region, contig_map)


# Cannot be at the top of the file because these parser functions need to be
# defined before adding them to the dictionary.
_REGION_FILE_READERS = {
    bed_parser: frozenset(['.bed']),
    bedpe_parser: frozenset(['.bedpe']),
}


def _get_parser_for_file(filename):
  for reader, exts in _REGION_FILE_READERS.iteritems():
    if any(filename.lower().endswith(ext) for ext in exts):
      return reader
  return None


def to_literal(range_pb):
  """Converts Range protobuf into string literal form.

  The string literal form looks like:

    reference_name:start+1-end

  since start and end are zero-based inclusive (start) and exclusive (end),
  while the literal form is one-based inclusive on both ends.

  Args:
    range_pb: A third_party.nucleus.protos.Range object.

  Returns:
    A string.
  """
  return '{}:{}-{}'.format(range_pb.reference_name, range_pb.start + 1,
                           range_pb.end)


def parse_literal(region_literal, contig_map=None):
  """Parses a Range from a string representation like chr:start-end.

  The region literal must conform to the following pattern:

    chromosome:start-end
    chromosome:position
    chromosome  [if contig_map is provided]

  chromosome can be any non-empty string without whitespace. start and end must
  both be positive integers. They can contain commas for readibility. start and
  end are positions not offsets, so start == 1 means an offset of zero. If only
  a single position is provided, this creates a 1 bp interval starting at
  position - 1 and ending at position.

  Inspired by the samtools region specification:
  http://www.htslib.org/doc/samtools.html

  Args:
    region_literal: string literal to parse.
    contig_map: An optional dictionary mapping from contig names to ContigInfo
      protobufs. If provided, allows literals of the format "contig_name", which
      will be parsed into a Range with reference_name=contig_name, start=0,
      end=n_bases where n_bases comes from the ContigInfo.

  Returns:
    third_party.nucleus.protos.Range.

  Raises:
    ValueError: if region_literal cannot be parsed.
  """

  def parse_position(pos_str):
    return int(pos_str.replace(',', ''))

  matched = _REGION_LITERAL_REGEXP.match(region_literal)
  if matched:
    chrom, start, end = matched.groups()
    return make_range(chrom, parse_position(start) - 1, parse_position(end))

  matched = _POSITION_LITERAL_REGEXP.match(region_literal)
  if matched:
    chrom, pos = matched.groups()
    pos = parse_position(pos)
    return make_range(chrom, pos - 1, pos)

  if contig_map and region_literal in contig_map:
    # If the region_literals is an exact contig name like chr1 or MT return a
    # range over the entire contig.
    return make_range(region_literal, 0, contig_map[region_literal].n_bases)
  else:
    raise ValueError('Could not parse: ', region_literal)


def parse_literals(region_literals, contig_map=None):
  """Parses each literal of region_literals in order."""
  return [parse_literal(literal, contig_map) for literal in region_literals]


def parse_lines(lines, file_format=None):
  """Creates a generator of Range objects from source lines in file_format.

  This is a genetic function for reading Range objects from a source. A parser
  will be created suitable for file_format, and each line of lines will be
  parsed.

  Args:
    lines: An iterable yielding single lines in file_format format, in order.
    file_format: A string specification of the format of the lines. Currently
      can be one of BED and BEDPE.

  Returns:
    A generator of Range objects.

  Raises:
    ValueError: If no suitable parser for file_format could be found.
  """
  readers = {
      'bedpe': bedpe_parser,
      'bed': bed_parser,
  }
  parser = readers.get(file_format.lower(), None)
  if not parser:
    raise ValueError('Unsupported file format', file_format)
  else:
    return parser(lines)


def contigs_n_bases(contigs):
  """Returns the sum of all n_bases of contigs."""
  return sum(c.n_bases for c in contigs)


def contigs_dict(contigs):
  """Creates a dictionary for contigs.

  Args:
    contigs: Iterable of ContigInfo protos.

  Returns:
    A dictionary mapping contig.name: contig for each contig in contigs.
  """
  return {contig.name: contig for contig in contigs}


def sorted_ranges(ranges, contigs=None):
  """Sorts ranges by reference_name, start, and end.

  Args:
    ranges: A sequence of google.v1.genomics.Range protos that we want to sort.
    contigs: None or an iterable of ConfigInfo protos. If not None, we will use
      the order of the contigs (as defined by their pos_in_fasta field values)
      to sort the Ranges on different contigs with respect to each other.

  Returns:
    A newly allocated list of google.v1.genomics.Range protos.
  """
  if contigs:
    contig_map = contigs_dict(contigs)

    def to_key(range_):
      pos = contig_map[range_.reference_name].pos_in_fasta
      return pos, range_.start, range_.end
  else:
    to_key = as_tuple

  return sorted(ranges, key=to_key)


def as_tuple(range_):
  """Returns a Python tuple (reference_name, start, end)."""
  return range_.reference_name, range_.start, range_.end


def overlap_len(range1, range2):
  """Computes the number of overlapping bases of range1 and range2.

  Args:
    range1: learning.genomics.genomics.Range.
    range2: learning.genomics.genomics.Range.

  Returns:
    int. The number of basepairs in common. 0 if the ranges are not on the same
    contig.
  """
  if range1.reference_name != range2.reference_name:
    return 0
  return max(0, (min(range1.end, range2.end) - max(range1.start, range2.start)))


def find_max_overlapping(query_range, search_ranges):
  """Gets the index of the element in search_ranges with max overlap with query.

  In case of ties, selects the lowest index range in search_ranges.

  Args:
    query_range: learning.genomics.genomics.core.Range, read genomic range.
    search_ranges: list[learning.genomics.genomics.core.Read]. Cannot be an
      iterable as we loop over the search_ranges multiple times. The list of
      regions we want to search for the maximal overlap with query_range.

  Returns:
    int, the search_ranges index with the maximum read overlap. Returns None
    when read has no overlap with any of the search_ranges or search_ranges is
    empty.
  """
  if not search_ranges:
    return None
  overlaps = [overlap_len(query_range, srange) for srange in search_ranges]
  argmax = max(range(len(search_ranges)), key=lambda i: overlaps[i])
  # We return None if the read doesn't overlap at all.
  return None if overlaps[argmax] == 0 else argmax


def expand(region, n_bp, contig_map=None):
  """Expands region by n_bp in both directions.

  Takes a Range(chrom, start, stop) and returns a new
  Range(chrom, new_start, new_stop), where:

  -- new_start is max(start - n_bp, 0)
  -- new_stop is stop + n_bp if contig_map is None, or min(stop + n_bp, max_bp)
     where max_bp is contig_map[chrom].n_bp.

  Args:
    region: A third_party.nucleus.protos.Range proto.
    n_bp: int >= 0; how many basepairs to increase region by.
    contig_map: None, or dict[string, ContigInfo]. If not None, used to get the
      maximum extent to increase stop by. Must have region.reference_name as a
      key.

  Returns:
    third_party.nucleus.protos.Range proto.

  Raises:
    ValueError: if n_bp is invalid.
    KeyError: contig_map is not None and region.reference_name isn't a key.
  """
  if n_bp < 0:
    raise ValueError('n_bp must be >= 0 but got {}'.format(n_bp))

  new_start = max(region.start - n_bp, 0)
  new_end = region.end + n_bp
  if contig_map is not None:
    new_end = min(new_end, contig_map[region.reference_name].n_bases)
  return make_range(region.reference_name, new_start, new_end)


def span(regions):
  """Returns a region that spans all of the bases in regions.

  This function returns a Range(chrom, start, stop), where start is the min
  of the starts in regions, and stop is the max end in regions. It may not be
  freshly allocated.

  Args:
    regions: list[Range]: an list of Range protos.

  Returns:
    Range proto.

  Raises:
    ValueError: if not all regions have the same reference_name.
    ValueError: if regions is empty.
  """
  if not regions:
    raise ValueError('regions is empty but must have at least one region')
  elif len(regions) == 1:
    return regions[0]
  elif any(r.reference_name != regions[0].reference_name for r in regions):
    raise ValueError('regions must be all on the same contig')
  else:
    start = min(r.start for r in regions)
    end = max(r.end for r in regions)
    return make_range(regions[0].reference_name, start, end)
