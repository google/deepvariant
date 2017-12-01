# Copyright 2017 Google Inc.
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
"""Variant utilities."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



import enum

from deepvariant.core import ranges
from deepvariant.core.genomics import struct_pb2
from deepvariant.core.genomics import variants_pb2

# The alternate allele string for reference (no alt).
NO_ALT_ALLELE = '.'
# The alternate allele string for the gVCF "any" alternate allele.
GVCF_ALT_ALLELE = '<*>'


def set_variantcall_gq(variant_call, gq):
  if 'GQ' in variant_call.info:
    del variant_call.info['GQ']
  variant_call.info['GQ'].values.extend([struct_pb2.Value(number_value=gq)])


def decode_variants(encoded_iter):
  """Yields a genomics.Variant from encoded_iter.

  Args:
    encoded_iter: An iterable that produces binary encoded
      learning.genomics.deepvariant.core.genomics.Variant strings.

  Yields:
    A parsed learning.genomics.deepvariant.core.genomics.Variant for each
    encoded element of encoded_iter
    in order.
  """
  for encoded in encoded_iter:
    yield variants_pb2.Variant.FromString(encoded)


def variant_position(variant):
  """Returns a new Range at the start position of variant.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A new Range with the same reference_name as variant and start but an end
    that is start + 1. This produces a range that is the single basepair of the
    start of variant, hence the name position.
  """
  return ranges.make_range(variant.reference_name, variant.start,
                           variant.start + 1)


def variant_range(variant):
  """Returns a new Range covering variant.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A new Range with the same reference_name, start, and end as variant.
  """
  return ranges.make_range(variant.reference_name, variant.start, variant.end)


@enum.unique
class GenotypeType(enum.Enum):
  """An enumeration of the types of genotypes."""
  hom_ref = ('homozygous reference', [0, 0], 0)
  het = ('heterozygous', [0, 1], 1)
  hom_var = ('homozygous non-reference', [1, 1], 2)
  no_call = ('no call', [-1, -1], -1)

  def __init__(self, full_name, example_gt, class_id):
    self.full_name = full_name
    self.example_gt = example_gt
    self.class_id = class_id


@enum.unique
class VariantType(enum.Enum):
  """An enumeration of the types of variants."""
  # a variant.proto where there is no alt allele
  ref = 0
  # a non-reference variant.proto where all ref and alt alleles
  # are single basepairs
  snp = 1
  # a non-reference variant.proto where at least one of ref or alt alleles
  # are longer than 1 bp
  indel = 2


def format_filters(variant):
  """Gets a human-readable string showing the filters applied to variant.

  Returns a string with the filter field values of variant separated by commas.
  If the filter field isn't set, returns '.'.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A string.
  """
  return ','.join(variant.filter) if variant.filter else '.'


def format_alleles(variant):
  """Gets a string representation of the variants alleles.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A string ref_bases/alt1,alt2 etc.
  """
  return '{}/{}'.format(variant.reference_bases, ','.join(
      variant.alternate_bases))


def format_position(variant):
  """Gets a string representation of the variants position.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A string chr:start + 1 (as start is zero-based).
  """
  return '{}:{}'.format(variant.reference_name, variant.start + 1)


def is_snp(variant):
  """Is variant a SNP?

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    True if all alleles of variant are 1 bp in length.
  """
  return (not is_ref(variant) and len(variant.reference_bases) == 1 and
          len(variant.alternate_bases) >= 1 and
          all(len(x) == 1 for x in variant.alternate_bases))


def is_indel(variant):
  """Is variant an indel?

  An indel event is simply one where the size of at least one of the alleles
  is > 1.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    True if the alleles in variant indicate an insertion/deletion event
    occurs at this site.
  """
  # redacted
  # redacted
  return (not is_ref(variant) and
          (len(variant.reference_bases) > 1 or
           any(len(alt) > 1 for alt in variant.alternate_bases)))


def is_biallelic(variant):
  """Returns True if variant has exactly one alternate allele."""
  return len(variant.alternate_bases) == 1


def is_multiallelic(variant):
  """Does variant have multiple alt alleles?

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    True if variant has more than one alt allele.
  """
  return len(variant.alternate_bases) > 1


def is_ref(variant):
  """Returns true if variant is a reference record.

  Variant protos can encode sites that aren't actually mutations in the
  sample.  For example, the record ref='A', alt='.' indicates that there is
  no mutation present (i.e., alt is the missing value).

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A boolean.
  """
  alts = variant.alternate_bases
  return not alts or (len(alts) == 1 and alts[0] == '.')


def variant_type(variant):
  """Gets the VariantType of variant.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    VariantType indicating the type of this variant.
  """
  if is_ref(variant):
    return VariantType.ref
  elif is_snp(variant):
    return VariantType.snp
  else:
    return VariantType.indel


def is_transition(allele1, allele2):
  """Is the pair of single bp alleles a transition?

  Args:
    allele1: A string of the first allele, must be 1 bp in length.
    allele2: A string of the second allele, must be 1 bp in length.

  Returns:
    True if allele1/allele2 are a transition SNP.

  Raises:
    ValueError: if allele1 and allele2 are equal or aren't 1 bp in length.
  """
  if allele1 == allele2:
    raise ValueError('Alleles must be unique:', allele1, allele2)
  if len(allele1) != 1:
    raise ValueError('Alleles must be 1 bp in length.', allele1)
  if len(allele2) != 1:
    raise ValueError('Alleles must be 1 bp in length.', allele2)

  alleles_set = {allele1, allele2}
  return any(alleles_set == x for x in [{'A', 'G'}, {'C', 'T'}])


def is_insertion(ref, alt):
  """Is alt an insertion w.r.t. ref?

  Args:
    ref: A string of the reference allele.
    alt: A string of the alternative allele.

  Returns:
    True if alt is an insertion w.r.t. ref.
  """
  return len(ref) < len(alt)


def is_deletion(ref, alt):
  """Is alt a deletion w.r.t. ref?

  Args:
    ref: A string of the reference allele.
    alt: A string of the alternative allele.

  Returns:
    True if alt is a deletion w.r.t. ref.
  """
  return len(ref) > len(alt)


def has_insertion(variant):
  """Does variant have an insertion?

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    True if the alleles in variant indicate an insertion event
    occurs at this site.
  """
  ref = variant.reference_bases
  return (is_indel(variant) and
          any(is_insertion(ref, alt) for alt in variant.alternate_bases))


def has_deletion(variant):
  """Does variant have a deletion?

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    True if the alleles in variant indicate an deletion event
    occurs at this site.
  """
  ref = variant.reference_bases
  return (is_indel(variant) and
          any(is_deletion(ref, alt) for alt in variant.alternate_bases))


@enum.unique
class AlleleMismatchType(enum.Enum):
  """An enumeration of the types of allele mismatches we detect."""
  # Duplicate alleles
  duplicate_eval_alleles = 1
  duplicate_true_alleles = 2
  # Truth has an allele that doesn't match any allele in eval.
  unmatched_true_alleles = 3
  # Eval has an allele that doesn't match any allele in truth.
  unmatched_eval_alleles = 4


def allele_mismatches(evalv, truev):
  """Determines the set of allele mismatch discordances between evalv and truev.

  Compares the alleles present in evalv and truev to determine if there are any
  disagreements between the set of called alleles in the two Variant protos. The
  type of differences basically boil down to:

  -- Are there duplicate alt alleles?
  -- Can we find a matching allele in the truev for each allele in evalv, and
    vice versa?

  Two alleles A and B match when they would produce the same sequence of bases
  in ref and alt haplotypes starting at the same position. So CA=>TA is the same
  as C=>T (position is the same, replacing A by A is a noop) but AC=>AT isn't
  the same as C=>T because the former event changes bases 1 bp further along in
  the reference genome than the C=>T allele.

  Args:
    evalv: A learning.genomics.deepvariant.core.genomics.Variant.
    truev: A learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A set of AlleleMismatchType values.
  """
  unmatched_eval_alleles = []
  # Use set removes duplicate alleles in truth and eval variants.
  allele_matches = {alt: [] for alt in set(truev.alternate_bases)}
  for eval_alt in set(evalv.alternate_bases):
    # Loop over each possible alt allele, adding eval_alt to each matching alt
    # allele.
    found_match = False
    for true_alt in allele_matches:
      if (simplify_alleles(evalv.reference_bases, eval_alt) == simplify_alleles(
          truev.reference_bases, true_alt)):
        # We are a match to true_alt, so record that fact in allele_matches
        allele_matches[true_alt].append(eval_alt)
        found_match = True
    if not found_match:
      # We never found a match for eval_alt.
      unmatched_eval_alleles.append(eval_alt)

  # At this point we've checked every alt against every eval allele, and are
  # ready to summarize the differences using our AlleleMismatchType enum.
  types = set()
  if len(set(evalv.alternate_bases)) != len(evalv.alternate_bases):
    types.add(AlleleMismatchType.duplicate_eval_alleles)
  if len(set(truev.alternate_bases)) != len(truev.alternate_bases):
    types.add(AlleleMismatchType.duplicate_true_alleles)
  if unmatched_eval_alleles:
    types.add(AlleleMismatchType.unmatched_eval_alleles)
  if any(len(match) != 1 for match in allele_matches.itervalues()):
    types.add(AlleleMismatchType.unmatched_true_alleles)
  return types


def simplify_alleles(*alleles):
  """Simplifies alleles by stripping off common postfix bases.

  For example, simplify("AC", "GC") would produce the tuple "A", "G" as the "C"
  base is a common postfix of both alleles. But simplify("AC", "GT") would
  produce "AC", "GT" as there is no common postfix.

  Note this function will never simplify any allele down to the empty string. So
  if alleles = ['CACA', 'CA'], the longest common postfix is 'CA' but we will
  not produce ['CA', ''] as this is an invalid Variant allele encoding. Instead
  we produce ['CAC', 'C'].

  Args:
    *alleles: A tuple of bases, each as a string, to simplify.

  Returns:
    A tuple, one for each allele in alleles in order, with any common postfix
    bases stripped off.
  """

  def all_the_same(items):
    first = next(items)
    return all(item == first for item in items)

  # Loop over the alleles to determine the length of the shared postfix. Start
  # at 1 so every allele, even after trimming the postfix, has at least len 1.
  # For example, alleles = ['ATT', 'TT'] reduces to ['AT', 'T'] not ['A', ''].
  shortest_allele_len = min(len(a) for a in alleles)
  common_postfix_len = 0
  for i in range(1, shortest_allele_len):
    if not all_the_same(a[-i] for a in alleles):
      break
    common_postfix_len = i

  if common_postfix_len:
    return tuple(a[0:-common_postfix_len] for a in alleles)
  else:
    # Fast path for the case where there's no shared postfix.
    return alleles


def is_filtered(variant):
  """Returns True if variant has a non-PASS filter field, or False otherwise."""
  return bool(variant.filter) and any(
      f not in {'PASS', '.'} for f in variant.filter)


def is_variant_call(variant,
                    require_non_ref_genotype=True,
                    no_calls_are_variant=False):
  """Is variant a non-reference call?

  A Variant proto doesn't always imply that there's a variant present in the
  genome. The call may not have alternate bases, may be filtered, may a have
  hom-ref genotype, etc. This function looks for all of those configurations
  and returns true iff the variant is asserting that a mutation is present
  in the same.

  Note that this code allows a variant without a calls field to be variant,
  but one with a genotype call must have a non-reference genotype to be
  considered variant (if require_non_ref_genotype is True, the default). If
  False, a variant that passes all fo the site-level requirements for being
  a variant_call will return a True value, regardless of the genotypes, which
  means that we'll consider a site with a sample with a hom-ref or no-call site
  a variant call.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.
    require_non_ref_genotype: Should we require a site with a genotype call to
      have a non-reference (het, hom-var) genotype for the site to be considered
      a variant call?
    no_calls_are_variant: If a site has genotypes, should we consider no_call
      genotypes as being variant or not?

  Returns:
    True if variant is really a mutation call.

  Raises:
    ValueError: If variant has more than one call (i.e., is multi-sample).
  """
  if not variant.alternate_bases:
    return False
  elif is_filtered(variant):
    return False
  elif not variant.calls or not require_non_ref_genotype:
    return True
  # All tests after this point should only look at genotype-based fields, as
  # we may have aborted out in the prev. line due to require_non_ref_genotype.
  elif len(variant.calls) > 1:
    raise ValueError('Unsupported: multiple genotypes found at', variant)
  elif any(g > 0 for g in variant.calls[0].genotype):
    return True
  elif no_calls_are_variant:
    return all(g == -1 for g in variant.calls[0].genotype)
  else:
    return False


def has_genotypes(variant):
  """Does variant have genotype calls?

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    True if variant has genotype calls.
  """
  # I don't want to return the actual data structure so I'm doing the
  # explicit True/False evaluation here.
  # pylint: disable=g-explicit-length-test
  return len(variant.calls) > 0


def genotype_type(variant):
  """Gets the GenotypeType for variant.

  If variant doesn't have genotypes, returns no_call. Otherwise
  returns one of no_call, hom_ref, het, or hom_var depending on the
  status of the genotypes in the call field of variant.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A GenotypeType.

  Raises:
    ValueError: If variant has more than one call (i.e., is multi-sample).
  """
  if not has_genotypes(variant):
    return GenotypeType.no_call
  elif len(variant.calls) > 1:
    raise ValueError('Unsupported: multiple genotypes found at', variant)
  else:
    gt = set(variant.calls[0].genotype)
    if gt == {-1}:
      return GenotypeType.no_call
    elif gt == {0}:
      return GenotypeType.hom_ref
    elif len(gt) > 1:
      return GenotypeType.het
    else:
      return GenotypeType.hom_var


def genotype_as_alleles(variant):
  """Gets genotype of the sample in variant as a list of actual alleles.

  Returns the alleles specified by the genotype indices of variant.calls[0].
  For example, if variant.reference_bases = 'A' and variant.alternative_bases
  = ['C'] and the genotypes are [0, 1], this function will return
  ['A', 'C'].

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    A list of allele (string) from variant, one for each genotype in
    variant.calls[0], in order.

  Raises:
    ValueError: If variant doesn't have genotypes.
    ValueError: If variant has more than one call (i.e., is multi-sample).
  """
  if not has_genotypes(variant):
    raise ValueError('Not genotypes present in', variant)
  elif len(variant.calls) > 1:
    raise ValueError('Unsupported: multiple genotypes found at', variant)
  else:
    # Genotypes are encoded as integers, where 0 is the reference allele,
    # indices > 0 refer to alt alleles, and the no-call genotypes is encoded
    # as -1 in the genotypes. This code relies on this encoding to quickly
    # reference into the alleles by adding 1 to the genotype index.
    alleles = ['.', variant.reference_bases] + list(variant.alternate_bases)
    return [alleles[i + 1] for i in variant.calls[0].genotype]


def genotype_quality(variant, default=None):
  """Gets the genotype quality (GQ) value the genotype call in variant.

  If variant doesn't have genotypes, returns default, otherwise tries
  to retrieve the GQ field of the call field, returning that value if
  present otherwise returning default if its absent.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.
    default: The value for GQ to return if variant has no genotypes or
      if GQ is present in the genotype call record.

  Returns:
    The GQ value (may be a string or whatever value default is).
  """
  if not has_genotypes(variant):
    return default
  call = variant.calls[0]
  if 'GQ' in call.info:
    return call.info['GQ'].values[0].number_value
  else:
    return default


def is_gvcf(variant):
  """Returns true if variant encodes a standard gVCF reference block.

  This means in practice that variant has a single alternate allele that is the
  canonical gVCF allele GVCF_ALT_ALLELE constant exported here.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Returns:
    Boolean. True if variant is a gVCF record, False otherwise.
  """
  return variant.alternate_bases == [GVCF_ALT_ALLELE]


def genotype_ordering_in_likelihoods(variant):
  """Yields (i, j, allele_i, allele_j) for the genotypes ordering in GLs.

  https://samtools.github.io/hts-specs/VCFv4.1.pdf
  "If A is the allele in REF and B,C,... are the alleles as ordered in ALT,
  the ordering of genotypes for the likelihoods is given by:
  F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the ordering is:
  AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
  The biallelic sites in our case are 0/0, 0/1, 1/1.
  The triallelic sites are 0/0, 0/1, 1/1, 0/2, 1/2, 2/2.
  This wiki page has more information that generalizes ot different ploidy.
  http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes

  Currently this function only implements for diploid cases.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant.

  Yields:
    allele indices and strings (i, j, allele_i, allele_j) in the correct order.
  """
  alleles = [variant.reference_bases] + list(variant.alternate_bases)
  for j in range(0, len(variant.alternate_bases) + 1):
    for i in range(0, j + 1):
      yield i, j, alleles[i], alleles[j]
