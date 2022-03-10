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
"""Variant utilities."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import enum
import itertools

import six

from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variantcall_utils
from third_party.nucleus.util import vcf_constants


def only_call(variant):
  """Ensures the Variant has exactly one VariantCall, and returns it.

  Args:
    variant: nucleus.genomics.v1.Variant. The variant of interest.

  Returns:
    The single nucleus.genomics.v1.VariantCall in the variant.

  Raises:
    ValueError: Not exactly one VariantCall is in the variant.
  """
  if len(variant.calls) != 1:
    raise ValueError('Expected exactly one VariantCall in {}'.format(variant))
  return variant.calls[0]


def decode_variants(encoded_iter):
  """Yields a genomics.Variant from encoded_iter.

  Args:
    encoded_iter: An iterable that produces binary encoded
      nucleus.genomics.v1.Variant strings.

  Yields:
    A parsed nucleus.genomics.v1.Variant for each encoded element of
    encoded_iter in order.
  """
  for encoded in encoded_iter:
    yield variants_pb2.Variant.FromString(encoded)


def variant_position(variant):
  """Returns a new Range at the start position of variant.

  Args:
    variant: nucleus.genomics.v1.Variant.

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
    variant: nucleus.genomics.v1.Variant.

  Returns:
    A new Range with the same reference_name, start, and end as variant.
  """
  return ranges.make_range(variant.reference_name, variant.start, variant.end)


def variant_range_tuple(variant):
  """Returns a new tuple of (reference_name, start, end) for the variant.

  A common use case for this function is to sort variants by chromosomal
  location, with usage like `sorted(variants, key=variant_range_tuple)`.

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    A three-tuple with the same reference_name, start, and end as variant.
  """
  return (variant.reference_name, variant.start, variant.end)


@enum.unique
class GenotypeType(enum.Enum):
  """An enumeration of the types of genotypes."""
  hom_ref = ('homozygous reference', [0, 0], 0)
  het = ('heterozygous', [0, 1], 1)
  hom_var = ('homozygous non-reference', [1, 1], 2)
  no_call = ('no call', [-1, -1], -1)

  def __init__(self, full_name, example_gt, class_id):
    """Create a GenotypeType with the given name, GT and class_id."""
    self.full_name = full_name
    self.example_gt = example_gt
    self.class_id = class_id


@enum.unique
class VariantType(enum.Enum):
  """An enumeration of the types of variants."""
  # A variant.proto where there is no alt allele.
  ref = 0
  # A non-reference variant.proto where all ref and alt alleles are single
  # basepairs.
  snp = 1
  # A non-reference variant.proto where at least one of ref or alt alleles are
  # longer than 1 bp.
  indel = 2


def format_filters(variant):
  """Returns a human-readable string showing the filters applied to variant.

  Returns a string with the filter field values of variant separated by commas.
  If the filter field isn't set, returns vcf_constants.MISSING_FIELD ('.').

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    A string.
  """
  if variant.filter:
    return ','.join(variant.filter)
  else:
    return vcf_constants.MISSING_FIELD


def format_alleles(variant):
  """Gets a string representation of the variant's alleles.

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    A string ref_bases/alt1,alt2 etc.
  """
  return '{}/{}'.format(variant.reference_bases, ','.join(
      variant.alternate_bases))


def format_position(variant):
  """Gets a string representation of the variant's position.

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    A string chr:start + 1 (as start is zero-based).
  """
  return '{}:{}'.format(variant.reference_name, variant.start + 1)


def _non_excluded_alts(alts, exclude_alleles=None):
  """Exclude any alts listed, by default: '<*>', '.', and '<NON_REF>'.

  These alleles are sometimes listed in ALT column but they shouldn't be
  analyzed and usually indicate reference blocks in formats like gVCF.

  E.g. 'A'->'<*>' is NOT an insertion, and 'A'->'.' is NOT a SNP.

  Args:
    alts: a list of strings representing the alternate alleles.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    alts alleles except those in exclude_alleles, by default excluding the GVCF
    '<*>' allele, the '<NON_REF>' symbolic allele, and '.' missing field by
    default.
  """
  if exclude_alleles is None:
    exclude_alleles = [
        vcf_constants.GVCF_ALT_ALLELE, vcf_constants.SYMBOLIC_ALT_ALLELE,
        vcf_constants.MISSING_FIELD
    ]
  return [a for a in alts if a not in exclude_alleles]


def is_snp(variant, exclude_alleles=None):
  """Is variant a SNP?

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if all alleles of variant are 1 bp in length.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  return (len(variant.reference_bases) == 1 and len(relevant_alts) >= 1 and
          all(len(x) == 1 for x in relevant_alts))


def is_indel(variant, exclude_alleles=None):
  """Is variant an indel?

  An indel event is simply one where the size of at least one of the alleles
  is > 1.

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if the alleles in variant indicate an insertion/deletion event
    occurs at this site.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  if not relevant_alts:
    return False
  return (len(variant.reference_bases) > 1 or
          any(len(alt) > 1 for alt in relevant_alts))


def is_biallelic(variant, exclude_alleles=None):
  """Returns True if variant has exactly one alternate allele.

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if the variant has exactly one alternate allele.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  return len(relevant_alts) == 1


def is_multiallelic(variant, exclude_alleles=None):
  """Does variant have multiple alt alleles?

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if variant has more than one alt allele.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  return len(relevant_alts) > 1


def variant_is_insertion(variant, exclude_alleles=None):
  """Are all the variant's alt alleles insertions?

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if variant has at least one alt allele and all alts are insertions.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  if not relevant_alts:
    return False
  return all(
      is_insertion(variant.reference_bases, alt) for alt in relevant_alts)


def variant_is_deletion(variant, exclude_alleles=None):
  """Are all the variant's alt alleles deletions?

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if variant has at least one alt allele and all alts are deletions.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  if not relevant_alts:
    return False
  return all(is_deletion(variant.reference_bases, alt) for alt in relevant_alts)


def is_ref(variant, exclude_alleles=None):
  """Returns true if variant is a reference record.

  Variant protos can encode sites that aren't actually mutations in the
  sample. For example, the record ref='A', alt='.' indicates that there is
  no mutation present (i.e., alt is the missing value).

  Args:
    variant: nucleus.genomics.v1.Variant.
    exclude_alleles: list(str). The alleles in this list will be ignored.

  Returns:
    True if there are no actual alternate alleles.
  """
  relevant_alts = _non_excluded_alts(variant.alternate_bases, exclude_alleles)
  return not relevant_alts


def variant_type(variant):
  """Gets the VariantType of variant.

  Args:
    variant: nucleus.genomics.v1.Variant.

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
    variant: nucleus.genomics.v1.Variant.

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
    variant: nucleus.genomics.v1.Variant.

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
  # Duplicate alleles.
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
    evalv: A nucleus.genomics.v1.Variant.
    truev: A nucleus.genomics.v1.Variant.

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
  if any(len(match) != 1 for match in six.itervalues(allele_matches)):
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


def simplify_variant_alleles(variant):
  """Replaces the alleles in variants with their simplified versions.

  This function takes a variant and replaces its ref and alt alleles with those
  produced by a call to variant_utils.simplify_alleles() to remove common
  postfix bases in the alleles that may be present due to pruning away alleles.

  Args:
    variant: learning.genomics.genomics.Variant proto we want to simplify.

  Returns:
    variant with its ref and alt alleles replaced with their simplified
      equivalents.
  """
  simplified_alleles = simplify_alleles(variant.reference_bases,
                                        *variant.alternate_bases)
  variant.reference_bases = simplified_alleles[0]
  variant.alternate_bases[:] = simplified_alleles[1:]
  variant.end = variant.start + len(variant.reference_bases)
  return variant


def is_filtered(variant):
  """Returns True if variant has a non-PASS filter field, or False otherwise."""
  return bool(variant.filter) and any(
      f not in {'PASS', vcf_constants.MISSING_FIELD} for f in variant.filter)


def is_variant_call(variant,
                    require_non_ref_genotype=True,
                    no_calls_are_variant=False,
                    call_indices=None,
                    apply_filter=True):
  """Is variant a non-reference call?

  A Variant proto doesn't always imply that there's a variant present in the
  genome. The call may not have alternate bases, may be filtered, may a have
  hom-ref genotype, etc. This function looks for all of those configurations
  and returns true iff the variant is asserting that a mutation is present
  in the same.

  Note that this code allows a variant without a calls field to be variant,
  but one with a genotype call must have a non-reference genotype to be
  considered variant (if require_non_ref_genotype is True, the default). If
  False, a variant that passes all of the site-level requirements for being
  a variant_call will return a True value, regardless of the genotypes, which
  means that we'll consider a site with a sample with a hom-ref or no-call site
  a variant call.

  Args:
    variant: nucleus.genomics.v1.Variant.
    require_non_ref_genotype: Should we require a site with a genotype call to
      have a non-reference (het, hom-var) genotype for the site to be considered
      a variant call?
    no_calls_are_variant: If a site has genotypes, should we consider no_call
      genotypes as being variant or not? e.g. -1/1 listed as ./. in VCF
    call_indices: A list of 0-based indices. If specified, only the calls
      at the given indices will be considered. The function will return
      True if any of those calls are variant.
    apply_filter: If set to True, will never treat this site as variant when
      any filter other than PASS or . is set.

  Returns:
    True if variant is really a mutation call.
  """
  if is_ref(variant):
    # No actual alt allele listed in ALT column
    return False
  elif apply_filter and is_filtered(variant):
    # Anything other than PASS or . in FILTER column
    return False
  elif not variant.calls or not require_non_ref_genotype:
    return True
  # All tests after this point should only look at genotype-based fields, as
  # we may have aborted out in the prev. line due to require_non_ref_genotype.
  else:
    # Check for non-ref genotypes and optionally no-call (-1) genotypes
    if call_indices is None:
      call_indices = range(len(variant.calls))
    for i in call_indices:
      for g in variant.calls[i].genotype:
        if g > 0 or (no_calls_are_variant and g < 0):
          return True
    return False


def has_calls(variant):
  """Does variant have any genotype calls?

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    True if variant has one or more VariantCalls.
  """
  return bool(variant.calls)


def genotype_type(variant):
  """Gets the GenotypeType for variant.

  If variant doesn't have genotypes, returns no_call. Otherwise
  returns one of no_call, hom_ref, het, or hom_var depending on the
  status of the genotypes in the call field of variant.

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    A GenotypeType.

  Raises:
    ValueError: If variant has more than one call (i.e., is multi-sample).
  """
  if not has_calls(variant):
    return GenotypeType.no_call
  elif len(variant.calls) > 1:
    raise ValueError('Unsupported: multiple genotypes found at', variant)
  else:
    gt = set(only_call(variant).genotype)
    if gt == {-1}:
      return GenotypeType.no_call
    elif gt == {0}:
      return GenotypeType.hom_ref
    elif len(gt) > 1:
      return GenotypeType.het
    else:
      return GenotypeType.hom_var


def genotype_as_alleles(variant, call_ix=0):
  """Gets genotype of the sample in variant as a list of actual alleles.

  Returns the alleles specified by the genotype indices of variant.calls[0].
  For example, if variant.reference_bases = 'A' and variant.alternative_bases
  = ['C'] and the genotypes are [0, 1], this function will return
  ['A', 'C'].

  Args:
    variant: nucleus.genomics.v1.Variant.
    call_ix: int. The index into the calls attribute indicating which
      VariantCall to return alleles for.

  Returns:
    A list of allele (string) from variant, one for each genotype in
    variant.calls[call_ix], in order.

  Raises:
    ValueError: If variant doesn't have a call at the specified index.
  """
  if not 0 <= call_ix < len(variant.calls):
    raise ValueError(
        'Unsupported: requesting call {} in variant with {} calls: {}'.format(
            call_ix, len(variant.calls), variant))
  else:
    # Genotypes are encoded as integers, where 0 is the reference allele,
    # indices > 0 refer to alt alleles, and the no-call genotypes is encoded
    # as -1 in the genotypes. This code relies on this encoding to quickly
    # reference into the alleles by adding 1 to the genotype index.
    alleles = ([vcf_constants.MISSING_FIELD, variant.reference_bases] +
               list(variant.alternate_bases))
    return [alleles[i + 1] for i in variant.calls[call_ix].genotype]


def unphase_all_genotypes(variant):
  """Sorts genotype and removes phasing bit of all calls in variant.

  This mutation is done in place rather than returning a different copy.

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    The variant with unphased calls.
  """
  for c in variant.calls:
    c.is_phased = False
    c.genotype.sort()
  return variant


def is_gvcf(variant):
  """Returns true if variant encodes a standard gVCF reference block.

  This means in practice that variant has a single alternate allele that is the
  canonical gVCF allele vcf_constants.GVCF_ALT_ALLELE.

  Args:
    variant: nucleus.genomics.v1.Variant.

  Returns:
    Boolean. True if variant is a gVCF record, False otherwise.
  """
  return variant.alternate_bases == [vcf_constants.GVCF_ALT_ALLELE]


def _genotype_order_in_likelihoods(num_alts, ploidy=2):
  """Yields tuples of `ploidy` ints for the given number of alt alleles.

  https://samtools.github.io/hts-specs/VCFv4.1.pdf
  "If A is the allele in REF and B,C,... are the alleles as ordered in ALT,
  the ordering of genotypes for the likelihoods is given by:
  F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the ordering is:
  AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
  The biallelic sites in our case are 0/0, 0/1, 1/1.
  The triallelic sites are 0/0, 0/1, 1/1, 0/2, 1/2, 2/2.
  This wiki page has more information that generalizes to different ploidy.
  http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes

  Args:
    num_alts: int. The number of alternate alleles at the site.
    ploidy: int. The ploidy for which to return genotypes.

  Yields:
    Tuples of `ploidy` ints representing allele indices in the order they appear
    in the corresponding genotype likelihood array.
  """
  if ploidy == 1:
    for i in range(num_alts + 1):
      yield (i,)
  elif ploidy == 2:
    for j in range(num_alts + 1):
      for i in range(j + 1):
        yield (i, j)
  else:
    raise NotImplementedError('Only haploid and diploid supported.')


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
    variant: nucleus.genomics.v1.Variant.

  Yields:
    allele indices and strings (i, j, allele_i, allele_j) in the correct order.
  """
  alleles = [variant.reference_bases] + list(variant.alternate_bases)
  for i, j in _genotype_order_in_likelihoods(
      len(variant.alternate_bases), ploidy=2):
    yield i, j, alleles[i], alleles[j]


def genotype_likelihood(variant_call, allele_indices):
  """Returns the genotype likelihood for the given allele indices.

  Args:
    variant_call: nucleus.genomics.v1.VariantCall. The VariantCall from
      which to extract the genotype likelihood of the allele indices.
    allele_indices: list(int). The list of allele indices for a given genotype.
      E.g. diploid heterozygous alternate can be represented as [0, 1].

  Returns:
    The float value of the genotype likelihood of this set of alleles.
  """
  return variant_call.genotype_likelihood[genotype_likelihood_index(
      allele_indices)]


def genotype_likelihood_index(allele_indices):
  """Returns the genotype likelihood index for the given allele indices.

  Args:
    allele_indices: list(int). The list of allele indices for a given genotype.
      E.g. diploid homozygous reference is represented as [0, 0].

  Returns:
    The index into the associated genotype likelihood array corresponding to
    the likelihood of this list of alleles.

  Raises:
    NotImplementedError: The allele_indices are more than diploid.
  """
  if len(allele_indices) == 1:
    # Haploid case.
    return allele_indices[0]
  elif len(allele_indices) == 2:
    # Diploid case.
    g1, g2 = sorted(allele_indices)
    return g1 + (g2 * (g2 + 1) // 2)
  else:
    raise NotImplementedError(
        'Genotype likelihood index only supports haploid and diploid: {}'.
        format(allele_indices))


def allele_indices_for_genotype_likelihood_index(gl_index, ploidy=2):
  """Returns a tuple of allele_indices corresponding to the given GL index.

  This is the inverse function to `genotype_likelihood_index`.

  Args:
    gl_index: int. The index within a genotype likelihood array for which to
      determine the associated alleles.
    ploidy: int. The ploidy of the result.

  Returns:
    A tuple of `ploidy` ints representing the allele indices at this GL index.

  Raises:
    NotImplementedError: The requested allele indices are more than diploid.
  """
  if ploidy == 1:
    return gl_index
  elif ploidy == 2:
    # TODO: Implement using algorithm described at
    # https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
    # rather than creating all genotypes explicitly.
    num_alts = 1
    while genotype_likelihood_index([num_alts, num_alts]) < gl_index:
      num_alts += 1
    genotypes = list(_genotype_order_in_likelihoods(num_alts, ploidy=ploidy))
    return genotypes[gl_index]
  else:
    raise NotImplementedError(
        'Allele calculations only supported for haploid and diploid.')


def allele_indices_with_num_alts(variant, num_alts, ploidy=2):
  """Returns a list of allele indices configurations with `num_alts` alternates.

  Args:
    variant: nucleus.genomics.v1.Variant. The variant of interest, which
      defines the candidate alternate alleles that can be used to generate
      allele indices configurations.
    num_alts: int in [0, `ploidy`]. The number of non-reference alleles for
      which to create the allele indices configurations.
    ploidy: int. The ploidy for which to return allele indices configurations.

  Returns: A list of tuples. Each tuple is of length `ploidy` and represents the
    allele indices of all `ploidy` genotypes that contain `num_alts`
    non-reference alleles.

  Raises:
    ValueError: The domain of `num_alts` is invalid.
    NotImplementedError: `ploidy` is not diploid.
  """
  if ploidy != 2:
    raise NotImplementedError(
        'allele_indices_with_num_alts only supports diploid.')
  if not 0 <= num_alts <= ploidy:
    raise ValueError(
        'Invalid number of alternate alleles requested: {} for ploidy {}'.
        format(num_alts, ploidy))

  max_candidate_alt_ix = len(variant.alternate_bases)
  if num_alts == 0:
    return [(0, 0)]
  elif num_alts == 1:
    return [(0, i) for i in range(1, max_candidate_alt_ix + 1)]
  else:
    return [(i, j)
            for i in range(1, max_candidate_alt_ix + 1)
            for j in range(i, max_candidate_alt_ix + 1)]


def variants_overlap(variant1, variant2):
  """Returns True if the range of variant1 and variant2 overlap.

  This is equivalent to:

    ranges_overlap(variant_range(variant1), variant_range(variant2))

  Args:
    variant1: nucleus.genomics.v1.Variant we want to compare for overlap.
    variant2: nucleus.genomics.v1.Variant we want to compare for overlap.

  Returns:
    True if the variants overlap, False otherwise.
  """
  return ranges.ranges_overlap(variant_range(variant1), variant_range(variant2))


def variant_key(variant, sort_alleles=True):
  """Gets a human-readable string key that is almost unique for Variant.

  Gets a string key that contains key information about the variant, formatted
  as:

    reference_name:start+1:reference_bases->alternative_bases

  where alternative bases is joined with a '/' for each entry in
  alternative_bases. The start+1 is so we display the position, which starts at
  1, and not the offset, which starts at 0.

  For example, a Variant(reference_name='20', start=10, reference_bases='AC',
  alternative_bases=['A', 'ACC']) would have a key of:

    20:11:AC->A/ACC

  The key is 'almost unique' in that the reference_name + start + alleles should
  generally occur once within a single VCF file, given the way the VCF
  specification works.

  Args:
    variant: nucleus.genomics.v1.Variant to make into a key.
    sort_alleles: bool. If True, the alternative_bases of variant will be sorted
      according to their lexicographic order. If False, the alternative_bases
      will be displayed in their order in the Variant.

  Returns:
    A str.
  """
  alts = variant.alternate_bases
  if sort_alleles:
    alts = sorted(alts)
  return '{}:{}:{}->{}'.format(variant.reference_name, variant.start + 1,
                               variant.reference_bases, '/'.join(alts))


def sorted_variants(variants):
  """Returns sorted(variants, key=variant_range_tuple)."""
  return sorted(variants, key=variant_range_tuple)


def variants_are_sorted(variants):
  """Returns True if variants are sorted w.r.t. variant_range.

  Args:
    variants: list[nucleus.genomics.v1.Variant]. A list of Variant
      protos that may or may not be sorted.

  Returns:
    True if variants are sorted, False otherwise.
  """
  def _pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return six.moves.zip(a, b)

  for r1, r2 in _pairwise(variant_range_tuple(v) for v in variants):
    if r2 < r1:
      return False
  return True


def set_info(variant, field_name, value, vcf_object=None):
  """Sets a field of the info map of the `Variant` to the given value(s).

  `variant.info` is analogous to the INFO field of a VCF record.

  Args:
    variant: Variant proto. The Variant to modify.
    field_name: str. The name of the field to set.
    value: A single value or list of values to update the Variant with. The type
      of the value is determined by the `vcf_object` if one is given, otherwise
      is looked up based on the reserved INFO fields in the VCF specification.
    vcf_object: (Optional) A VcfReader or VcfWriter object. If not None, the
      type of the field is inferred from the associated VcfReader or VcfWriter
      based on its name. Otherwise, the type is inferred if it is a reserved
      field.
  """
  if vcf_object is None:
    set_field_fn = vcf_constants.reserved_info_field_set_fn(field_name)
  else:
    set_field_fn = vcf_object.field_access_cache.info_field_set_fn(field_name)
  set_field_fn(variant.info, field_name, value)


def get_info(variant, field_name, vcf_object=None):
  """Returns the value of the `field_name` INFO field.

  The `vcf_object` is used to determine the type of the resulting value. If it
  is a single value or a Flag, that single value will be returned. Otherwise,
  the list of values is returned.

  Args:
    variant: Variant proto. The Variant of interest.
    field_name: str. The name of the field to retrieve values from.
    vcf_object: (Optional) A VcfReader or VcfWriter object. If not None, the
      type of the field is inferred from the associated VcfReader or VcfWriter
      based on its name. Otherwise, the type is inferred if it is a reserved
      field.
  """
  if vcf_object is None:
    get_field_fn = vcf_constants.reserved_info_field_get_fn(field_name)
  else:
    get_field_fn = vcf_object.field_access_cache.info_field_get_fn(field_name)
  return get_field_fn(variant.info, field_name)


def calc_ac(variant):
  """Returns a list of alt counts based on variant.calls."""
  counts = [0] * len(variant.alternate_bases)
  for call in variant.calls:
    for gt in call.genotype:
      if gt > 0:
        counts[gt - 1] += 1
  return counts


def calc_an(variant):
  """Returns the total number of alleles in called genotypes in variant."""
  return sum(
      len([1 for gt in call.genotype if gt > -1]) for call in variant.calls)


def is_singleton(variant):
  """Returns True iff the variant has exactly one non-ref VariantCall."""
  non_ref_count = 0
  for c in variant.calls:
    if variantcall_utils.has_variation(c):
      non_ref_count += 1
      if non_ref_count > 1:
        return False
  return non_ref_count == 1


def major_allele_frequency(variant):
  """Returns the frequency of the most common allele in the variant."""
  count_dict = collections.defaultdict(int)
  for call in variant.calls:
    for geno in call.genotype:
      count_dict[geno] += 1
  denom = sum(cnt for geno, cnt in count_dict.items() if geno >= 0)
  if denom > 0:
    numer = max(cnt for geno, cnt in count_dict.items() if geno >= 0)
    return float(numer) / denom
  else:
    return 0
