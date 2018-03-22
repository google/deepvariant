"""Validates that a VCF file and a FASTA reference file correspond.

They correspond if:
a) they cover the same contigs,
b) the reference covers every variant in the vcf file, and
c) they agree on the reference bases covered by the variants.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

from absl import app

from deepvariant.util.io import fasta
from deepvariant.util.io import vcf
from deepvariant.util import variant_utils


def validate_contigs(ref_contigs, vcf_contigs):
  """Validate that the two lists of ContigInfos have the same set of names."""
  ref_names = {ci.name for ci in ref_contigs}
  vcf_names = {ci.name for ci in vcf_contigs}
  if ref_names != vcf_names:
    print('Contig names differ!')
    print('The following contigs are in one but not both: ')
    print(ref_names ^ vcf_names)
    sys.exit(-1)


def validate_variant(ref_reader, variant):
  """Validate that variant is covered by the reference and agrees with it."""
  var_range = variant_utils.variant_range(variant)
  ref_bases = ref_reader.query(var_range)
  if ref_bases != variant.reference_bases:
    print('In range {}:{}-{} '.format(
        var_range.reference_name, var_range.start, var_range.end))
    print('Reference says ', ref_bases)
    print('But variant says ', variant.reference_bases)
    sys.exit(-1)


def main(argv):
  if len(argv) != 3:
    print('Usage: {} <input_ref> <input_vcf>'.format(argv[0]))
    sys.exit(-1)
  in_ref = argv[1]
  in_vcf = argv[2]

  with fasta.RefFastaReader(in_ref) as ref_reader:
    with vcf.VcfReader(in_vcf, use_index=False) as vcf_reader:
      validate_contigs(ref_reader.header.contigs, vcf_reader.header.contigs)
      for variant in vcf_reader:
        validate_variant(ref_reader, variant)

  # VCF is valid!
  sys.exit(0)


if __name__ == '__main__':
  app.run(main)
