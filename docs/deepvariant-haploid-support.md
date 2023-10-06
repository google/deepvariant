# DeepVariant support for variant calling in chromosome X and Y

## Case study

A case study on how to use the parameters mentioned here are described in
[DeepVariant haploid support](deepvariant-xy-calling-case-study.md).

## Haploid calling support

As DeepVariant is a diploid variant caller, it assigns genotypes as {Hom-ref,
Het, Hom-alt} for each candidate allele it observes. For samples with karyotype
XY, the chromosome X and Y are effectively haploid. So, we are introducing two
flags to re-adjust the genotypes in regions that are considered to be haploid
for those samples.

You can use `--haploid_contigs` and `--par_regions_bed` parameters to readjust
the genotypes in haploid regions. For samples with XY karyotype, it is expected
that users will set `--haploid_contigs="chrX,chrY"` for
[GRCh38](https://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa)
and `--haploid_contigs="X,Y"` for
[GRCh37](https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa).
You can also provide a PAR region bed file with
`--par_regions_bed="/input/GRCh3X_par.bed"` parameter. The regions in the PAR
bed file will be skipped from genotype readjustment. You can download the PAR
bed files from here:
[GRCh38_par.bed](https://storage.googleapis.com/deepvariant/case-study-testdata/GRCh38_PAR.bed),
[GRCh37_par.bed](https://storage.googleapis.com/deepvariant/case-study-testdata/GRCh37_PAR.bed).

## How it works

The genotype re-adjustment is implemented in the `postprocess_variants` stage of
DeepVariant. For any variant, that is in the`--haploid_contigs` regions and
**not** in the `--par_regions_bed` regions, the genotype likelihoods of
heterozygous variants are set as 0 and the genotypes are normalized again after
re-adjusting the likelihoods. After that the most-likely genotype is assigned to
the allele which excludes any heterozygous calls.

For example, suppose we observe an alternate allele `ALT1` at a position that we
consider to be haploid. So the observed alleles at that position are:
`Candidates: {REF, ALT1}` The neural network generates likelihoods for the
genotypes for this candidate as such:

```
Homozygous reference:   likelihood(REF,REF)
Heterozygous alternate: likelihood(REF,ALT1)
Homozygous alternaate:  likelihood(ALT1,ALT1)
```

So the likelihood vector looks like: `L={L[(REF, REF)], L[(REF, ALT1)], L[(ALT1,
ALT1)]}` In the post processing step we set `L[(REF, ALT1)] = 0` as that is a
likelihood associated with heterozygous genotype and heterozygous calling is
excluded in haploid regions. The likelihood vector becomes: `L={L[(REF, REF), 0,
L(ALT1, ALT1)]}`. Then we normalize the likelihood vector and assign the
genotype based on the adjusted values from the vector.
