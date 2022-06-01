# DeepTrio

## Overview

DeepTrio is built on top of DeepVariant. It is intended for variant calling of
trios or duos. The main advantage of DeepTrio is that genetic inheritance is
considered by a neural network for calling variants in trio samples. Also,
variant candidates are generated from all samples at once, which ensures a
genotype call is made for any position in the trio with a variant. Since
DeepTrio is built on top of DeepVariant,
[general information](deepvariant-details.md) for DeepVariant also applies to
DeepTrio. At the highest level, a user needs to provide the following:

1.  A reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
    format and its corresponding
    [.fai index file](http://www.htslib.org/doc/faidx.html) generated using the
    `samtools faidx` command.

1.  An aligned reads files for child and one or two parents in
    [BAM](http://genome.sph.umich.edu/wiki/BAM) format and its corresponding
    index file (.bai). The reads must be aligned to the reference genome
    described above.

The output of DeepTrio is a set of variants in
[VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format representing the
child and one or two parents.

Similar to DeepVariant, DeepTrio is composed of three stages: `make_examples`,
`call_variants`, and `postprocess_variants`. Some of the components (
`call_variants`, `postprocess_variants`) are shared with DeepVariant, and
`make_examples` is specialized for DeepTrio. More details about each program are
described in detail in the
[Inputs and outputs](deepvariant-details.md#inputs-and-outputs) section of the
DeepVariant documentation.

DeepTrio comes with three models for different types of input data:

*   Illumina whole genome data (WGS).
*   Illumina whole exome data (WES).
*   PacBio HiFi whole genome data (PacBio WGS).

## Running DeepTrio

The easiest and recommended way to run DeepTrio is using
`google/deepvariant:deeptrio-latest` docker image. Please refer to the
[quick start guide](deeptrio-quick-start.md) for more details on how to run
DeepTrio using docker.

Merging VCFs can be done using
[GLnexus](https://github.com/dnanexus-rnd/GLnexus) which has been optimized for
use with DeepVariant gVCFs. The process is described in the DeepTrio case
studies
([DeepTrio whole genome sequencing case study](deeptrio-wgs-case-study.md) and
[Using DeepTrio for small variant calling from the trio sequenced with PacBio
HiFi](deeptrio-pacbio-case-study.md)), and in the manuscript,
["Accurate, scalable cohort variant calls using DeepVariant and GLnexus"](https://www.biorxiv.org/content/10.1101/2020.02.10.942086v2).

Please note that DeepTrio can be run with a `run_deeptrio.py` script that
automates all DeepTrio steps and thus greatly simplifies the inference pipeline.
The details of using this script can be found in the section below as well as in
the DeepTrio case studies.

Also please note: for the non-PAR regions of the sex chromosomes (X and Y), we
recommend running these providing only the parent who contributed the child's
chromosome (e.g. for chromosomeX, only the mother and son samples and for
chromosomeY only the father and son samples).

If needed, DeepTrio can be built from source. For more details please refer to
[Building DeeepTrio](deeptrio-build-test.md).

## DeepTrio Input assumptions

The reference genome FASTA, passed in using the `--ref` flag, must be indexed
and can either be uncompressed or compressed with `bgzip`.

All BAM files should be aligned to a "compatible" version of the genome
reference provided as the `--ref`. DeepTrio will only process contigs shared by
both the BAM and reference. BAM files must be also sorted and indexed. They must
exist on disk, so you cannot pipe them into DeepTrio. Duplicate marking may be
performed. In our analyses, there is almost no difference in accuracy with and
without duplicate marking except at lower (<20x) coverages. Finally, we
recommend that you do not perform BQSR. Running BQSR has a small decrease on
accuracy.

If you are providing `--regions` or other similar arguments, these should refer
to contigs present in the reference genome. These arguments accept
space-separated lists, so all of the follow examples are valid arguments for
`--regions` or similar arguments:

*   `--regions chr20` => only process all of chromosome 20
*   `--regions chr20:10,000,000-11,000,000` => only process 10-11mb of chr20
*   `--regions "chr20 chr21"` => only process chromosomes 20 and 21

## Training data

DeepTrio models are trained using the latest publicly avavilable GIAB
benchmarks. You can find more details about the training data for each DeepTrio
model in the
[DeepTrio Training Data document](deeptrio-details-training-data.md).

## DeepVariant dependency

DeepTrio is built on top of DeepVariant and they share most of the components.
Please see [DeepVariant usage guide](deepvariant-details.md) for a full
description of DeepVariant components as well as other consideration for running
DeepVariant pipeline.
