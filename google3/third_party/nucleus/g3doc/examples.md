# Example programs

After running ``source install.sh``, the executable versions of the
example programs can be located in ``bazel-bin/nucleus/examples/``.
For example, to run ``ascii_pileup``, you would actually run a command
like

```shell
bazel-bin/nucleus/examples/ascii_pileup input.sam chr3:99393
```

If you would like to rebuild just the example programs -- after modifying
one of them, perhaps -- this can be done with

```shell
bazel build -c opt $COPT_FLAGS nucleus/examples:all
```

Here is a summary of the [example
programs](https://github.com/google/nucleus/blob/master/nucleus/examples)
included with Nucleus:

* ``add_ad_to_vcf input.vcf output.vcf``

  If ``input.vcf`` is a VCF file in which the variant calls have
  allele depths (in their 'AD' FORMAT fields), then ``output.vcf`` will
  be the same but with the allele depths summed into the 'AD' INFO field
  of the variants.  This program is a good example of writing out a file
  with a modified header, as well as demonstrating the
  ``variant_utils`` and ``variantcall_utils`` routines for setting and
  getting INFO and FORMAT information.

* ``ascii_pileup input.sam chrX:3029``

  This will print the reads from ``input.sam`` that overlap with the
  given location.  The location is highlighted in each read and the
  reads are sorted and visually aligned by location, so the end
  result is surprisingly pretty, particularly if you have a wide monitor.
  Demonstrates the SAM reader and using range queries.

* ``count_variants input.vcf``

  Prints the number of variants in input.vcf, broken down by type
  (ref/SNP/indel) and by chromosome.  If you want to just read in a VCF
  file and print out some information about it, this is a good one
  to copy.

* ``filter_vcf input.vcf output.vcf``

  ``output.vcf`` will contain all of the variants in ``input.vcf`` that
  have a quality score greater than 3.01.  This is the example program
  used in the [overview](overview.md).

* ``ngs_errors --ref /tmp/hs37d5.fa.gz --vcf /tmp/NA12878_calls.vcf.gz --bam /tmp/NA12878_sliced.bam  --examples_out /tmp/examples.tfrecord``

  This program creates examples which are suitable (but not ideal)
  for training a TensorFlow model to predict the true sequence for a read.
  The examples are taken from the reads in the BAM file, but only those that
  don't overlap with any variants in the input VCF file.  For each read,
  a tf.Example protocol buffer is made containing
    + the read sequence,
    + the true sequence from the reference FASTA input,
    + the read name,
    + the read alignment cigar, and
    + the alignment qualities.
  In addition, the program comments explain how to download the required
  input files and transform them into the right formats.

  This is the longest example, but it really displays the power of Nucleus in
  taking genomics data and turning it into machine learning inputs.

* ``validate_vcf ref.fasta input.vcf``

  This will print a warning if the input FASTA reference file and the
  input VCF file are mismatched, which can happen if
    a) they mention a different set of contigs,
    b) the VCF files contains variants not covered by the reference,
    c) the VCF file contains a variant covered by the reference, but
       the VCF and FASTA file disagree about the correct reference bases
       for the range.
  This program shows how to use the ``RefFastaReader``, and to query it
  using ranges taken from a VCF file.
