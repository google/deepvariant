# Overview

Here is an example Nucleus program:

```python
from nucleus.io import vcf

with vcf.VcfReader('/tmp/example.vcf.gz') as reader:
  print('Sample names in VCF: ', ' '.join(reader.header.sample_names))
  with vcf.VcfWriter('/tmp/filtered.tfrecord', header=reader.header) as writer:
    for variant in reader:
      if variant.quality > 3.01:
        writer.write(variant)
```

Let's go through it line by line.

```python
from nucleus.io import vcf
```

`nucleus/io` is the directory containing Nucleus's Python classes for
doing reading and writing.  There are modules for each of the genomics
formats that Nucleus supports, including FASTA, FASTQ, SAM/BAM and VCF.
That directory also contains the base classes that define the APIs for all
of the various readers and writers.  In brief, each reader can either be
iterated over or queried by a region, and every writer has just the single
method `write`.

```python
with vcf.VcfReader('/tmp/example.vcf.gz') as reader:
```

We use Python's `with` pattern because the reader object returned by
`vcf.VcfReader` is smart enough to close the file at the end of the block
when we do.  `vcf.VcfReader` is also smart enough to realize that this
particular input file is compressed, and deal with that transparently.
Of course, vcf files don't have to be compressed; `VcfReader` would have
accepted an uncompressed `/tmp/example2.vcf` just as well.

By default, we *do* assume that VCF files come along with a
[TABIX](http://www.htslib.org/doc/tabix.html) index that allows us to
efficiently do random access queries on genomic ranges.  For
`/tmp/example.vcf.gz` we would assume the index is in the file
`/tmp/example.vcf.gz.tbi`.  To turn off this assumption, pass the
`use_index=False` keyword parameter to `VcfReader`.

```python
  print('Sample names in VCF: ', ' '.join(reader.header.sample_names))
```

This does what you would expect; it prints the names of all the samples
present in the input VCF file.  Every Nucleus reader object contains a
`.header` field with information about the file's header and metadata.
(For some file types, like reference FASTA files, the header information
also comes from the file's index's header.)

The names of the subfields of the `.header` of a VCF file can be found in
the definition of the protocol buffer message `VcfHeader` in
[variants.proto](https://github.com/google/nucleus/blob/master/nucleus/protos/variants.proto).
Using protocol buffers is mostly intuitive -- you can treat messages as
named tuples, repeated fields as lists, map fields as dictionaries, etc. --
and for more information you can refer to the [protocol buffer
documentation](https://developers.google.com/protocol-buffers/docs/pythontutorial).

```python
  with vcf.VcfWriter('/tmp/filtered.tfrecord', header=reader.header) as writer:
```

You might be surprised by the file name here.  `vcf.VcfWriter` can certainly
write to VCF files, and if called with say `vcf.VcfWriter('/tmp/out.vcf')`
or `vcf.VcfWriter('/tmp/out.vcf.gz')` it would do just that.

But Nucleus reader and writers can also read and write from the
[TFRecords](https://www.tensorflow.org/api_guides/python/python_io) file
format used by [TensorFlow](https://www.tensorflow.org).  As the above
example suggests, the desired file format is indicated by the file extensions:
'.tfrecord' (or '.tfrecord.gz') means TFRecords, while '.vcf' means VCF.

It has been our experience that this ability to easily switch between native
genomics files and TFRecords files makes it easier to gradually introduce
machine learning into genomics workflows.  We hope you find it convenient as
well!

```python
    for variant in reader:
```

As mentioned above, every Nucleus reader supports iteration.  You can do it
this way, directly using the reader object, or if you need a Python iterable
to pass to another routine, you can get one from calling `reader.iterate()`.

Each `variant` in this loop will be a `Variant` protocol buffer object,
defined in
[variants.proto](https://github.com/google/nucleus/blob/master/nucleus/protos/variants.proto).
Most of the `Variant` fields are easy to work with:  `variant.reference_name`
is the CHROM field, `variant.start` is the starting position, etc.  But the
`variant.info` and `variant_call.info` map fields are a little bit trickier
because they contain typed information, and those types depend on the header
definitions.  To make things easier, Nucleus defines some convenience routines
for setting and getting these fields in
[util/variant_utils.py](https://github.com/google/nucleus/blob/master/nucleus/util/variant_utils.py)
and 
[util/variantcall_utils.py](https://github.com/google/nucleus/blob/master/nucleus/util/variantcall_utils.py).

```python
      if variant.quality > 3.01:
        writer.write(variant)
```

These lines should now be easy to understand -- we test the `quality` field
of the `Variant`, and if it is high enough, we write it to the output.
Unsurprisingly, the `write` method of Nucleus writers takes a protocol buffer
of the appropriate type as input.

(If you are wondering where the 3.01 came from, quality is measured on a
Phred-scale, quality = -10 * log_10(probability), and
-10 * log_10(0.5) is approximately 3.01.)

Here's a table to summarize the file types currently supported by Nucleus,
and their associated protocol buffer types:

Format    | Record Type | Header Type    | Reader? | Writer?
--------- | ----------- | -------------- | ------- | -------
BED       | BedRecord   | BedHeader      | Y       | Y
FASTA Ref | string      | RefFastaHeader | Y       | N
FastQ     | FastqRecord | none           | Y       | Y
SAM/BAM   | Read        | SamHeader      | Y       | N
VCF       | Variant     | VcfHeader      | Y       | Y

Don't despair if your favorite genomics format isn't listed, though, as we
hope to add more soon.  (And
[contributions](https://github.com/google/nucleus/blob/master/CONTRIBUTING.md)
are welcomed!)
