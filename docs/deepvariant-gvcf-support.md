# DeepVariant Genomic VCF (gVCF) support

Beginning with the 0.5.0 release, DeepVariant supports the creation of Genomic
VCF (gVCF) output. This has the same underlying format specification as the
[VCF format] but also includes additional records that distinguish regions that
have sequence coverage that appears to match the reference genome from regions
without sequence coverage, in which the genotype is unknown.

gVCF files are required as input for analyses that create a set of variants in
a cohort of individuals, such as cohort merging or joint genotyping.

## Description of gVCF format

When run with gVCF output enabled, DeepVariant generates both the VCF output
containing only variant calls as well as an additional gVCF output file that
contains both variants and non-variant sites. The gVCF file includes both
variant calls and regions that are confidently called as matching the reference
genome. The non-variant sites compare the reference allele to an "unspecified
alternate" allele, represented by `<*>`. To minimize output file size, adjacent
records with equal (or similar, see discussion below) genotype qualities are
merged into a single record.

Section 5.5 of the [VCF format] specification gives a description of the gVCF
format and example output, partially reproduced below. The gVCF output of
DeepVariant is syntactically and semantically equivalent to this example.

```bash
#CHROM  POS ID  REF  ALT     QUAL  FILTER  INFO       FORMAT   Sample
1      4370  .   G   <*>       .     .     END=4383   GT:GQ    0/0:37
1      4384  .   C   <*>       .     .     END=4388   GT:GQ    0/0:41
1      4389  .   T   TC,<*>   50     .     .          GT:GQ    0/1:50
1      4390  .   C   <*>       .     .     END=4390   GT:GQ    0/0:3
```

## Creating gVCF output with DeepVariant

The exact same three programs (`make_examples`, `call_variants`, and
`postprocess_variants`) are used when creating gVCF output as in the [WGS case
study]. However, additional flags must be passed to the `make_examples` and
`postprocess_variants` steps.

### `make_examples`

The `make_examples` program is where the gVCF records are computed.

One additional flag is required in `make_examples`, the `--gvcf <filename>`
flag. This specifies an additional output, which is a TFRecord file of Variant
protocol buffers. If running with multiple processes, the sharding applied to
this output filename must be the same as that applied to the `--examples`
output.

A concrete example call, using variables defined in the [WGS case study]:

```bash
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"

( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    python "${BIN_DIR}"/make_examples.zip \
    --mode calling \
    --ref "${REF}" \
    --reads "${BAM}" \
    --examples "${EXAMPLES}" \
    --gvcf "${GVCF_TFRECORDS}" \
    --task {}
) >"${LOG_DIR}/make_examples.log" 2>&1`
```

NOTE: gVCF outputs are only valid when `make_examples` is run in "calling" mode;
if attempted to run in "training" mode the program will exit and notify the user
of the error.

### `postprocess_variants`

When run in gVCF mode, the `postprocess_variants` program handles the creation
of the final gVCF file that incorporates both the non-variant records and the
true variants discovered by the previous programs.

Two additional flags are required in `postprocess_variants`, the input
`--nonvariant_site_tfrecord_path <filename>` which corresponds to the TFRecord
of Variant protocol buffers created in `make_examples`, and the output
`--gvcf_outfile <filename>` which is the final gVCF output.

A concrete example call, using variables defined in the [WGS case study] and in
the above `make_examples` example:

```bash
OUTPUT_GVCF="${OUTPUT_DIR}/HG002.output.g.vcf.gz"

( time python "${BIN_DIR}"/postprocess_variants.zip \
    --ref "${REF}" \
    --infile "${CALL_VARIANTS_OUTPUT}" \
    --outfile "${OUTPUT_VCF}" \
    --nonvariant_site_tfrecord_path "${GVCF_TFRECORDS}" \
    --gvcf_outfile "${OUTPUT_GVCF}"
) >"${LOG_DIR}/postprocess_variants.log" 2>&1
```

## Storage and runtime considerations

The number of non-variant records created when running DeepVariant in gVCF
depends highly on the sequencing depth of the input sample. This is because the
gVCF records at adjacent sites are merged when the genotype qualities are equal,
and we limit the possible genotype quality seen to be at most 50. For
deeply-sequenced individuals (e.g. 30-50x coverage), many sites hit the GQ=50
cap and are merged into few records. Samples with lower sequencing depth have
more sites within the dynamic range of the binomial model used to estimate
non-variant site genotype quality, and thus more records are created.

To mitigate this effect, the `make_examples` program has a flag
`--gvcf_gq_binsize <int>`. This flag allows the merging of adjacent records that
all have GQ values within a bin of the given size, and for each record emits the
minimum GQ value seen within the bin. PL value is taken from the same record
containing the minimum GQ value.

For example, setting `--gvcf_gq_binsize 5` has the effect that adjacent records
with GQ=0; GQ in [1, 5]; GQ in [6, 10]; GQ in [11, 15]; etc. are binned
together.

A concrete example shown below has non-variant sites at each of positions 1-9 on
a hypothetical chromosome:

```bash
Example input records:
Genome position   |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |
GQ of position    |  8 | 10 |  9 | 27 | 47 | 50 | 50 | 45 | 33 |
```

They would create five resultant gVCF record values with `--gvcf_gq_binsize 5`,
with relevant values of:

```bash
start | INFO  | GQ
------------------
    1 | END=3 |  8
    4 | END=4 | 27
    5 | END=7 | 47
    8 | END=8 | 45
    9 | END=9 | 33
```

By synthetically downsampling a 50x coverage whole genome and applying different
GQ binning strategies, we see how the size of the resultant data varies as the
two factors change. The below figure shows the size of output (measured as the
number of records generated relative to the baseline of a 50x whole genome with
`--gvcf_gq_binsize 1`) at different coverage levels, for GQ bins of size 1, 3,
5, and 10. The value of each bar is written in blue font above it for clarity.

![gVCF size](images/DeepVariant-gvcf-sizes-figure.png?raw=true "DeepVariant gVCF sizes")

### Runtime

Despite the creation of many additional records, the running time of
`make_examples` increases minimally when gVCF support is enabled. The
single-threaded `postprocess_variants` program is more adversely affected, with
observed runtimes increasing on the [WGS case study] from ~25 minutes to 5-7
hours depending on genome coverage.

### New option to include MED_DP

Starting in v1.2.0, we added a flag to enable adding MED_DP (median read
coverage seen in the block) in addition to the default MIN_DP (minimum read
coverage seen in the block).

To test it, you can follow the steps in [Quick Start], and in the step where
you run the one-step script `/opt/deepvariant/bin/run_deepvariant`, add this
flag:

```bash
--make_examples_extra_args="include_med_dp=true"
```

Then, if you look at your output gVCF, you'll see the additional MED_DP
information, like:

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
chr20   10000000        .       T       <*>     0       .       END=10000116    GT:GQ:MIN_DP:MED_DP:PL  0/0:50:45:58:0,135,1349
```

[VCF format]: https://samtools.github.io/hts-specs/VCFv4.3.pdf
[WGS case study]: deepvariant-case-study.md
[Quick Start]: deepvariant-quick-start.md
