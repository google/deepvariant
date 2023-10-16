# DeepVariant FAQ

## How does DeepVariant work?

See this
[overview](https://github.com/google/deepvariant#how-deepvariant-works).

## Why does DeepVariant not call a specific variant in my data?

**Missing variants due to a candidate not being generated:**

There are multiple reasons that DeepVariant may not call a variant. It is
important to first determine whether a candidate variant was proposed by
DeepVariant. A potential variant requires at least 2 reads to support a variant
and a minimum fraction of reads supporting the variant (0.12 for SNPs and PacBio
Indels, 0.06 for Illumina Indels). All sites that have been generated as
candidates are written in the VCF file, so if you do not see a row in the VCF
file for the variant in question, it means that a candidate was not made.
However, within these sites certain possible alleles may have been pruned in
reporting. To see all alleles, you may add: `--debug_output_all_candidates=ALT`
in the postprocess_variants step.

To increase the sensitivity of DeepVariant to these sites, you may add the
following parameters, here shown with their defaults:

```
--make_examples_extra_args="vsc_min_count_snps=2,vsc_min_fraction_snps=0.12,vsc_min_count_indels=2,vsc_min_fraction_indels=0.06"
```

It is sometimes also the case that realignment of the reads within DeepVariant
changes or reduces the evidence supporting the variant. To check for this, try
using the `--norealign_reads` flag to turn off realignment temporarily. Note
that we don't recommend turning off the realigner for Illumina data in general
cases because the realigner improves accuracy overall.

There is also the option to output the realigned reads, e.g. to inspect the new
alignments in IGV. See the "What is the realigner and how does it work?" section
for instructions.

**Missing variants where a candidate is generated:**

If a candidate is made, but is called as reference (either 0/0 or ./.) it means
that the neural network processed the genomic region, but based on all of its
learned experience from training data, it decided the highest probability for
the position was as non-variant. Some of the reasons that DeepVariant may
suspect a false positive are: strand-bias in reads, low mapping quality in
reads, low base quality in reads, and overall low coverage.

In addition, there is another pattern that causes DeepVariant to suspect variant
positions which can initially seem counterintuitive to human observers. This
occurs when a dense set of variants appears on one haplotype while the other
haplotype is fully reference, and humans often perceive this as missing a
clearly heterozygous position. DeepVariant seems to have learned that this
signature often indicates a region which is a segmental duplication, copy number
variant, or structural variant where multiple copies of similar genomic regions
are mapping to the same reference location. In this case, it may be worthwhile
to inspect the region to see if it has elevated coverage, and whether you can
identify more than 2 haplotypes present by overlapping the reads. If you can, it
suggests that the region may have a copy number variation. Some analysis of this
was presented at AGBT as a poster
“[Uncaptured segmental duplication creates artifacts in workflows using GRCh37](https://pbs.twimg.com/media/ERe2bSyWsAcE00h?format=jpg&name=4096x4096)”.

This pattern of undercalling positions at high variant density may affect
variant-dense non-human species (those with a variant density of >1 in 40
positions). For an analysis of this, please see our blog
“[Improved non-human variant calling using species-specific DeepVariant models](https://google.github.io/deepvariant/posts/2018-12-05-improved-non-human-variant-calling-using-species-specific-deepvariant-models/)”.

If these reasons seem applicable, there could be some other reason DeepVariant
determined the position is not variant. You can catalog the variant position and
its support. The way to improve variant calling for these positions is to train
new models, but be aware that training is already a balance between reducing
false negatives and positives, and it may not be possible to call variants like
the one you are seeing without increasing overall false positives by a greater
amount.

## How does DeepVariant use pileup images to call variants?

See this
[blog post](https://google.github.io/deepvariant/posts/2020-02-20-looking-through-deepvariants-eyes/).

## What happens if I change the pileup_image_height?

If the actual depth in a particular region is greater than the pileup image
height, DeepVariant randomly downsamples reads until the image has been filled
up. For the default DeepVariant models (height 100), an image can accommodate at
most 95 reads in a given region (5 rows are reserved for the reference
sequence).

You may be able to successfully run our pretrained models with a different
pileup image height (via `--pileup_image_height` in `make_examples.py`),
depending on the new height, but we generally do not recommend using different
image heights at training and inference time. If you wish to use a different
pileup image height, we recommend retraining a new model with images of that
height.

If you are working with extremely high coverage sequencing data for applications
such as somatic sequencing, we recommend using a somatic caller instead of
DeepVariant, which is a germline caller.

## Can I use DeepVariant for somatic (non-germline) calling?

We do not recommend using DeepVariant for somatic calling. We do have a
prototype implementation for somatic calling, which can take a tumor and normal
BAM and call subclonal variants. However, we don't yet have enough confidence in
the available truth sets, and that they come from a diverse enough sampling of
cancers with mutational profiles, for us to be certain in releasing something of
high quality. We're watching developments in the area of these truth sets and
hope to be able to further develop the somatic caller in the future.

## Can I use DeepVariant on plant genomes?

DeepVariant has previously been applied to plant species. In the case of rice,
there was good evidence of high accuracy. You can see
[some results in this blog post](https://cloud.google.com/blog/products/data-analytics/analyzing-3024-rice-genomes-characterized-by-deepvariant).
However, these rice genomes were diploid and with a similar variant density of
humans.

DeepVariant is currently written to be a diploid variant caller. So if the plant
species you are working with is polyploid, it is not yet clear how DeepVariant
will perform. That is because even with re-training, DeepVariant can only
produce variant calls that are homozygous alternate, heterozygous, or homozygous
reference, which don't have much meaning in a tetraploid genome, for example.

## Can I use DeepVariant on other non-human species?

See this
[blog post](https://google.github.io/deepvariant/posts/2018-12-05-improved-non-human-variant-calling-using-species-specific-deepvariant-models/).

## How do I build/run DeepVariant?

In general, we recommend running DeepVariant using Docker for the simplest
setup. If you are building from source because you want to experiment with
changes to the codebase, we still recommend Docker. You can clone the
DeepVariant repo, modify the source code, and build a Docker image with your
changes using the provided Dockerfile.

## Why can't it find one of the input files? E.g., "Could not open"

This often happens because the way Docker works, input and output directories
have to be mounted and then files are referred to by their mounted location,
which can be confusing. To check that files are visible inside the Docker
container, you can `ls` inside the container. For example, using the setup shown
in the README and looking inside the `/input` volume:

```
BIN_VERSION="1.6.0"
docker run \
  -v "YOUR_INPUT_DIR":"/input" \
  -v "YOUR_OUTPUT_DIR:/output" \
  google/deepvariant:"${BIN_VERSION}" \
  ls /input
```

Mounting directories with Docker can be confusing. One trick to make this
simpler is to set both sides as your `$HOME`, so the paths are the same inside
and outside the Docker container.

```
echo $HOME # see what your home directory is first.
ls $HOME
BIN_VERSION="1.6.0"
sudo docker run \
  -v "${HOME}":"${HOME}" \
  google/deepvariant:"${BIN_VERSION}" \
  ls $HOME
```

## How do I run multi-sample calling?

Since the DeepVariant v0.9 release, we recommend
"[Best practices for multi-sample variant calling with DeepVariant](https://github.com/google/deepvariant/blob/r0.9/docs/trio-merge-case-study.md)".

For specifically calling on duos or trios, we introduced
[DeepTrio](https://github.com/google/deepvariant/blob/r1.6/docs/deeptrio-details.md)
in v1.1.

## Can call_variants be run on multiple GPUs?

No. TensorFlow's Estimator API does not support running predictions on multiple
GPUs, so call_variants can only use 1 GPU at prediction time.

## Can model_train be run on multiple GPUs?

No. TensorFlow's Estimator API does provide support for running training on
multiple GPUs through the use of a DistributionStrategy. However,
DistributionStrategy cannot be used with exponential moving average (EMA), which
is present in the DeepVariant codebase.

## What is the realigner and how does it work?

From the
[DeepVariant 2018 manuscript](https://www.nature.com/articles/nbt.4235.epdf?author_access_token=q4ZmzqvvcGBqTuKyKgYrQ9RgN0jAjWel9jnR3ZoTv0NuM3saQzpZk8yexjfPUhdFj4zyaA4Yvq0LWBoCYQ4B9vqPuv8e2HHy4vShDgEs8YxI_hLs9ov6Y1f_4fyS7kGZ):

> Mapped reads are preprocessed using an error-tolerant, local
> De-Bruijn-graph-based read assembly procedure that realigns them according to
> their most likely derived haplotype. Candidate windows across the genome are
> selected for reassembly by looking for any evidence of possible genetic
> variation, such as mismatching or soft clipped bases. The selection criteria
> for a candidate window are very permissive so that true variation is unlikely
> to be missed. All candidate windows across the genome are considered
> independently. De Bruijn graphs are constructed using multiple fixed k-mer
> sizes (from 20 to 75, inclusive, with increments of 5) out of the reference
> genome bases for the candidate window, as well as all overlapping reads. Edges
> are given a weight determined by how many times they are observed in the
> reads. We trim any edges with weight less than three, except that edges found
> in the reference are never trimmed. Candidate haplotypes are generated by
> traversing the assembly graphs and the top two most likely haplotypes are
> selected that best explain the read evidence. The likelihood function used to
> score haplotypes is a traditional pair HMM with fixed parameters that do not
> depend on base quality scores. This likelihood function assumes that each read
> is independent. Finally, each read is then realigned to its most likely
> haplotype. This procedure updates both the position and the CIGAR string for
> each read.

Local realignment is not performed for long reads (PacBio, and other similar
technologies). The realigner step can optionally be switched off using
`--norealign_reads`.

There is also the option to output the realigned reads, e.g. to inspect the new
alignments in IGV. This can be done by passing the following parameters:
`--make_examples_extra_args="emit_realigned_reads=true,realigner_diagnostics=/output/realigned_reads"`

Note that this is meant for debugging and produces a bam file for every
candidate variant, which can result in millions of tiny bam files, so when using
this, narrow down the DeepVariant run using `--regions` to just the variants you
want to inspect more closely.

## How are `AD` and `DP` values calculated?

In order to efficiently perform variant calling, DeepVariant partitions the
genome into chunks (set by `--partition_size`), and will read in a max number of
reads into each partition (set by `--max_reads_per_partition`). By default,
`--partition_size` is set to 1000 and `--max_reads_per_partition` is set to
1500. The `AD` and `DP` values are based on the read depths constrained by
`--max_reads_per_partition`.

For example, if you have a depth of 2000x at a given site, DeepVariant will
subsample 1500 reads, and `DP` or `AD` will be capped at 1500. If you want to
calculate the true `AD` and `DP` values at high-depth regions, you can set
`--max_reads_per_partition=0` to calculate `AD` and `DP` using all reads. In
practice, capping reads per partition reduces runtimes with little/no impact on
accuracy.

## Missing variant calls near the edge of a contig

This is a known issue that we don't currently address. Please see:
https://github.com/google/deepvariant/issues/505 for more context.

## Singularity related questions:

### `TMPDIR`

If you have issues with `TMPDIR` when running with Singularity, try adding this
to your command:

```bash
export TMPDIR="$PWD/tmp_dir"
```

See https://github.com/google/deepvariant/issues/524#issuecomment-1067597987.

### Issues with `/mnt/`

User reported that sometimes their setup uses `/mnt/`, which exists in our
Docker image, and it has caused an issue in Singularity.

You can use `-B` in Singularity to avoid this issue. See:
https://github.com/google/deepvariant/issues/530#issuecomment-1076923302 for
more details.

