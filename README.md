# DeepVariant

DeepVariant is an analysis pipeline that uses a deep neural network to call
genetic variants from next-generation DNA sequencing data.

## Availability

<!-- mdlint off(URL_BAD_G3DOC_PATH) -->

DeepVariant is a suite of Python/C++ programs that run on any Unix-like
operating system. For convenience the documentation refers to building and
running DeepVariant on [Google Cloud Platform](https://cloud.google.com/), but
the tools themselves can be built and run on any standard Linux computer,
including on-premise machines. Note that DeepVariant currently requires
Python 2.7 and does not yet work with Python 3.

Pre-built binaries are available at
[gs://deepvariant/](https://console.cloud.google.com/storage/browser/deepvariant).
These are compiled to use SSE4 and AVX instructions, so you'll need a CPU (such
as Intel Sandy Bridge) that supports them. (The file /proc/cpuinfo lists these
features under "flags".)

Alternatively, see [Building and testing
DeepVariant](docs/deepvariant-build-test.md) for more information on building
DeepVariant from sources for your platform.

For managed pipeline execution of DeepVariant see the [cost- and
speed-optimized, Docker-based
pipelines](https://cloud.google.com/genomics/deepvariant) created for Google
Cloud Platform.

## Documentation

*   [DeepVariant release notes](https://github.com/google/deepvariant/releases)
*   [Building and testing DeepVariant](docs/deepvariant-build-test.md)
*   [DeepVariant quick start](docs/deepvariant-quick-start.md)
*   [DeepVariant via Docker](docs/deepvariant-docker.md)
*   [DeepVariant whole genome case study](docs/deepvariant-case-study.md)
*   [DeepVariant exome case study](docs/deepvariant-exome-case-study.md)
*   [DeepVariant Genomic VCF (gVCF) support](docs/deepvariant-gvcf-support.md)
*   [DeepVariant usage guide](docs/deepvariant-details.md)
*   [DeepVariant model training](docs/deepvariant-model-training.md)
*   [Datalab example: visualizing pileup
    images/tensors](docs/visualizing_examples.ipynb)
*   [Getting Started with GCP](deepvariant-gcp-info.md) (It is not required to
    run DeepVariant on GCP.)

## Other resources

*   [Google Developer Codelab: Variant Calling on a Rice genome with
    DeepVariant](https://codelabs.developers.google.com/codelabs/genomics-deepvariant)
*   [Improve DeepVariant for BGISEQ germline variant
    calling](http://bit.ly/train-deepvariant) |
    [slides](https://github.com/SVAI/RecausalNucleotideNetworks/blob/master/ReCausalNucleotideNetwork.pdf)

<!-- mdlint on -->

<a name="about"></a>
## About DeepVariant

For technical details describing how DeepVariant works please see our
[preprint](https://doi.org/10.1101/092890).

![DeepVariant workflow](docs/DeepVariant-workflow-figure.png?raw=true "DeepVariant workflow")

Briefly, we started with some of the reference genomes from [Genome in a
Bottle](http://jimb.stanford.edu/giab/), for which there is high-quality ground
truth available (or the closest approximation currently possible). Using
multiple replicates of these genomes, we produced approximately one hundred
million training examples in the form of multi-channel tensors encoding the
sequencing instrument data, and then trained a TensorFlow-based image
classification model ([inception-v3](https://arxiv.org/abs/1512.00567)) to
assign genotype likelihoods from the experimental data produced by the
instrument. Read additional information on the [Google Research
blog](https://research.googleblog.com/2017/12/deepvariant-highly-accurate-genomes.html).

Under the hood, DeepVariant relies on
[Nucleus](https://github.com/google/nucleus), a library of Python and C++ code
for reading and writing data in common genomics file formats (like SAM and VCF)
designed for painless integration with the
[TensorFlow](https://www.tensorflow.org/) machine learning framework.

## Evaluating DeepVariant

We are delighted to see several external evaluations of the DeepVariant method.

The 2016 PrecisionFDA Truth Challenge, administered by the FDA, assessed several
community-submitted variant callsets on the (at the time) blinded evaluation
sample, HG002. DeepVariant won the [Highest SNP
Performance](https://precision.fda.gov/challenges/truth/results) award in the
challenge.

DNAnexus [posted an extensive
evaluation](https://blog.dnanexus.com/2017-12-05-evaluating-deepvariant-googles-machine-learning-variant-caller/)
of several variant calling methods, including DeepVariant, using a variety of
read sets from HG001, HG002, and HG005. They have also evaluated DeepVariant
under a variety of [noisy sequencing
conditions](https://blog.dnanexus.com/2018-01-16-evaluating-the-performance-of-ngs-pipelines-on-noisy-wgs-data/).

Independent evaluations of DeepVariant v0.6 from both
[DNAnexus](https://blog.dnanexus.com/2018-04-18-deepvariant-amplified/) and
[bcbio](https://github.com/bcbio/bcbio_validations/tree/master/deepvariant#deepvariant-v06-release-strelka2-stratification-and-initial-gatk-cnn)
are also available. Their analyses support our findings of improved indel
accuracy, and also include comparisons to other variant calling tools.

## Support

The [Genomics team in Google Brain](https://research.google.com/teams/brain/genomics/)
actively supports DeepVariant and are always interested in improving the quality
of DeepVariant. If you run into an issue, please report the problem on our [Issue
tracker](https://github.com/google/deepvariant/issues). Make sure to add enough
detail to your report that we can reproduce the problem and fix it. We encourage
including links to snippets of BAM/VCF/etc. files that provoke the bug, if
possible. Depending on the severity of the issue we may patch DeepVariant
immediately with the fix or roll it into the next release.

If you have questions about next-generation sequencing, bioinformatics, or other
general topics not specific to DeepVariant we recommend you post your question
to a community discussion forum such as [BioStars](https://www.biostars.org/).

## Contributing

Interested in contributing? See [CONTRIBUTING](CONTRIBUTING.md).

## License

DeepVariant is licensed under the terms of the [BSD-3-Clause license](LICENSE).

## Acknowledgements

DeepVariant happily makes use of many open source packages.  We'd like to
specifically call out a few key ones:

*   [Boost Graph
    Library](http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html)

*   [abseil-cpp](https://github.com/abseil/abseil-cpp) and
    [abseil-py](https://github.com/abseil/abseil-py)

*   [CLIF](https://github.com/google/clif)

*   [GNU Parallel](https://www.gnu.org/software/parallel/)

*   [htslib & samtools](http://www.htslib.org/)

*   [Nucleus](https://github.com/google/nucleus)

*   [numpy](http://www.numpy.org/)

*   [scipy](https://www.scipy.org/)

*   [SSW
    Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

*   [TensorFlow and Slim](https://www.tensorflow.org/)

We thank all of the developers and contributors to these packages for their
work.


## Disclaimer

*   This is not an official Google product.
