# DeepVariant

DeepVariant is an analysis pipeline that uses a deep neural network to call
genetic variants from next-generation DNA sequencing data.

# Documentation

<!-- mdlint off(URL_BAD_G3DOC_PATH) -->

*   [DeepVariant release notes](docs/deepvariant-release-notes.md)
*   [Building and testing DeepVariant](docs/deepvariant-build-test.md)
*   [DeepVariant quick start](docs/deepvariant-quick-start.md)
*   [DeepVariant via Docker](docs/deepvariant-docker.md)
*   [DeepVariant whole genome case study](docs/deepvariant-case-study.md)
*   [DeepVariant exome case study](docs/deepvariant-exome-case-study.md)
*   [DeepVariant details](docs/deepvariant-details.md)
*   [DeepVariant model training](docs/deepvariant-model-training.md)

<!-- mdlint on -->

<a name="about"></a>
## About DeepVariant

![DeepVariant workflow](docs/DeepVariant-workflow-figure.png?raw=true "DeepVariant workflow")

For technical details describing how DeepVariant works please see our
[preprint](https://www.biorxiv.org/content/early/2016/12/21/092890). Briefly,
we started with GIAB reference genomes, for which there is high-quality ground
truth available. Using multiple replicates of these genomes, we produced tens of
millions of training examples in the form of multi-channel tensors encoding the
sequencing instrument data, and then trained a TensorFlow-based image
classification model ([inception-v3](https://arxiv.org/abs/1512.00567)) to
assign genotype likelihoods from the experimental data produced by the
instrument. Read additional information on the [Google Research
blog](https://research.googleblog.com/2017/12/deepvariant-highly-accurate-genomes.html).

## Support

The [Genomics team in Google Brain](https://research.google.com/teams/brain/genomics/)
actively supports DeepVariant and are always interested in improving the quality
of DeepVariant. If you run into an issue, we recommend you follow one of two
approaches to getting the issue resolved.

If you have found a bug in DeepVariant - i.e., the code itself needs to be
fixed - please report the problem on our [Issue
tracker](https://github.com/google/deepvariant/issues). Make sure to add enough
detail to your report that we can reproduce the problem and fix it. We encourage
including links to snippets of BAM/VCF/etc. files that provoke the bug, if
possible. Depending on the severity of the issue we may patch DeepVariant
immediately with the fix or roll it into the next release.

If you have general questions about DeepVariant usage, please post your question
to [BioStars](https://www.biostars.org/), adding the tag 'deepvariant'. We
monitor [BioStars posts tagged with
DeepVariant](https://www.biostars.org/t/deepvariant/) and will respond as needed
there.

## Contributing

Interested in contributing? See [CONTRIBUTING](CONTRIBUTING.md).

## License

DeepVariant is licensed under the terms of the [BSD-3-Clause license](LICENSE).

## Acknowledgements

DeepVariant happily makes use of many open source packages.  We'd like to
specifically call out a few key ones:

* [Boost Graph Library](http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html)

* [CLIF](https://github.com/google/clif)

* [GNU Parallel](https://www.gnu.org/software/parallel/)

* [htslib & samtools](http://www.htslib.org/)

* [numpy](http://www.numpy.org/)

* [scipy](https://www.scipy.org/)

* [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

* [TensorFlow and Slim](https://www.tensorflow.org/)

We thank all of the developers and contributors to these packages for their
work.


## Disclaimer

*   This is not an official Google product.
