# DeepVariant release notes

## Binary Releases

### 0.4.1

This fixes a problem with htslib_gcp_oauth when network access is unavailable.

### 0.4.0

This is the initial open source release of the DeepVariant binaries!

We've backed off on the previous aggressive compiler optimizations. We now build
for the lowest common denominator CPU chipset on Google Compute Platform (which
is Sandy Bridge) to ensure compatibility. This means we build specifically as:

```shell
bazel --batch build -c opt --copt=-msse4.1 --copt=-msse4.2 --copt=-mavx --copt=-O3 :binaries --build_python_zip
```

Additionally, we now use (and make available) an optimized, pre-built pip wheel
of Tensorflow 1.4 (built using the same options shown above). In our tests this
reduces the runtime of call_variants by more than 20%.

We now use a default HTS_BLOCK_SIZE of 128 MB in make_examples. In our tests
this enables near-local-file-system-level performance when accessing BAM files
in Google Cloud buckets.

### 0.3.0

This release adds GPU-acceleration-capable executables, compiled on Ubuntu
16.04. Users can take advantage of GPUs during the `call_variants` step of the
DeepVariant pipeline. GPU support is only available with Ubuntu 16.

Improved input data pipeline for `call_variants`.

We apply an exponential moving average to model parameters during inference.
This helps reduce checkpoint to checkpoint variability in the model.

### 0.2.3

This release builds with Tensorflow outside the hermetic zip files.
It makes use of Abseil, instead of getting those libraries from
Tensorflow.  It builds against the open-source version of Clif, supplied
as precompiled binaries.  Native python logging is used where possible.
The scikit-bio dependency has been removed.  Sundry improvements to the
rest.

Because we are now fetching the tf_nightly whl in run-prereq.sh, instead
of building from source, the Tensorflow runtime that will be installed
isn't compiled with aggressive optimization.  (If this is an issue, tf
could be compiled locally from recent source with custom options and
that whl installed instead.)

### 0.2.2

Improvement on the realigner implementation.

### 0.2.1

This release is built with more aggressive optimization, specifically

```shell
bazel --batch build -c opt --copt=-msse4.1 --copt=-msse4.2 --copt=-mavx --copt=-mavx2 --copt=-mfma --copt=-O3 :binaries --build_python_zip
```

Previous releases used just `-c opt`.

This means DeepVariant will require a relatively modern CPU,
one that supports SSE4.2 and AVX2.


### 0.2.0

This release provides CPU-only executables, compiled on Ubuntu 14.04.

Error handling is currently less graceful than we intend, and is a work in
progress. Some operations "check fail" and print a stack trace or core dump when
they should simply print a message and exit.

Reading from Google Cloud Storage is possible, but is a work in progress and may
have poor performance.

OSS binaries with GPU support are in development, and will be released soon.

We are working on support for Pip.

Many runtime improvements are coming.

More documentation on advanced ways to run DeepVariant are coming.

## Model Releases

### 0.4.0

This is the initial open source release of the DeepVariant model!

The model is identical to the 0.3.0 model.

### 0.3.0

For training we added in copies of the 9 replicates of HG001 each downsampled at
50% coverage. In our tests this additional training data means DeepVariant
can generalize to a wider variety of input sequencing data.

### 0.2.0

The model is trained on 9 replicates of HG001. We use the truth set [v.3.3.2
from Genome in a
Bottle](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/)
for training. The underlying model is Inception V3.
