# Linux x86 reference capture

One-time dev capture of upstream DeepVariant's outputs on a small chr20 fixture,
used as the parity reference for our Apple Silicon native rewrite.

This **uses Docker** (`google/deepvariant:1.10.0`) under qemu emulation on the
Mac. Docker appears only here, in a one-time dev-time pipeline. It is **not**
in the user product, which the user-facing constraints already guarantee.

## What gets captured

For each model variant we want to bench (wgs / wes / pacbio / ont_r104 / pangenome):

```text
tools/reference/output/<variant>/
├── examples_chr20.tfrecord       # output of make_examples
├── call_variants_chr20.tfrecord  # output of call_variants (the parity reference)
├── output.vcf.gz                 # output of postprocess_variants
└── manifest.json                 # exact upstream version, args, sha256s
```

The 1000-example slice in `cache/wgs_chr20_1000.tfrecord` is what `bench.py`
reads as input.

## Pipeline

```sh
# one-time downloads
./fetch_chr20_fixture.sh                    # ~120 MB BAM, ~80 MB ref FASTA

# per-variant
./capture_linux_x86.sh wgs
```

`capture_linux_x86.sh` runs upstream DeepVariant's three stages under qemu
linux/amd64 emulation, capturing the inputs and outputs at each stage
boundary.

## Caveats

- qemu emulation is slow — expect 30-60 minutes per variant on the M4 Max
  (vs 10 minutes natively). One-time cost.
- Storage: ~5 GB total per variant after capture.
- The captured TFRecords are committed under `testdata/reference/` via Git LFS
  (set up with `git lfs track "testdata/reference/**"` before the first commit).
