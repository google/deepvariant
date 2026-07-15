# DeepVariant — native arm64 Apple Silicon port

[![status](https://img.shields.io/badge/status-Phase%204%20PASS-brightgreen)](docs/validation.md)
[![build](https://img.shields.io/badge/build-CMake-blue)](CMakeLists.txt)
[![license](https://img.shields.io/badge/license-BSD--3--Clause-orange)](LICENSE)

A fully native arm64 macOS port of Google's
[DeepVariant 1.10.0](README_UPSTREAM.md) for Apple Silicon. Single
statically-linked Mach-O binary, **no Python interpreter at runtime**,
**no Docker**, **no Rosetta 2**. Inference runs on Apple Metal
Performance Shaders Graph (MPSGraph) in FP32 across all 188
Inception-v3 conv layers; the final dense + softmax falls back to
BNNS-CPU FP32 single-thread for threshold-flip determinism.

> **Status (2026-05-02)**: Phase 4 spec gates met across all four tools
> (WGS, DeepTrio, DeepSomatic, Pangenome) — all at 100 % Docker FILTER
> parity on chr20 fixtures. HG002 whole-genome F1 SNP/INDEL
> **bit-identical to Docker** at 6 decimal places (SNP 0.996440,
> INDEL 0.995766); 99.9935 % PASS-set agreement; 1.84× wall-time vs
> Docker on the same M4 Max. gVCF, alt-aligned pileup, methylation, and
> DirectPhasing flags implemented. Phase 5 packaging + Phase 6 Homebrew
> tap pending.

## Why this port

| Metric | Linux x86 Docker (Rosetta 2) | This port (native arm64) | Speedup |
|--------|------------------------------|--------------------------|---------|
| chr20 wall-time on M4 Max | ~17 min | **6 m 27 s** | **2.6×** |
| HG002 WG (whole genome) on M4 Max | ~6 h | **3 h 16 min** | **1.84×** |
| GPU residency | 0 (CPU-only emulation) | ≥ 40 % during inference | — |
| Python interpreter | required | **none at runtime** | — |
| Docker daemon | required | **none** | — |

Equivalence with upstream `google/deepvariant:1.10.0` Docker is
**clinical-grade** (not bit-exact — fundamentally unachievable on
Apple GPU due to FP32 non-associativity in any parallel reduction).
We define equivalence by four criteria, in order:

1. Site set identical (CHROM/POS/REF/ALT)
2. FILTER class identical (PASS / RefCall / NoCall / LowQual)
3. GT identical
4. PASS variant set identical

See [`docs/scientific_report.md`](docs/scientific_report.md) for the
full mathematical framework, methods, biological-impact analysis of
FILTER mismatches, and rare-variant impact assessment.

## Validation summary — chr20 trio WGS (vs GIAB v4.2.1)

| Sample | Type  | F1      | Δ vs upstream Docker    | Phase 4 gate |
|--------|-------|---------|-------------------------|--------------|
| HG002  | SNP   | 0.99740 | 0.00000 (bit-identical) | **PASS** ✓ |
| HG002  | INDEL | 0.99598 | 0.00000 (bit-identical) | **PASS** ✓ |
| HG003  | SNP   | 0.99777 | within FP-drift residue | **PASS** ✓ |
| HG003  | INDEL | 0.99688 | within FP-drift residue | **PASS** ✓ |
| HG004  | SNP   | 0.99767 | within FP-drift residue | **PASS** ✓ |
| HG004  | INDEL | 0.99636 | within FP-drift residue | **PASS** ✓ |

NovaSeq 35× PCR-free Illumina chr20, evaluated against GIAB v4.2.1
high-confidence regions.

### Docker FILTER parity — all four tools (chr20:10M-10.1M)

| Tool | Shared | FM | PASS identical | Result |
|------|-------:|---:|---------------:|--------|
| WGS (HG002) | 313 | **0** | 261 / 261 | **PASS** ✓ |
| DeepTrio child (HG002) | 262 | **0** | 262 / 262 | **PASS** ✓ |
| DeepTrio parent1 (HG003) | 265 | **0** | 265 / 265 | **PASS** ✓ |
| DeepTrio parent2 (HG004) | 222 | **0** | 222 / 222 | **PASS** ✓ |
| DeepSomatic (HG002 tumor + HG003 normal) | 693 | **0** | 34 PASS + 92 GERMLINE | **PASS** ✓ |
| Pangenome-aware WGS (HG002 + GBZ BAM) | 322 | **0** | 247 / 247 | **PASS** ✓ |

### HG002 whole-genome WGS (vs GIAB v4.2.1)

| Type  | F1      | Δ vs Docker           | PASS-set Δ             | GT-disagree PASS-PASS |
|-------|---------|-----------------------|------------------------|-----------------------|
| SNP   | 0.99644 | **0** (bit-identical) | 317 / 4.84 M (0.007%) | 136 / 4.84 M (0.003%) |
| INDEL | 0.99577 | **0** (bit-identical) | —                      | —                     |

Wall-time: 3 h 16 min native vs 5 h 59 min Docker → **1.84× faster** on
the same M4 Max machine with identical inputs and `--num_shards=14`.

Full benchmark: [`validation/output/HG002_wg_benchmark.md`](validation/output/HG002_wg_benchmark.md)

## Quick start

### Build

```bash
git clone https://github.com/IPNP-BIPN/deepvariant && cd deepvariant
git checkout feature/apple-silicon-native-v2
./scripts/build-prereq-macos.sh                 # Homebrew deps
cmake -S . -B build-macos -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-macos --target deepvariant
```

Build prerequisites:

- macOS ≥ 14, Apple Silicon (M1/M2/M3/M4)
- Apple Xcode Command Line Tools
- Homebrew: `cmake`, `ninja`, `htslib`, `abseil`, `protobuf`,
  `samtools`, `bcftools`, `tabix`, `bgzip`

### Run a chr20 trio benchmark

```bash
# Pre-extracted chr20 fixture (HG002/HG003/HG004 NovaSeq 35× BAMs +
# GIAB v4.2.1 truth + GRCh38 no_alt chr20 reference, ~3 GB)
./tools/reference/fetch_chr20_fixture.sh

# Run the full pipeline + hap.py per sample (~30 min total)
./validation/run_giab_chr20_trio.sh

# Inspect results
column -t -s, validation/output/HG00*_chr20/happy.summary.csv | less -S
cat validation/output/chr20_trio_summary.tsv
```

### Run a whole-genome benchmark (~10 h, ~120 GB download)

```bash
./validation/download_giab_full_genome.sh       # one-time, ~120 GB
./validation/run_giab_wg_chunked.sh HG002       # ~3 h 16 min on M4 Max
```

The chunked runner processes one chromosome at a time, freeing
intermediate files between chunks to stay within ~90 GB peak disk.
See [`docs/wg_benchmark_audit.md`](docs/wg_benchmark_audit.md).

### Run a one-shot WGS variant call

```bash
./build-macos/bin/deepvariant run \
  --reads=/path/to/sample.bam \
  --ref=/path/to/GRCh38.fa \
  --regions=chr20 \
  --output_vcf=/tmp/out.vcf.gz \
  --inference_backend=metal \
  --model_type=WGS \
  --checkpoint=validation/work/wgs.dvw \
  --small_model_path=validation/work/wgs_small_weights \
  --num_shards=14
```

### Run DeepTrio (child + 2 parents)

```bash
./build-macos/bin/deepvariant trio \
  --reads_child=HG002.bam \
  --reads_parent1=HG003.bam \
  --reads_parent2=HG004.bam \
  --ref=GRCh38.fa --regions=chr20 \
  --output_vcf_child=/tmp/child.vcf.gz \
  --output_vcf_parent1=/tmp/parent1.vcf.gz \
  --output_vcf_parent2=/tmp/parent2.vcf.gz \
  --checkpoint_child=validation/work/deeptrio.wgs_child.dvw \
  --checkpoint_parent=validation/work/deeptrio.wgs_parent.dvw \
  --num_shards=14
```

### Run DeepSomatic (tumor + normal)

```bash
./build-macos/bin/deepvariant somatic \
  --reads_tumor=tumor.bam \
  --reads_normal=normal.bam \
  --ref=GRCh38.fa --regions=chr20 \
  --output_vcf=/tmp/somatic.vcf.gz \
  --checkpoint=validation/work/deepsomatic.wgs.dvw \
  --num_shards=14
```

### Run with gVCF output

```bash
./build-macos/bin/deepvariant run \
  --reads=sample.bam --ref=GRCh38.fa \
  --output_vcf=/tmp/out.vcf.gz \
  --output_gvcf=/tmp/out.g.vcf.gz \
  --checkpoint=validation/work/wgs.dvw \
  --num_shards=14
```

Subcommands: `run`, `make_examples`, `call_variants`,
`postprocess_variants`, `trio`, `somatic`, `pangenome`.

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│ deepvariant run (single-binary, native arm64)               │
├──────────────────┬──────────────────────┬───────────────────┤
│  make_examples   │   call_variants      │ postprocess_      │
│  (CPU, N threads)│   (GPU + BNNS-CPU)   │  variants (CPU)   │
├──────────────────┼──────────────────────┼───────────────────┤
│ - SamReader      │ - MPSGraph FP32      │ - CombineLikeli- │
│   (htslib mmap)  │   (Inception-v3,     │   hoods           │
│ - AlleleCounter  │    188 conv layers)  │ - simplify_       │
│ - DBG realigner  │ - BNNS-CPU FP32      │   alleles         │
│ - PileupImage    │   single-thread      │ - haplotype res   │
│ - NEON encoding  │   (2048→3 dense +    │   (Boost-graph)   │
│ - libstdc++      │    softmax)          │ - VCF + gVCF      │
│   shuffle        │                      │   emission        │
│ - NumPy MT19937  │ Optional backend:    │ - DirectPhasing   │
│   reservoir      │ - ANE FP16 first     │                   │
│   sampling       │   + GPU FP32 rerun   │                   │
│                  │   (ane_speculate)    │                   │
└──────────────────┴──────────────────────┴───────────────────┘
       ↓ examples.tfrecord    ↓ cvo.tfrecord    ↓ output.vcf.gz
                                                  output.g.vcf.gz
```

Seven Phase 5.5d root-cause fixes close 1.13 % FILTER drift (pre-fix)
to 0 FM (post-fix on HG002 chr20 full):

1. libstdc++-compatible `std::shuffle` (vs libc++ default)
2. NumPy MT19937 + Algorithm-R reservoir sampling
3. Multi-allelic CombineLikelihoods CVO-prune
4. Haplotype resolution port (Boost-graph max-weight)
5. `simplify_variant_alleles` postfix strip
6. BNNS-CPU FP32 small-model + AltAlleleQual rounding
7. PL log-space subtract + truncation

## Supported models

| Model type | `--model_type` | Pileup shape | Tool |
|------------|----------------|--------------|------|
| WGS Illumina | `WGS` | 100×221×7 | `run` |
| WES Illumina | `WES` | 100×221×7 | `run` |
| PacBio HiFi | `PACBIO` | 100×147×10 | `run` |
| Oxford Nanopore | `ONT` | 100×199×10 | `run` |
| Hybrid PacBio+Illumina | `HYBRID_PACBIO_ILLUMINA` | 100×221×6 | `run` |
| MaSeq | `MASSEQ` | 100×199×9 | `run` |
| RNA-seq | `RNASEQ` | 100×221×6 | `run` |
| DeepTrio WGS | `WGS` | 140×221×7 | `trio` |
| DeepTrio WES | `WES` | 140×221×7 | `trio` |
| DeepSomatic WGS | `WGS` | 200×221×7 | `somatic` |
| DeepSomatic WES | `WES` | 200×221×7 | `somatic` |
| DeepSomatic PacBio | `PACBIO` | 200×147×9 | `somatic` |
| DeepSomatic ONT | `ONT` | 200×99×9 | `somatic` |
| DeepSomatic FFPE WGS | `WGS --ffpe` | 200×221×7 | `somatic` |
| Pangenome-aware WGS | — | 200×221×7 | `pangenome` |

## Inference backends

| Backend | Flag | Speed | Docker FILTER parity |
|---------|------|-------|----------------------|
| `metal` (default) | `--inference_backend=metal` | 1.84× vs Docker | 100 % (gate met) |
| `ane_speculate` | `--inference_backend=ane_speculate` | ~2.5–3× vs Docker | empirical (in progress) |
| `coreml` | `--inference_backend=coreml` | debug only | — |

## Documentation

| Document | Audience |
|----------|----------|
| [`docs/scientific_report.md`](docs/scientific_report.md) | Publication-grade report: math, methods, results, FM analysis, rare-variant impact |
| [`docs/validation.md`](docs/validation.md) | Methods + all-mode F1 results + WG benchmark |
| [`docs/wg_benchmark_audit.md`](docs/wg_benchmark_audit.md) | Whole-genome benchmark: measured results, disk budget |
| [`CLAUDE.md`](CLAUDE.md) | Project memory: phase status, root-cause fix log, constraints |
| [`README_UPSTREAM.md`](README_UPSTREAM.md) | Original Google DeepVariant 1.10.0 README (attribution) |

## Test fixtures + reference data

- `validation/work/wgs.dvw` — WGS weights (387 tensors, ~83 MB)
  SHA-256: `57fcefeaf230e7a795bb1fdbc275e5f02039f010de2ebcf8a9fde0cb9f006479`
- `validation/work/wgs_small_weights/` — WGS BNNS-CPU small-model weights
- `validation/work/deeptrio.wgs_{child,parent}.dvw` — DeepTrio WGS weights
- `validation/work/deepsomatic.wgs.dvw` — DeepSomatic WGS weights
- `validation/work/pangenome.wgs.dvw` — Pangenome-aware WGS weights
- `validation/output/chr20_trio_summary.tsv` — chr20 trio F1 numbers
- `testdata/reference/per_layer/*.npy` — per-tap TF reference outputs (Git LFS)

## Performance

Measured on Apple M4 Max (16 cores, 128 GB unified memory,
macOS 26.4.1) with `--num_shards=14`:

| Stage | chr20 wall-time |
|-------|-----------------|
| make_examples | ~3 min (210 390 candidates, 14 threads) |
| call_variants | ~30 s (441 batches × MPSGraph FP32) |
| postprocess_variants | ~5 s (haplotype res + VCF emit) |
| **Total `deepvariant run`** | **~3 min** |
| **HG002 whole genome** | **3 h 16 min** (vs Docker 5 h 59 min → **1.84×**) |

GPU residency confirmed via `powermetrics --samplers gpu_power -i 500`
(GPU ≥ 40 % active during inference).

## Repository layout

```
deepvariant/
├── deepvariant/              # upstream C++ sources (BSD-3, Google)
│   └── native/               # this port (BSD-3, Demaille)
│       ├── make_examples_main.cc      # Stage 1 orchestrator
│       ├── call_variants_main.cc      # Stage 2 (Metal + BNNS)
│       ├── postprocess_main.cc        # Stage 3 + gVCF merge
│       ├── cli.cc                     # `deepvariant run` dispatcher
│       ├── metal_inference.{h,mm}     # MPSGraph Inception-v3 build
│       ├── bnns_finalize.{h,mm}       # BNNS-CPU FP32 final dense
│       ├── neon_base_color.h          # NEON pileup encoding (A2.1)
│       ├── neon_cigar_classify.h      # NEON CIGAR walk (A2.2)
│       ├── numpy_mt19937.h            # NumPy-compat reservoir sampling
│       ├── libstdcxx_shuffle.h        # libstdc++-compat std::shuffle
│       ├── haplotypes.{h,cc}          # haplotype resolution port
│       └── gvcf_emit.{h,cc}           # gVCF block emitter
├── third_party/nucleus/      # nucleus io (sam/vcf/fasta) — upstream
├── docs/                     # validation + scientific report
├── validation/               # benchmark scripts + reference outputs
├── tools/conversion/         # weight extraction + per-layer dumps
├── tools/reference/          # Docker reference capture scripts
├── release/                  # codesign + notarize scripts (Phase 5)
├── scripts/build-prereq-macos.sh
├── CMakeLists.txt
├── CLAUDE.md                 # project memory
└── README.md                 # this file
```

## Hard constraints (from project plan)

- macOS ≥ 14, arm64 only
- No Docker / Rosetta 2 / CUDA at runtime
- No Python interpreter in the runtime artefact
- SNP F1 ≥ upstream − 0.05 %, INDEL F1 ≥ upstream − 0.10 %
- GPU residency verified via `powermetrics`
- Speedup ≥ 2.5× vs published Linux x86 reference
- 100 % FILTER-class parity vs `google/deepvariant:1.10.0` Docker
  on chr20 full — **met for WGS, DeepTrio, DeepSomatic, Pangenome**

## Reproducibility

Each `deepvariant run` invocation is deterministic on the same
hardware (verified by repeated runs producing byte-identical CVOs).

Cross-chip determinism (M1 vs M2 vs M3 vs M4) preserves FILTER class
by construction (FP32 cumulative drift bounded by the threshold-flip
sensitivity analysis in
[`docs/scientific_report.md`](docs/scientific_report.md) §2.4).

Build provenance:

| Component | Version |
|-----------|---------|
| Apple clang | 21.0.0 (`clang-2100.0.123.102`) |
| CMake | 4.3.2 |
| macOS | 26.4.1 (build 25E253) |
| Docker (validation only) | 29.2.1 (Docker Desktop 4.63.0) |
| `jmcdani20/hap.py` | v0.3.12 |

## Citation

If you use this port in academic work, please cite both:

1. The original DeepVariant paper:

   Poplin R., Chang P-C., Alexander D., et al. (2018).
   *A universal SNP and small-indel variant caller using deep neural
   networks*. Nature Biotechnology **36**, 983-987.

2. This port (preprint forthcoming on bioRxiv):

   Demaille B. (2026). *Native Apple Silicon port of DeepVariant
   1.10.0: characterising FP32 non-associativity in clinical-grade
   variant calling on heterogeneous hardware*. (Preprint URL TBD)

## License + attribution

This port is BSD-3-Clause licensed (see [`LICENSE`](LICENSE)).
Original DeepVariant code is © 2020 Google LLC, BSD-3-Clause.
Pre-trained model weights distributed by Google at
`gs://deepvariant/models/DeepVariant/1.10.0/` are used under the
same BSD-3-Clause license.

This is a derivative work. Google is not affiliated with this port
and provides no endorsement of it. The "DeepVariant" name is a
Google trademark used here for nominative reference to the
underlying open-source project.

## Contact

Benjamin Demaille — benjamin.demaille@icloud.com

Repository: [IPNP-BIPN/deepvariant](https://github.com/IPNP-BIPN/deepvariant)

## Acknowledgements

- Google DeepVariant team for the original method, codebase,
  pre-trained models, and the public Linux x86 Docker reference
  used as our parity baseline.
- NIST Genome in a Bottle (GIAB) consortium for the v4.2.1 truth
  sets used in F1 evaluation.
- Apple for the Metal Performance Shaders Graph framework + BNNS
  Accelerate library.
- htslib / nucleus / abseil / protobuf maintainers for the
  underlying open-source dependencies.
