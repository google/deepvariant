# Whole-Genome GIAB Benchmark — Audit (publication-ready)

Audit conducted: 2026-05-01
Branch / commit: `feature/apple-silicon-native-v2` @ `a3d7247b`
Hardware: M4 Max, 14 cores, 64 GB unified memory
Backend: Apple Metal MPSGraph FP32 + BNNS-CPU finalize

## Goal

Produce publication-ready F1 numbers (SNP, INDEL; recall, precision,
F1, with stratification) for our native arm64 DeepVariant on the GIAB
HG002/HG003/HG004 trio, whole-genome (~3.1 Gb), against GIAB v4.2.1
truth sets, and benchmark them against Google's published v1.10.0
Linux x86 numbers.

## Spec gates (from master plan)

| Gate | Threshold | Per |
|------|-----------|-----|
| SNP F1 | ≥ Linux x86 F1 − 0.05 % | sample |
| INDEL F1 | ≥ Linux x86 F1 − 0.10 % | sample |
| 100 % Docker FILTER parity | already met chr20 | — |

## Current state of evidence

### Already established (chr20 only)

| Sample | Region | SNP F1 | INDEL F1 | Source |
|--------|--------|--------|----------|--------|
| HG002 | chr20 | **0.997402** | **0.995942** | `validation/output/HG002_chr20_full/happy.summary.csv` |

This number is **bit-identical to Google's `google/deepvariant:1.10.0`
Docker baseline** (Phase 5.5d/10, 2026-04-29). Every digit matches.

### Missing for publication

- HG003 chr20 + WG F1
- HG004 chr20 + WG F1
- HG002 WG F1
- HG002/3/4 stratified F1 (lowcomplexity, segdup, MHC, GC bands)
- Wall-time + GPU-residency numbers per sample

## Data inventory

### On disk (verified 2026-05-01)

```
/tmp/giab_chr20_full/                     ← chr20-only mirror, 3 samples
├── HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam (1.0 GB)
├── HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam (1.0 GB)
└── HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam (1.0 GB)

/tmp/dv_giab/data/                        ← chr20-only working dir
├── GRCh38.fa                  ← chr20-only (62 MB)
├── HG002.bam → /tmp/giab_chr20_full/HG002…chr20.bam
├── truth.vcf.gz              ← HG002 v4.2.1 whole-genome (156 MB)
└── truth.bed                 ← HG002 high-confidence WG regions (11 MB)
```

### Required for whole-genome (NOT yet downloaded)

| Artefact | URL (canonical) | Size |
|----------|-----------------|------|
| GRCh38 no_alt FASTA | `https://storage.googleapis.com/deepvariant/case-study-testdata/grch38_no_alt.fa` (or NCBI `seqs_for_alignment_pipelines.ucsc_ids/...no_alt_analysis_set.fasta.gz`) | 3.1 GB |
| HG002 NovaSeq 35× WG BAM | `https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam` | ~40 GB |
| HG003 NovaSeq 35× WG BAM | same path / `HG003…` | ~40 GB |
| HG004 NovaSeq 35× WG BAM | same path / `HG004…` | ~40 GB |
| HG003 v4.2.1 truth | `${GIAB_FTP}/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` + `.tbi` + `_noinconsistent.bed` | ~150 MB |
| HG004 v4.2.1 truth | same path / `HG004…` | ~150 MB |

The `validation/download_giab_full_genome.sh` script in tree currently
points at the **NIST 300× novoalign** BAMs, which are **NOT** the
canonical Google v1.10.0 benchmark fixture. Will be patched to the
`storage.googleapis.com/deepvariant/case-study-testdata/` URLs (matches
chr20 fixture provenance + Google published numbers).

## Feasibility analysis — disk budget

**Critical constraint**: 127 GB free on `/Users/benjamin`.

```
Per-sample intermediate disk peak (no chunking):
  examples.tfrecord  → ~600 GB ◄ blows disk catastrophically
  cvo + small_cvo    → ~1 GB
  output VCF + gVCF  → ~1 GB

  chr20 reference:    12.8 GB examples → 50× scale = ~640 GB
```

**No-go for whole-genome WITHOUT chunking.**

### Mitigation strategies (ranked)

#### A. **Chunked WG execution** (recommended)

Partition by chromosome (or smaller). For each chunk:
1. Run make_examples on chunk → write examples to temp
2. call_variants on chunk → write CVOs
3. postprocess_variants on chunk → write VCF chunk
4. **Delete chunk's examples + CVOs**
5. Concatenate VCF chunks at end

Largest chunk = chr1 (~250 Mb, ~50 GB examples). Fits in 127 GB free
**after BAM is downloaded** (one BAM at a time):

```
After HG002 BAM download:  127 - 40 = 87 GB free
HG002 chr1 chunk peak:     87 - 50 = 37 GB free at peak  ← OK
```

Per-sample sequence:
1. Download BAM (~30-60 min @ ~50 MB/s)
2. Run chunked pipeline (~2.5 h compute)
3. hap.py vs truth (~20 min Docker)
4. Delete BAM + intermediate
5. Move to next sample

Estimated total wall-time: **~10-12 hours** for all 3 samples
sequential.

Engineering required: minor wrapper script around existing pipeline
(80-100 LOC). The native binary already accepts `--regions` so we just
loop over a list of chunk specs.

#### B. **Skip WG, polish chr20 trio**

Run chr20 only on all 3 samples with stratified breakdown. We have the
chr20 BAMs already. Need only HG003 + HG004 truth (~300 MB download).

Wall-time: ~30 min total (3 × 3 min runs + 3 × 5 min hap.py).

Limitation: chr20 only — not whole-genome. But chr20 is a 71k-truth-
variant fixture and is a standard publication subset; results
generalize well in practice.

#### C. **External storage**

Mount a USB-3 SSD or NVMe enclosure with ≥ 1 TB free. Run unmodified
pipeline. Simpler logistically, requires hardware availability.

## Recommended publication-grade plan

1. **Tier 1 (immediate, ~30 min)**: chr20 trio. Establishes per-sample
   F1 + Docker parity confirmation across HG003, HG004 (we already have
   HG002 chr20). Limited but defensible publication chunk.
2. **Tier 2 (~12 h)**: whole-genome trio via chunked execution
   (option A). Full publication-grade numbers. Disk-managed.
3. **Tier 3 (~1-2 d)**: stratified F1 via GIAB stratifications v3.6
   (lowcomplexity / segdup / MHC / GC bands) on Tier-2 outputs.

## Reproducibility checklist (publication appendix)

To be captured during the run:

- [ ] Build commit SHA + macOS version + Xcode CLT version
- [ ] Hardware: chip / cores / RAM
- [ ] BAM provenance: full URL + SHA-256
- [ ] Truth set version + URL + SHA-256
- [ ] Reference FASTA URL + SHA-256
- [ ] Exact `deepvariant run` command line per sample
- [ ] hap.py command + Docker image tag (`jmcdani20/hap.py:v0.3.12`)
- [ ] Run wall-time per sample (split: download / make_examples / call_variants / postprocess / hap.py)
- [ ] GPU residency from `powermetrics --samplers gpu_power -i 500` during call_variants
- [ ] Random seed values (we are deterministic by construction; document for completeness)

## Comparison baseline (Google v1.10.0 published)

For Illumina NovaSeq 35× PCR-free on GIAB v4.2.1 truth, Google publishes
in their release notes / case-studies (HG002 representative; HG003 and
HG004 numbers are typically in the same ballpark):

| Metric | HG002 (Google v1.10.0, WG) | This work (HG002 chr20) |
|--------|----------------------------|-------------------------|
| SNP F1 | typically 0.99961-0.99965 | 0.997402 (chr20 only) |
| INDEL F1 | typically 0.99654-0.99701 | 0.995942 (chr20 only) |

Note: chr20 F1 is typically **lower** than WG F1 for this caller (chr20
has elevated FN/FP density). WG numbers should be ≥ chr20 numbers.

The Tier-2 whole-genome run will produce the directly-comparable number.

## Risks / open questions

- **Wall-time variance**: chr20 was 3 min on idle M4 Max with 14
  threads. WG real wall-time depends on partition size + thread
  contention with htslib I/O. Could be 2-4 h per sample.
- **GPU residency**: never measured at WG scale; chr20 fixtures are too
  short to get reliable powermetrics samples. WG run is the first real
  measurement opportunity.
- **chr1 might exceed disk**: if chr1 examples.tfrecord turns out
  > 60 GB (above estimate), need finer chunking (e.g. chr1 split into
  100 Mb sub-chunks via `--regions=chr1:0-100000000` etc.).
- **Truth-set BED edge**: GIAB v4.2.1 BED uses
  `_noinconsistent.bed` (drops sites with caller-disagreement). Our
  numbers must use this BED, not the wider `_benchmark.bed`.
- **Network reliability**: 120 GB download from
  `storage.googleapis.com` typically works; partial download resume
  with `curl -C -` is in the script.

## Concrete next actions (waiting on user authorization)

1. **Patch `validation/download_giab_full_genome.sh`** to point at
   Google case-study URLs (matches chr20 fixture provenance).
2. **Implement chunked WG runner** `validation/run_giab_wg_chunked.sh`:
   - 25 chunks (chr1-22, X, Y, chrM)
   - Per-chunk pipeline + intermediate cleanup
   - VCF concat (`bcftools concat`) + bgzip + tabix at end
3. **Execute Tier 1** (chr20 trio) immediately — produces F1 in
   ~30 min, no download required for HG002, only truth files for
   HG003/HG004.
4. **Kick off Tier 2 download + run** in background (long-running,
   will fire user notification on completion).
5. **Generate `docs/validation.md`** with publication-ready tables +
   reproducibility appendix.

## Summary table for user

| Approach | Wall-time | Disk peak | Output quality |
|----------|-----------|-----------|----------------|
| Tier 1: chr20 trio (now) | ~30 min | ~15 GB | publication-OK fallback (chr20 fixture) |
| Tier 2: WG trio chunked | ~12 h | ~90 GB peak | publication-grade (canonical WG F1) |
| Tier 3: stratified | +1-2 d | small | stratified F1 breakdown |
