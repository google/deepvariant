# Whole-Genome Validation Plan — pre-production hardening

Date: 2026-05-01

This document tracks the three blockers identified before declaring
the port "production-ready for whole-genome cohort runs":

1. macOS disk-write quota / Jetsam crash on chr1
2. Stability validation across ≥ 10 consecutive WGS runs
3. Whole-genome Docker FILTER-parity not measured (only chr20)

---

## 1. macOS disk-write quota crash — **FIX SHIPPED**

**Symptom**: deepvariant process killed silently mid-run on chr1 of a
WG benchmark. `time` reports "Invalid argument" signal. macOS unified
log (`Library/Logs/DiagnosticReports/deepvariant_*.diag`) shows:

```
Writes:        137.44 GB of file backed memory dirtied over 10720 s
               (12.82 MB/s avg), exceeding limit of 1590.73 KB/s over 86400 s
Action taken:  none (warning) → eventually SIGKILL via Jetsam
```

**Root cause**: `std::ofstream` (the previous TFRecordWriter backend)
keeps writes in the userspace buffer and lets the kernel page-cache
absorb them. Pages stay dirty for seconds until macOS flushes them to
disk asynchronously. At our sustained ~120 MB/s write rate (1 GB
shard × 14 shards in ~80 s), dirty pages accumulate faster than the
flusher can clear, hitting macOS's per-process Jetsam quota.

**Fix** (commit pending, post-trio): refactor TFRecordWriter to use
raw POSIX fd with `fcntl(F_NOCACHE, 1)`. F_NOCACHE bypasses the
unified buffer cache; writes go straight to the SSD device with no
kernel-side dirty-page accounting. A 1 MiB userspace coalescing
buffer is preserved so the SSD can still batch writes efficiently;
chr20:10M-10.1M smoke test confirms 0 perf regression and 100 %
Docker FILTER parity preserved.

**Verification needed**: re-run a full chr1 chunk after the fix lands
to confirm no Jetsam crash. ETA: ~50 min compute.

## 2. Stability across ≥ 10 consecutive WGS runs — **NOT DONE**

**Risk**: any rare deterministic bug (assertion, OOM, htslib edge
case, Metal driver hiccup) that fires once per N WGS will derail a
400-patient cohort. The current "Phase 4 PASS" verdict is on chr20
trio (3 samples × 63 Mb each). No WGS has yet completed end-to-end on
this build.

**Validation plan**:

```bash
# Pilot batch: 10 WGS samples, sequential, full pipeline
for i in $(seq 1 10); do
  rm -rf /tmp/dv_pilot_${i}
  ./build-macos/bin/deepvariant run \
    --reads=/path/to/sample_${i}.bam \
    --ref=/tmp/dv_giab/full/GRCh38.fa \
    --output_vcf=/tmp/dv_pilot_${i}/out.vcf.gz \
    --intermediate_results_dir=/tmp/dv_pilot_${i} \
    --inference_backend=metal \
    --model_type=WGS \
    --checkpoint=validation/work/wgs.dvw \
    --num_shards=14 \
    --batch_size=512 \
    > /tmp/dv_pilot_${i}.log 2>&1
  ec=$?
  echo "Sample $i exit code: $ec"
done
```

**Pass criteria**:
- 10/10 runs complete with exit code 0
- Each VCF has ≥ 4 M variants (sanity check on yield)
- No Jetsam events in `Library/Logs/DiagnosticReports/`
- No detectable monotonic memory leak (Activity Monitor:
  RSS at end of run #10 ≤ 1.5× RSS at end of run #1)
- Wall-time variance ≤ 10 % across the 10 runs

**Compute budget**: 10 × ~5h25 = ~54 h on 1 M4 Max. Practical: launch
overnight × 5 nights, 2 samples per night.

**Status**: NOT STARTED. Pre-requisite: fix #1 must be in.

## 3. Whole-genome Docker FILTER-parity — **NOT MEASURED**

**Risk**: chr20-only validation extrapolates the FP-drift residue
linearly across the genome. At ~10⁻⁵ softmax drift per call and ~5 M
calls genome-wide, residue FILTER mismatches could be 1000-10000 sites
(0.02-0.2 % of PASS-set drift). We have not measured this.

**Validation plan**:

```bash
# 1. Run our pipeline on HG002 WGS
./build-macos/bin/deepvariant run \
  --reads=/tmp/dv_giab/full/HG002.bam \
  --ref=/tmp/dv_giab/full/GRCh38.fa \
  --output_vcf=/tmp/dv_wg_ours/HG002.vcf.gz \
  --intermediate_results_dir=/tmp/dv_wg_ours \
  --inference_backend=metal \
  --model_type=WGS \
  --checkpoint=validation/work/wgs.dvw \
  --num_shards=14 \
  --batch_size=512

# 2. Run google/deepvariant:1.10.0 Docker on the SAME inputs at the
#    SAME shard count (8) — important for shard-count-conditional
#    parity.
docker run --rm \
  -v /tmp/dv_giab/full:/data:ro \
  -v /tmp/dv_wg_docker:/work \
  google/deepvariant:1.10.0 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/GRCh38.fa \
    --reads=/data/HG002.bam \
    --output_vcf=/work/output.vcf.gz \
    --num_shards=14    # IMPORTANT: match our num_shards

# 3. Diff
bash validation/diff_filter_classes.sh \
  /tmp/dv_wg_ours/HG002.vcf.gz \
  /tmp/dv_wg_docker/output.vcf.gz
```

**Pass criteria**:
- ≥ 99.95 % FILTER-class parity on shared sites (i.e. ≤ 0.05 % FM)
- 100 % CHROM/POS/REF/ALT identity on shared
- 100 % GT identity on shared
- PASS-set asymmetric difference ≤ 1 000 sites of ~5 M (≤ 0.02 %)

**Compute budget**:
- Ours WGS: ~5h25 (post-optims, M4 Max)
- Docker WGS (Rosetta 2 emulation): ~22 h
- Combined: ~28 h on the same Mac, sequential
- Practical: run ours during day, Docker overnight × 1 day

**Status**: NOT STARTED. Pre-requisite: fix #1, ideally fix #2.

---

## Decision tree before "production for 400 WGS"

```
                  Fix #1 deployed
                        ↓
              Pilot 10 WGS (54h)
                        ↓
              All 10 pass? ──no──→ Investigate failure pattern; loop
                        ↓ yes
              WG Docker parity ≥ 99.95 %?
                        ↓
                no ─→ Document residue magnitude; decide whether
                      to ship anyway (clinical context dependent)
                ↓ yes
            ✓ Production-ready for 400 WGS cohort
```

**Estimated time to production-ready**: 1-2 weeks of mostly-overnight
runs once Fix #1 is in tree.

---

## Why Fix #1 alone is not enough

A single passing chr1 chunk after F_NOCACHE doesn't prove WGS
stability. macOS Jetsam has multiple triggers (RSS, CPU time, file
descriptors, mach ports, …) and we've only addressed the dirty-page
one. Other failure modes that could surface at WGS scale but not
chr20:

- htslib mmap pressure: 46 GB BAM × num_shards parallel readers can
  trigger VM map exhaustion
- MPSGraph executable cache: per-batch_size compile is cached, but
  we instantiate a new MetalInception per stage — 25 chunks × 1 inst
  could exhaust some Metal pool
- TCC kernel auth events for read access (we've seen these in
  `.ips` reports as `EXC_CRASH` during `libsystem_info.dylib::User by ID`)

These need empirical exposure → 10 WGS pilot is the test.

---

## Post-validation: production runbook for 400 WGS

Once #1, #2, #3 are GREEN:

1. Set up 5-Mac fleet (Mac Studio M4 Max recommended for unified
   memory + GPU power, 64 GB+ RAM each)
2. Per-Mac script (sequential, retry-on-failure):
   ```bash
   for sample in $(cat samples_for_mac_$n.txt); do
     ./scripts/run_one_wgs.sh ${sample} || \
       (echo "RETRY $sample" >> retries.log && \
        ./scripts/run_one_wgs.sh ${sample})
   done
   ```
3. Monitor (per Mac):
   - `top -o cpu` for CPU saturation
   - `powermetrics --samplers gpu_power -i 30000` for GPU residency
   - `df -h` for disk
   - `vmstat 60` for swap pressure
4. Output: 1 VCF per sample to a shared NFS / SMB / etc.
5. Failure handling: any sample that fails twice → flag for manual
   review, do NOT block the cohort progression
6. Expected wall-time: 80 samples per Mac × 5h25 = ~18 days; 5 Macs
   parallel → cohort done in ~18 days

If the cohort is time-critical (< 1 week), use cloud GPU instead.

---

## Open question: the residue may be cohort-sample-dependent

Different patient BAMs have different read distributions, coverage
patterns, and structural variation. The FP-drift residue at any given
GQ-borderline site depends on the input vector through 188 conv
layers. We have measured it on HG002/HG003/HG004 GIAB Ashkenazi trio
fixtures only. A Krebs / cancer / FFPE BAM could produce different
residue magnitudes.

For a clinical cohort, this means: validate the residue on **a
representative sub-sample of the cohort** (e.g. 5 patients spanning
the expected BAM characteristics) before committing to the full 400.

Compute budget: 5 patients × 5h25 = 27 h ours + 5 × 22 h Docker = 137 h
total. ~6 days of sequential validation.
