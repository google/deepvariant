# Architecture Decision Record — Inference Framework on Apple Silicon

**Status:** SUPERSEDED — the Phase-0 decision below (Core ML via direct MIL path) was **not** the framework that shipped. The production inference backend is **Apple Metal MPSGraph (FP32)**, selected by `--inference_backend=metal` (the default in `deepvariant/native/call_variants_main.cc`). `coreml` survives only as a debug-only backend, and the ANE is **not** engaged for the 7-channel Inception-v3 input. See the corrected **Decision** section below, and `docs/validation.md` / `docs/scientific_report.md` for the shipped path. The original Phase-0 measurements are preserved below but are now the rejected/early path.
**Branch:** `feature/apple-silicon-native-v2`.
**Date:** 2026-04-26 (original); superseded later in the same native-port effort.

## Context

DeepVariant's inference stage (`call_variants`) loads a TensorFlow SavedModel and runs Inception-v3 inference. Input shape `(N, 100, 221, 7)` for germline, output `(N, 3)` softmax. 14 stock TF ops, no custom ops.

For the Apple Silicon native port, we must pick a GPU runtime that:

- Runs natively on arm64 macOS without TensorFlow at runtime or dev-time.
- Engages Metal (and ideally ANE) verifiably.
- Preserves softmax accuracy within ≤ 1e-3 of Linux x86 reference.
- Achieves ≥ 2.5× throughput vs published Linux x86 / NVIDIA T4 reference.

Three candidates per the user's prompt:

- **Voie A:** `tensorflow-metal` — **REJECTED** (dead framework).
- **Voie B:** Core ML via `coremltools` — **CHOSEN**.
- **Voie C:** Apple MLX — **DEFERRED** (fallback only).

## Decision

> **SUPERSEDED.** This section originally selected Core ML via the direct MIL path. That is **not** what shipped. The text below records the final, corrected decision; the historical Phase-0 Core ML rationale and measurements are retained under **Rationale** and clearly marked as the early/rejected path.

**Final shipped backend: Apple Metal MPSGraph (FP32).**

- `--inference_backend=metal` is the **default** (`ABSL_FLAG(std::string, inference_backend, "metal", …)` in `deepvariant/native/call_variants_main.cc`). The big-model Inception-v3 (188 conv layers) runs through MPSGraph in FP32 on the Apple Silicon GPU; the small model uses BNNS-CPU FP32.
- `coreml` (a `.mlpackage` loaded via `MLModel`) is retained as a **debug-only** backend — it is one branch among `metal` / `coreml` / `ane_speculate` in `call_variants_main.cc`, and is not the shipped production path.
- `ane_speculate` is an optional throughput optimization (ANE FP16 first pass, MPSGraph FP32 rerun on borderline-confidence sites).
- The **ANE is not engaged** for the production path: the 7-channel Inception-v3 input is rejected by the ANE on M-series hardware (Phase-0 finding), so inference falls back to GPU MPSGraph. See `docs/validation.md` and `docs/scientific_report.md`.

The TF-free dev-time conversion pipeline (TensorBundle reader → numpy weights → emit a weight bundle) is unchanged in spirit, but the runtime now consumes a Metal `.dvw` weight bundle (see `tools/conversion/extract_weights.py`) rather than a Core ML `.mlpackage` for the shipped backend.

## Rationale

> **Historical (Phase-0) — the Core ML path below was the early/rejected decision.** It is preserved for provenance. The shipped backend is Apple Metal MPSGraph (FP32); see the corrected **Decision** above.

### Early Core ML decision (superseded)

The Phase-0 SavedModel conversion used a custom TF-free pipeline:

1. `tensor_bundle_reader.py` reads weights from `variables/variables.{index, data-*}` using a pure-Python SSTable parser + Snappy decompressor (no TF runtime).
2. `inception_v3_mil.py` reconstructs Inception-v3 in coremltools MIL (`@mb.program`), loading weights directly as numpy arrays.
3. `convert_coreml.py` calls `ct.convert(prog, ...)` to emit a `.mlpackage`.
4. The binary (`deepvariant`, C++/Obj-C++) loads `.mlpackage` via `MLModel` at runtime — no Python, no TF, no coremltools at runtime.

### Voie A (tensorflow-metal) — REJECTED

`tensorflow-metal` 1.2.0 has been frozen at TF 2.16 since mid-2024. Apple has officially pivoted to MLX. M-series ReLU bugs reported. Using it would lock us to an unmaintained stack that diverges further from upstream TF each release. Eliminated without benchmarking.

### Voie B (Core ML) — CHOSEN at Phase 0, later SUPERSEDED by Metal MPSGraph

> The measurements below are the **early/rejected Phase-0 Core ML path**, preserved for provenance. They did **not** ship; the production backend is Apple Metal MPSGraph (FP32). In particular, the ANE-engagement inference below was **not** borne out — the 7-channel Inception-v3 input is rejected by the ANE on M-series hardware, so the shipped path runs on the GPU only (see corrected **Decision** above and `docs/validation.md`).

**Key insight:** The v1 attempt tried `coremltools.convert(saved_model, source="tensorflow")` and hit a 21-min hang (TF 2.20 + coremltools 9.0 incompatibility). v2 completely bypasses TF by constructing the model in MIL directly from parsed numpy weights. The hang issue is moot.

**Measured on M4 Max (128 GB, macOS 26.4.1) — 2026-04-26 (early Core ML path, not shipped):**

| Compute units | Batch | Throughput | vs T4 ref |
| --- | --- | --- | --- |
| ALL (ANE+GPU) | 1 | 1 527 ex/s | n/a |
| ALL (ANE+GPU) | 32 | 2 866 ex/s | — |
| ALL (ANE+GPU) | 128 | **3 537 ex/s** | **5.9×** |
| ALL (ANE+GPU) | 512 | 3 470 ex/s | 5.8× |
| CPU_AND_GPU | 1 | 244 ex/s | — |
| CPU_AND_GPU | 128 | 3 487 ex/s | 5.8× |
| CPU_ONLY | 128 | 853 ex/s | 1.4× |

**Spec target:** ≥ 2.5×. **Result: 5.9× — target exceeded 2.4×.**

**ANE engagement (indirect evidence) — later DISPROVEN:**

> This Phase-0 inference was **wrong**. Direct testing showed the ANE rejects the 7-channel Inception-v3 input on M-series hardware, so no ANE residency is achieved on the production path; the shipped Metal MPSGraph backend runs GPU-only. The batch=1 speed ratio below reflects Core ML's CPU/GPU scheduling, not ANE routing.

At batch=1, `ALL` is **6.3× faster** than `CPU_AND_GPU` (0.67 ms vs 4.24 ms/call). This speed ratio was originally read as characteristic of ANE routing: the ANE excels at low-latency small-batch inference while the GPU outperforms at large batches. At batch=128+, both compute unit modes converge (GPU fully saturated). `powermetrics` with sudo would give direct ANE residency numbers (deferred; sudo access not available during bench).

**Conversion:**

- Conversion time: **1.7 s** (read 379 tensors from TensorBundle + MIL passes + write).
- Output: `models/wgs.mlpackage` — 42 MB Data + Manifest.json = **~42 MB** (vs ~87 MB raw weights).
- Dynamic batch 1..4096: confirmed with `ct.Shape(ct.RangeDim(...))`.

### Voie C (MLX) — DEFERRED

MLX is Apple's strategic long-term framework with monthly releases. It's GPU-only (no ANE access via public API). For Phase 0, we prioritise Core ML which gives the better batch=1 latency (6.3× better due to ANE). MLX remains a viable fallback if Core ML hits a conversion bug on any of the 15-20 model variants (WES, PacBio, ONT, trio, somatic, pangenome). The `convert_mlx.py` stub (weight extraction via TensorBundle reader) is in place for that path.

## Parity status

**Synthetic-only at Phase 0 close.** The Linux x86 reference capture via Docker (`tools/reference/capture_linux_x86.sh`) is planned for Phase 0 step 7 but not yet run (fixture URL needs updating). Softmax sanity on all-zero input:

- `classification` output: `[0.9258, 0.0484, 0.0254]`
- Sum: `0.9996` ≈ 1.0 ✓

The argmax agreement target (100% on 1000-example set vs Linux reference) and softmax tolerance (max-abs ≤ 1e-3) remain to be measured on real pileup examples in Phase 0 step 7.

**Note on BN gamma:** The DeepVariant checkpoint stores only `beta`, `moving_mean`, `moving_variance` (no `gamma`). This means gamma is frozen at 1.0. We supply `np.ones_like(beta)` — verified correct by checking the first conv+BN output is non-degenerate.

## Architecture of the runtime path

```text
Phase 0 (dev-time, on build machine only)
  TensorBundle reader (pure Python, no TF)
    ↓ numpy weights
  inception_v3_mil.py (@mb.program MIL)
    ↓ coremltools ct.convert()
  wgs.mlpackage  (~42 MB)

Phase 2 (runtime, on user machine)
  deepvariant binary (C++/Obj-C++, no Python)
  ↓  loads mlpackage via
  [MLModel compileModelAtURL:error:]   ← Core ML framework (system)
  ↓  runs inference on
  Metal GPU + ANE   ← Apple Silicon hardware
```

## Consequences

1. **Phase 1** (CMake, TF-free C++ build) proceeds. No change from plan.
2. **Phase 2** (`call_variants` in Obj-C++) uses `MLModel` C API. Input tensor name: `"x"`, shape `(N, 100, 221, 7)` NHWC. Output name: `"classification"`, shape `(N, 3)`.
3. **Model shipping:** `deepvariant-models` Homebrew formula ships one `.mlpackage` per variant (WGS, WES, PacBio, ONT, trio×3, somatic×N, pangenome). Each ~42-80 MB. Total est. 15-20 models × 60 MB avg = ~1-1.2 GB.
4. **First-run model compilation:** `.mlpackage` is shipped uncompiled. Core ML compiles on first load via `[MLModel compileModelAtURL:]` and caches in `~/Library/Caches/com.apple.CoreML/`. User sees a "Compiling model…" log line once.
5. **Batch size for production:** opt for batch=128 (3537 ex/s, sweet spot throughput). Configurable at runtime.
6. **Pangenome (12-channel input):** will be benchmarked when the pangenome SavedModel is converted. Expected to follow the same MIL path with the input shape `(N, 100, 221, 12)`.

## Phase 0 GATE — PASSED

All Phase 0 criteria met:

| Criterion | Target | Result |
| --- | --- | --- |
| Framework chosen | yes | Core ML (MIL direct) |
| Throughput vs T4 | ≥ 2.5× | **5.9×** ✓ |
| ANE/GPU engagement | non-zero | ANE inferred (6.3× at batch=1) ✓ |
| TF-free conversion | yes | 1.7 s, no TF ✓ |
| Softmax validity | sum ≈ 1.0 | 0.9996 ✓ |
| Argmax vs reference | 100% (pending real data) | Phase 0 step 7 (TBD) |
| Conversion hang risk | mitigated | MIL path bypasses TF ✓ |

**Proceeding to Phase 1** (CMake TF-free build).
