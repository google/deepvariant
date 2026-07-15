# Packaging — Single-Binary Distribution on Homebrew

**Status:** Draft (will be filled in during Phase 5).
**Branch:** `feature/apple-silicon-native-v2`.

## Goal

One signed/notarized arm64 Mach-O at ~150-300 MB, plus a separate ~8.5 GB `deepvariant-models` formula. Both installed via:

```sh
brew tap benjamindemaille/deepvariant
brew install deepvariant deepvariant-models
deepvariant run --model_type=WGS --reads=in.bam --ref=ref.fa --output_vcf=out.vcf
```

No compilation on the user's machine. Cold-cache `brew install deepvariant` < 60 s.

## Binary layout (planned)

```text
$HOMEBREW_PREFIX/
├── Cellar/deepvariant/<version>/
│   └── bin/deepvariant            (single signed Mach-O, all deps static)
├── share/deepvariant-models/<version>/
│   ├── wgs.mlpackage
│   ├── wes.mlpackage
│   ├── pacbio.mlpackage
│   ├── ont.mlpackage
│   ├── trio_parent.mlpackage
│   ├── trio_child.mlpackage
│   ├── ...
│   ├── somatic_*.mlpackage
│   └── pangenome_*.mlpackage      (~15-20 mlpackages total, ~8.5 GB)
```

## Static linking inventory

| Lib | Source | Static? |
| --- | --- | --- |
| htslib 1.18 | FetchContent / submodule | Yes |
| libssw 1.2.5 | submodule | Yes |
| abseil-cpp 20240722 | FetchContent | Yes |
| protobuf 21.9 | FetchContent | Yes |
| gbwt / gbwtgraph / sdsl-lite / libdivsufsort / libhandlegraph | submodules | Yes |
| Core ML.framework | system | Dynamic (system) |
| Foundation / Metal | system | Dynamic (system) |

Verification: `otool -L bin/deepvariant` should show only `/usr/lib/*` and `/System/*` paths.

## Code signing & notarization

- Sign with Apple Developer ID Application certificate via `codesign --options=runtime --timestamp`.
- Notarize via `xcrun notarytool submit ... --wait`.
- Staple ticket with `xcrun stapler staple`.
- Verify with `spctl --assess --verbose ./deepvariant` (must pass).

All four tools are in Xcode CLT — **no full Xcode required** on the build/release machine.

## Core ML model compilation strategy

We ship `.mlpackage` files **uncompiled**. The binary calls `[MLModel compileModelAtURL:url error:&err]` at first load; Core ML caches the resulting `.mlmodelc` in `~/Library/Caches/com.apple.CoreML/`. Subsequent runs are unaffected.

- Avoids requiring full Xcode (which bundles `xcrun coremlcompiler` for ahead-of-time compilation).
- Cost: first run adds a few seconds per model used. Logged as `Compiling Core ML model for first run…`.
- Cache invalidation is handled by Core ML (it re-compiles if the `.mlpackage` mtime changes).

## Bottle build flow (CI)

A self-hosted M-series GitHub Actions runner triggered on tag:

1. Build `deepvariant` static-linked.
2. Sign + notarize + staple.
3. Run conversion pipeline (for models bottle): produces all `.mlpackage`s, signs them, packs.
4. Upload bottles to GitHub Release.
5. Update tap formula sha256s.

Reproducibility: every dep pinned with sha256 in CMake `FetchContent_Declare`.

## Open questions (deferred to Phase 5)

- Bottle hosting beyond GitHub Releases? (Cloudflare R2 mirror if downloads scale.)
- Hardened runtime entitlements: do we need any? (Probably none — no JIT, no Metal capture.)
- Per-macOS-version bottle tags: `arm64_sequoia` (macOS 15) and `arm64_sonoma` (macOS 14) at minimum.
