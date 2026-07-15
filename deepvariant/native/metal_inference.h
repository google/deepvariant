// MPSGraph + Metal builder for the DeepVariant Inception-v3 big-model
// inference path. Phase 5.5 — replaces coreml_inference for the shipped
// binary. Reads weights from a `.dvw` bundle (see dv_weights.h).
//
// Architecture (mirrors tools/conversion/inception_v3_mil.py):
//
//     input (N, 100, 221, 7) NHWC FP32
//     ↓ NHWC → NCHW transpose
//     ↓ stem: 5× conv-bn-relu + 2× maxpool
//     ↓ 3× InceptionA (Mixed_5b, 5c, 5d)
//     ↓ Reduction-A (Mixed_6a)
//     ↓ 4× InceptionB (Mixed_6b, 6c, 6d, 6e)
//     ↓ Reduction-B (Mixed_7a)
//     ↓ 2× InceptionC (Mixed_7b, 7c)
//     ↓ global avg pool → (N, 2048)
//     output: (N, 2048) FP32 features (pre-dense, pre-softmax)
//
// The final dense (2048→3) + softmax goes through BnnsFinalize for
// deterministic CPU reduction, NOT through this MPSGraph (see
// bnns_finalize.h). That split is what gets us bit-parity with TF on
// the final per-class probabilities.
//
// Threadsafe for Predict() once Create() succeeds; the graph is
// immutable after build.
#pragma once

#include <cstddef>
#include <memory>
#include <string>

namespace deepvariant {

class MetalInception {
 public:
  // Open the `.dvw` weight bundle and build the MPSGraph.  Returns
  // nullptr on any error (file missing, weight tensor missing, MPSGraph
  // failure).
  //
  // input_height/input_width/input_channels parameterize the placeholder
  // input shape. WGS: (100,221,7). DeepTrio WGS: (140,221,7).
  // PacBio germline: (100,147,10). ONT: (100,199,10). MASSEQ: (100,199,9).
  // Somatic PacBio TN: (200,147,9). Somatic ONT TN: (200,99,9).
  static std::unique_ptr<MetalInception> Create(
      const std::string& dvw_path,
      int input_height = 100,
      int input_channels = 7,
      int input_width = 221);

  ~MetalInception();

  // Run inference on a batch of pileup images.
  //
  //   input  : (batch_size, 100, 221, 7) FP32 NHWC, contiguous
  //   output : (batch_size, 2048) FP32 features
  //
  // Returns false on dispatch error.
  bool Predict(const float* input, int batch_size, float* output);

  // Debug-only: run the graph but stop at one of the named tap points
  // and return that tensor's output instead of the global-avg-pool
  // features. Used by tools/debug_metal_layer.cc to localise where
  // Metal output diverges from the Core ML / TF reference.
  //
  // Tap names (in order of execution):
  //   "stem_s1a"   — output of CBR(conv=0, bn=1) — shape (B, 32, 49, 110)
  //   "stem_s2a"   — output of CBR(conv=2, bn=3) — (B, 32, 47, 108)
  //   "stem_s2b"   — CBR(4,5) — (B, 64, 47, 108)
  //   "stem_mp3a"  — maxpool — (B, 64, 23, 53)
  //   "stem_s3b"   — CBR(6,7) — (B, 80, 21, 51)
  //   "stem_s4a"   — CBR(8,9) — (B, 192, 19, 49)
  //   "stem_mp5a"  — maxpool — (B, 192, 9, 24)
  //   "5b" / "5c" / "5d" / "6a" / "6b" / "6c" / "6d" / "6e" / "7a" / "7b" / "7c"
  //   "gap"        — global avg pool — (B, 2048)  (default Predict tap)
  //
  // The output buffer must be sized for the requested tap. Returns
  // false on unknown tap name or dispatch error.
  bool PredictAtTap(const std::string& tap_name,
                    const float* input, int batch_size,
                    float* output, int* out_total_elems_per_image);

  // Number of per-example floats Predict() writes:
  //   - default: 2048 (post-GAP feature vector, BnnsFinalize follows)
  //   - DV_METAL_GPU_FINALIZE=1: 3 (post-softmax probabilities; bypass
  //     BnnsFinalize)
  int FeatureDim() const;

  // True if DV_METAL_GPU_FINALIZE=1 selected at Create() — Predict()
  // emits softmax probabilities directly. Callers should skip
  // BnnsFinalize::ApplyBatch when this returns true.
  bool IsGpuFinalize() const;

  MetalInception(const MetalInception&) = delete;
  MetalInception& operator=(const MetalInception&) = delete;

 private:
  MetalInception();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
