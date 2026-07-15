// Smoke test for the Phase 5.5 MPSGraph Inception-v3 backend.
// Verifies that:
//   1) MetalInception::Create() opens a .dvw bundle, builds the graph
//      without exception, and returns a valid object.
//   2) Predict() dispatches a zero-input batch without crashing.
//   3) Output is shape (B, 2048) FP32 and finite.
//
// Does NOT verify numerical correctness against TF here — that's
// parity_check_metal.py's job. This test just exercises the build +
// dispatch path on the real WGS weights.

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <vector>

#include "deepvariant/native/metal_inference.h"
#include "gtest/gtest.h"

namespace deepvariant {
namespace {

const char* DvwPath() {
  if (const char* p = std::getenv("DV_WGS_DVW")) return p;
  return "validation/work/wgs.dvw";
}

TEST(MetalInferenceSmoke, BuildAndDispatch) {
  const std::string path = DvwPath();
  if (!std::filesystem::exists(path)) {
    GTEST_SKIP() << "DVW file not available at " << path
                 << ". Build via tools/conversion/extract_weights.py.";
  }

  auto inf = MetalInception::Create(path);
  ASSERT_NE(inf, nullptr) << "MetalInception::Create failed";
  EXPECT_EQ(inf->FeatureDim(), 2048);

  // One-image batch of all zeros — Inception-v3 should produce some
  // deterministic feature vector; we just check shape + finiteness.
  constexpr int B = 1;
  constexpr int H = 100;
  constexpr int W = 221;
  constexpr int C = 7;
  std::vector<float> input((size_t)B * H * W * C, 0.0f);
  std::vector<float> output((size_t)B * 2048, 0.0f);
  ASSERT_TRUE(inf->Predict(input.data(), B, output.data()));

  // No NaNs / Infs.
  size_t n_nonzero = 0;
  for (float v : output) {
    ASSERT_TRUE(std::isfinite(v)) << "non-finite output element";
    if (v != 0.0f) ++n_nonzero;
  }
  // With BN biases learned from real data, a zero pileup should produce
  // many non-zero activations after 188 conv layers.
  EXPECT_GT(n_nonzero, 100u);
}

}  // namespace
}  // namespace deepvariant
