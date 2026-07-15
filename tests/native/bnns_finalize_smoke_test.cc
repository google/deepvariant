// Smoke test for the BNNS finalize layer (dense + softmax).
// Verifies that:
//   1) Create() loads layer-188 weights from a .dvw bundle.
//   2) ApplyBatch() produces (B, 3) probability vectors that sum to 1.
//   3) The implementation is deterministic across runs.

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <vector>

#include "deepvariant/native/bnns_finalize.h"
#include "gtest/gtest.h"

namespace deepvariant {
namespace {

const char* DvwPath() {
  if (const char* p = std::getenv("DV_WGS_DVW")) return p;
  return "validation/work/wgs.dvw";
}

TEST(BnnsFinalizeSmoke, LoadAndApply) {
  const std::string path = DvwPath();
  if (!std::filesystem::exists(path)) {
    GTEST_SKIP() << "DVW file not available at " << path;
  }
  auto fz = BnnsFinalize::Create(path);
  ASSERT_NE(fz, nullptr);
  EXPECT_EQ(fz->InputDim(), 2048);
  EXPECT_EQ(fz->OutputDim(), 3);

  // 4 dummy feature vectors, each filled with i / 2048.0 (i.e. 0..1).
  constexpr int B = 4;
  std::vector<float> features((size_t)B * 2048);
  for (int n = 0; n < B; ++n) {
    for (int i = 0; i < 2048; ++i) {
      features[(size_t)n * 2048 + i] =
          static_cast<float>(i) / 2048.0f * (n + 1);
    }
  }
  std::vector<float> probs((size_t)B * 3, 0.0f);
  ASSERT_TRUE(fz->ApplyBatch(features.data(), B, probs.data()));

  // Each row sums to 1 (within FP32 epsilon) and all entries are
  // non-negative.
  for (int n = 0; n < B; ++n) {
    float total = 0.0f;
    for (int o = 0; o < 3; ++o) {
      const float p = probs[(size_t)n * 3 + o];
      EXPECT_GE(p, 0.0f) << "negative prob at row " << n << " col " << o;
      EXPECT_LE(p, 1.0f);
      total += p;
    }
    EXPECT_NEAR(total, 1.0f, 1e-5f) << "row " << n << " does not sum to 1";
  }
}

TEST(BnnsFinalizeSmoke, Deterministic) {
  const std::string path = DvwPath();
  if (!std::filesystem::exists(path)) GTEST_SKIP();
  auto fz = BnnsFinalize::Create(path);
  ASSERT_NE(fz, nullptr);

  std::vector<float> features(2048);
  for (int i = 0; i < 2048; ++i) {
    features[i] = std::sin(0.01f * static_cast<float>(i));
  }
  std::vector<float> p1(3), p2(3);
  ASSERT_TRUE(fz->ApplyBatch(features.data(), 1, p1.data()));
  ASSERT_TRUE(fz->ApplyBatch(features.data(), 1, p2.data()));
  for (int o = 0; o < 3; ++o) {
    EXPECT_EQ(p1[o], p2[o]) << "non-deterministic output at " << o;
  }
}

}  // namespace
}  // namespace deepvariant
