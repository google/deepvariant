// Smoke test: load the small_model .mlpackage and run a forward pass.
// Does NOT verify outputs (those depend on the model weights); only
// confirms the wrapper loads, compiles, and produces 3-class softmax
// that sums to ~1.0 on a vector of 70 zeros.

#include "gtest/gtest.h"
#include "deepvariant/native/small_model_inference.h"

#include <cstdlib>
#include <filesystem>
#include <vector>

namespace deepvariant {
namespace {

const char* MlpackagePath() {
  if (const char* p = std::getenv("DV_SMALL_MODEL_MLPACKAGE")) return p;
  return "tools/conversion/models/wgs_small.mlpackage";
}

TEST(SmallModelSmoke, LoadAndPredict) {
  const std::string path = MlpackagePath();
  if (!std::filesystem::exists(path)) {
    GTEST_SKIP() << "Small model not available at " << path
                 << ". Build it via tools/conversion/convert_small_model.sh.";
  }
  auto m = SmallModel::Load(path);
  ASSERT_NE(m, nullptr) << "Failed to load " << path;

  // Predict on a vector of 70 zeros — verify shape and softmax sum.
  std::vector<float> features(70, 0.0f);
  std::vector<float> probs(3, 0.0f);
  ASSERT_TRUE(m->Predict(features.data(), 1, probs.data()));
  const float sum = probs[0] + probs[1] + probs[2];
  EXPECT_GT(sum, 0.99f);
  EXPECT_LT(sum, 1.01f);
}

}  // namespace
}  // namespace deepvariant
