// Smoke test: open a .dvw file, walk its tensor table, sanity-check
// the first conv kernel against the known TF SavedModel shape.

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "deepvariant/native/dv_weights.h"
#include "gtest/gtest.h"

namespace deepvariant {
namespace {

const char* DvwPath() {
  if (const char* p = std::getenv("DV_WGS_DVW")) return p;
  return "validation/work/wgs.dvw";
}

TEST(DvWeightsSmoke, LoadAndLookup) {
  const std::string path = DvwPath();
  if (!std::filesystem::exists(path)) {
    GTEST_SKIP() << "DVW file not available at " << path
                 << ". Build via tools/conversion/extract_weights.py.";
  }

  auto w = DvwWeights::Open(path);
  ASSERT_NE(w, nullptr) << "DvwWeights::Open(" << path << ") failed";
  EXPECT_EQ(w->Version(), 1u);
  EXPECT_GT(w->Names().size(), 100u)
      << "expected at least 100 tensors in the WGS bundle";

  // First conv kernel of WGS Inception-v3: 3×3 conv, 7→32 channels in
  // HWIO order (TF's stored layout).
  const std::string first_conv =
      "layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE";
  const auto* k0 = w->Get(first_conv);
  ASSERT_NE(k0, nullptr) << first_conv << " not found";
  ASSERT_EQ(k0->shape.size(), 4u);
  EXPECT_EQ(k0->shape[0], 3u);
  EXPECT_EQ(k0->shape[1], 3u);
  EXPECT_EQ(k0->shape[2], 7u);
  EXPECT_EQ(k0->shape[3], 32u);
  EXPECT_EQ(k0->n_elements, 3u * 3u * 7u * 32u);
  EXPECT_EQ(k0->n_bytes, k0->n_elements * sizeof(float));

  // Sanity-check the data is reachable and looks like real weights
  // (not all zero). At least one element should be non-zero.
  bool any_nonzero = false;
  for (size_t i = 0; i < k0->n_elements && !any_nonzero; ++i) {
    if (k0->data[i] != 0.0f) any_nonzero = true;
  }
  EXPECT_TRUE(any_nonzero);
}

TEST(DvWeightsSmoke, MissingNameReturnsNull) {
  const std::string path = DvwPath();
  if (!std::filesystem::exists(path)) {
    GTEST_SKIP() << "DVW file not available";
  }
  auto w = DvwWeights::Open(path);
  ASSERT_NE(w, nullptr);
  EXPECT_EQ(w->Get("does/not/exist"), nullptr);
}

}  // namespace
}  // namespace deepvariant
