// Lightweight wrapper around the small_model 3-layer MLP. Input
// dimension is detected from the .npy weight files at load time:
//   WGS:           70 → 750 → 750 → 3 (single sample)
//   DeepTrio WGS:  106 → 750 → 750 → 3 (3 samples × 12 base feats + 70)
//   DeepSomatic:   94 → 750 → 750 → 3 (2 samples × 12 base feats + 70)
// One file per layer / role: layer_{0,1,2}_{kernel,bias}.npy.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace deepvariant {

class SmallModel {
 public:
  // Returns nullptr on load failure.
  static std::unique_ptr<SmallModel> Load(const std::string& mlpackage_path);
  ~SmallModel();

  // features: flat vector of N * input_dim() floats, row-major (one row
  // per candidate). probs: caller-allocated, size N * 3.
  bool Predict(const float* features, int N, float* probs);

  // Number of input features the loaded model expects (70 for WGS,
  // 106 for DeepTrio, 94 for DeepSomatic).
  int input_dim() const;

  SmallModel(const SmallModel&) = delete;
  SmallModel& operator=(const SmallModel&) = delete;

 private:
  SmallModel();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace deepvariant
