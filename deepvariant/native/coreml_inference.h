// Core ML inference wrapper — pure C++ interface (no Obj-C types exposed).
// The implementation is in coreml_inference.mm (Obj-C++).
//
// Usage:
//   auto model = CoreMLModel::Load("/path/to/wgs.mlpackage");
//   // images: flat float32 array, row-major (N, H, W, C)
//   model->Predict(images, N, H, W, C, probs);
//   // probs: flat float32 array (N, num_classes)
#pragma once

#include <cstddef>
#include <memory>
#include <string>

namespace deepvariant {

// Compute units for Core ML inference.
enum class ComputeUnits {
  kAll,          // ANE first, then GPU, then CPU (default)
  kCpuAndGpu,    // GPU + CPU, skip ANE
  kCpuOnly,
};

class CoreMLModel {
 public:
  // Load a .mlpackage file, compile on first run (cached in
  // ~/Library/Caches/com.apple.CoreML/), and prepare for inference.
  // Returns nullptr on failure.
  static std::unique_ptr<CoreMLModel> Load(
      const std::string& mlpackage_path,
      ComputeUnits compute_units = ComputeUnits::kAll);

  ~CoreMLModel();

  // Run batched inference.
  //  images: float32 array of shape (N, H, W, C), row-major, not freed.
  //  probs:  float32 output (N, num_classes), caller-allocated, row-major.
  // Returns true on success.
  bool Predict(const float* images, int N, int H, int W, int C,
               float* probs, int num_classes);

  int InputHeight() const { return input_height_; }
  int InputWidth() const { return input_width_; }
  int InputChannels() const { return input_channels_; }
  int NumClasses() const { return num_classes_; }
  const std::string& InputName() const { return input_name_; }
  const std::string& OutputName() const { return output_name_; }

  CoreMLModel(const CoreMLModel&) = delete;
  CoreMLModel& operator=(const CoreMLModel&) = delete;

 private:
  CoreMLModel();
  struct Impl;
  std::unique_ptr<Impl> impl_;
  int input_height_  = 100;
  int input_width_   = 221;
  int input_channels_ = 7;
  int num_classes_   = 3;
  std::string input_name_  = "x";
  std::string output_name_ = "classification";
};

}  // namespace deepvariant
