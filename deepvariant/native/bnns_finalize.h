// Deterministic CPU dense (2048→3) + softmax for the Inception-v3
// classifier head. Phase 5.5 — designed to be bit-identical to TF's
// CPU `tf.nn.softmax(tf.matmul(x, W) + b)` output.
//
// The MPSGraph Inception backbone (`metal_inference.{h,mm}`) emits a
// (B, 2048) feature vector. We finalize on CPU with a sequential
// reduction (no SIMD, no parallel sum tree) to guarantee a single
// well-defined FP32 ordering, which is the only way to match TF's
// reference output reproducibly across M-series chip generations.
//
// Despite the "BNNS" name we currently use a hand-rolled sequential
// matmul (3 outputs × 2048 inputs = 6144 FMA operations per example
// — well under 10 µs even single-threaded). The BNNS framework is
// kept as a future optimization if we ever need to push throughput
// higher; the *deterministic* path stays the hand-rolled one.
//
// Weights are read from a .dvw bundle:
//     layer_with_weights-188/kernel  shape (2048, 3)  HWIO-style
//     layer_with_weights-188/bias    shape (3,)
//
// Threadsafe for ApplyBatch() once the constructor returns.
#pragma once

#include <memory>
#include <string>

namespace deepvariant {

class DvwWeights;  // forward-declared

class BnnsFinalize {
 public:
  // Open the .dvw and pull layer_with_weights-188's kernel + bias.
  // Returns nullptr if the bundle doesn't have a matching dense layer
  // (e.g. wrong model variant).
  static std::unique_ptr<BnnsFinalize> Create(const std::string& dvw_path);

  // As Create() but consumes a pre-opened DvwWeights (sharing the
  // mmap with metal_inference).  Does NOT take ownership.
  static std::unique_ptr<BnnsFinalize> CreateFromWeights(
      const DvwWeights& weights);

  ~BnnsFinalize();

  // Apply dense + softmax to a batch of feature vectors.
  //   features : (batch_size, 2048) FP32, row-major
  //   probs    : (batch_size, 3)    FP32, row-major
  // Returns false on size mismatch.
  bool ApplyBatch(const float* features, int batch_size,
                  float* probs) const;

  int InputDim() const { return in_dim_; }
  int OutputDim() const { return out_dim_; }

  BnnsFinalize(const BnnsFinalize&) = delete;
  BnnsFinalize& operator=(const BnnsFinalize&) = delete;

 private:
  BnnsFinalize();
  // Owns: the kernel matrix in row-major (out_dim, in_dim) layout
  // (transposed from the .dvw's (in_dim, out_dim) so the inner loop
  // strides 1 along the input axis — same as TF's MatMul kernel
  // when transpose_b=False) and the bias.
  int in_dim_ = 0;
  int out_dim_ = 0;
  std::unique_ptr<float[]> kernel_;  // [out_dim_ * in_dim_]
  std::unique_ptr<float[]> bias_;    // [out_dim_]
};

}  // namespace deepvariant
