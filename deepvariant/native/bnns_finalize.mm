#include "deepvariant/native/bnns_finalize.h"

#include <cmath>
#include <cstring>
#include <utility>

#include "absl/log/log.h"
#include "deepvariant/native/dv_weights.h"

namespace deepvariant {

namespace {

constexpr const char* kDenseKernel =
    "layer_with_weights-188/kernel/.ATTRIBUTES/VARIABLE_VALUE";
constexpr const char* kDenseBias =
    "layer_with_weights-188/bias/.ATTRIBUTES/VARIABLE_VALUE";

bool LoadDense(const DvwWeights& weights,
               int* in_dim, int* out_dim,
               std::unique_ptr<float[]>* kernel,
               std::unique_ptr<float[]>* bias) {
  const auto* k = weights.Get(kDenseKernel);
  const auto* b = weights.Get(kDenseBias);
  if (!k || !b) {
    LOG(ERROR) << "BnnsFinalize: missing layer-188 kernel/bias";
    return false;
  }
  if (k->shape.size() != 2u || b->shape.size() != 1u) {
    LOG(ERROR) << "BnnsFinalize: bad shape for dense layer";
    return false;
  }
  // Source kernel is (in_dim, out_dim) — TF stores Dense as (input, output).
  const int in = static_cast<int>(k->shape[0]);
  const int out = static_cast<int>(k->shape[1]);
  if ((int)b->shape[0] != out) {
    LOG(ERROR) << "BnnsFinalize: bias size mismatch";
    return false;
  }
  *in_dim = in;
  *out_dim = out;

  // Transpose to (out, in) row-major so the inner loop is a contiguous
  // dot product over `in_dim` — same memory access pattern TF's matmul
  // uses on x86 (transpose_b=False).
  kernel->reset(new float[(size_t)out * in]);
  for (int o = 0; o < out; ++o) {
    for (int i = 0; i < in; ++i) {
      (*kernel)[(size_t)o * in + i] = k->data[(size_t)i * out + o];
    }
  }
  bias->reset(new float[out]);
  std::memcpy(bias->get(), b->data, (size_t)out * sizeof(float));
  return true;
}

}  // namespace

BnnsFinalize::BnnsFinalize() = default;
BnnsFinalize::~BnnsFinalize() = default;

std::unique_ptr<BnnsFinalize> BnnsFinalize::Create(
    const std::string& dvw_path) {
  auto w = DvwWeights::Open(dvw_path);
  if (!w) {
    LOG(ERROR) << "BnnsFinalize::Create: cannot open " << dvw_path;
    return nullptr;
  }
  return CreateFromWeights(*w);
}

std::unique_ptr<BnnsFinalize> BnnsFinalize::CreateFromWeights(
    const DvwWeights& w) {
  auto self = std::unique_ptr<BnnsFinalize>(new BnnsFinalize());
  if (!LoadDense(w, &self->in_dim_, &self->out_dim_,
                 &self->kernel_, &self->bias_)) {
    return nullptr;
  }
  return self;
}

bool BnnsFinalize::ApplyBatch(const float* features, int batch_size,
                              float* probs) const {
  if (!features || !probs || batch_size <= 0 ||
      !kernel_ || !bias_ || in_dim_ <= 0 || out_dim_ <= 0) {
    LOG(ERROR) << "BnnsFinalize::ApplyBatch: bad args";
    return false;
  }
  for (int n = 0; n < batch_size; ++n) {
    const float* x = features + (size_t)n * in_dim_;
    float* p = probs + (size_t)n * out_dim_;

    // Dense: logits[o] = sum_i x[i] * W[o, i]  + bias[o]
    //                    -- inner loop is sequential, single-threaded.
    // Each accumulator is a fresh FP32 register, so the order is
    // strictly i = 0, 1, …, in_dim_-1 with no parallel reduction.
    for (int o = 0; o < out_dim_; ++o) {
      const float* row = kernel_.get() + (size_t)o * in_dim_;
      float acc = 0.0f;
      for (int i = 0; i < in_dim_; ++i) {
        acc += x[i] * row[i];
      }
      p[o] = acc + bias_[o];
    }

    // Softmax with max-shift for numeric stability:
    //   m       = max_o logits[o]
    //   exp_o   = expf(logits[o] - m)
    //   probs_o = exp_o / sum(exp)
    float m = p[0];
    for (int o = 1; o < out_dim_; ++o) {
      if (p[o] > m) m = p[o];
    }
    float total = 0.0f;
    for (int o = 0; o < out_dim_; ++o) {
      const float e = std::exp(p[o] - m);
      p[o] = e;
      total += e;
    }
    const float inv = 1.0f / total;
    for (int o = 0; o < out_dim_; ++o) {
      p[o] *= inv;
    }
  }
  return true;
}

}  // namespace deepvariant
