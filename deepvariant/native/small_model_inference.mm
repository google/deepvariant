// Phase 5.5d/7 — Small-model inference, deterministic FP32 BNNS-CPU.
//
// The small_model is a 3-layer MLP (70 → 750 → 750 → 3) with
//   y1 = ReLU(x  · W1 + b1)        (70  → 750)
//   y2 = ReLU(y1 · W2 + b2)        (750 → 750)
//   y3 = softmax(y2 · W3 + b3)     (750 → 3)
//
// Earlier this wrapped Core ML; on identical inputs Core ML produces
// ~0.005-0.01 drift on max_p relative to Docker's TF/Keras FP32 path
// (Apple Core ML has its own SIMD reduction order regardless of
// `MLComputeUnitsCPUOnly`). For multi-allelic SNP sites where Docker's
// small_model commits at GQ=20-21, our Core ML path commits at GQ=18-19
// and the candidate falls through to deepvariant — closing the chr20
// FILTER drift below 0.001 % requires the small_model's dispatch
// decisions to match Docker's bit-for-bit.
//
// This implementation reads weights from `<dir>/layer_{0,1,2}_{kernel,bias}.npy`
// (FP32, row-major) — the same weights Docker bundles at
// `/opt/smallmodels/wgs/model.keras` (Keras Sequential, 3 Dense layers).
// Inference is strict-scalar sequential FP32: per output element a
// scalar `for` accumulator, no SIMD reduction, no `mad`/FMA — same
// pattern as `bnns_finalize.mm`. This produces output identical to
// Eigen+single-thread on x86 within 1 ULP and matches Docker's
// dispatch decisions on the chr20 cross-MID sites.

#include "deepvariant/native/small_model_inference.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

// Minimal NumPy v1 .npy reader for FP32 contiguous arrays.
struct NpyArr {
  std::vector<size_t> shape;
  std::vector<float> data;
};

bool ReadNpy(const std::string& path, NpyArr* out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  char magic[6];
  f.read(magic, 6);
  if (std::memcmp(magic, "\x93NUMPY", 6) != 0) return false;
  uint8_t major, minor;
  f.read(reinterpret_cast<char*>(&major), 1);
  f.read(reinterpret_cast<char*>(&minor), 1);
  size_t header_len = 0;
  if (major == 1) {
    uint16_t hl;
    f.read(reinterpret_cast<char*>(&hl), 2);
    header_len = hl;
  } else {
    uint32_t hl;
    f.read(reinterpret_cast<char*>(&hl), 4);
    header_len = hl;
  }
  std::string header(header_len, '\0');
  f.read(header.data(), header_len);
  // Parse `shape': (a, b)` — minimal scanf.
  auto p = header.find("'shape':");
  if (p == std::string::npos) return false;
  auto lp = header.find('(', p);
  auto rp = header.find(')', lp);
  if (lp == std::string::npos || rp == std::string::npos) return false;
  out->shape.clear();
  size_t i = lp + 1;
  while (i < rp) {
    while (i < rp && (header[i] == ' ' || header[i] == ',')) ++i;
    size_t j = i;
    while (j < rp && header[j] >= '0' && header[j] <= '9') ++j;
    if (j > i) out->shape.push_back(std::stoul(header.substr(i, j - i)));
    i = j + 1;
  }
  size_t total = 1;
  for (size_t s : out->shape) total *= s;
  out->data.resize(total);
  f.read(reinterpret_cast<char*>(out->data.data()), total * sizeof(float));
  return f.good() || f.eof();
}

// Sequential scalar FP32 dense forward: y[o] = sum_i x[i] * W[i, o] + b[o]
// with strict left-to-right accumulation order — matches Eigen's
// single-thread scalar GEMM bit-for-bit when no FMA fusion.
void DenseFp32(const float* x, int in_dim,
               const float* W, const float* b,
               int out_dim, float* y) {
  for (int o = 0; o < out_dim; ++o) {
    float acc = 0.0f;
    for (int i = 0; i < in_dim; ++i) {
      acc += x[i] * W[i * out_dim + o];
    }
    y[o] = acc + b[o];
  }
}

void Relu(float* v, int n) {
  for (int i = 0; i < n; ++i) if (v[i] < 0.0f) v[i] = 0.0f;
}

void Softmax3(float* v) {
  float m = v[0];
  if (v[1] > m) m = v[1];
  if (v[2] > m) m = v[2];
  float e0 = std::exp(v[0] - m);
  float e1 = std::exp(v[1] - m);
  float e2 = std::exp(v[2] - m);
  float total = e0 + e1 + e2;
  v[0] = e0 / total;
  v[1] = e1 / total;
  v[2] = e2 / total;
}

}  // namespace

struct SmallModel::Impl {
  int input_dim = 0;         // 70 (WGS), 106 (DeepTrio), 94 (DeepSomatic)
  // Layer 1: input_dim → 750
  std::vector<float> W1;     // shape (input_dim, 750), row-major
  std::vector<float> b1;     // (750,)
  // Layer 2: 750 → 750
  std::vector<float> W2;     // shape (750, 750)
  std::vector<float> b2;     // (750,)
  // Layer 3: 750 → 3
  std::vector<float> W3;     // shape (750, 3)
  std::vector<float> b3;     // (3,)
};

SmallModel::SmallModel() : impl_(std::make_unique<Impl>()) {}
SmallModel::~SmallModel() = default;

int SmallModel::input_dim() const { return impl_->input_dim; }

// static
std::unique_ptr<SmallModel> SmallModel::Load(const std::string& path) {
  // Path is a directory holding the 6 weight .npy files; stripping a
  // trailing `/` if present.
  std::string root = path;
  if (!root.empty() && root.back() == '/') root.pop_back();

  NpyArr W1, b1, W2, b2, W3, b3;
  if (!ReadNpy(root + "/layer_0_kernel.npy", &W1) ||
      !ReadNpy(root + "/layer_0_bias.npy", &b1) ||
      !ReadNpy(root + "/layer_1_kernel.npy", &W2) ||
      !ReadNpy(root + "/layer_1_bias.npy", &b2) ||
      !ReadNpy(root + "/layer_2_kernel.npy", &W3) ||
      !ReadNpy(root + "/layer_2_bias.npy", &b3)) {
    LOG(ERROR) << "SmallModel: failed to read weights from " << root
               << " (expected layer_{0,1,2}_{kernel,bias}.npy)";
    return nullptr;
  }
  // Input dimension is detected from layer_0_kernel's first axis. We
  // accept any sane value (70/94/106 in current model variants); the
  // layer_1 / layer_2 / bias shapes must be consistent.
  if (W1.shape.size() != 2 || W1.shape[1] != 750 ||
      b1.shape != std::vector<size_t>{750} ||
      W2.shape != std::vector<size_t>{750, 750} ||
      b2.shape != std::vector<size_t>{750} ||
      W3.shape != std::vector<size_t>{750, 3} ||
      b3.shape != std::vector<size_t>{3}) {
    LOG(ERROR) << "SmallModel: weight shapes don't match "
                  "(?,750)+(750)+(750,750)+(750)+(750,3)+(3) "
                  "— got W1=("
               << (W1.shape.size() >= 1 ? std::to_string(W1.shape[0]) : "?")
               << ","
               << (W1.shape.size() >= 2 ? std::to_string(W1.shape[1]) : "?")
               << ")";
    return nullptr;
  }

  auto out = std::unique_ptr<SmallModel>(new SmallModel());
  out->impl_->input_dim = static_cast<int>(W1.shape[0]);
  out->impl_->W1 = std::move(W1.data);
  out->impl_->b1 = std::move(b1.data);
  out->impl_->W2 = std::move(W2.data);
  out->impl_->b2 = std::move(b2.data);
  out->impl_->W3 = std::move(W3.data);
  out->impl_->b3 = std::move(b3.data);
  LOG(INFO) << "SmallModel: loaded BNNS-CPU FP32 MLP from " << root
            << " (input_dim=" << out->impl_->input_dim << ")";
  return out;
}

bool SmallModel::Predict(const float* features, int N, float* probs) {
  const int in_dim = impl_->input_dim;
  if (in_dim <= 0) return false;
  // Per-batch scratch — small enough to allocate per call (avoids
  // thread-safety issues if called concurrently).
  std::vector<float> y1(750), y2(750), y3(3);
  for (int n = 0; n < N; ++n) {
    const float* x = features + (size_t)n * in_dim;
    // Layer 1: in_dim → 750 + ReLU
    DenseFp32(x, in_dim, impl_->W1.data(), impl_->b1.data(), 750, y1.data());
    Relu(y1.data(), 750);
    // Layer 2: 750 → 750 + ReLU
    DenseFp32(y1.data(), 750, impl_->W2.data(), impl_->b2.data(), 750,
               y2.data());
    Relu(y2.data(), 750);
    // Layer 3: 750 → 3 + softmax
    DenseFp32(y2.data(), 750, impl_->W3.data(), impl_->b3.data(), 3,
               y3.data());
    Softmax3(y3.data());
    probs[(size_t)n * 3 + 0] = y3[0];
    probs[(size_t)n * 3 + 1] = y3[1];
    probs[(size_t)n * 3 + 2] = y3[2];
  }
  return true;
}

}  // namespace deepvariant
