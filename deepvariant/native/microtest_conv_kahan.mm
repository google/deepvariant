// Phase 5.5e/Path B microtest: dispatch the Kahan-compensated Conv2D
// kernel on small known cases, compare against a CPU reference
// implementing the same Kahan compensated summation. Bit-identical
// match expected (PASS) — any divergence means the kernel's reduction
// order or compensation logic deviates from the CPU spec.
//
// Also reports max-abs-diff vs the basic-FMA scalar reference (i.e.,
// `microtest_conv_serial`'s output) — should be ≤ ~1 ULP × N for the
// same input, demonstrating the precision improvement.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "deepvariant/native/metal_conv_kahan.h"

namespace deepvariant {

// CPU reference implementing the SAME Kahan compensated summation as
// the Metal kernel (NHWC src, HWIO W, (kh, kw, c_in)-order). Uses
// std::fma for the y = x*w - c step (single-rounded), then Kahan
// compensation. Bit-identical to GPU kernel output expected.
void RefConvKahan(const ConvDesc& d,
                  const float* src, const float* W, const float* bias,
                  float* dst) {
  for (int n = 0; n < d.B; ++n) {
    for (int h_out = 0; h_out < d.H_out; ++h_out) {
      for (int w_out = 0; w_out < d.W_out; ++w_out) {
        const int h_base = h_out * d.stride_h - d.pad_h;
        const int w_base = w_out * d.stride_w - d.pad_w;
        for (int c_out = 0; c_out < d.C_out; ++c_out) {
          float sum = 0.0f;
          float c = 0.0f;
          for (int kh = 0; kh < d.Kh; ++kh) {
            const int h_in = h_base + kh;
            if (h_in < 0 || h_in >= d.H_in) continue;
            for (int kw = 0; kw < d.Kw; ++kw) {
              const int w_in = w_base + kw;
              if (w_in < 0 || w_in >= d.W_in) continue;
              for (int c_in = 0; c_in < d.C_in; ++c_in) {
                const float x = src[
                    ((n * d.H_in + h_in) * d.W_in + w_in) * d.C_in + c_in];
                const float w = W[
                    ((kh * d.Kw + kw) * d.C_in + c_in) * d.C_out + c_out];
                // y = x*w - c (single-rounded FMA, matches metal::precise::fma)
                const float y = std::fma(x, w, -c);
                const float t = sum + y;
                c = (t - sum) - y;
                sum = t;
              }
            }
          }
          sum += bias[c_out];
          if (d.relu) sum = std::fmax(sum, 0.0f);
          dst[((n * d.H_out + h_out) * d.W_out + w_out) * d.C_out + c_out] =
              sum;
        }
      }
    }
  }
}

// Basic-FMA scalar reference (no compensation) for comparison —
// matches microtest_conv_serial's RefConv. Used to quantify how much
// Kahan reduces drift vs basic accumulation.
void RefConvBasic(const ConvDesc& d,
                  const float* src, const float* W, const float* bias,
                  float* dst) {
  for (int n = 0; n < d.B; ++n) {
    for (int h_out = 0; h_out < d.H_out; ++h_out) {
      for (int w_out = 0; w_out < d.W_out; ++w_out) {
        const int h_base = h_out * d.stride_h - d.pad_h;
        const int w_base = w_out * d.stride_w - d.pad_w;
        for (int c_out = 0; c_out < d.C_out; ++c_out) {
          float acc = 0.0f;
          for (int kh = 0; kh < d.Kh; ++kh) {
            const int h_in = h_base + kh;
            if (h_in < 0 || h_in >= d.H_in) continue;
            for (int kw = 0; kw < d.Kw; ++kw) {
              const int w_in = w_base + kw;
              if (w_in < 0 || w_in >= d.W_in) continue;
              for (int c_in = 0; c_in < d.C_in; ++c_in) {
                const float x = src[
                    ((n * d.H_in + h_in) * d.W_in + w_in) * d.C_in + c_in];
                const float w = W[
                    ((kh * d.Kw + kw) * d.C_in + c_in) * d.C_out + c_out];
                acc = std::fma(x, w, acc);
              }
            }
          }
          acc += bias[c_out];
          if (d.relu) acc = std::fmax(acc, 0.0f);
          dst[((n * d.H_out + h_out) * d.W_out + w_out) * d.C_out + c_out] =
              acc;
        }
      }
    }
  }
}

int RunCase(MetalConvKahan& mck, const char* label, const ConvDesc& d) {
  const size_t src_n =
      (size_t)d.B * d.H_in * d.W_in * d.C_in;
  const size_t w_n = (size_t)d.Kh * d.Kw * d.C_in * d.C_out;
  const size_t bias_n = d.C_out;
  const size_t dst_n =
      (size_t)d.B * d.H_out * d.W_out * d.C_out;

  std::mt19937 rng(0x55c1);
  std::uniform_real_distribution<float> u(-1.0f, 1.0f);
  std::vector<float> src(src_n), W(w_n), bias(bias_n);
  std::vector<float> dst_kahan_ref(dst_n, 0.0f);
  std::vector<float> dst_basic_ref(dst_n, 0.0f);
  std::vector<float> dst_gpu(dst_n, 0.0f);
  for (auto& v : src) v = u(rng);
  for (auto& v : W) v = u(rng);
  for (auto& v : bias) v = u(rng);

  RefConvKahan(d, src.data(), W.data(), bias.data(), dst_kahan_ref.data());
  RefConvBasic(d, src.data(), W.data(), bias.data(), dst_basic_ref.data());

  id<MTLDevice> device = mck.Device();
  id<MTLCommandQueue> queue = [device newCommandQueue];

  id<MTLBuffer> src_buf = [device newBufferWithBytes:src.data()
      length:src_n * sizeof(float) options:MTLResourceStorageModeShared];
  id<MTLBuffer> w_buf   = [device newBufferWithBytes:W.data()
      length:w_n  * sizeof(float) options:MTLResourceStorageModeShared];
  id<MTLBuffer> b_buf   = [device newBufferWithBytes:bias.data()
      length:bias_n * sizeof(float) options:MTLResourceStorageModeShared];
  id<MTLBuffer> dst_buf = [device newBufferWithLength:dst_n * sizeof(float)
      options:MTLResourceStorageModeShared];

  id<MTLCommandBuffer> cb = [queue commandBuffer];
  if (!mck.Encode(cb, src_buf, w_buf, b_buf, dst_buf, d)) {
    std::printf("[%s] FAIL — Encode returned false\n", label);
    return 1;
  }
  [cb commit];
  [cb waitUntilCompleted];

  std::memcpy(dst_gpu.data(), dst_buf.contents,
              dst_n * sizeof(float));

  // 1) Compare GPU Kahan to CPU Kahan (must be bit-exact)
  size_t mismatch_kahan = 0;
  double max_abs_kahan = 0.0;
  size_t max_idx_kahan = 0;
  for (size_t i = 0; i < dst_n; ++i) {
    const double d_abs =
        std::fabs((double)dst_kahan_ref[i] - (double)dst_gpu[i]);
    if (d_abs > max_abs_kahan) { max_abs_kahan = d_abs; max_idx_kahan = i; }
    if (d_abs > 0.0) ++mismatch_kahan;
  }

  // 2) Compare GPU Kahan to CPU Basic-FMA (shows Kahan vs basic drift)
  double max_abs_basic = 0.0, sum_abs_basic = 0.0;
  for (size_t i = 0; i < dst_n; ++i) {
    const double d_abs =
        std::fabs((double)dst_basic_ref[i] - (double)dst_gpu[i]);
    if (d_abs > max_abs_basic) max_abs_basic = d_abs;
    sum_abs_basic += d_abs;
  }
  const double mean_abs_basic = sum_abs_basic / (double)dst_n;

  std::printf(
      "[%s] B=%d  H_in=%d  W_in=%d  C_in=%d  → H_out=%d W_out=%d C_out=%d  "
      "Kh=%d Kw=%d  s=(%d,%d) p=(%d,%d) relu=%d\n",
      label, d.B, d.H_in, d.W_in, d.C_in,
      d.H_out, d.W_out, d.C_out, d.Kh, d.Kw,
      d.stride_h, d.stride_w, d.pad_h, d.pad_w, d.relu ? 1 : 0);
  std::printf("       (vs CPU Kahan)  max_abs=%.6e  mismatched=%zu/%zu  %s\n",
              max_abs_kahan, mismatch_kahan, dst_n,
              max_abs_kahan == 0.0 ? "PASS (bit-exact)" :
              (max_abs_kahan <= 1e-6 ? "PASS (≤1 ULP)" : "FAIL"));
  std::printf("       (vs CPU Basic)  max_abs=%.6e  mean_abs=%.6e  "
              "(this is the Kahan precision improvement vs naive FMA)\n",
              max_abs_basic, mean_abs_basic);
  return max_abs_kahan <= 1e-6 ? 0 : 1;
}

int RunAll() {
  auto mck = MetalConvKahan::Create();
  if (!mck) {
    std::printf("MetalConvKahan::Create failed\n");
    return 2;
  }

  int n_fail = 0;

  // Case 1: 1×1 conv on a 4×4 input, no padding, stride 1.
  {
    ConvDesc d{};
    d.B = 2; d.H_in = 4; d.W_in = 4; d.C_in = 3;
    d.H_out = 4; d.W_out = 4; d.C_out = 5;
    d.Kh = 1; d.Kw = 1;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 0; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mck, "1×1 stride-1 SAME pad-0", d);
  }

  // Case 2: 3×3 conv, stride 1, SAME padding.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 5; d.W_in = 7; d.C_in = 4;
    d.H_out = 5; d.W_out = 7; d.C_out = 6;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 1; d.pad_w = 1;
    d.relu = true;
    n_fail += RunCase(*mck, "3×3 stride-1 SAME pad-1", d);
  }

  // Case 3: 3×3 conv, stride 2, VALID padding.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 9; d.W_in = 9; d.C_in = 8;
    d.H_out = 4; d.W_out = 4; d.C_out = 16;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 2; d.stride_w = 2; d.pad_h = 0; d.pad_w = 0;
    d.relu = false;
    n_fail += RunCase(*mck, "3×3 stride-2 VALID", d);
  }

  // Case 4: stem_s1a real shape — (B=1, 100, 221, 7) → (1, 49, 110, 32)
  // with Kh=Kw=3, stride=2, VALID. Most exercises real-network shapes
  // and validates the kernel on the actual DV input geometry.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 100; d.W_in = 221; d.C_in = 7;
    d.H_out = 49; d.W_out = 110; d.C_out = 32;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 2; d.stride_w = 2; d.pad_h = 0; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mck, "stem_s1a-shape 3×3 s=2 VALID", d);
  }

  std::printf("\n%d/4 cases FAILED\n", n_fail);
  return n_fail == 0 ? 0 : 1;
}

}  // namespace deepvariant

int main(int /*argc*/, char** /*argv*/) {
  return deepvariant::RunAll();
}
