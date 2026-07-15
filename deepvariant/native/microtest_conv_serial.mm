// Phase 5.5c microtest: dispatch the deterministic Conv2D kernel on a
// small known case, compare against a CPU reference implementing the
// same scalar (kh, kw, c_in)-order accumulation. Bit-identical match
// expected (PASS) — any divergence means the kernel's reduction order
// or padding handling deviates from the CPU spec.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "deepvariant/native/metal_conv_serial.h"

namespace deepvariant {

// Reference scalar conv (matches the kernel exactly: NHWC src, HWIO W,
// (kh, kw, c_in)-order, FMA via std::fma, optional ReLU).
void RefConv(const ConvDesc& d,
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

int RunCase(MetalConvSerial& mcs, const char* label, const ConvDesc& d) {
  const size_t src_n =
      (size_t)d.B * d.H_in * d.W_in * d.C_in;
  const size_t w_n = (size_t)d.Kh * d.Kw * d.C_in * d.C_out;
  const size_t bias_n = d.C_out;
  const size_t dst_n =
      (size_t)d.B * d.H_out * d.W_out * d.C_out;

  std::mt19937 rng(0x55c1);
  std::uniform_real_distribution<float> u(-1.0f, 1.0f);
  std::vector<float> src(src_n), W(w_n), bias(bias_n);
  std::vector<float> dst_ref(dst_n, 0.0f), dst_gpu(dst_n, 0.0f);
  for (auto& v : src) v = u(rng);
  for (auto& v : W) v = u(rng);
  for (auto& v : bias) v = u(rng);

  RefConv(d, src.data(), W.data(), bias.data(), dst_ref.data());

  id<MTLDevice> device = mcs.Device();
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
  if (!mcs.Encode(cb, src_buf, w_buf, b_buf, dst_buf, d)) {
    std::printf("[%s] FAIL — Encode returned false\n", label);
    return 1;
  }
  [cb commit];
  [cb waitUntilCompleted];

  std::memcpy(dst_gpu.data(), dst_buf.contents,
              dst_n * sizeof(float));

  size_t mismatch = 0;
  double max_abs = 0.0;
  size_t max_idx = 0;
  for (size_t i = 0; i < dst_n; ++i) {
    const double d_abs =
        std::fabs((double)dst_ref[i] - (double)dst_gpu[i]);
    if (d_abs > max_abs) { max_abs = d_abs; max_idx = i; }
    if (d_abs > 0.0) ++mismatch;
  }

  std::printf(
      "[%s] B=%d  H_in=%d  W_in=%d  C_in=%d  → H_out=%d W_out=%d C_out=%d  "
      "Kh=%d Kw=%d  s=(%d,%d) p=(%d,%d) relu=%d\n",
      label, d.B, d.H_in, d.W_in, d.C_in,
      d.H_out, d.W_out, d.C_out, d.Kh, d.Kw,
      d.stride_h, d.stride_w, d.pad_h, d.pad_w, d.relu ? 1 : 0);
  std::printf("       n_elems=%zu  max_abs=%.6e  mismatched=%zu/%zu  "
              "first_diff_idx=%zu  ref=%.6e  gpu=%.6e  %s\n",
              dst_n, max_abs, mismatch, dst_n, max_idx,
              dst_ref[max_idx], dst_gpu[max_idx],
              max_abs == 0.0 ? "PASS (bit-exact)" :
              (max_abs <= 1e-6 ? "PASS (≤1 ULP)" : "FAIL"));
  return max_abs <= 1e-6 ? 0 : 1;
}

int RunAll() {
  auto mcs = MetalConvSerial::Create();
  if (!mcs) {
    std::printf("MetalConvSerial::Create failed\n");
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
    n_fail += RunCase(*mcs, "1×1 stride-1 SAME pad-0", d);
  }

  // Case 2: 3×3 conv, stride 1, SAME padding.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 5; d.W_in = 7; d.C_in = 4;
    d.H_out = 5; d.W_out = 7; d.C_out = 6;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 1; d.pad_w = 1;
    d.relu = true;
    n_fail += RunCase(*mcs, "3×3 stride-1 SAME pad-1", d);
  }

  // Case 3: 3×3 conv, stride 2, VALID padding.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 9; d.W_in = 9; d.C_in = 8;
    d.H_out = 4; d.W_out = 4; d.C_out = 16;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 2; d.stride_w = 2; d.pad_h = 0; d.pad_w = 0;
    d.relu = false;
    n_fail += RunCase(*mcs, "3×3 stride-2 VALID", d);
  }

  // Case 4: stem_s1a real shape — (B=1, 100, 221, 7) → (1, 49, 110, 32)
  // with Kh=Kw=3, stride=2, VALID. Most exercises real-network shapes.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 100; d.W_in = 221; d.C_in = 7;
    d.H_out = 49; d.W_out = 110; d.C_out = 32;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 2; d.stride_w = 2; d.pad_h = 0; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mcs, "stem_s1a-shape 3×3 s=2 VALID", d);
  }

  // ===== Tier 6.0 shape coverage: Inception blocks 5b–7c =====
  // Validate MetalConvSerial supports every conv shape used by the
  // 11 Mixed_X blocks. If any case FAILs, the conv_serial-full-network
  // refactor cannot proceed without a kernel-side fix.

  // Case 5: 5×5 stride-1 SAME (Mixed_5b/5c/5d branch5x5).
  // E.g. Mixed_5b: 48 → 64 channels, padded to maintain spatial size.
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 23; d.W_in = 53; d.C_in = 48;
    d.H_out = 23; d.W_out = 53; d.C_out = 64;
    d.Kh = 5; d.Kw = 5;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 2; d.pad_w = 2;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_5b-shape 5×5 s=1 SAME 48→64", d);
  }

  // Case 6: 7×1 stride-1 SAME (Mixed_6b–6e branch7x7 asymmetric).
  // E.g. Mixed_6b: 128 → 128 channels with kernel (7,1).
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 11; d.W_in = 26; d.C_in = 128;
    d.H_out = 11; d.W_out = 26; d.C_out = 128;
    d.Kh = 7; d.Kw = 1;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 3; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_6b-shape 7×1 s=1 SAME 128→128", d);
  }

  // Case 7: 1×7 stride-1 SAME (Mixed_6b–6e branch7x7 asymmetric).
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 11; d.W_in = 26; d.C_in = 128;
    d.H_out = 11; d.W_out = 26; d.C_out = 192;
    d.Kh = 1; d.Kw = 7;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 0; d.pad_w = 3;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_6b-shape 1×7 s=1 SAME 128→192", d);
  }

  // Case 8: 1×3 stride-1 SAME (Mixed_7b/7c branch3x3 asymmetric split).
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 5; d.W_in = 12; d.C_in = 384;
    d.H_out = 5; d.W_out = 12; d.C_out = 384;
    d.Kh = 1; d.Kw = 3;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 0; d.pad_w = 1;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_7b-shape 1×3 s=1 SAME 384→384", d);
  }

  // Case 9: 3×1 stride-1 SAME (Mixed_7b/7c branch3x3 asymmetric split).
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 5; d.W_in = 12; d.C_in = 384;
    d.H_out = 5; d.W_out = 12; d.C_out = 384;
    d.Kh = 3; d.Kw = 1;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 1; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_7b-shape 3×1 s=1 SAME 384→384", d);
  }

  // Case 10: 1×1 stride-1 SAME, large channels (Mixed_7c branch1x1 320→320).
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 5; d.W_in = 12; d.C_in = 1280;
    d.H_out = 5; d.W_out = 12; d.C_out = 320;
    d.Kh = 1; d.Kw = 1;
    d.stride_h = 1; d.stride_w = 1; d.pad_h = 0; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_7c-shape 1×1 s=1 1280→320", d);
  }

  // Case 11: 3×3 stride-2 VALID (Mixed_6a/7a reduction blocks).
  {
    ConvDesc d{};
    d.B = 1; d.H_in = 23; d.W_in = 53; d.C_in = 96;
    d.H_out = 11; d.W_out = 26; d.C_out = 96;
    d.Kh = 3; d.Kw = 3;
    d.stride_h = 2; d.stride_w = 2; d.pad_h = 0; d.pad_w = 0;
    d.relu = true;
    n_fail += RunCase(*mcs, "Mixed_6a-shape 3×3 s=2 VALID 96→96", d);
  }

  std::printf("\n%d/11 cases FAILED\n", n_fail);
  return n_fail == 0 ? 0 : 1;
}

}  // namespace deepvariant

int main(int /*argc*/, char** /*argv*/) {
  return deepvariant::RunAll();
}
