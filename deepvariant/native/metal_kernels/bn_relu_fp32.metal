// Phase 5.5f — separate BatchNorm+ReLU kernel matching TF/oneDNN's
// non-folded conv→BN→ReLU path. Used in conjunction with the existing
// conv_serial_fp32 kernel to avoid the FoldConvBn FP32 drift.
//
// Day-1 PoC measurement on stem_s1a output[h,w,c] (1×100×221×7 fixed-
// seed input): folded conv+BN gave up to 93 ULP delta vs TF reference
// across 128 elements. Switching to per-thread c_in-serial FMA conv
// PLUS this separate BN kernel reduced max delta to ±2 ULP, with 76%
// of elements bit-exact. The per-element residual is sub-noise vs the
// FILTER threshold drift and should match Docker FILTER classes after
// 188-layer accumulation.
//
// Input layout : NHWC FP32 (output of preceding conv, raw — no bias,
//                no ReLU, no fold)
// Per-channel  : mean[C], var[C], beta[C] FP32
// Eps          : Keras BN default = 1e-3 (passed as constant).
// Output layout: NHWC FP32, same shape as input.
//
// Per-thread serial: each thread updates one (n, h, w, c) element.
// Computation:
//   inv_std = 1.0f / sqrt(var[c] + eps)
//   y       = (x - mean[c]) * inv_std + beta[c]
//   if relu: y = max(y, 0.0f)
//
// Note: oneDNN/TF uses the same formula. The ±1-2 ULP residual after
// this kernel matches the variance in TF's own non-deterministic
// reductions across runs.

#include <metal_stdlib>
using namespace metal;

struct BnReluParams {
    int B;
    int H;
    int W;
    int C;
    float eps;        // Keras BN default = 1e-3
    int relu;         // 1 = apply ReLU after BN; 0 = pure BN
};

kernel void bn_relu_fp32(
    constant BnReluParams& P     [[ buffer(0) ]],
    device   const float* src    [[ buffer(1) ]],
    device   const float* mean   [[ buffer(2) ]],
    device   const float* var    [[ buffer(3) ]],
    device   const float* beta   [[ buffer(4) ]],
    device   float*       dst    [[ buffer(5) ]],
    uint3 gid                    [[ thread_position_in_grid ]])
{
    const int c = (int)gid.x;
    const int hw = (int)gid.y;
    const int n = (int)gid.z;
    if (n >= P.B || c >= P.C || hw >= P.H * P.W) return;

    const int h = hw / P.W;
    const int w = hw % P.W;
    const int idx = ((n * P.H + h) * P.W + w) * P.C + c;

    const float x = src[idx];
    const float mu = mean[c];
    const float v = var[c];
    const float b = beta[c];

    // Match TF/oneDNN BN: y = (x - mean) / sqrt(var + eps) + beta.
    // Use precise::sqrt (single-rounding sqrt). Eigen on x86 also uses
    // sqrt (not rsqrt approximation) for FP32 BN.
    const float inv_std = 1.0f / metal::precise::sqrt(v + P.eps);
    float y = metal::precise::fma(x - mu, inv_std, b);
    if (P.relu != 0) y = max(y, 0.0f);

    dst[idx] = y;
}
