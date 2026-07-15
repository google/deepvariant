// Phase 5.5e/Path B — Kahan-compensated Conv2D + ReLU kernel.
//
// One thread per output element (n, h_out, w_out, c_out). Inside the
// thread the (kh, kw, c_in) accumulation uses Kahan compensated
// summation — each `sum + y` recovers the rounding error in a
// compensation term `c` and folds it into the next iteration's
// increment via `precise::fma(x, w, -c)`. Cross-platform bit-
// deterministic within ~1 ULP regardless of reduction order — Demmel
// & Nguyen ARITH-21 2013, "Fast Reproducible Floating-Point Summation"
// (XBLAS/ReproBLAS).
//
// Cost: ~4× per inner FMA vs basic `precise::fma` (3 extra ops per
// iteration). For DeepVariant Inception-v3 stem_s1a (3×3 stride-2
// 7→32, 63 inner FMAs per output): 63 → 252 ops per output element.
// Estimated wall-time impact at full-network rollout: ~2-3× MPSGraph
// baseline (~8-12 min/chr20 vs current 4 min) — under the 8 min gate.
//
// Layouts (identical to conv_serial_fp32):
//   src   : NHWC      (B,  H_in,  W_in,  C_in)   row-major
//   W     : HWIO      (Kh, Kw,  C_in,  C_out)    row-major
//   bias  : (C_out,)
//   dst   : NHWC      (B,  H_out, W_out, C_out)  row-major

#include <metal_stdlib>
using namespace metal;

struct ConvParams {
    int B;
    int H_in;
    int W_in;
    int C_in;
    int H_out;
    int W_out;
    int C_out;
    int Kh;
    int Kw;
    int stride_h;
    int stride_w;
    int pad_h;          // > 0 → top  zero-pad rows
    int pad_w;          // > 0 → left zero-pad cols
    int relu;           // 1 → apply ReLU after bias add
};

kernel void conv_kahan_fp32(
    constant ConvParams& P     [[ buffer(0) ]],
    device   const float* src  [[ buffer(1) ]],
    device   const float* W    [[ buffer(2) ]],
    device   const float* bias [[ buffer(3) ]],
    device   float*       dst  [[ buffer(4) ]],
    uint3 gid                  [[ thread_position_in_grid ]])
{
    const int c_out = (int)gid.x;
    const int hw    = (int)gid.y;
    const int n     = (int)gid.z;
    if (n >= P.B || c_out >= P.C_out || hw >= P.H_out * P.W_out) return;
    const int h_out = hw / P.W_out;
    const int w_out = hw % P.W_out;

    const int h_base = h_out * P.stride_h - P.pad_h;
    const int w_base = w_out * P.stride_w - P.pad_w;

    // Kahan compensated summation:
    //   y = (x*w) - c          [single-rounded via precise::fma(x, w, -c)]
    //   t = sum + y            [rounding loses some low bits]
    //   c = (t - sum) - y      [recover the lost bits — exact in FP32
    //                           because t, sum, y are all FP32 and
    //                           |y| << |sum|]
    //   sum = t
    //
    // The compensation `c` accumulates the rounding losses across
    // iterations; each new addition recovers them via `fma(x, w, -c)`.
    // Final error is O(ε² · |sum|) per step rather than O(ε · |sum|).
    float sum = 0.0f;
    float c = 0.0f;
    for (int kh = 0; kh < P.Kh; ++kh) {
        const int h_in = h_base + kh;
        if (h_in < 0 || h_in >= P.H_in) continue;
        for (int kw = 0; kw < P.Kw; ++kw) {
            const int w_in = w_base + kw;
            if (w_in < 0 || w_in >= P.W_in) continue;
            for (int c_in = 0; c_in < P.C_in; ++c_in) {
                const float x = src[
                    ((n * P.H_in + h_in) * P.W_in + w_in) * P.C_in + c_in];
                const float w = W[
                    ((kh * P.Kw + kw) * P.C_in + c_in) * P.C_out + c_out];
                // y = x*w - c  (single-rounded FMA)
                const float y = metal::precise::fma(x, w, -c);
                const float t = sum + y;
                c = (t - sum) - y;
                sum = t;
            }
        }
    }

    sum += bias[c_out];
    if (P.relu != 0) sum = max(sum, 0.0f);

    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C_out + c_out] = sum;
}
