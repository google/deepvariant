// Phase 5.5c — deterministic-reduction-order Conv2D + ReLU kernel.
//
// One thread per output element (n, h_out, w_out, c_out). Inside the
// thread the (kh, kw, c_in) accumulation is a strict scalar `for`
// loop using `metal::precise::fma(x, w, acc)` — single-rounding,
// IEEE 754 fused multiply-add, identical bit pattern to Eigen+AVX-512
// FMA on x86. No SIMD-group reduction, no atomics, no `mad`-style
// fast-math contraction.
//
// Layouts:
//   src   : NHWC      (B,  H_in,  W_in,  C_in)   row-major
//   W     : HWIO      (Kh, Kw,  C_in,  C_out)    row-major
//   bias  : (C_out,)
//   dst   : NHWC      (B,  H_out, W_out, C_out)  row-major
//
// Padding follows the standard "explicit pad" model: pad_h / pad_w
// pre-computed by the host (host emulates SAME or VALID).
//
// Built with the embedded Metal toolchain at host build time; loaded
// from the dv_metal_kernels.metallib produced by CMake.

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

kernel void conv_serial_fp32(
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

    // Strict (kh, kw, c_in)-order scalar accumulation. metal::precise::fma
    // emits IEEE 754 single-rounded fused multiply-add, matching Eigen's
    // fma() path on x86-AVX-512 bit-for-bit.
    float acc = 0.0f;
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
                acc = metal::precise::fma(x, w, acc);
            }
        }
    }

    acc += bias[c_out];
    if (P.relu != 0) acc = max(acc, 0.0f);

    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C_out + c_out] = acc;
}
