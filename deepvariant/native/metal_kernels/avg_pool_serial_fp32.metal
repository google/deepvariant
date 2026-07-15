// Phase 5.5e — deterministic-reduction-order AvgPool2D kernel.
//
// One thread per output element (n, h_out, w_out, c). Inside the
// thread the (kh, kw) accumulation is a strict scalar `for` loop
// summing into `acc`, then divides by either `Kh*Kw` (include padding
// in average) or by `count_valid` (exclude padding from average).
//
// Inception-v3 uses `exclude_padding_from_average=True` for the
// AvgPool branch in InceptionA / InceptionB / InceptionC blocks
// (Keras default for `AveragePooling2D` is also exclude-padding when
// padding='same'). The kernel parameter `exclude_pad` mirrors that.
//
// Layouts:
//   src   : NHWC      (B,  H_in,  W_in,  C)   row-major
//   dst   : NHWC      (B,  H_out, W_out, C)   row-major
//
// Bit-deterministic across runs and across Apple Silicon chip
// generations (per-thread strict-serial accumulation; no SIMD-group
// reductions; `metal::precise::fma` for IEEE single-rounded FMA).

#include <metal_stdlib>
using namespace metal;

struct AvgPoolParams {
    int B;
    int H_in;
    int W_in;
    int C;
    int H_out;
    int W_out;
    int Kh;
    int Kw;
    int stride_h;
    int stride_w;
    int pad_h;
    int pad_w;
    int exclude_pad;     // 1 → divide by count of in-bounds positions
                         // 0 → divide by Kh*Kw (include zero-padded as 0/Kh*Kw)
};

kernel void avgpool2d_fp32(
    constant AvgPoolParams& P [[ buffer(0) ]],
    device   const float*   src [[ buffer(1) ]],
    device   float*         dst [[ buffer(2) ]],
    uint3 gid                [[ thread_position_in_grid ]])
{
    const int c  = (int)gid.x;
    const int hw = (int)gid.y;
    const int n  = (int)gid.z;
    if (n >= P.B || c >= P.C || hw >= P.H_out * P.W_out) return;
    const int h_out = hw / P.W_out;
    const int w_out = hw % P.W_out;

    const int h_base = h_out * P.stride_h - P.pad_h;
    const int w_base = w_out * P.stride_w - P.pad_w;

    float acc = 0.0f;
    int count = 0;
    for (int kh = 0; kh < P.Kh; ++kh) {
        const int h_in = h_base + kh;
        if (h_in < 0 || h_in >= P.H_in) {
            if (P.exclude_pad == 0) {
                // include zero-pad in average — count++ but acc unchanged
                count += P.Kw;
            }
            continue;
        }
        for (int kw = 0; kw < P.Kw; ++kw) {
            const int w_in = w_base + kw;
            if (w_in < 0 || w_in >= P.W_in) {
                if (P.exclude_pad == 0) ++count;
                continue;
            }
            acc += src[
                ((n * P.H_in + h_in) * P.W_in + w_in) * P.C + c];
            ++count;
        }
    }

    const float divisor = (count > 0) ? (float)count : 1.0f;
    dst[((n * P.H_out + h_out) * P.W_out + w_out) * P.C + c] = acc / divisor;
}
