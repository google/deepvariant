// Phase 5.5e — deterministic global-avg-pool kernel.
//
// Inception-v3 ends with a global avg pool over the last spatial
// volume (8×8 = 64 elements per channel for the 100-row input;
// dimensions parameterised here). One thread per output element
// (n, c). Inside: scalar `for (h, w)` summing into `acc`, divide by
// `H_in * W_in` at the end.
//
// Layout:
//   src   : NHWC (B, H_in, W_in, C)   row-major
//   dst   : (B, C)                     row-major
//
// Bit-deterministic: per-thread strict-serial accumulation, no SIMD
// reductions. Output matches `np.mean(x, axis=(1, 2))` bit-for-bit
// when summed in the same order.

#include <metal_stdlib>
using namespace metal;

struct GlobalAvgPoolParams {
    int B;
    int H_in;
    int W_in;
    int C;
};

kernel void global_avg_pool_fp32(
    constant GlobalAvgPoolParams& P [[ buffer(0) ]],
    device   const float*         src [[ buffer(1) ]],
    device   float*               dst [[ buffer(2) ]],
    uint2 gid                       [[ thread_position_in_grid ]])
{
    const int c = (int)gid.x;
    const int n = (int)gid.y;
    if (n >= P.B || c >= P.C) return;

    float acc = 0.0f;
    for (int h = 0; h < P.H_in; ++h) {
        for (int w = 0; w < P.W_in; ++w) {
            acc += src[((n * P.H_in + h) * P.W_in + w) * P.C + c];
        }
    }
    const int n_elems = P.H_in * P.W_in;
    dst[n * P.C + c] = acc / (float)n_elems;
}
