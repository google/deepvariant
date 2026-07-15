// Phase 5.5e — channel-axis concat for NHWC FP32 tensors.
//
// One thread per output element (n, h, w, c_out). Pure data movement:
// each thread reads from the appropriate input branch and writes to
// dst. NO reductions, so deterministic by construction.
//
// Up to 4 input branches (matches Inception-v3 Mixed_5b/c/d/6b/c/d/e/
// 7b/c which all have ≤ 4 branches). For Reduction-A/B (3 branches)
// the last input is unused (set offset to negative).
//
// Layouts:
//   src_i : NHWC (B, H, W, C_i)   row-major
//   dst   : NHWC (B, H, W, sum(C_i))  row-major
//
// `c_offset[i]` is the starting channel index for branch i in the
// output. `c_size[i]` is the channel count of branch i. The thread
// for output (n, h, w, c_out) finds which branch owns c_out by
// linear search across the 4 offsets — O(4) at most, unrolled.

#include <metal_stdlib>
using namespace metal;

struct ConcatParams {
    int B;
    int H;
    int W;
    int n_branches;     // 1 .. 4
    int c_size_0;
    int c_size_1;
    int c_size_2;
    int c_size_3;
    int c_total;        // sum of c_size_*
};

kernel void concat_channels_fp32(
    constant ConcatParams& P [[ buffer(0) ]],
    device   const float*  src0 [[ buffer(1) ]],
    device   const float*  src1 [[ buffer(2) ]],
    device   const float*  src2 [[ buffer(3) ]],
    device   const float*  src3 [[ buffer(4) ]],
    device   float*        dst  [[ buffer(5) ]],
    uint3 gid                 [[ thread_position_in_grid ]])
{
    const int c_out = (int)gid.x;
    const int hw    = (int)gid.y;
    const int n     = (int)gid.z;
    if (n >= P.B || c_out >= P.c_total || hw >= P.H * P.W) return;
    const int h = hw / P.W;
    const int w = hw % P.W;

    // Find which branch owns c_out + the local channel index.
    int b = 0;
    int c_local = c_out;
    int c_size = P.c_size_0;
    if (c_local < c_size) {
        b = 0;
    } else {
        c_local -= c_size;
        c_size = P.c_size_1;
        if (c_local < c_size) {
            b = 1;
        } else {
            c_local -= c_size;
            c_size = P.c_size_2;
            if (c_local < c_size) {
                b = 2;
            } else {
                c_local -= c_size;
                b = 3;
            }
        }
    }

    float v;
    const int hw_off = (n * P.H + h) * P.W + w;
    switch (b) {
        case 0: v = src0[hw_off * P.c_size_0 + c_local]; break;
        case 1: v = src1[hw_off * P.c_size_1 + c_local]; break;
        case 2: v = src2[hw_off * P.c_size_2 + c_local]; break;
        default: v = src3[hw_off * P.c_size_3 + c_local]; break;
    }
    dst[hw_off * P.c_total + c_out] = v;
}
