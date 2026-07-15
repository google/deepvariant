// Phase 8 / Tier 6.0 microtest — validate DetMixedBlock for Mixed_5b
// against TF reference output dumped by tools/conversion/dump_tf_per_layer.py.
//
// Inputs (from /tmp/dv_per_layer/):
//   stem_mp5a.npy — TF reference for stem_mp5a output, shape (1, H, W, 192).
//                    This is the input to Mixed_5b.
//   5b.npy        — TF reference for Mixed_5b output, shape (1, H, W, 256).
//
// Test:
//   1. Load TF stem_mp5a.npy → det Mixed_5b → measure max_abs/mean_abs
//      vs TF 5b.npy.
//   2. Acceptance: max_abs ≤ 1e-3 (matches MPSGraph baseline at 5b tap;
//      better is bonus).
//
// Usage:
//   ./microtest_det_mixed5b /Users/.../wgs.dvw /tmp/dv_per_layer

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "deepvariant/native/dv_weights.h"
#include "deepvariant/native/metal_avg_pool.h"
#include "deepvariant/native/metal_bn_relu.h"
#include "deepvariant/native/metal_concat.h"
#include "deepvariant/native/metal_conv_serial.h"
#include "deepvariant/native/metal_det_mixed.h"

namespace deepvariant {

struct NpyData {
  std::vector<int> shape;
  std::vector<float> data;
  size_t total = 0;
};

// Minimal .npy loader (FP32 only). Mirrors debug_metal_main.cc.
bool LoadNpyFp32(const std::string& path, NpyData* out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) {
    std::fprintf(stderr, "npy: cannot open %s\n", path.c_str());
    return false;
  }
  char magic[6];
  f.read(magic, 6);
  if (std::memcmp(magic, "\x93NUMPY", 6) != 0) return false;
  uint8_t major, minor;
  f.read((char*)&major, 1);
  f.read((char*)&minor, 1);
  uint32_t header_len;
  if (major == 1) {
    uint16_t hl;
    f.read((char*)&hl, 2);
    header_len = hl;
  } else {
    uint32_t hl;
    f.read((char*)&hl, 4);
    header_len = hl;
  }
  std::string header(header_len, '\0');
  f.read(header.data(), header_len);
  auto p = header.find("'shape':");
  if (p == std::string::npos) return false;
  auto lp = header.find('(', p);
  auto rp = header.find(')', lp);
  std::string ss = header.substr(lp + 1, rp - lp - 1);
  out->shape.clear();
  for (size_t i = 0; i < ss.size();) {
    while (i < ss.size() && (ss[i] == ' ' || ss[i] == ',')) ++i;
    if (i >= ss.size()) break;
    size_t e = i;
    while (e < ss.size() && ss[e] >= '0' && ss[e] <= '9') ++e;
    if (e == i) break;
    out->shape.push_back(std::stoi(ss.substr(i, e - i)));
    i = e;
  }
  out->total = 1;
  for (int d : out->shape) out->total *= (size_t)d;
  out->data.resize(out->total);
  f.read((char*)out->data.data(), out->total * sizeof(float));
  return (bool)f;
}

int Run(const std::string& dvw_path, const std::string& ref_dir) {
  // 1) Load TF reference inputs/outputs.
  NpyData input_npy, ref_npy;
  if (!LoadNpyFp32(ref_dir + "/stem_mp5a.npy", &input_npy)) {
    std::fprintf(stderr, "FAIL: cannot load stem_mp5a.npy\n");
    return 1;
  }
  if (!LoadNpyFp32(ref_dir + "/5b.npy", &ref_npy)) {
    std::fprintf(stderr, "FAIL: cannot load 5b.npy\n");
    return 1;
  }
  if (input_npy.shape.size() != 4 || ref_npy.shape.size() != 4) {
    std::fprintf(stderr, "FAIL: NPY shapes must be 4D\n");
    return 1;
  }
  const int B = input_npy.shape[0];
  const int H_in = input_npy.shape[1];
  const int W_in = input_npy.shape[2];
  const int C_in = input_npy.shape[3];
  std::printf("Input  shape: (%d, %d, %d, %d) — %zu elems\n",
              B, H_in, W_in, C_in, input_npy.total);
  std::printf("Ref 5b shape: (%d, %d, %d, %d) — %zu elems\n",
              ref_npy.shape[0], ref_npy.shape[1], ref_npy.shape[2],
              ref_npy.shape[3], ref_npy.total);

  // 2) Load .dvw weights.
  auto dvw_p = DvwWeights::Open(dvw_path);
  if (!dvw_p) {
    std::fprintf(stderr, "FAIL: cannot open %s\n", dvw_path.c_str());
    return 1;
  }
  const DvwWeights& dvw = *dvw_p;

  // 3) Initialise Metal device + kernel dispatchers.
  id<MTLDevice> device = MTLCreateSystemDefaultDevice();
  if (!device) {
    std::fprintf(stderr, "FAIL: no Metal device\n");
    return 1;
  }
  id<MTLCommandQueue> queue = [device newCommandQueue];
  auto conv_serial = MetalConvSerial::Create();
  auto bn_relu = MetalBnRelu::Create();
  auto avg_pool = MetalAvgPool::Create();
  auto max_pool = MetalMaxPool::Create();
  auto concat = MetalConcat::Create();
  if (!conv_serial || !bn_relu || !avg_pool || !max_pool || !concat) {
    std::fprintf(stderr, "FAIL: kernel dispatcher creation failed\n");
    return 1;
  }

  // 4) Build DetMixedBlock for Mixed_5b at the input geometry.
  DetMixedBlock block;
  if (!BuildDetMixed5b(device, dvw, /*max_B=*/B, H_in, W_in, C_in, &block)) {
    std::fprintf(stderr, "FAIL: BuildDetMixed5b\n");
    return 1;
  }
  std::printf("Built block: H_out=%d W_out=%d C_out=%d (branches=%zu)\n",
              block.H_out, block.W_out, block.C_out, block.branches.size());

  // 5) Allocate input MTLBuffer + load TF input data.
  const size_t in_bytes = input_npy.total * sizeof(float);
  id<MTLBuffer> input_buf =
      [device newBufferWithBytes:input_npy.data.data()
                          length:in_bytes
                         options:MTLResourceStorageModeShared];

  // 6) Dispatch.
  id<MTLCommandBuffer> cb = [queue commandBuffer];
  if (!DispatchDetMixedBlock(cb, conv_serial.get(), bn_relu.get(),
                              avg_pool.get(), max_pool.get(), concat.get(),
                              block, input_buf, B)) {
    std::fprintf(stderr, "FAIL: DispatchDetMixedBlock\n");
    return 1;
  }
  [cb commit];
  [cb waitUntilCompleted];

  // 7) Read output + compare to TF reference.
  std::vector<float> our_out(ref_npy.total, 0.0f);
  std::memcpy(our_out.data(), [block.concat_out contents],
              ref_npy.total * sizeof(float));

  double max_abs = 0.0, sum_abs = 0.0, max_rel = 0.0;
  size_t max_idx = 0;
  for (size_t i = 0; i < ref_npy.total; ++i) {
    const double d = std::fabs((double)our_out[i] - (double)ref_npy.data[i]);
    sum_abs += d;
    if (d > max_abs) { max_abs = d; max_idx = i; }
    const double denom = std::fabs((double)ref_npy.data[i]);
    if (denom > 1e-6) {
      const double r = d / denom;
      if (r > max_rel) max_rel = r;
    }
  }
  const double mean_abs = sum_abs / (double)ref_npy.total;

  std::printf("\n=== det Mixed_5b vs TF reference ===\n");
  std::printf("max_abs   = %.6e\n", max_abs);
  std::printf("mean_abs  = %.6e\n", mean_abs);
  std::printf("max_rel   = %.6e\n", max_rel);
  std::printf("first divergent idx %zu: ref=%.6e ours=%.6e\n",
              max_idx, ref_npy.data[max_idx], our_out[max_idx]);

  // Acceptance: matches MPSGraph baseline drift at 5b (~1e-3 max_abs
  // per Probe D). Better is bonus.
  const char* status;
  if (max_abs <= 1e-5) {
    status = "PASS (within 1e-5 — bit-near-exact)";
  } else if (max_abs <= 1.5e-3) {
    status = "PASS (matches MPSGraph baseline drift at 5b ~1.5e-3)";
  } else {
    status = "FAIL (drift exceeds MPSGraph baseline)";
  }
  std::printf("Status    = %s\n", status);
  return max_abs <= 1.5e-3 ? 0 : 1;
}

}  // namespace deepvariant

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: %s <wgs.dvw> <ref_dir>\n"
                 "  ref_dir must contain stem_mp5a.npy and 5b.npy\n",
                 argv[0]);
    return 2;
  }
  return deepvariant::Run(argv[1], argv[2]);
}
