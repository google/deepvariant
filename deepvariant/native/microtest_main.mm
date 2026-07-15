// Phase 5.5a investigation: hand-verifiable MPSGraph conv micro-tests.
//
// Builds tiny self-contained MPSGraphs (no .dvw, no full network) with
// inputs and weights small enough to compute the expected output by
// pencil-and-paper. If MPSGraph fails the trivial 1×1 case, the bug is
// in matmul/conv itself. If 1×1 passes but 3×3 fails, the bug is in
// spatial / imToCol handling.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShadersGraph/MetalPerformanceShadersGraph.h>
#import <MetalPerformanceShadersGraph/MPSGraphImToColOps.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

// Print a float vector vs expected, return PASS/FAIL on max-abs ≤ tol.
bool CompareVec(const char* label, const std::vector<float>& got,
                const std::vector<float>& expected, float tol) {
  if (got.size() != expected.size()) {
    std::printf("  %s: SIZE mismatch (got %zu, expected %zu)\n",
                label, got.size(), expected.size());
    return false;
  }
  float max_abs = 0.0f;
  for (size_t i = 0; i < got.size(); ++i) {
    max_abs = std::max(max_abs, std::fabs(got[i] - expected[i]));
  }
  std::printf("  got      :");
  for (size_t i = 0; i < got.size() && i < 16; ++i) {
    std::printf(" %9.3f", got[i]);
  }
  if (got.size() > 16) std::printf(" ...");
  std::printf("\n  expected :");
  for (size_t i = 0; i < expected.size() && i < 16; ++i) {
    std::printf(" %9.3f", expected[i]);
  }
  if (expected.size() > 16) std::printf(" ...");
  std::printf("\n  max-abs  : %.6e   verdict: %s\n",
              max_abs, max_abs <= tol ? "PASS" : "FAIL");
  std::fflush(stdout);
  return max_abs <= tol;
}

// Run a graph that takes one input and produces one output. Returns
// the output as a flat float vector.
bool RunGraph(MPSGraph* g, MPSGraphTensor* input, MPSGraphTensor* output,
              NSArray<NSNumber*>* in_shape, const float* in_data,
              std::vector<float>* out_buf,
              NSArray<NSNumber*>** out_shape_ret) {
  id<MTLDevice> device = MTLCreateSystemDefaultDevice();
  if (!device) return false;
  id<MTLCommandQueue> queue = [device newCommandQueue];

  MPSGraphCompilationDescriptor* desc = [MPSGraphCompilationDescriptor new];
  desc.optimizationLevel = MPSGraphOptimizationLevel0;
  desc.waitForCompilationCompletion = YES;
  MPSGraphShapedType* in_st = [[MPSGraphShapedType alloc]
      initWithShape:in_shape dataType:MPSDataTypeFloat32];
  MPSGraphExecutable* exe =
      [g compileWithDevice:[MPSGraphDevice deviceWithMTLDevice:device]
                     feeds:@{input: in_st}
             targetTensors:@[output]
          targetOperations:nil
     compilationDescriptor:desc];
  if (!exe) {
    std::printf("  compile FAILED\n");
    std::fflush(stdout);
    return false;
  }

  NSUInteger n_in = 1;
  for (NSNumber* d in in_shape) n_in *= [d unsignedIntegerValue];
  NSData* in_nsdata = [NSData dataWithBytes:in_data
                                     length:n_in * sizeof(float)];
  MPSGraphTensorData* in_td = [[MPSGraphTensorData alloc]
      initWithDevice:[MPSGraphDevice deviceWithMTLDevice:device]
                data:in_nsdata
               shape:in_shape
            dataType:MPSDataTypeFloat32];
  MPSGraphExecutableExecutionDescriptor* runDesc =
      [MPSGraphExecutableExecutionDescriptor new];
  runDesc.waitUntilCompleted = YES;
  NSArray<MPSGraphTensorData*>* outs =
      [exe runWithMTLCommandQueue:queue
                      inputsArray:@[in_td]
                     resultsArray:nil
              executionDescriptor:runDesc];
  if (!outs || outs.count != 1) return false;
  *out_shape_ret = outs[0].shape;
  NSUInteger total = 1;
  for (NSNumber* d in outs[0].shape) total *= [d unsignedIntegerValue];
  out_buf->resize(total);
  [outs[0].mpsndarray readBytes:out_buf->data() strideBytes:nil];
  return true;
}

void PrintShape(NSArray<NSNumber*>* shape) {
  std::printf("  out shape: ");
  for (NSNumber* d in shape) {
    std::printf("%lu ", (unsigned long)[d unsignedIntegerValue]);
  }
  std::printf("\n");
}

// ===========================================================================
// Test 1: trivial 1×1 conv (PASSED — kept as smoke test)
// ===========================================================================

void Test1_Conv1x1() {
  std::printf("\n=== Test 1: 1×1 conv via convolution2DWithSourceTensor ===\n");
  std::fflush(stdout);
  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape = @[@1, @1, @1, @2];
    MPSGraphTensor* x = [g placeholderWithShape:in_shape
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    float w[6] = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f};
    NSData* w_data = [NSData dataWithBytes:w length:sizeof(w)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@1, @1, @2, @3]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:1 strideInY:1
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    float in_data[2] = {3.f, 5.f};
    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y, in_shape, in_data, &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    CompareVec("conv1x1", got, {23.f, 31.f, 39.f}, 1e-4f);
  }
}

// ===========================================================================
// Test 2: 3×3 single-channel conv on 3×3 input.
// Input:  1..9 NHWC
// Weight: 1..9 HWIO
// VALID stride 1 → out (1,1,1,1) value = sum(i*i for i in 1..9) = 285
// ===========================================================================

void Test2_Conv3x3SingleCh() {
  std::printf("\n=== Test 2: 3×3 conv 1→1 ===\n");
  std::fflush(stdout);
  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape = @[@1, @3, @3, @1];
    MPSGraphTensor* x = [g placeholderWithShape:in_shape
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    float w[9];
    for (int i = 0; i < 9; ++i) w[i] = (float)(i + 1);
    NSData* w_data = [NSData dataWithBytes:w length:sizeof(w)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@3, @3, @1, @1]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:1 strideInY:1
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    float in_data[9];
    for (int i = 0; i < 9; ++i) in_data[i] = (float)(i + 1);
    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y, in_shape, in_data, &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    CompareVec("conv3x3_1to1", got, {285.f}, 1e-2f);
  }
}

// ===========================================================================
// Test 3: 3×3 conv 7→1, single output position (kernel matches stem_s1a's
// shape but output channel = 1).
//
// Input  (1, 3, 3, 7) — values input[h, w, c] = (h*3 + w)*7 + c + 1
// Weight (3, 3, 7, 1) HWIO — same numbering 1..63
// VALID stride 1 → output (1, 1, 1, 1) = sum_{i=1..63} i*i = 85344
// ===========================================================================

void Test3_Conv3x3_7to1() {
  std::printf("\n=== Test 3: 3×3 conv 7→1 ===\n");
  std::fflush(stdout);
  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape = @[@1, @3, @3, @7];
    MPSGraphTensor* x = [g placeholderWithShape:in_shape
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    float w[3 * 3 * 7];
    for (int i = 0; i < 63; ++i) w[i] = (float)(i + 1);
    NSData* w_data = [NSData dataWithBytes:w length:sizeof(w)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@3, @3, @7, @1]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:1 strideInY:1
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    float in_data[3 * 3 * 7];
    for (int i = 0; i < 63; ++i) in_data[i] = (float)(i + 1);
    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y, in_shape, in_data, &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    CompareVec("conv3x3_7to1", got, {85344.f}, 1.f);
  }
}

// ===========================================================================
// Test 4: 3×3 conv 7→32 (full stem_s1a kernel shape, single position).
//
// Weight[h, w, c, o] = (h*3 + w)*7 + c + 1 + o*0.001
// Expected[o] = sum_{i=1..63} i * (i + o*0.001) = 85344 + o * 0.001 * 2016
// ===========================================================================

void Test4_Conv3x3_7to32() {
  std::printf("\n=== Test 4: 3×3 conv 7→32 (matches stem_s1a shape) ===\n");
  std::fflush(stdout);
  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape = @[@1, @3, @3, @7];
    MPSGraphTensor* x = [g placeholderWithShape:in_shape
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    std::vector<float> w(3 * 3 * 7 * 32);
    for (int h = 0; h < 3; ++h)
      for (int wj = 0; wj < 3; ++wj)
        for (int c = 0; c < 7; ++c)
          for (int o = 0; o < 32; ++o) {
            float v = (float)((h * 3 + wj) * 7 + c + 1) + (float)o * 0.001f;
            w[((h * 3 + wj) * 7 + c) * 32 + o] = v;
          }
    NSData* w_data = [NSData dataWithBytes:w.data()
                                     length:w.size() * sizeof(float)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@3, @3, @7, @32]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:1 strideInY:1
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    float in_data[3 * 3 * 7];
    for (int i = 0; i < 63; ++i) in_data[i] = (float)(i + 1);
    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y, in_shape, in_data, &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    std::vector<float> expected(32);
    for (int o = 0; o < 32; ++o) {
      expected[o] = 85344.f + (float)o * 0.001f * 2016.f;
    }
    CompareVec("conv3x3_7to32", got, expected, 1.f);
  }
}

// ===========================================================================
// Test 5: same as Test 4 but at stride 2, 4×4 input → output (1,1,1,32).
// This is the closest analogue of stem_s1a (3×3 stride 2 valid).
// ===========================================================================

void Test5_Conv3x3_S2_7to32() {
  std::printf("\n=== Test 5: 3×3 stride-2 valid conv 7→32 (stem_s1a-like) ===\n");
  std::fflush(stdout);
  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape = @[@1, @3, @3, @7];  // produces 1x1 output at stride 2 valid
    MPSGraphTensor* x = [g placeholderWithShape:in_shape
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    std::vector<float> w(3 * 3 * 7 * 32);
    for (int h = 0; h < 3; ++h)
      for (int wj = 0; wj < 3; ++wj)
        for (int c = 0; c < 7; ++c)
          for (int o = 0; o < 32; ++o) {
            float v = (float)((h * 3 + wj) * 7 + c + 1) + (float)o * 0.001f;
            w[((h * 3 + wj) * 7 + c) * 32 + o] = v;
          }
    NSData* w_data = [NSData dataWithBytes:w.data()
                                     length:w.size() * sizeof(float)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@3, @3, @7, @32]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:2 strideInY:2
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    float in_data[3 * 3 * 7];
    for (int i = 0; i < 63; ++i) in_data[i] = (float)(i + 1);
    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y, in_shape, in_data, &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    // Same expected as Test 4 (only 1 output position, stride doesn't matter)
    std::vector<float> expected(32);
    for (int o = 0; o < 32; ++o) {
      expected[o] = 85344.f + (float)o * 0.001f * 2016.f;
    }
    CompareVec("conv3x3_s2_7to32", got, expected, 1.f);
  }
}

// ===========================================================================
// Test 6: same as Test 5 but on a LARGE input (100×221) — exactly the
// stem_s1a input size. Stride 2 valid 3×3, 7→32 channels.
//
// We compute expected output[h_out=0, w_out=0, o=0] by hand from a
// known-pattern input: input[h, w, c] = 1.0 if (h, w, c) == (1, 1, 0),
// else 0. Then output[0, 0, 0, 0] = weight[1, 1, 0, 0] * 1.0 = 1+0.001*0
// = 1.0. All other outputs at (h_out=0, w_out=0) = weight[1, 1, 0, o]
// = 1 + o*0.001.
// ===========================================================================

void Test6_Conv3x3_S2_7to32_LargeInput() {
  std::printf("\n=== Test 6: 3×3 s=2 conv 7→32 on (100, 221) input ===\n");
  std::fflush(stdout);
  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape = @[@1, @100, @221, @7];
    MPSGraphTensor* x = [g placeholderWithShape:in_shape
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    // Same weight pattern as Test 4/5
    std::vector<float> w(3 * 3 * 7 * 32);
    for (int h = 0; h < 3; ++h)
      for (int wj = 0; wj < 3; ++wj)
        for (int c = 0; c < 7; ++c)
          for (int o = 0; o < 32; ++o) {
            float v = (float)((h * 3 + wj) * 7 + c + 1) + (float)o * 0.001f;
            w[((h * 3 + wj) * 7 + c) * 32 + o] = v;
          }
    NSData* w_data = [NSData dataWithBytes:w.data()
                                     length:w.size() * sizeof(float)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@3, @3, @7, @32]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:2 strideInY:2
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    // Sparse input: only (h=1, w=1, c=0) = 1.0
    std::vector<float> in_data(100 * 221 * 7, 0.0f);
    in_data[(1 * 221 + 1) * 7 + 0] = 1.0f;
    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y, in_shape, in_data.data(), &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    // Output (1, 49, 110, 32). Position (h_out=0, w_out=0) covers input
    // window (h=0..2, w=0..2). Only (h=1, w=1, c=0) is non-zero (=1.0).
    // So out[0, 0, 0, o] = weight[h=1, w=1, c=0, o] = (1*3+1)*7+0+1 = 29 + o*0.001
    std::vector<float> expected(32);
    for (int o = 0; o < 32; ++o) {
      expected[o] = 29.0f + (float)o * 0.001f;
    }
    // We extract just out[0, 0, 0, *] = first 32 of the flat output.
    std::vector<float> out_first32(got.begin(), got.begin() + 32);
    CompareVec("conv3x3_s2_large_first_pixel", out_first32, expected, 1e-2f);

    // Also check a middle position to verify the kernel applies correctly
    // away from the corner. Set input[h=20, w=30, c=3] = 7.0, rerun.
    std::printf("\n  -- mid-position test (input[h=20, w=30, c=3] = 7.0) --\n");
    std::fflush(stdout);
    std::vector<float> in_data2(100 * 221 * 7, 0.0f);
    in_data2[(20 * 221 + 30) * 7 + 3] = 7.0f;
    std::vector<float> got2;
    NSArray<NSNumber*>* out_shape2 = nil;
    // Re-run on the same compiled exe — but we don't have a handle here.
    // Just rebuild: same graph, different input.
    if (!RunGraph(g, x, y, in_shape, in_data2.data(), &got2, &out_shape2)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    // Output position (h_out, w_out) = (10, 15) covers input window
    // (h=20..22, w=30..32). The non-zero is at (20, 30, 3) = kernel
    // position (kh=0, kw=0, c=3). So:
    //   out[0, 10, 15, o] = 7.0 * weight[0, 0, 3, o] = 7.0 * (4 + o*0.001)
    std::vector<float> expected2(32);
    for (int o = 0; o < 32; ++o) {
      expected2[o] = 7.0f * (4.0f + (float)o * 0.001f);
    }
    // Output is (1, 49, 110, 32). Position (h_out=10, w_out=15) flat =
    // (10 * 110 + 15) * 32 = 36800.
    const size_t off = (10 * 110 + 15) * 32;
    std::vector<float> out_mid32(got2.begin() + off, got2.begin() + off + 32);
    CompareVec("conv3x3_s2_large_mid_pixel", out_mid32, expected2, 1e-2f);
  }
}

}  // namespace

// ===========================================================================
// Test 7: same conv as Test 6 but with the REAL stem_s1a folded weights
// and the REAL seed-0 input.npy. Output is compared to TF reference
// stem_s1a.npy. This isolates whether the bug is in MPSGraph proper or
// in our wrapper code in metal_inference.mm.
// ===========================================================================

#include <cstdint>
#include <fstream>

namespace {

bool LoadNpyFp32_Mini(const std::string& path, std::vector<float>* out,
                     std::vector<int>* shape) {
  std::ifstream f(path, std::ios::binary);
  if (!f) { std::fprintf(stderr, "  cannot open %s\n", path.c_str()); return false; }
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
  auto lp = header.find('(', p);
  auto rp = header.find(')', lp);
  shape->clear();
  std::string ss = header.substr(lp + 1, rp - lp - 1);
  for (size_t i = 0; i < ss.size();) {
    while (i < ss.size() && (ss[i] == ' ' || ss[i] == ',')) ++i;
    if (i >= ss.size()) break;
    size_t e = i;
    while (e < ss.size() && ss[e] >= '0' && ss[e] <= '9') ++e;
    if (e == i) break;
    shape->push_back(std::stoi(ss.substr(i, e - i)));
    i = e;
  }
  size_t total = 1;
  for (int d : *shape) total *= (size_t)d;
  out->resize(total);
  f.read((char*)out->data(), total * sizeof(float));
  return (bool)f;
}

}  // namespace (unnamed continuation)

void Test7_RealStemS1a(const std::string& ref_dir) {
  std::printf("\n=== Test 7: real stem_s1a weights + real input vs TF ref ===\n");
  std::fflush(stdout);

  // Load TF reference stem_s1a output (the gold)
  std::vector<float> tf_out;
  std::vector<int> tf_shape;
  if (!LoadNpyFp32_Mini(ref_dir + "/stem_s1a.npy", &tf_out, &tf_shape)) {
    std::printf("  FAILED to load TF reference\n");
    return;
  }
  std::printf("  TF stem_s1a shape: ");
  for (int d : tf_shape) std::printf("%d ", d);
  std::printf("\n");

  // Load real input
  std::vector<float> in_data;
  std::vector<int> in_shape;
  if (!LoadNpyFp32_Mini(ref_dir + "/_input.npy", &in_data, &in_shape)) {
    std::printf("  FAILED to load input\n");
    return;
  }

  // Hand-fold layer-0 + layer-1 weights from the bundle. We can't reuse
  // the dvw_weights C++ class here without the dependency graph, so we
  // expect the user to pass `<ref_dir>/_handroll_W_hwio.bin` and
  // `_handroll_bias.bin` as raw FP32 dumps produced by Python (see
  // tools/conversion/dump_stem_s1a_weights.py).
  std::vector<float> w_hwio(3 * 3 * 7 * 32);
  std::vector<float> bias(32);
  std::ifstream wf(ref_dir + "/_handroll_W_hwio.bin", std::ios::binary);
  std::ifstream bf(ref_dir + "/_handroll_bias.bin", std::ios::binary);
  if (!wf || !bf) {
    std::printf("  FAILED to load _handroll_W_hwio.bin / _handroll_bias.bin\n");
    return;
  }
  wf.read((char*)w_hwio.data(), w_hwio.size() * sizeof(float));
  bf.read((char*)bias.data(), bias.size() * sizeof(float));

  @autoreleasepool {
    MPSGraph* g = [MPSGraph new];
    NSArray<NSNumber*>* in_shape_ns = @[@1, @100, @221, @7];
    MPSGraphTensor* x = [g placeholderWithShape:in_shape_ns
                                        dataType:MPSDataTypeFloat32
                                            name:@"x"];
    NSData* w_data = [NSData dataWithBytes:w_hwio.data()
                                     length:w_hwio.size() * sizeof(float)];
    MPSGraphTensor* W = [g constantWithData:w_data
                                       shape:@[@3, @3, @7, @32]
                                    dataType:MPSDataTypeFloat32];
    NSData* b_data = [NSData dataWithBytes:bias.data()
                                     length:bias.size() * sizeof(float)];
    MPSGraphTensor* B = [g constantWithData:b_data
                                       shape:@[@32]
                                    dataType:MPSDataTypeFloat32];
    MPSGraphConvolution2DOpDescriptor* d =
        [MPSGraphConvolution2DOpDescriptor
            descriptorWithStrideInX:2 strideInY:2
                    dilationRateInX:1 dilationRateInY:1
                             groups:1
                       paddingStyle:MPSGraphPaddingStyleTF_VALID
                         dataLayout:MPSGraphTensorNamedDataLayoutNHWC
                      weightsLayout:MPSGraphTensorNamedDataLayoutHWIO];
    MPSGraphTensor* y = [g convolution2DWithSourceTensor:x weightsTensor:W
                                                descriptor:d name:@"conv"];
    // Bias broadcast (1, 1, 1, 32) + ReLU
    MPSGraphTensor* B4 = [g reshapeTensor:B
                                 withShape:@[@1, @1, @1, @32]
                                      name:@"b4"];
    MPSGraphTensor* y_bias = [g additionWithPrimaryTensor:y
                                          secondaryTensor:B4
                                                     name:@"add_bias"];
    MPSGraphTensor* y_relu = [g reLUWithTensor:y_bias name:@"relu"];

    std::vector<float> got;
    NSArray<NSNumber*>* out_shape = nil;
    if (!RunGraph(g, x, y_relu, in_shape_ns, in_data.data(),
                  &got, &out_shape)) {
      std::printf("  RUN FAILED\n");
      return;
    }
    PrintShape(out_shape);
    std::printf("  Metal[0..8]:");
    for (int i = 0; i < 8; ++i) std::printf(" %9.3f", got[i]);
    std::printf("\n  TF   [0..8]:");
    for (int i = 0; i < 8; ++i) std::printf(" %9.3f", tf_out[i]);
    std::printf("\n");
    float max_abs = 0.0f;
    int n_close = 0;
    for (size_t i = 0; i < got.size(); ++i) {
      float d = std::fabs(got[i] - tf_out[i]);
      if (d <= 1e-3f) ++n_close;
      if (d > max_abs) max_abs = d;
    }
    std::printf("  max-abs : %.6e  close (≤1e-3): %d / %zu  verdict: %s\n",
                max_abs, n_close, got.size(),
                max_abs <= 1e-2f ? "PASS" : "FAIL");
    std::fflush(stdout);
  }
}

int main(int argc, char** argv) {
  std::printf("microtest start\n");
  std::fflush(stdout);
  Test1_Conv1x1();
  Test2_Conv3x3SingleCh();
  Test3_Conv3x3_7to1();
  Test4_Conv3x3_7to32();
  Test5_Conv3x3_S2_7to32();
  Test6_Conv3x3_S2_7to32_LargeInput();
  if (argc > 1) {
    Test7_RealStemS1a(argv[1]);
  } else {
    std::printf("\n(skipping Test 7 — pass <ref_dir> as argv[1] to enable)\n");
  }
  std::printf("\nmicrotest done\n");
  return 0;
}
