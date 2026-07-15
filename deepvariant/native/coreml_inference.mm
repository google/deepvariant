// Core ML inference implementation (Obj-C++).
// Loads a .mlpackage, compiles on first run (Core ML caches the
// .mlmodelc in ~/Library/Caches/com.apple.CoreML/), and runs
// batched prediction via MLModel.predictionsFromBatch:error:.

#include "deepvariant/native/coreml_inference.h"

#import <CoreML/CoreML.h>
#import <Foundation/Foundation.h>

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace deepvariant {

namespace {

// True when `arr`'s memory is densely packed C-contiguous (row-major) for its
// shape, i.e. the stride of each dimension equals the product of the extents of
// all dimensions inside it. Only then is a flat memcpy against a C-contiguous
// host buffer valid. MLMultiArray makes no layout guarantee — `initWithShape:`
// allocations and especially prediction-result arrays can carry padded strides.
bool IsPackedRowMajor(MLMultiArray* arr) {
  NSArray<NSNumber*>* shape = arr.shape;
  NSArray<NSNumber*>* strides = arr.strides;
  NSInteger expected = 1;
  for (NSInteger d = (NSInteger)shape.count - 1; d >= 0; --d) {
    if (strides[d].integerValue != expected) return false;
    expected *= shape[d].integerValue;
  }
  return true;
}

// Copy `count` floats between a C-contiguous row-major host buffer and `arr`,
// honoring `arr.strides` (element units). When `to_array` is true the host
// buffer is the source (scatter into `arr`); otherwise `arr` is the source
// (gather into the host buffer). `count` must equal the product of the shape.
void StridedCopyFloat(MLMultiArray* arr, float* host, size_t count,
                      bool to_array) {
  float* data = (float*)arr.dataPointer;
  const NSUInteger nd = arr.shape.count;
  std::vector<NSInteger> dims(nd), strides(nd), idx(nd, 0);
  for (NSUInteger d = 0; d < nd; ++d) {
    dims[d] = arr.shape[d].integerValue;
    strides[d] = arr.strides[d].integerValue;
  }
  for (size_t flat = 0; flat < count; ++flat) {
    NSInteger off = 0;
    for (NSUInteger d = 0; d < nd; ++d) off += idx[d] * strides[d];
    if (to_array) {
      data[off] = host[flat];
    } else {
      host[flat] = data[off];
    }
    // Increment the multi-index, last dimension varying fastest (row-major).
    for (NSInteger d = (NSInteger)nd - 1; d >= 0; --d) {
      if (++idx[d] < dims[d]) break;
      idx[d] = 0;
    }
  }
}

}  // namespace

struct CoreMLModel::Impl {
  MLModel* model = nil;
  NSString* input_name  = @"x";
  NSString* output_name = @"classification";
};

CoreMLModel::CoreMLModel() : impl_(std::make_unique<Impl>()) {}
CoreMLModel::~CoreMLModel() = default;

// static
std::unique_ptr<CoreMLModel> CoreMLModel::Load(
    const std::string& path, ComputeUnits compute_units) {
  @autoreleasepool {
    NSError* error = nil;
    NSURL* url = [NSURL fileURLWithPath:
        [NSString stringWithUTF8String:path.c_str()]];

    // Compile the .mlpackage to .mlmodelc (cached by Core ML).
    NSURL* compiled = [MLModel compileModelAtURL:url error:&error];
    if (!compiled) {
      NSLog(@"CoreML compile failed: %@", error.localizedDescription);
      return nullptr;
    }

    MLModelConfiguration* cfg = [[MLModelConfiguration alloc] init];
    switch (compute_units) {
      case ComputeUnits::kAll:
        cfg.computeUnits = MLComputeUnitsAll;
        break;
      case ComputeUnits::kCpuAndGpu:
        cfg.computeUnits = MLComputeUnitsCPUAndGPU;
        break;
      case ComputeUnits::kCpuOnly:
        cfg.computeUnits = MLComputeUnitsCPUOnly;
        break;
    }

    MLModel* model = [MLModel modelWithContentsOfURL:compiled
                                        configuration:cfg
                                                error:&error];
    if (!model) {
      NSLog(@"CoreML load failed: %@", error.localizedDescription);
      return nullptr;
    }

    // Inspect input/output names + shapes from the model description.
    auto out = std::unique_ptr<CoreMLModel>(new CoreMLModel());
    out->impl_->model = model;

    MLModelDescription* desc = model.modelDescription;
    if (desc.inputDescriptionsByName.count > 0) {
      NSString* name = desc.inputDescriptionsByName.allKeys.firstObject;
      out->impl_->input_name = name;
      out->input_name_ = name.UTF8String;
      MLFeatureDescription* fd = desc.inputDescriptionsByName[name];
      if (fd.type == MLFeatureTypeMultiArray) {
        NSArray<NSNumber*>* shape = fd.multiArrayConstraint.shape;
        if (shape.count >= 4) {
          // shape = (N, H, W, C) or (N, C, H, W); our model uses NHWC.
          out->input_height_   = shape[1].intValue;
          out->input_width_    = shape[2].intValue;
          out->input_channels_ = shape[3].intValue;
        }
      }
    }
    if (desc.outputDescriptionsByName.count > 0) {
      NSString* name = desc.outputDescriptionsByName.allKeys.firstObject;
      out->impl_->output_name = name;
      out->output_name_ = name.UTF8String;
      MLFeatureDescription* fd = desc.outputDescriptionsByName[name];
      if (fd.type == MLFeatureTypeMultiArray) {
        NSArray<NSNumber*>* shape = fd.multiArrayConstraint.shape;
        if (shape.count >= 2) {
          out->num_classes_ = shape[1].intValue;
        }
      }
    }

    return out;
  }
}

bool CoreMLModel::Predict(const float* images, int N, int H, int W, int C,
                          float* probs, int num_classes) {
  @autoreleasepool {
    NSError* error = nil;
    MLModel* model = impl_->model;
    NSString* in_name  = impl_->input_name;
    NSString* out_name = impl_->output_name;

    const NSInteger elemPerImage = H * W * C;

    // Single (N,H,W,C) MLMultiArray covering the whole batch — lets Core ML
    // route the whole batch through GPU/ANE in one shot instead of N
    // separate predictionFromFeatures: calls (which dominate runtime when
    // GPU dispatch overhead > inference time).
    NSArray<NSNumber*>* shape = @[@(N), @(H), @(W), @(C)];
    MLMultiArray* arr = [[MLMultiArray alloc]
        initWithShape:shape
            dataType:MLMultiArrayDataTypeFloat32
               error:&error];
    if (!arr) {
      NSLog(@"MLMultiArray alloc failed: %@", error.localizedDescription);
      return false;
    }
    const size_t in_count = (size_t)N * (size_t)elemPerImage;
    if (IsPackedRowMajor(arr)) {
      std::memcpy(arr.dataPointer, images, in_count * sizeof(float));
    } else {
      // Padded strides: scatter element-by-element so we write the right cells.
      StridedCopyFloat(arr, const_cast<float*>(images), in_count,
                       /*to_array=*/true);
    }

    MLDictionaryFeatureProvider* fp =
        [[MLDictionaryFeatureProvider alloc]
            initWithDictionary:@{in_name: arr}
                         error:&error];
    if (!fp) {
      NSLog(@"Feature provider failed: %@", error.localizedDescription);
      return false;
    }

    id<MLFeatureProvider> result =
        [model predictionFromFeatures:fp error:&error];
    if (!result) {
      NSLog(@"Batch prediction failed: %@", error.localizedDescription);
      return false;
    }

    MLMultiArray* out_arr =
        [result featureValueForName:out_name].multiArrayValue;
    if (!out_arr) {
      NSLog(@"Output '%@' missing in batch result", out_name);
      return false;
    }
    // Output is FP32 (we requested it at conversion time) and shape (N, K).
    const size_t out_count = (size_t)N * (size_t)num_classes;
    // Require an exact logical-element match: a larger array means a different
    // shape than (N, num_classes), and copying/gathering out_count elements
    // against its real dims would silently scramble the per-genotype
    // probabilities (the strided gather walks out_arr's own dims).
    if ((size_t)out_arr.count != out_count) {
      NSLog(@"Output array has %ld elements, expected %zu",
            (long)out_arr.count, out_count);
      return false;
    }
    if (IsPackedRowMajor(out_arr)) {
      std::memcpy(probs, out_arr.dataPointer, out_count * sizeof(float));
    } else {
      // Core ML may hand back a strided/padded array; gather honoring strides.
      StridedCopyFloat(out_arr, probs, out_count, /*to_array=*/false);
    }
    return true;
  }
}

}  // namespace deepvariant
