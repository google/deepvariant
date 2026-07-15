// tf_compat.h — umbrella header pulled into every nucleus/deepvariant
// compilation unit via -include (CMakeLists.txt target_compile_options).
// Maps TF platform macros to abseil equivalents; also provides commonly
// used abseil includes that were transitively pulled in by TF in Bazel.
#pragma once
#include "tensorflow/core/platform/types.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/core/platform/macros.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/core/errors.h"
// Common abseil headers that TF code always provided transitively.
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/strings/str_format.h"
#include "absl/memory/memory.h"
#include "absl/types/optional.h"
