// TF platform macros → abseil / compiler builtins.
#pragma once
#include "absl/base/optimization.h"

#ifndef TF_PREDICT_FALSE
#  define TF_PREDICT_FALSE(x) ABSL_PREDICT_FALSE(x)
#  define TF_PREDICT_TRUE(x)  ABSL_PREDICT_TRUE(x)
#endif

#ifndef TF_MUST_USE_RESULT
#  define TF_MUST_USE_RESULT [[nodiscard]]
#endif

#ifndef TF_DISALLOW_COPY_AND_ASSIGN
#  define TF_DISALLOW_COPY_AND_ASSIGN(T) \
    T(const T&) = delete;                \
    void operator=(const T&) = delete
#endif

#ifndef TF_ATTRIBUTE_NOINLINE
#  define TF_ATTRIBUTE_NOINLINE __attribute__((noinline))
#endif
