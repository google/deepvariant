// Maps TF logging macros to abseil equivalents.
// absl/log/log.h already defines LOG(severity) with INFO/WARNING/ERROR/FATAL.
// absl/log/check.h already defines CHECK, DCHECK, CHECK_EQ, etc.
// We just expose these without redefining any token names.
#pragma once
#include "absl/log/check.h"
#include "absl/log/log.h"
