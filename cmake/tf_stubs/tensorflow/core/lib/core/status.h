// tensorflow::Status → absl::Status (same gRPC code semantics).
#pragma once
#include <string>
#include "absl/status/status.h"
#include "absl/status/statusor.h"

namespace tensorflow {

using Status = absl::Status;

namespace error {
using Code = absl::StatusCode;
}  // namespace error

inline Status OkStatus() { return absl::OkStatus(); }
inline bool IsOk(const Status& s) { return s.ok(); }

}  // namespace tensorflow
