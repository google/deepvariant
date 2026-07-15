// tensorflow::errors::* → absl::*Error factory functions.
#pragma once
#include "tensorflow/core/lib/core/status.h"
#include "absl/strings/string_view.h"

namespace tensorflow {
namespace errors {

inline Status InvalidArgument(absl::string_view msg) {
  return absl::InvalidArgumentError(msg);
}
inline Status NotFound(absl::string_view msg) {
  return absl::NotFoundError(msg);
}
inline Status AlreadyExists(absl::string_view msg) {
  return absl::AlreadyExistsError(msg);
}
inline Status Internal(absl::string_view msg) {
  return absl::InternalError(msg);
}
inline Status Unimplemented(absl::string_view msg) {
  return absl::UnimplementedError(msg);
}
inline Status FailedPrecondition(absl::string_view msg) {
  return absl::FailedPreconditionError(msg);
}
inline Status OutOfRange(absl::string_view msg) {
  return absl::OutOfRangeError(msg);
}
inline Status DataLoss(absl::string_view msg) {
  return absl::DataLossError(msg);
}
inline Status Aborted(absl::string_view msg) {
  return absl::AbortedError(msg);
}

inline bool IsNotFound(const Status& s) { return absl::IsNotFound(s); }
inline bool IsInvalidArgument(const Status& s) { return absl::IsInvalidArgument(s); }
inline bool IsInternal(const Status& s) { return absl::IsInternal(s); }

}  // namespace errors
}  // namespace tensorflow
