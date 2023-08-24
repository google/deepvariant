/*
 * Copyright 2018 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef THIRD_PARTY_NUCLEUS_VENDOR_STATUS_H_
#define THIRD_PARTY_NUCLEUS_VENDOR_STATUS_H_

#include "absl/status/status.h"
#include "absl/types/optional.h"
#include "tensorflow/tsl/platform/macros.h"

namespace nucleus {

// This class is a fork of the Tensorflow tsl::Status on March 2023.
// It's internal to nucleus. We forked it to keep a distinct type compared
// to absl::Status, so we can have a dedicated CLIF converter without ODR
// violations with other absl::Status CLIF converters.
class Status {
 public:
  /// Create a success status.
  Status() {}
  ~Status();  // Not inlined to save code space

  /// \brief Create a status with the specified error code and msg as a
  /// human-readable string containing more detailed information.
  Status(absl::StatusCode code, absl::string_view msg);

  /// Copy the specified status.
  Status(const Status& s);
  Status& operator=(const Status& s);
#ifndef SWIG
  Status(Status&& s) noexcept;
  Status& operator=(Status&& s) noexcept;
#endif  // SWIG

  /// Returns true iff the status indicates success.
  bool ok() const { return (state_ == nullptr); }

  absl::StatusCode code() const {
    return ok() ? absl::StatusCode::kOk : state_->code;
  }
  int raw_code() const { return static_cast<int>(code()); }

  const std::string& error_message() const {
    return ok() ? empty_string() : state_->msg;
  }

  bool operator==(const Status& x) const;
  bool operator!=(const Status& x) const;

  /// \brief If `ok()`, stores `new_status` into `*this`.  If `!ok()`,
  /// preserves the current status, but may augment with additional
  /// information about `new_status`.
  ///
  /// Convenient way of keeping track of the first error encountered.
  /// Instead of:
  ///   `if (overall_status.ok()) overall_status = new_status`
  /// Use:
  ///   `overall_status.Update(new_status);`
  void Update(const Status& new_status);

  /// \brief Return a string representation of this status suitable for
  /// printing. Returns the string `"OK"` for success.
  ///
  /// By default, it returns combination of the error code name, the message and
  /// any associated payload messages. This string is designed simply to be
  /// human readable and its exact format should not be load bearing. Do not
  /// depend on the exact format of the result of `ToString()` which is subject
  /// to change.
  std::string ToString() const;

  // Ignores any errors. This method does nothing except potentially suppress
  // complaints from any tools that are checking that errors are not dropped on
  // the floor.
  void IgnoreError() const;

  //----------------------------------------------------------------------------
  // Payload Management APIs (Cloned from absl::Status)
  //----------------------------------------------------------------------------
  // A payload may be attached to a status to provide additional context to an
  // error that may not be satisfied by an existing `tsl::error::Code`.
  // Typically, this payload serves one of several purposes:
  //
  //   * It may provide more fine-grained semantic information about the error
  //     to facilitate actionable remedies.
  //   * It may provide human-readable contexual information that is more
  //     appropriate to display to an end user.
  //
  // A payload consists of a [key,value] pair, where the key is a string
  // referring to a unique "type URL" and the value is an object of type
  // `absl::Cord` to hold the contextual data.
  //
  // The "type URL" should be unique and follow the format of a URL
  // (https://en.wikipedia.org/wiki/URL) and, ideally, provide some
  // documentation or schema on how to interpret its associated data. For
  // example, the default type URL for a protobuf message type is
  // "type.googleapis.com/packagename.messagename". Other custom wire formats
  // should define the format of type URL in a similar practice so as to
  // minimize the chance of conflict between type URLs.
  // Users should ensure that the type URL can be mapped to a concrete
  // C++ type if they want to deserialize the payload and read it effectively.
  //
  // To attach a payload to a status object, call `Status::SetPayload()`,
  // passing it the type URL and an `absl::Cord` of associated data. Similarly,
  // to extract the payload from a status, call `Status::GetPayload()`. You
  // may attach multiple payloads (with differing type URLs) to any given
  // status object, provided that the status is currently exhibiting an error
  // code (i.e. is not OK).
  // TODO: Use absl::Cord for payload value type.

  // The Payload-related APIs are cloned from absl::Status.
  //
  // Returns the payload of a status given its unique `type_url` key, if
  // present.
  absl::optional<absl::Cord> GetPayload(absl::string_view type_url) const;

  // Sets the payload for a non-ok status using a `type_url` key, overwriting
  // any existing payload for that `type_url`.
  //
  // This function does nothing if the Status is ok.
  void SetPayload(absl::string_view type_url, absl::Cord payload);

  // Erases the payload corresponding to the `type_url` key.  Returns `true` if
  // the payload was present.
  bool ErasePayload(absl::string_view type_url);

  // Iterates over the stored payloads and calls the
  // `visitor(type_key, payload)` callable for each one.
  //
  // The order of calls to `visitor()` is not specified and may change at
  // any time and any mutation on the same Status object during visitation is
  // forbidden and could result in undefined behavior.
  void ForEachPayload(
      absl::FunctionRef<void(absl::string_view, const absl::Cord&)> visitor)
      const;

 private:
  friend Status FromAbslStatus(const absl::Status& s);

  static const std::string& empty_string();
  struct State {
    State() TF_ATTRIBUTE_NOINLINE = default;
    ~State() TF_ATTRIBUTE_NOINLINE = default;
    State(const State&) TF_ATTRIBUTE_NOINLINE = default;
    State& operator=(const State&) TF_ATTRIBUTE_NOINLINE = default;

    absl::StatusCode code;
    std::string msg;
    std::unordered_map<std::string, absl::Cord> payloads;
  };

  // OK status has a `NULL` state_.  Otherwise, `state_` points to
  // a `State` structure containing the error code and message(s)
  std::unique_ptr<State> state_;

  void SlowCopyFrom(const State* src);
  State* NewStateFromNonOKStatus(const Status& s);
};

inline Status::Status(const Status& s)
    : state_((s.state_ == nullptr) ? nullptr : NewStateFromNonOKStatus(s)) {}

inline Status& Status::operator=(const Status& s) {
  // The following condition catches both aliasing (when this == &s),
  // and the common case where both s and *this are ok.
  if (state_ != s.state_) {
    SlowCopyFrom(s.state_.get());
  }
  return *this;
}

#ifndef SWIG
inline Status::Status(Status&& s) noexcept : state_(std::move(s.state_)) {}

inline Status& Status::operator=(Status&& s) noexcept {
  if (state_ != s.state_) {
    state_ = std::move(s.state_);
  }
  return *this;
}
#endif  // SWIG

inline bool Status::operator==(const Status& x) const {
  return (this->state_ == x.state_) || (ToString() == x.ToString());
}

inline bool Status::operator!=(const Status& x) const { return !(*this == x); }

/// @ingroup core
std::ostream& operator<<(std::ostream& os, const Status& x);

std::string* CheckOpHelperOutOfLine(const ::nucleus::Status& v,
                                    const char* msg);

inline std::string* TfCheckOpHelper(::nucleus::Status v, const char* msg) {
  if (v.ok()) return nullptr;
  return CheckOpHelperOutOfLine(v, msg);
}

#define NUCLEUS_DO_CHECK_OK(val, level)                         \
  while (auto* _result = ::nucleus::TfCheckOpHelper(val, #val)) \
  LOG(level) << *(_result)

#define NUCLEUS_CHECK_OK(val) NUCLEUS_DO_CHECK_OK(val, FATAL)
#define NUCLEUS_QCHECK_OK(val) NUCLEUS_DO_CHECK_OK(val, QFATAL)

inline Status FailedPrecondition(absl::string_view message) {
  return Status(absl::StatusCode::kFailedPrecondition, message);
}

inline Status NotFound(absl::string_view message) {
  return Status(absl::StatusCode::kNotFound, message);
}

inline Status Unknown(absl::string_view message) {
  return Status(absl::StatusCode::kUnknown, message);
}

inline Status Internal(absl::string_view message) {
  return Status(absl::StatusCode::kInternal, message);
}

inline Status OutOfRange(absl::string_view message) {
  return Status(absl::StatusCode::kOutOfRange, message);
}

inline Status InvalidArgument(absl::string_view message) {
  return Status(absl::StatusCode::kInvalidArgument, message);
}

inline Status DataLoss(absl::string_view message) {
  return Status(absl::StatusCode::kDataLoss, message);
}

inline Status Unimplemented(absl::string_view message) {
  return Status(absl::StatusCode::kUnimplemented, message);
}

inline bool IsOutOfRange(const Status& status) {
  return status.code() == absl::StatusCode::kOutOfRange;
}
}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_VENDOR_STATUS_H_
