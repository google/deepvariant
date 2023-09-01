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

// StatusOr<T> is the union of a Status object and a T
// object. StatusOr models the concept of an object that is either a
// usable value, or an error Status explaining why such a value is
// not present. To this end, StatusOr<T> does not allow its Status
// value to be Status::OK. Further, StatusOr<T*> does not allow the
// contained pointer to be NULL.
//
// The primary use-case for StatusOr<T> is as the return value of a
// function which may fail.
//
// Example client usage for a StatusOr<T>, where T is not a pointer:
//
//  StatusOr<float> result = DoBigCalculationThatCouldFail();
//  if (result.ok()) {
//    float answer = result.ValueOrDie();
//    printf("Big calculation yielded: %f", answer);
//  } else {
//    LOG(ERROR) << result.status();
//  }
//
// Example client usage for a StatusOr<T*>:
//
//  StatusOr<Foo*> result = FooFactory::MakeNewFoo(arg);
//  if (result.ok()) {
//    std::unique_ptr<Foo> foo(result.ValueOrDie());
//    foo->DoSomethingCool();
//  } else {
//    LOG(ERROR) << result.status();
//  }
//
// Example client usage for a StatusOr<std::unique_ptr<T>>:
//
//  StatusOr<std::unique_ptr<Foo>> result = FooFactory::MakeNewFoo(arg);
//  if (result.ok()) {
//    std::unique_ptr<Foo> foo = std::move(result.ValueOrDie());
//    foo->DoSomethingCool();
//  } else {
//    LOG(ERROR) << result.status();
//  }
//
// Example factory implementation returning StatusOr<T*>:
//
//  StatusOr<Foo*> FooFactory::MakeNewFoo(int arg) {
//    if (arg <= 0) {
//      return Status(port::error::INVALID_ARGUMENT,
//                            "Arg must be positive");
//    } else {
//      return new Foo(arg);
//    }
//  }
//
//
// This is a copy of StatusOr from
// tensorflow/compiler/xla/stream_executor/lib/statusor.h based on the original
// internal google sources with the key changes to statusor.h reapplied by hand.
#ifndef THIRD_PARTY_NUCLEUS_VENDOR_STATUSOR_H_
#define THIRD_PARTY_NUCLEUS_VENDOR_STATUSOR_H_

#include <new>
#include <type_traits>
#include <utility>

#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/core/status.h"
#include "tensorflow/core/lib/core/errors.h"
#include "tensorflow/core/lib/core/status.h"

namespace nucleus {

template <typename T>
class StatusOr {
  template <typename U>
  friend class StatusOr;

 public:
  // Construct a new StatusOr with Status::UNKNOWN status
  StatusOr() : status_(absl::StatusCode::kUnknown, "") {}

  // Construct a new StatusOr with the given non-ok status. After calling
  // this constructor, calls to ValueOrDie() is invalid.
  //
  // NOTE: Not explicit - we want to use StatusOr<T> as a return
  // value, so it is convenient and sensible to be able to do 'return
  // Status()' when the return type is StatusOr<T>.
  //
  // REQUIRES: status != Status::OK.
  // In optimized builds, passing Status::OK here will have the effect
  // of passing PosixErrorSpace::EINVAL as a fallback.
  StatusOr(const nucleus::Status& status);  // NOLINT

  // Construct a new StatusOr with the given value. If T is a plain pointer,
  // value must not be NULL. After calling this constructor, calls to
  // ValueOrDie() will succeed, and calls to status() will return OK.
  //
  // NOTE: Not explicit - we want to use StatusOr<T> as a return type
  // so it is convenient and sensible to be able to do 'return T()'
  // when the return type is StatusOr<T>.
  //
  // REQUIRES: if T is a plain pointer, value != NULL.
  // In optimized builds, passing a NULL pointer here will have
  // the effect of passing PosixErrorSpace::EINVAL as a fallback.
  StatusOr(const T& value);  // NOLINT

  // Conversion copy constructor, T must be copy constructible from U
  template <typename U>
  StatusOr(const StatusOr<U>& other)  // NOLINT
      : status_(other.status_), value_(other.value_) {}

  // Conversion assignment operator, T must be assignable from U
  template <typename U>
  StatusOr& operator=(const StatusOr<U>& other) {
    status_ = other.status_;
    value_ = other.value_;
    return *this;
  }

  // Rvalue-reference overloads of the other constructors and assignment
  // operators, to support move-only types and avoid unnecessary copying.
  StatusOr(T&& value);  // NOLINT

  // Move conversion operator to avoid unnecessary copy.
  // T must be assignable from U.
  // Not marked with explicit so the implicit conversion can happen.
  template <typename U>
  StatusOr(StatusOr<U>&& other)  // NOLINT
      : status_(std::move(other.status_)), value_(std::move(other.value_)) {}

  // Move assignment operator to avoid unnecessary copy.
  // T must be assignable from U
  template <typename U>
  StatusOr& operator=(StatusOr<U>&& other) {
    status_ = std::move(other.status_);
    value_ = std::move(other.value_);
    return *this;
  }

  // Returns a reference to our status. If this contains a T, then
  // returns Status::OK.
  const nucleus::Status& status() const { return status_; }

  // Returns this->status().ok()
  bool ok() const { return status_.ok(); }
  // Returns this->status().error_message()
  const string& error_message() const { return status_.error_message(); }
  // Returns this->status().code()
  tensorflow::error::Code code() const {
    return static_cast<tensorflow::error::Code>(status_.code());
  }

  // Returns a reference to our current value, requires that this->ok().
  // If you need to initialize a T object from the stored value,
  // ConsumeValueOrDie() may be more efficient.
  const T& ValueOrDie() const;
  T& ValueOrDie();

  // Returns our current value, requires this->ok(). Use this if
  // you would otherwise want to say std::move(s.ValueOrDie()), for example
  // if you need to initialize a T object from the stored value and you don't
  // need subsequent access to the stored value. It uses T's move constructor,
  // if it has one, so it will work with move-only types, and will often be
  // more efficient than ValueOrDie, but may leave the stored value
  // in an arbitrary valid state.
  T ConsumeValueOrDie();

 private:
  nucleus::Status status_;
  T value_;

  void CheckValueNotNull(const T& value);

  template <typename U>
  struct IsNull {
    // For non-pointer U, a reference can never be NULL.
    static inline bool IsValueNull(const U& t) { return false; }
  };

  template <typename U>
  struct IsNull<U*> {
    static inline bool IsValueNull(const U* t) { return t == NULL; }
  };
};

////////////////////////////////////////////////////////////////////////////////
// Implementation details for StatusOr<T>

template <typename T>
StatusOr<T>::StatusOr(const T& value) : status_(), value_(value) {
  CheckValueNotNull(value);
}

template <typename T>
const T& StatusOr<T>::ValueOrDie() const {
  NUCLEUS_CHECK_OK(status_);
  return value_;
}

template <typename T>
T& StatusOr<T>::ValueOrDie() {
  NUCLEUS_CHECK_OK(status_);
  return value_;
}

template <typename T>
T StatusOr<T>::ConsumeValueOrDie() {
  NUCLEUS_CHECK_OK(status_);
  return std::move(value_);
}

template <typename T>
StatusOr<T>::StatusOr(const nucleus::Status& status) : status_(status) {
  assert(!status.ok());
  if (status.ok()) {
    status_ = nucleus::Status(
        absl::StatusCode::kInternal,
        "Status::OK is not a valid constructor argument to StatusOr<T>");
  }
}

template <typename T>
StatusOr<T>::StatusOr(T&& value) : status_() {
  CheckValueNotNull(value);
  value_ = std::move(value);
}

template <typename T>
void StatusOr<T>::CheckValueNotNull(const T& value) {
  assert(!IsNull<T>::IsValueNull(value));
  if (IsNull<T>::IsValueNull(value)) {
    status_ = nucleus::Status(
        absl::StatusCode::kInternal,
        "NULL is not a valid constructor argument to StatusOr<T*>");
  }
}

#define NUCLEUS_RETURN_IF_ERROR(...)                     \
  do {                                                   \
    ::nucleus::Status _status = (__VA_ARGS__);           \
    if (TF_PREDICT_FALSE(!_status.ok())) return _status; \
  } while (0)

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_VENDOR_STATUSOR_H_
