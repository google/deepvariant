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

#include "third_party/nucleus/core/status.h"

#include <memory>
#include <string>

#include "absl/strings/escaping.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/logging.h"
#include "tensorflow/tsl/platform/stacktrace.h"

namespace nucleus {

Status::~Status() {}

Status::Status(absl::StatusCode code, absl::string_view msg) {
  assert(code != absl::StatusCode::kOk);
  state_ = std::make_unique<State>();
  state_->code = code;
  state_->msg = std::string(msg);
  VLOG(5) << "Generated non-OK status: \"" << *this << "\". "
          << tsl::CurrentStackTrace();
}

void Status::Update(const Status& new_status) {
  if (ok()) {
    *this = new_status;
  }
}

void Status::SlowCopyFrom(const State* src) {
  if (src == nullptr) {
    state_ = nullptr;
  } else {
    state_ = std::make_unique<State>(*src);
  }
}

Status::State* Status::NewStateFromNonOKStatus(const Status& s) {
  return new State(*s.state_);
}

const std::string& Status::empty_string() {
  static std::string* empty = new std::string;
  return *empty;
}

std::string error_name(absl::StatusCode code) {
  switch (code) {
    case absl::StatusCode::kOk:
      return "OK";
      break;
    case absl::StatusCode::kCancelled:
      return "CANCELLED";
      break;
    case absl::StatusCode::kUnknown:
      return "UNKNOWN";
      break;
    case absl::StatusCode::kInvalidArgument:
      return "INVALID_ARGUMENT";
      break;
    case absl::StatusCode::kDeadlineExceeded:
      return "DEADLINE_EXCEEDED";
      break;
    case absl::StatusCode::kNotFound:
      return "NOT_FOUND";
      break;
    case absl::StatusCode::kAlreadyExists:
      return "ALREADY_EXISTS";
      break;
    case absl::StatusCode::kPermissionDenied:
      return "PERMISSION_DENIED";
      break;
    case absl::StatusCode::kUnauthenticated:
      return "UNAUTHENTICATED";
      break;
    case absl::StatusCode::kResourceExhausted:
      return "RESOURCE_EXHAUSTED";
      break;
    case absl::StatusCode::kFailedPrecondition:
      return "FAILED_PRECONDITION";
      break;
    case absl::StatusCode::kAborted:
      return "ABORTED";
      break;
    case absl::StatusCode::kOutOfRange:
      return "OUT_OF_RANGE";
      break;
    case absl::StatusCode::kUnimplemented:
      return "UNIMPLEMENTED";
      break;
    case absl::StatusCode::kInternal:
      return "INTERNAL";
      break;
    case absl::StatusCode::kUnavailable:
      return "UNAVAILABLE";
      break;
    case absl::StatusCode::kDataLoss:
      return "DATA_LOSS";
      break;
    default:
      char tmp[30];
      snprintf(tmp, sizeof(tmp), "UNKNOWN_CODE(%d)", static_cast<int>(code));
      return tmp;
      break;
  }
}

std::string Status::ToString() const {
  if (state_ == nullptr) {
    return "OK";
  } else {
    std::string result(error_name(state_->code));
    result += ": ";
    result += state_->msg;

    for (const std::pair<const std::string, absl::Cord>& element :
         state_->payloads) {
      absl::StrAppend(&result, " [", element.first, "='",
                      absl::CHexEscape(std::string(element.second)), "']");
    }

    return result;
  }
}

void Status::IgnoreError() const {
  // no-op
}

void Status::SetPayload(absl::string_view type_url, absl::Cord payload) {
  if (ok()) return;
  state_->payloads[std::string(type_url)] = payload;
}

absl::optional<absl::Cord> Status::GetPayload(
    absl::string_view type_url) const {
  if (ok()) return absl::nullopt;
  auto payload_iter = state_->payloads.find(std::string(type_url));
  if (payload_iter == state_->payloads.end()) return absl::nullopt;
  return payload_iter->second;
}

bool Status::ErasePayload(absl::string_view type_url) {
  if (ok()) return false;
  auto payload_iter = state_->payloads.find(std::string(type_url));
  if (payload_iter == state_->payloads.end()) return false;
  state_->payloads.erase(payload_iter);
  return true;
}

void Status::ForEachPayload(
    absl::FunctionRef<void(absl::string_view, const absl::Cord&)> visitor)
    const {
  if (ok()) return;
  for (const auto& payload : state_->payloads) {
    visitor(payload.first, payload.second);
  }
}

std::ostream& operator<<(std::ostream& os, const Status& x) {
  os << x.ToString();
  return os;
}

std::string* CheckOpHelperOutOfLine(const ::nucleus::Status& v,
                                    const char* msg) {
  std::string r("Non-OK-status: ");
  r += msg;
  r += " status: ";
  r += v.ToString();
  // Leaks string but this is only to be used in a fatal error message
  return new std::string(r);
}

}  // namespace nucleus
