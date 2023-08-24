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
 *
 */

#include "third_party/nucleus/io/reader_base.h"

#include "third_party/nucleus/core/status.h"

namespace nucleus {

// Reader class methods

Reader::~Reader() {
  // If there is an outstanding iterable, we need to tell it that
  // the reader is dead so it doesn't still try to use it.
  absl::MutexLock lock(&mutex_);
  if (live_iterable_ != nullptr) {
    live_iterable_->reader_ = nullptr;
  }
}

// IterableBase class methods

IterableBase::IterableBase(const Reader* reader) : reader_(reader) {}

IterableBase::~IterableBase() {
  // We cannot return a Status from our destructor, so the best we can do
  // if we need to release resources and cannot is CHECK-fail.
  NUCLEUS_CHECK_OK(Release());
}

nucleus::Status IterableBase::Release() {
  if (IsAlive()) {
    absl::MutexLock lock(&reader_->mutex_);
    if (reader_->live_iterable_ == nullptr) {
      return ::nucleus::FailedPrecondition("reader_->live_iterable_ is null");
    }
    reader_->live_iterable_ = nullptr;
    reader_ = nullptr;
  }
  return ::nucleus::Status();
}

bool IterableBase::IsAlive() const { return reader_ != nullptr; }

::nucleus::Status IterableBase::CheckIsAlive() const {
  if (!IsAlive()) return ::nucleus::FailedPrecondition("Reader is not alive");
  return nucleus::Status();
}

::nucleus::Status IterableBase::PythonEnter() { return nucleus::Status(); }

::nucleus::Status IterableBase::PythonExit() { return Release(); }

}  // namespace nucleus
