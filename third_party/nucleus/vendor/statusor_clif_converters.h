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

#ifndef THIRD_PARTY_NUCLEUS_VENDOR_STATUSOR_CLIF_CONVERTERS_H_
#define THIRD_PARTY_NUCLEUS_VENDOR_STATUSOR_CLIF_CONVERTERS_H_

#include "clif/python/postconv.h"
#include "clif/python/types.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/vendor/statusor.h"

// Note: comment below is an instruction to CLIF.
// NOLINTNEXTLINE
// CLIF use `::nucleus::StatusOr` as StatusOr, NumTemplateParameter:1, HasPyObjFromOnly
// CLIF use `::nucleus::Status` as Status, HasPyObjFromOnly

namespace nucleus {

PyObject* Clif_PyObjFrom(const nucleus::Status& c, const ::clif::py::PostConv&);

}  // namespace nucleus

namespace nucleus {
namespace internal {

void ErrorFromStatus(const ::nucleus::Status& status);

}  // namespace internal

template <typename T>
PyObject* Clif_PyObjFrom(const StatusOr<T>& c, const ::clif::py::PostConv& pc) {
  if (!c.ok()) {
    internal::ErrorFromStatus(c.status());
    return nullptr;
  } else {
    using ::clif::Clif_PyObjFrom;
    return Clif_PyObjFrom(c.ValueOrDie(), pc.Get(0));
  }
}

template <typename T>
PyObject* Clif_PyObjFrom(StatusOr<T>&& c, const ::clif::py::PostConv& pc) {
  if (!c.ok()) {
    internal::ErrorFromStatus(c.status());
    return nullptr;
  } else {
    using ::clif::Clif_PyObjFrom;
    return Clif_PyObjFrom(c.ConsumeValueOrDie(), pc.Get(0));
  }
}

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_VENDOR_STATUSOR_CLIF_CONVERTERS_H_
