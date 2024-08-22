/*
 * Copyright 2024 Google LLC.
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

#pragma once

#if true  // Trick to stop tooling from moving the #include around.
// MUST appear before any standard headers are included.
#include <pybind11/pybind11.h>
#endif

// IWYU pragma: always_keep // See pybind11/docs/type_caster_iwyu.rst

#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/core/status.h"
#include "third_party/nucleus/core/statusor.h"

namespace pybind11 {
namespace detail {

template <typename PayloadType>
struct type_caster<nucleus::StatusOr<PayloadType>> {
 public:
  static constexpr auto name = const_name("StatusOr");

  using PayloadCaster = make_caster<PayloadType>;
  using StatusCaster = make_caster<nucleus::Status>;

  // Convert C++ -> Python.
  // Adapted from:
  // https://github.com/pybind/pybind11_abseil/blob/baadf890e80e72969f03976489a5294bd0c6efa7/pybind11_abseil/statusor_caster.h#L84-L102
  static handle cast(const nucleus::StatusOr<PayloadType>* src,
                     return_value_policy policy, handle parent) {
    if (!src) return none().release();
    return cast_impl(*src, policy, parent);
  }

  static handle cast(const nucleus::StatusOr<PayloadType>& src,
                     return_value_policy policy, handle parent) {
    return cast_impl(src, policy, parent);
  }

  static handle cast(nucleus::StatusOr<PayloadType>&& src,
                     return_value_policy policy, handle parent) {
    return cast_impl(std::move(src), policy, parent);
  }

 private:
  template <typename CType>
  static handle cast_impl(CType&& src, return_value_policy policy, handle parent) {
    if (src.ok()) {
      return PayloadCaster::cast(std::forward<CType>(src).ConsumeValueOrDie(),
                                 policy, parent);
    }
    return StatusCaster::cast(std::forward<CType>(src).status(), return_value_policy::move, parent);
  }

  // Python to C++ conversion is not enabled:
  // bool load(handle src, bool);
};

}  // namespace detail
}  // namespace pybind11
