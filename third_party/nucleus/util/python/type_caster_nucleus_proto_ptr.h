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

#if true  // Trick to stop tooling from moving the #include around.
// MUST appear before any standard headers are included.
#include <pybind11/pybind11.h>
#endif

// IWYU pragma: always_keep // See pybind11/docs/type_caster_iwyu.rst

#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "google/protobuf/message.h"
#include "third_party/py/google/protobuf/proto_api.h"
#include "third_party/pybind11/include/pybind11/gil_safe_call_once.h"

namespace nucleus {
namespace internal {

inline const google::protobuf::python::PyProto_API* GetPyProtoApi() {
  PYBIND11_CONSTINIT static pybind11::gil_safe_call_once_and_store<
      const google::protobuf::python::PyProto_API*>
      storage;
  return storage
      .call_once_and_store_result([]() {
        auto* ptr = static_cast<const google::protobuf::python::PyProto_API*>(
            PyCapsule_Import(google::protobuf::python::PyProtoAPICapsuleName(), 0));
        if (ptr == nullptr) {
          pybind11::raise_from(PyExc_RuntimeError,
                               "nucleus::internal::GetPyProtoApi() failed.");
          throw pybind11::error_already_set();
        }
        return ptr;
      })
      .get_stored();
}

inline bool handle_get_message_pointer_errors(const google::protobuf::Message* cpb) {
  if (cpb == nullptr) {
    ABSL_CHECK(PyErr_Occurred());
    if (PyErr_ExceptionMatches(PyExc_TypeError)) {
      PyErr_Clear();
      return false;
    }
    throw pybind11::error_already_set();
  }
  if (PyErr_Occurred()) {
    ABSL_LOG(ERROR) << "UNEXPECTED Python error!";
    throw pybind11::error_already_set();
  }
  return true;
}

template <typename ProtoPtrT>
struct type_caster_proto_ptr_base {
  template <typename T_>
  using cast_op_type = pybind11::detail::cast_op_type<T_>;

  explicit operator ProtoPtrT*() { return &proto_ptr; }
  explicit operator ProtoPtrT&() { return proto_ptr; }

  ProtoPtrT proto_ptr;
};

template <typename T>
class type_caster_empty_proto_ptr
    : public type_caster_proto_ptr_base<nucleus::EmptyProtoPtr<T>> {
 public:
  static constexpr auto name = pybind11::detail::const_name("EmptyProtoPtr");

  // Python to C++ conversion.
  bool load(pybind11::handle src, bool /*convert*/) {
    auto* py_proto_api = GetPyProtoApi();
    google::protobuf::Message* cpb = py_proto_api->GetMutableMessagePointer(src.ptr());
    if (!handle_get_message_pointer_errors(cpb)) {
      return false;
    }
    this->proto_ptr.p_ = dynamic_cast<T*>(cpb);
    return (this->proto_ptr.p_ != nullptr);
  }
};

template <typename T>
class type_caster_const_proto_ptr
    : public type_caster_proto_ptr_base<nucleus::ConstProtoPtr<T>> {
 public:
  static constexpr auto name = pybind11::detail::const_name("ConstProtoPtr");

  // Python to C++ conversion.
  bool load(pybind11::handle src, bool /*convert*/) {
    auto* py_proto_api = GetPyProtoApi();
    const google::protobuf::Message* cpb = py_proto_api->GetMessagePointer(src.ptr());
    if (!handle_get_message_pointer_errors(cpb)) {
      return false;
    }
    this->proto_ptr.p_ = dynamic_cast<const T*>(cpb);
    return (this->proto_ptr.p_ != nullptr);
  }
};

}  // namespace internal
}  // namespace nucleus

namespace pybind11 {
namespace detail {

template <typename T>
struct type_caster<nucleus::EmptyProtoPtr<T>>
    : nucleus::internal::type_caster_empty_proto_ptr<T> {};

template <typename T>
struct type_caster<nucleus::ConstProtoPtr<T>>
    : nucleus::internal::type_caster_const_proto_ptr<T> {};

}  // namespace detail
}  // namespace pybind11
