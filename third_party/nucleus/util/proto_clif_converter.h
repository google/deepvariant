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
#ifndef THIRD_PARTY_NUCLEUS_UTIL_PROTO_CLIF_CONVERTER_H_
#define THIRD_PARTY_NUCLEUS_UTIL_PROTO_CLIF_CONVERTER_H_

#include "google/protobuf/message.h"
#include "python/google/protobuf/proto_api.h"
#include "clif/python/types.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

// Note: the comments below are instructions to CLIF.
// CLIF use `::nucleus::EmptyProtoPtr` as EmptyProtoPtr
// CLIF use `::nucleus::ConstProtoPtr` as ConstProtoPtr

const ::google::protobuf::python::PyProto_API* GetPyProtoApi(PyObject* py);

// Convert from Python protocol buffer object py to a C++ pointer.
// Unlike the conversions that CLIF automatically generates for protocol
// buffers, this one does no copying if the Python protocol buffer uses
// the C++ memory layout.
template <typename T>
bool Clif_PyObjAs(PyObject* py, EmptyProtoPtr<T>* c) {
  CHECK(c != nullptr);

  auto* py_proto_api = GetPyProtoApi(py);
  if (py_proto_api == nullptr) {
    PyErr_SetString(PyExc_RuntimeError, "Could not load PyProto API");
    return false;
  }

  ::google::protobuf::Message* cpb = py_proto_api->GetMutableMessagePointer(py);
  if (cpb == nullptr) {
    // Clients might depend on our non-copying semantics, so we can't fall
    // back on CLIF here but instead must fail loudly.
    PyErr_SetString(PyExc_RuntimeError,
                    "Python protobuf did not contain a mutable C++ protobuf");
    return false;
  } else {
    c->p_ = dynamic_cast<T*>(cpb);
    if (c->p_ == nullptr) {
      // DO NOT DELETE THIS WARNING!  Without it, the above dynamic_cast
      // will fail when running from a Python 3 pip package.
      LOG(WARNING) << "Failed to cast type " << typeid(*cpb).name();
      PyErr_SetString(PyExc_RuntimeError, "Dynamic cast failed");
      return false;
    }
    return true;
  }
}

// Convert from Python protocol buffer object py to a C++ pointer.
// Unlike the conversions that CLIF automatically generates for protocol
// buffers, this one does no copying if the Python protocol buffer uses
// the C++ memory layout.
template <typename T>
bool Clif_PyObjAs(PyObject* py, ConstProtoPtr<T>* c) {
  CHECK(c != nullptr);

  auto* py_proto_api = GetPyProtoApi(py);
  if (py_proto_api == nullptr) {
    PyErr_SetString(PyExc_RuntimeError, "Could not load PyProto API");
    return false;
  }

  const ::google::protobuf::Message* cpb = py_proto_api->GetMessagePointer(py);
  if (cpb == nullptr) {
    // Clients might depend on our non-copying semantics, so we can't fall
    // back on CLIF here but instead must fail loudly.
    PyErr_SetString(PyExc_RuntimeError,
                    "Python protobuf did not contain a C++ protobuf");
    return false;
  } else {
    c->p_ = dynamic_cast<const T*>(cpb);
    if (c->p_ == nullptr) {
      // DO NOT DELETE THIS WARNING!  Without it, the above dynamic_cast
      // will fail when running from a Python 3 pip package.
      LOG(WARNING) << "Failed to cast type " << typeid(*cpb).name();
      PyErr_SetString(PyExc_RuntimeError, "Dynamic cast failed");
      return false;
    }
    return true;
  }
}

}  // namespace nucleus

#endif  // THIRD_PARTY_NUCLEUS_UTIL_PROTO_CLIF_CONVERTER_H_
