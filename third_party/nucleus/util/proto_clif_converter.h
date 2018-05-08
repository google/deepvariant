/*
 * Copyright 2018 Google Inc.
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

#include "clif/python/types.h"
#include "google/protobuf/message.h"
#include "net/proto2/python/public/python_protobuf.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace clif {

// Note: comment below is an instruction to CLIF.
// CLIF use `nucleus::ProtoPtr` as ProtoPtr

// Convert from Python protocol buffer object py to a C++ pointer.
// Unlike the conversions that CLIF automatically generates for protocol
// buffers, this one does no copying if the Python protocol buffer uses
// the C++ memory layout.
template <typename T>
bool Clif_PyObjAs(PyObject* py, nucleus::ProtoPtr<T>* c) {
  CHECK(c != nullptr);

  ::google::protobuf::Message* cpb = google::protobuf::python::MutableCProtoInsidePyProto(py);
  if (cpb == nullptr) {
    // Clients might depend on our non-copying semantics, so we can't fall
    // back on CLIF here but instead must fail loudly.
    return false;
  }

  c->p_ = dynamic_cast<T*>(cpb);

  return true;
}

}  // namespace clif
#endif  // THIRD_PARTY_NUCLEUS_UTIL_PROTO_CLIF_CONVERTER_H_
