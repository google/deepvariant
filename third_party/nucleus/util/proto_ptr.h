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
#ifndef THIRD_PARTY_NUCLEUS_UTIL_PROTO_PTR_H_
#define THIRD_PARTY_NUCLEUS_UTIL_PROTO_PTR_H_

// CLIF normally will serialize/deserialize protocol
// buffers when passing them from C++ to/from Python.
// These wrappers disable that default handling.
namespace nucleus {

// Use this wrapper when the C++ code fills in an EMPTY
// protocol buffer.  DO NOT use this to pass a non-empty
// protocol buffer from Python to C++; it will fail at
// runtime.
template <class T>
class EmptyProtoPtr {
 public:
  EmptyProtoPtr(T* p) : p_(p) {}
  EmptyProtoPtr() : p_(nullptr) {}

  T* p_;
};

// Use this wrapper when the C++ code reads, but does
// not modify, the Python protocol buffer.
template <class T>
class ConstProtoPtr {
 public:
  ConstProtoPtr(T* p) : p_(p) {}
  ConstProtoPtr() : p_(nullptr) {}

  const T* p_;
};

}  // namespace nucleus
#endif  // THIRD_PARTY_NUCLEUS_UTIL_PROTO_PTR_H_
