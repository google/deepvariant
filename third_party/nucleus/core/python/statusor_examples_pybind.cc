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

#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/core/statusor_examples.h"

PYBIND11_MODULE(statusor_examples, m) {
  using namespace ::nucleus;  // NOLINT

  namespace py = pybind11;

  // Intended usage: `hasattr(statusor_examples, 'USING_PYBIND')`
  m.attr("USING_PYBIND") = py::none();

  m.def("MakeIntOK", MakeIntOK);
  m.def("MakeIntFail", MakeIntFail);
  m.def("MakeStrOK", MakeStrOK);
  m.def("MakeStrFail", MakeStrFail);
  // These may not be relevant to production code situations:
  // m.def("MakeIntUniquePtrOK", MakeIntUniquePtrOK);
  // m.def("MakeIntUniquePtrFail", MakeIntUniquePtrFail);
  // m.def("MakeIntVectorOK", MakeIntVectorOK);
  // m.def("MakeIntVectorFail", MakeIntVectorFail);
  m.def("FuncReturningStatusOK", FuncReturningStatusOK);
  m.def("FuncReturningStatusFail", FuncReturningStatusFail);

  py::classh<StringOwner>(m, "StringOwner")
      .def_static("Factory", &StringOwner::Factory)
      .def("GetText", &StringOwner::GetText);
}