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

#include "third_party/nucleus/io/gff_reader.h"
#include "third_party/pybind11/include/pybind11/chrono.h"
#include "third_party/pybind11/include/pybind11/complex.h"
#include "third_party/pybind11/include/pybind11/functional.h"
#include "third_party/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(gff_reader, m) {
  using namespace ::nucleus;
  py::class_<GffReader>(m, "GffReader")
      .def_static("from_file", &GffReader::FromFile, py::arg("gffPath"),
                  py::arg("options"))
      .def("iterate", &GffReader::Iterate)
      .def("__enter__", &GffReader::PythonEnter)
      .def("__exit__", &GffReader::Close)
      .def_property_readonly("header", &GffReader::Header);
  // Do I need this?
  // GffReader_class.def(
  //     "as_nucleus_Reader",
  //     [](::nucleus::GffReader* self) {
  //         return pybind11::capsule(static_cast<void *>(self));
  //     }
  // );
  py::class_<GffIterable>(m, "GffIterable")
      .def("PythonNext", &GffIterable::PythonNext, py::arg("gff"))
      .def("Release", &GffIterable::Release)
      .def("__enter__", &GffIterable::PythonEnter)
      .def("__exit__", &GffIterable::PythonExit);
  // Do I need this?
  // GffIterable_class.def(
  //     "as_nucleus_IterableBase",
  //     [](::nucleus::Iterable< ::nucleus::genomics::v1::GffRecord>* self) {
  //         return pybind11::capsule(static_cast<void *>(self));
  //     }
  // );
}