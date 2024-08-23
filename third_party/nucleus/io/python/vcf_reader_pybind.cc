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

#include "third_party/pybind11/include/pybind11/cast.h"
#if true  // Trick to stop tooling from moving the #include around.
// MUST appear before any standard headers are included.
#include <pybind11/pybind11.h>
#endif

#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

namespace py = pybind11;

PYBIND11_MODULE(vcf_reader, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::nucleus;  // NOLINT
  py::classh<VcfReader>(m, "VcfReader")
      .def_static("from_file", &VcfReader::FromFile, py::arg("vcfPath"),
                  py::arg("options"))
      .def_static("from_file_with_header", &VcfReader::FromFileWithHeader,
                  py::arg("variantsPath"), py::arg("options"),
                  py::arg("header"))
      .def("iterate",
           [](VcfReader& self) {
             auto cpp_result = self.Iterate();
             auto ret0 = py::cast(std::move(cpp_result));
             auto postproc = py::module_::import(
                 "third_party.nucleus.io.clif_postproc");
             return postproc.attr("WrappedVariantIterable")(ret0);
           })
      .def(
          "query",
          [](VcfReader& self, const nucleus::genomics::v1::Range& region) {
            auto cpp_result = self.Query(region);
            auto ret0 = py::cast(std::move(cpp_result));
            auto postproc = py::module_::import(
                "third_party.nucleus.io.clif_postproc");
            return postproc.attr("WrappedVariantIterable")(ret0);
          },
          py::arg("region"))
      .def(
          "from_string",
          [](::nucleus::VcfReader* self, const absl::string_view& vcf_line) {
            nucleus::genomics::v1::Variant v;
            StatusOr<bool> cpp_result = self->FromStringPython(vcf_line, &v);
            py::object ret0 = py::cast(std::move(cpp_result));
            py::object ret1 = py::cast(std::move(v));
            auto postproc = py::module_::import(
                "third_party.nucleus.io.clif_postproc");
            return postproc.attr("ValueErrorOnFalse")(ret0, ret1);
          },
          py::arg("vcf_line"))
      .def("__enter__", [](py::object self) { return self; })
      .def("__exit__", [](VcfReader& self, py::args) { return self.Close(); })
      .def_property_readonly("header", &VcfReader::Header);
  py::classh<VariantIterable>(m, "VariantIterable")
      .def("PythonNext", &VariantIterable::PythonNext, py::arg("variant"))
      .def("Release", &VariantIterable::Release)
      .def("__enter__", [](py::object self) { return self; })
      .def("__exit__", &VariantIterable::PythonExit);
}
