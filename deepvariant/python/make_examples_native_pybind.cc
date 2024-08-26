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
 */

#if true  // Trick to stop tooling from moving the #include around.
// MUST appear before any standard headers are included.
#include <pybind11/pybind11.h>
#endif

#include <pybind11/stl.h>

#include "deepvariant/make_examples_native.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

namespace py = pybind11;

PYBIND11_MODULE(make_examples_native, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::learning::genomics::deepvariant;  // NOLINT

  py::classh<ExamplesGenerator>(m, "ExamplesGenerator")
      .def(
          py::init<const MakeExamplesOptions&,
                   const std::unordered_map<std::string, std::string>&, bool>(),
          py::arg("options"), py::arg("example_filenames"),
          py::arg("test_mode") = false)
      // TODO Do I need to actually implement or call the constructor??

      .def("write_examples_in_region",
           [](ExamplesGenerator* self,
              const std::vector<nucleus::ConstProtoPtr<DeepVariantCall>>&
                  candidates,
              const std::vector<std::vector<nucleus::ConstProtoPtr<
                  ::nucleus::genomics::v1::Read>>>& reads_per_sample,
              const std::vector<int>& sample_order, const std::string& role,
              const std::vector<float>& mean_coverage_per_sample) {
             std::vector<int> image_shape;
             std::unordered_map<std::string, int> cpp_result =
                 self->WriteExamplesInRegion(
                     candidates, reads_per_sample, sample_order, role,
                     mean_coverage_per_sample, &image_shape);
             py::object ret0 = py::cast(std::move(cpp_result));
             py::object ret1 = py::cast(std::move(image_shape));
             return py::make_tuple(ret0, ret1);
           })
      .def("append_label", &ExamplesGenerator::AppendLabel, py::arg("label"))
      .def("signal_shard_finished", &ExamplesGenerator::SignalShardFinished);

    py::classh<VariantLabel>(m, "VariantLabel")
        .def(py::init<bool, const nucleus::genomics::v1::Variant&,
                     const std::vector<int>&, bool>())
      .def_readwrite("is_confident", &VariantLabel::is_confident)
      .def_readwrite("variant", &VariantLabel::variant)
      .def_readwrite("genotype", &VariantLabel::genotype)
      .def_readwrite("is_denovo", &VariantLabel::is_denovo);

    py::classh<CustomizedClassesLabel, VariantLabel>(m,
                                                     "CustomizedClassesLabel")
        .def(py::init<>())
        .def(py::init<bool, const nucleus::genomics::v1::Variant&,
                      const nucleus::genomics::v1::Variant&,
                      const std::unordered_map<std::string, int>&,
                      const std::string&>())
        .def_readwrite("is_confident", &CustomizedClassesLabel::is_confident)
        .def_readwrite("variant", &CustomizedClassesLabel::variant)
        .def_readwrite("truth_variant", &CustomizedClassesLabel::truth_variant)
        .def_readwrite("classes_dict", &CustomizedClassesLabel::classes_dict)
        .def_readwrite("info_field_name",
                       &CustomizedClassesLabel::info_field_name);
}
