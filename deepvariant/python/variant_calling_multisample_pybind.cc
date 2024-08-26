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

#include "deepvariant/variant_calling_multisample.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

namespace py = pybind11;

PYBIND11_MODULE(variant_calling_multisample, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::learning::genomics::deepvariant;  // NOLINT

  py::classh<multi_sample::VariantCaller>(m, "VariantCaller")
      .def(py::init<const VariantCallerOptions&>(), py::arg("options"))
      .def("calls_from_allele_counts",
           [](const multi_sample::VariantCaller& self,
              const std::unordered_map<std::string, AlleleCounter*>&
                  allele_counts_wrapper,
              const std::string& target_sample) {
             auto cpp_result = self.CallsFromAlleleCounts(allele_counts_wrapper,
                                                          target_sample);
             return cpp_result;
           })
      .def("call_positions_from_allele_counts",
           [](const multi_sample::VariantCaller& self,
              const std::unordered_map<std::string, AlleleCounter*>&
                  allele_counters,
              const std::string& target_sample) {
             auto cpp_result = self.CallPositionsFromAlleleCounts(
                 allele_counters, target_sample);
             return cpp_result;
           });
}
