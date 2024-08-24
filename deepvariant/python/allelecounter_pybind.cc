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

#include "deepvariant/allelecounter.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

namespace py = pybind11;

PYBIND11_MODULE(allelecounter, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::learning::genomics::deepvariant;  // NOLINT

  py::classh<AlleleCounter>(m, "AlleleCounter")
      // I have no idea whether these are right. They seem to build at least.
      .def_static(
          "Default",
          [](const nucleus::GenomeReference* arg0,
             const nucleus::genomics::v1::Range arg1,
             const nucleus::genomics::v1::Range arg2,
             const std::vector<int>& arg3, const AlleleCounterOptions& arg4) {
            return std::make_unique<AlleleCounter>(
                arg0, std::move(arg1), std::move(arg2), std::move(arg3),
                std::move(arg4));
          })
      .def(py::init<const nucleus::GenomeReference*,
                    const nucleus::genomics::v1::Range&,
                    const std::vector<int>&, const AlleleCounterOptions&>())
      .def("add", &AlleleCounter::AddPython, py::arg("read"), py::arg("sample"))
      .def("NormalizeAndAddPython",
           [](AlleleCounter* self,
              ::nucleus::ConstProtoPtr<const ::nucleus::genomics::v1::Read>
                  wrapped,
              const string& sample) {
             int read_shift;
             std::unique_ptr<std::vector<nucleus::genomics::v1::CigarUnit>>
                 cpp_result =
                     self->NormalizeAndAddPython(wrapped, sample, &read_shift);
             py::object ret0 = py::cast(std::move(*cpp_result));
             py::object ret1 = py::cast(std::move(read_shift));
             return py::make_tuple(ret0, ret1);
           })
      .def("summary_counts", &AlleleCounter::SummaryCounts)
      .def("counts", &AlleleCounter::Counts);
}
