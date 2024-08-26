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

#include "deepvariant/realigner/fast_pass_aligner.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

namespace py = pybind11;
PYBIND11_MODULE(fast_pass_aligner, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::learning::genomics::deepvariant;  // NOLINT

  py::classh<FastPassAligner>(m, "FastPassAligner")
      .def(py::init<>())
      .def("set_reference", &FastPassAligner::set_reference, py::arg("ref"))
      .def("set_ref_start", &FastPassAligner::set_ref_start, py::arg("chr"),
           py::arg("pos"))
      .def("set_haplotypes", &FastPassAligner::set_haplotypes,
           py::arg("hyplotypes"))
      .def("set_options", &FastPassAligner::set_options, py::arg("options"))
      .def("set_is_debug", &FastPassAligner::set_is_debug, py::arg("isDebug"))
      .def("set_normalize_reads", &FastPassAligner::set_normalize_reads,
           py::arg("normalizeReads"))
      .def("set_debug_read_id", &FastPassAligner::set_debug_read_id,
           py::arg("readId"))
      .def("set_ref_prefix_len", &FastPassAligner::set_ref_prefix_len,
           py::arg("ref_prefix_len"))
      .def("set_ref_suffix_len", &FastPassAligner::set_ref_suffix_len,
           py::arg("set_ref_suffix_len"))
      .def("realign_reads",
           [](FastPassAligner* self,
              const std::vector<nucleus::genomics::v1::Read>& reads_param) {
             auto cpp_result = self->AlignReads(reads_param);
             return std::move(*cpp_result);
             // TODO): possible that `return *cpp_result` could work.
             // bug doing a move to be safe.
             // https://sigcpp.github.io/2020/06/08/return-value-optimization
           });
}
