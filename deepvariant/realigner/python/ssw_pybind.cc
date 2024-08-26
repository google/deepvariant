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

#include "deepvariant/realigner/ssw.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

namespace py = pybind11;

PYBIND11_MODULE(ssw, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::learning::genomics::deepvariant;  // NOLINT

  py::classh<Aligner>(m, "Aligner")
      .def(py::init<>())
      .def(py::init<uint8_t, uint8_t, uint8_t, uint8_t>())
      .def_static(
          "construct",
          [](uint8_t match_score, uint8_t mismatch_penalty,
             uint8_t gap_opening_penalty, uint8_t gap_extending_penalty) {
            return std::make_unique<Aligner>(std::move(match_score),
                                             std::move(mismatch_penalty),
                                             std::move(gap_opening_penalty),
                                             std::move(gap_extending_penalty));
          },
          py::arg("match_score"), py::arg("mismatch_penalty"),
          py::arg("gap_opening_penalty"), py::arg("gap_extending_penalty"))
      .def("set_reference_sequence", &Aligner::SetReferenceSequence,
           py::arg("seq"))
      .def("align", [](Aligner* self, const string& query, const Filter& filter,
                       int maskLen) {
        Alignment alignment;
        int cpp_result = self->Align(query, filter, maskLen, &alignment);
        py::object ret1 = py::cast(std::move(alignment));
        auto postproc =
            py::module_::import("third_party.nucleus.io.clif_postproc");
        py::object ret0 =
            postproc.attr("ValueErrorOnInaccurate")(cpp_result, ret1);
        return ret0;
      });

  py::classh<Filter>(m, "Filter")
      .def(py::init())
      .def(py::init<const bool&, const bool&, uint16_t, uint16_t>())
      .def_static(
          "construct",
          [](bool pos, bool cigar, uint16_t score, uint16_t dis) {
            return std::make_unique<Filter>(std::move(pos), std::move(cigar),
                                            std::move(score), std::move(dis));
          },
          py::arg("pos"), py::arg("cigar"), py::arg("int"), py::arg("dis"))
      .def_readwrite("report_begin_position", &Filter::report_begin_position)
      .def_readwrite("report_cigar", &Filter::report_cigar)
      .def_readwrite("score_filter", &Filter::score_filter)
      .def_readwrite("distance_filter", &Filter::distance_filter);

  py::classh<Alignment>(m, "Alignment")
      .def_readwrite("sw_score", &Alignment::sw_score)
      .def_readwrite("sw_score_next_best", &Alignment::sw_score_next_best)
      .def_readwrite("ref_begin", &Alignment::ref_begin)
      .def_readwrite("ref_end", &Alignment::ref_end)
      .def_readwrite("query_begin", &Alignment::query_begin)
      .def_readwrite("query_end", &Alignment::query_end)
      .def_readwrite("ref_end_next_best", &Alignment::ref_end_next_best)
      .def_readwrite("mismatches", &Alignment::mismatches)
      .def_property(
          "cigar_string",
          [](const Alignment* self) { return py::bytes(self->cigar_string); },
          [](Alignment* self, const std::string& value) {
            self->cigar_string = value;
          });
}
