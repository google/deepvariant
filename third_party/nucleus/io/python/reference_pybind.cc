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

#include "third_party/nucleus/io/reference.h"
#include "third_party/pybind11/include/pybind11/chrono.h"
#include "third_party/pybind11/include/pybind11/complex.h"
#include "third_party/pybind11/include/pybind11/functional.h"
#include "third_party/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(reference, m) {
  using namespace ::nucleus;  // NOLINT
  py::class_<GenomeReference>(m, "GenomeReference")
      .def_property_readonly("contig_names", &GenomeReference::ContigNames)
      .def("contig", &GenomeReference::Contig, py::arg("chrom"))
      .def("bases", &GenomeReference::GetBases, py::arg("region"))
      .def("iterate", &GenomeReference::Iterate)
      .def("has_contig", &GenomeReference::HasContig, py::arg("contig_name"))
      .def("is_valid_interval", &GenomeReference::IsValidInterval,
           py::arg("region"))
      .def("__enter__", [](py::object self) { return self; })
      .def("__exit__", &GenomeReference::Close);

  py::class_<InMemoryFastaReader>(m, "InMemoryFastaReader")
      .def_static("create", &InMemoryFastaReader::Create, py::arg("contigs"),
                  py::arg("seqs"))
      .def_property_readonly("reference_sequences",
                             &InMemoryFastaReader::ReferenceSequences);

  py::class_<IndexedFastaReader>(m, "IndexedFastaReader")
      .def_static(
          "from_file",
          py::overload_cast<const string&, const string&,
                            const nucleus::genomics::v1::FastaReaderOptions&,
                            int>(&IndexedFastaReader::FromFile),
          py::arg("fasta_path"), py::arg("fai_path"), py::arg("options"),
          py::arg("cache_size_bases") =
              INDEXED_FASTA_READER_DEFAULT_CACHE_SIZE);

  py::class_<UnindexedFastaReader>(m, "UnindexedFastaReader")
      .def_static("from_file", &UnindexedFastaReader::FromFile,
                  py::arg("fasta_path"));

  py::class_<GenomeReferenceRecordIterable>(m, "GenomeReferenceRecordIterable")
      .def("Next", &GenomeReferenceRecordIterable::PythonNext)
      .def("Release", &GenomeReferenceRecordIterable::Release)
      .def("__enter__", [](py::object self) { return self; })
      .def("__exit__", &GenomeReferenceRecordIterable::PythonExit);
}
