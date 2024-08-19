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

namespace py = pybind11;

PYBIND11_MODULE(reference, m) {
  using namespace ::nucleus;
  py::class_<GenomeReference>(m, "GenomeReference")
      .def_property_readonly("contig_names", &GenomeReference::ContigNames)
      .def("contig", &GenomeReference::Contig, py::arg("chrom"))
      .def("bases", &GenomeReference::GetBases, py::arg("region"))
      .def("iterate", &GenomeReference::Iterate)
      .def("has_contig", &GenomeReference::HasContig, py::arg("contig_name"))
      .def("is_valid_interval", &GenomeReference::IsValidInterval,
           py::arg("region"))
      .def("__enter__", &GenomeReference::PythonEnter)
      .def("__exit__", &GenomeReference::Close);
  // I don't know what to do with this one:
  //  GenomeReference_class.def(
  //      "as_nucleus_Reader",
  //      [](::nucleus::GenomeReference* self) {
  //          return pybind11::capsule(static_cast<void *>(self));
  //      }
  //  );

  py::class_<InMemoryFastaReader>(m, "InMemoryFastaReader")
      .def_static("create", &InMemoryFastaReader::Create, py::arg("contigs"),
                  py::arg("seqs"))
      .def_property_readonly("reference_sequences",
                             &InMemoryFastaReader::ReferenceSequences);
  // I don't know what to do with this one:
  // InMemoryFastaReader_class.def(
  //     "as_nucleus_GenomeReference",
  //     [](::nucleus::InMemoryFastaReader* self) {
  //         return pybind11::capsule(static_cast<void *>(self));
  //     }
  // );

  // // I couldn't get this one to compile yet:
  // py::class_<IndexedFastaReader>(m, "IndexedFastaReader")
  //     .def_static("from_file", &IndexedFastaReader::FromFile,
  //                 py::arg("fasta_path"),
  //                 py::arg("fai_path"),
  //                 py::arg("options"),
  //                 py::arg("cache_size_bases") = 65536);

  // // And I'm not sure what to with these:
  //     IndexedFastaReader_class.def(
  //         "as_nucleus_GenomeReference",
  //         [](::nucleus::IndexedFastaReader* self) {
  //             return pybind11::capsule(static_cast<void *>(self));
  //         }
  //     );
  //     IndexedFastaReader_class.def(
  //         "as_nucleus_Reader",
  //         [](::nucleus::IndexedFastaReader* self) {
  //             return pybind11::capsule(static_cast<void *>(self));
  //         }
  //     );

  py::class_<UnindexedFastaReader>(m, "UnindexedFastaReader")
      .def_static("from_file", &UnindexedFastaReader::FromFile,
                  py::arg("fasta_path"));

  // I don't know what to do with these:
  // UnindexedFastaReader_class.def(
  //     "as_nucleus_GenomeReference",
  //     [](::nucleus::UnindexedFastaReader* self) {
  //         return pybind11::capsule(static_cast<void *>(self));
  //     }
  // );
  // UnindexedFastaReader_class.def(
  //     "as_nucleus_Reader",
  //     [](::nucleus::UnindexedFastaReader* self) {
  //         return pybind11::capsule(static_cast<void *>(self));
  //     }
  // );
}

//       using namespace ::nucleus;
//       pybind11::classh<::nucleus::Iterable< ::std::pair<
//       ::std::basic_string<char, ::std::char_traits<char>,
//       ::std::allocator<char>>, ::std::basic_string<char,
//       ::std::char_traits<char>, ::std::allocator<char>>>>>
//       GenomeReferenceRecordIterable_class(m, "GenomeReferenceRecordIterable",
//       pybind11::metaclass((PyObject*) &PyType_Type),
//       pybind11::release_gil_before_calling_cpp_dtor());
//       GenomeReferenceRecordIterable_class.def("Next", [](pybind11::object
//       self_py) {
//         auto self = pybind11::cast<::nucleus::Iterable< ::std::pair<
//         ::std::basic_string<char, ::std::char_traits<char>,
//         ::std::allocator<char>>, ::std::basic_string<char,
//         ::std::char_traits<char>, ::std::allocator<char>>>>*>(self_py);
//         pybind11::object ret0;
//         ::std::pair< ::std::basic_string<char, ::std::char_traits<char>,
//         ::std::allocator<char>>, ::std::basic_string<char,
//         ::std::char_traits<char>, ::std::allocator<char>>> ret1{};
//         {
//           pybind11::gil_scoped_release gil_release;
//           ::nucleus::StatusOr<bool> ret0_ = self->Next(&ret1);
//           {
//             pybind11::gil_scoped_acquire gil_acquire;
//             ret0 = pybind11::cast(std::move(ret0_),
//             pybind11::return_value_policy_pack(pybind11::return_value_policy_pack(std::vector<pybind11::return_value_policy_pack>({pybind11::return_value_policy::_clif_automatic}),
//             pybind11::return_value_policy::_clif_automatic)), self_py);
//             ::clif::ThrowErrorAlreadySetIfPythonErrorOccurred();
//           }
//         }
//         return std::make_tuple(ret0, pybind11::cast(std::move(ret1),
//         pybind11::return_value_policy_pack({pybind11::return_value_policy::_clif_automatic,
//         pybind11::return_value_policy::_clif_automatic}), self_py));
//       }, pybind11::return_value_policy::_clif_automatic);
//       GenomeReferenceRecordIterable_class.def("Release", [](pybind11::object
//       self_py) {
//         auto self = pybind11::cast<::nucleus::Iterable< ::std::pair<
//         ::std::basic_string<char, ::std::char_traits<char>,
//         ::std::allocator<char>>, ::std::basic_string<char,
//         ::std::char_traits<char>, ::std::allocator<char>>>>*>(self_py);
//         pybind11::object ret0;
//         {
//           pybind11::gil_scoped_release gil_release;
//           ::nucleus::Status ret0_ = self->Release();
//           {
//             pybind11::gil_scoped_acquire gil_acquire;
//             ret0 = pybind11::cast(std::move(ret0_),
//             pybind11::return_value_policy_pack(pybind11::return_value_policy::_clif_automatic),
//             self_py);
//             ::clif::ThrowErrorAlreadySetIfPythonErrorOccurred();
//           }
//         }
//         return ret0;
//       }, pybind11::return_value_policy::_clif_automatic);
//       GenomeReferenceRecordIterable_class.def("__enter__",
//       [](pybind11::object self_py) {
//         auto self = pybind11::cast<::nucleus::Iterable< ::std::pair<
//         ::std::basic_string<char, ::std::char_traits<char>,
//         ::std::allocator<char>>, ::std::basic_string<char,
//         ::std::char_traits<char>, ::std::allocator<char>>>>*>(self_py);
//         pybind11::object ret0;
//         {
//           pybind11::gil_scoped_release gil_release;
//           ::nucleus::Status ret0_ = self->PythonEnter();
//           {
//             pybind11::gil_scoped_acquire gil_acquire;
//             ret0 = pybind11::cast(std::move(ret0_),
//             pybind11::return_value_policy_pack(pybind11::return_value_policy::_clif_automatic),
//             self_py);
//             ::clif::ThrowErrorAlreadySetIfPythonErrorOccurred();
//           }
//         }
//         return self;
//       }, pybind11::return_value_policy::_clif_automatic);
//       GenomeReferenceRecordIterable_class.def("__exit__", [](pybind11::object
//       self_py, pybind11::args) {
//         auto self = pybind11::cast<::nucleus::Iterable< ::std::pair<
//         ::std::basic_string<char, ::std::char_traits<char>,
//         ::std::allocator<char>>, ::std::basic_string<char,
//         ::std::char_traits<char>, ::std::allocator<char>>>>*>(self_py);
//         pybind11::object ret0;
//         {
//           pybind11::gil_scoped_release gil_release;
//           ::nucleus::Status ret0_ = self->PythonExit();
//           {
//             pybind11::gil_scoped_acquire gil_acquire;
//             ret0 = pybind11::cast(std::move(ret0_),
//             pybind11::return_value_policy_pack(pybind11::return_value_policy::_clif_automatic),
//             self_py);
//             ::clif::ThrowErrorAlreadySetIfPythonErrorOccurred();
//           }
//         }
//         return pybind11::none();
//       }, pybind11::return_value_policy::_clif_automatic);
//       GenomeReferenceRecordIterable_class.def(
//           "as_nucleus_IterableBase",
//           [](::nucleus::Iterable< ::std::pair< ::std::basic_string<char,
//           ::std::char_traits<char>, ::std::allocator<char>>,
//           ::std::basic_string<char, ::std::char_traits<char>,
//           ::std::allocator<char>>>>* self) {
//               return pybind11::capsule(static_cast<void *>(self));
//           }
//       );
//       GenomeReferenceRecordIterable_class.def("__reduce_ex__",
//       ::clif_pybind11::ReduceExImpl, pybind11::arg("protocol")=-1);
//     }
// }
