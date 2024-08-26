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

#include "third_party/pybind11/include/pybind11/cast.h"
#if true  // Trick to stop tooling from moving the #include around.
// MUST appear before any standard headers are included.
#include <pybind11/pybind11.h>
#endif

#include <pybind11/stl.h>

#include "deepvariant/pileup_image_native.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_status.h"
#include "third_party/nucleus/core/python/type_caster_nucleus_statusor.h"
#include "third_party/nucleus/util/python/type_caster_nucleus_proto_ptr.h"
#include "third_party/pybind11_protobuf/native_proto_caster.h"

// TODO: Use <pybind11/numpy.h> instead.
// These includes (probably) need to stay here (do not move them up).
#include "numpy/arrayobject.h"
#include "numpy/ndarrayobject.h"

namespace py = pybind11;

namespace {

py::object ConvertImageRowToNumPyArray(
    std::unique_ptr<learning::genomics::deepvariant::ImageRow> img_row) {
  if (!img_row) {
    return py::none();
  }
  npy_intp dims[]{1, img_row->Width(), img_row->num_channels};
  PyArrayObject* res =
      reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(3, dims, NPY_UBYTE));
  if (res == nullptr) {
    throw py::error_already_set();
  }

  unsigned char* data = reinterpret_cast<unsigned char*>(PyArray_DATA(res));
  if (data == nullptr) {
    Py_DECREF(res);
    throw py::error_already_set();
  }

  unsigned char* cur = data;
  for (int i = 0; i < img_row->Width(); i++) {
    if (!img_row->channel_data.empty()) {
      // Iterate over channels here and fill data...
      for (int j = 0; j < img_row->channel_data.size(); j++) {
        *cur++ = img_row->channel_data[j][i];
      }
    }
  }

  py::object retval = py::reinterpret_steal<py::object>(PyArray_Return(res));
  if (!retval) {
    throw py::error_already_set();
  }
  return retval;
}

}  // namespace

PYBIND11_MODULE(pileup_image_native, m) {
  pybind11_protobuf::ImportNativeProtoCasters();
  using namespace ::learning::genomics::deepvariant;  // NOLINT

  py::classh<PileupImageEncoderNative>(m, "PileupImageEncoderNative")
      .def(py::init<const PileupImageOptions&>(), py::arg("options"))
      .def("all_channels_enum", &PileupImageEncoderNative::AllChannelsEnum,
           py::arg("alt_aligned_pileup"))
      .def("build_pileup_for_one_sample",
           &PileupImageEncoderNative::BuildPileupForOneSamplePython,
           py::arg("dv_call"), py::arg("ref_bases"), py::arg("reads"),
           py::arg("image_start_pos"), py::arg("alt_alleles"),
           py::arg("sample_options"))
      .def(
          "encode_read",
          [](PileupImageEncoderNative* self,
             const nucleus::ConstProtoPtr<const DeepVariantCall>&
                 wrapped_dv_call,
             const string& ref_bases,
             const nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>&
                 wrapped_read,
             int image_start_pos, const std::vector<std::string>& alt_alleles) {
            std::unique_ptr<ImageRow> cpp_result =
                self->EncodeReadPython(wrapped_dv_call, ref_bases, wrapped_read,
                                       image_start_pos, alt_alleles);
            return ConvertImageRowToNumPyArray(std::move(cpp_result));
          },
          py::arg("dv_call"), py::arg("ref_bases"), py::arg("read"),
          py::arg("image_start_pos"), py::arg("alt_alleles"))
      .def("encode_reference", [](PileupImageEncoderNative* self,
                                  const string& ref_bases) {
        std::unique_ptr<ImageRow> cpp_result = self->EncodeReference(ref_bases);
        return ConvertImageRowToNumPyArray(std::move(cpp_result));
      });
}
