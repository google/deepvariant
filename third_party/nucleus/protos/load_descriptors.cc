/*
 * Copyright 2018 Google Inc.
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

// This is a Python extension module that exists soley to load the
// C++ generated messages for Nucleus's protocol buffers BEFORE the
// Python code for those protocol buffers is imported.

#include <Python.h>

#include "third_party/nucleus/protos/bed.pb.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/fasta.pb.h"
#include "third_party/nucleus/protos/fastq.pb.h"
#include "third_party/nucleus/protos/gff.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

static PyMethodDef load_descriptors_methods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

#if PY_MAJOR_VERSION > 2
static struct PyModuleDef load_descriptors_module = {
  PyModuleDef_HEAD_INIT,
  "lib_load_descriptors",
  NULL,
  -1,
  load_descriptors_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif  /* Python 3 */

PyMODINIT_FUNC
#if PY_MAJOR_VERSION > 2
PyInit_lib_load_descriptors(void)
#else
initlib_load_descriptors(void)
#endif
{
  nucleus::genomics::v1::BedRecord().descriptor();
  nucleus::genomics::v1::BedHeader().descriptor();
  nucleus::genomics::v1::BedReaderOptions().descriptor();
  nucleus::genomics::v1::BedWriterOptions().descriptor();
  nucleus::genomics::v1::CigarUnit().descriptor();
  nucleus::genomics::v1::FastaRecord().descriptor();
  nucleus::genomics::v1::FastaReaderOptions().descriptor();
  nucleus::genomics::v1::FastaWriterOptions().descriptor();
  nucleus::genomics::v1::FastqRecord().descriptor();
  nucleus::genomics::v1::FastqReaderOptions().descriptor();
  nucleus::genomics::v1::FastqWriterOptions().descriptor();
  nucleus::genomics::v1::GffRecord().descriptor();
  nucleus::genomics::v1::GffHeader().descriptor();
  nucleus::genomics::v1::GffReaderOptions().descriptor();
  nucleus::genomics::v1::GffWriterOptions().descriptor();
  nucleus::genomics::v1::Position().descriptor();
  nucleus::genomics::v1::Range().descriptor();
  nucleus::genomics::v1::LinearAlignment().descriptor();
  nucleus::genomics::v1::Read().descriptor();
  nucleus::genomics::v1::SamHeader().descriptor();
  nucleus::genomics::v1::ReadGroup().descriptor();
  nucleus::genomics::v1::Program().descriptor();
  nucleus::genomics::v1::SamReaderOptions().descriptor();
  nucleus::genomics::v1::ReadRequirements().descriptor();
  nucleus::genomics::v1::ContigInfo().descriptor();
  nucleus::genomics::v1::ReferenceSequence().descriptor();
  nucleus::genomics::v1::Struct().descriptor();
  nucleus::genomics::v1::Value().descriptor();
  nucleus::genomics::v1::ListValue().descriptor();
  nucleus::genomics::v1::Variant().descriptor();
  nucleus::genomics::v1::VariantCall().descriptor();
  nucleus::genomics::v1::VcfHeader().descriptor();
  nucleus::genomics::v1::VcfFilterInfo().descriptor();
  nucleus::genomics::v1::VcfInfo().descriptor();
  nucleus::genomics::v1::VcfFormatInfo().descriptor();
  nucleus::genomics::v1::VcfStructuredExtra().descriptor();
  nucleus::genomics::v1::VcfExtra().descriptor();
  nucleus::genomics::v1::VcfReaderOptions().descriptor();
  nucleus::genomics::v1::VcfWriterOptions().descriptor();

  PyObject *m;

#if PY_MAJOR_VERSION > 2
  m = PyModule_Create(&load_descriptors_module);
  if (m == NULL) {
    return NULL;
  }
#else
  m = Py_InitModule("lib_load_descriptors", load_descriptors_methods);
  if (m == NULL) {
      return;
  }
#endif

#if PY_MAJOR_VERSION > 2
  return m;
#endif
}
