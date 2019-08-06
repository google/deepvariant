/*
 * Copyright 2019 Google LLC.
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

#include "third_party/nucleus/io/tabix_indexer.h"

#include "third_party/nucleus/io/hts_path.h"
#include "third_party/nucleus/platform/types.h"
#include "tensorflow/core/lib/core/errors.h"

namespace nucleus {

namespace tf = tensorflow;

tf::Status TbxIndexBuild(const string& path) {
  int val = tbx_index_build_x(path, 0, &tbx_conf_vcf);
  if (val < 0) {
    LOG(WARNING) << "Return code: " << val << "\nFile path: " << path;
    return tf::errors::Internal("Failure to write tabix index.");
  }
  return tf::Status::OK();
}

tf::Status CSIIndexBuild(string path, int min_shift) {
  // Create a index file in CSI format by setting min_shift as a non-zero value.
  int val = tbx_index_build_x(path, min_shift, &tbx_conf_vcf);
  if (val < 0) {
    LOG(WARNING) << "Return code: " << val << "\nFile path: " << path;
    return tf::errors::Internal("Failure to write CSI index.");
  }
  return tf::Status::OK();
}

}  // namespace nucleus
