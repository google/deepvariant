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

#include "tensorflow/core/framework/op.h"  // copybara: quotes ""
#include "tensorflow/core/framework/shape_inference.h"  // copybara: quotes ""

namespace learning {
namespace genomics {
namespace deepvariant {

REGISTER_OP("StreamExamplesInit")
    .Input("shm_dir: string")
    .Input("num_shards: int32")
    .Output("resource: resource")
    .Attr("container: string = ''")
    .Attr("shared_name: string = ''")
    .SetShapeFn([](tensorflow::shape_inference::InferenceContext* c) {
      c->set_output(0, c->Scalar());
      return tensorflow::Status();
    });

REGISTER_OP("StreamExamplesNext")
    .Input("input: resource")
    .Input("index: int64")
    .Output("image: string")
    .Output("variant: string")
    .Output("alt_allele_idx: string")
    .SetShapeFn([](tensorflow::shape_inference::InferenceContext* c) {
      c->set_output(0, c->MakeShape({c->UnknownDim()}));
      c->set_output(1, c->MakeShape({c->UnknownDim()}));
      c->set_output(2, c->MakeShape({c->UnknownDim()}));
      return tensorflow::Status();
    });

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
