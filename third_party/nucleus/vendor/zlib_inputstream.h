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
 */

#ifndef THIRD_PARTY_NUCLEUS_VENDOR_ZLIB_INPUTSTREAM_H_
#define THIRD_PARTY_NUCLEUS_VENDOR_ZLIB_INPUTSTREAM_H_

#ifndef TENSORFLOW_LIB_IO_ZLIB_INPUTSTREAM_H_
#define TENSORFLOW_LIB_IO_ZLIB_INPUTSTREAM_H_

#include <zlib.h>

#include <string>

#include "third_party/nucleus/vendor/zlib_compression_options.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/io/inputstream_interface.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/platform/macros.h"
#include "tensorflow/core/platform/types.h"

namespace tensorflow {
namespace io {

// An ZlibInputStream provides support for reading from a stream compressed
// using zlib (http://www.zlib.net/). Buffers the contents of the file.
//
// A given instance of an ZlibInputStream is NOT safe for concurrent use
// by multiple threads
class ZlibInputStream : public InputStreamInterface {
 public:
  // Create a ZlibInputStream for `input_stream` with a buffer of size
  // `input_buffer_bytes` bytes for reading contents from `input_stream` and
  // another buffer with size `output_buffer_bytes` for caching decompressed
  // contents. Does *not* take ownership of "input_stream".
  ZlibInputStream(InputStreamInterface* input_stream, size_t input_buffer_bytes,
                  size_t output_buffer_bytes,
                  const ZlibCompressionOptions& zlib_options);

  ~ZlibInputStream();

  // Reads bytes_to_read bytes into *result, overwriting *result.
  //
  // Return Status codes:
  // OK:           If successful.
  // OUT_OF_RANGE: If there are not enough bytes to read before
  //               the end of the stream.
  // ABORTED:      If inflate() fails, we return the error code with the
  //               error message in `z_stream_->msg`.
  // others:       If reading from stream failed.
  Status ReadNBytes(int64 bytes_to_read, string* result) override;

  int64 Tell() const override;

  Status Reset() override;

 private:
  void InitZlibBuffer();

  InputStreamInterface* input_stream_;  // Not owned
  size_t input_buffer_capacity_;        // Size of z_stream_input_
  size_t output_buffer_capacity_;       // Size of z_stream_output_
  char* next_unread_byte_;              // Next unread byte in z_stream_output_

  // Buffer for storing contents read from compressed stream.
  // redacted
  // the implementation.
  std::unique_ptr<Bytef[]> z_stream_input_;

  // Buffer for storing inflated contents of `input_stream_`.
  std::unique_ptr<Bytef[]> z_stream_output_;

  ZlibCompressionOptions const zlib_options_;

  // Configuration passed to `inflate`.
  //
  // z_stream_->next_in:
  //   Next byte to de-compress. Points to some byte in z_stream_input_ buffer.
  // z_stream_->avail_in:
  //   Number of bytes available to be decompressed at this time.
  // z_stream_->next_out:
  //   Next byte to write de-compressed data to. Points to some byte in
  //   z_stream_output_ buffer.
  // z_stream_->avail_out:
  //   Number of free bytes available at write location.
  std::unique_ptr<z_stream> z_stream_;

  // Reads data from `input_stream_` and tries to fill up `z_stream_input_` if
  // enough unread data is left in `input_stream_`.
  //
  // Looks up z_stream_->next_in to check how much data in z_stream_input_
  // has already been read. The used data is removed and new data is added to
  // after any unread data in z_stream_input_.
  // After this call z_stream_->next_in points to the start of z_stream_input_
  // and z_stream_->avail_in stores the number of readable bytes in
  // z_stream_input_.
  //
  // Returns OutOfRange error if NO data could be read from stream. Note that
  // this won't return an OutOfRange if there wasn't sufficient data in stream
  // to completely fill up z_stream_input_.
  Status ReadFromStream();

  // Calls `inflate()` and returns DataLoss Status if it failed.
  Status Inflate();

  // Starts reading bytes at `next_unread_byte_` till either `bytes_to_read`
  // bytes have been read or `z_stream_->next_out` is reached.
  // Returns the number of bytes read and advances the `next_unread_byte_`
  // pointer to the next location to read from.
  size_t ReadBytesFromCache(size_t bytes_to_read, string* result);

  // The number of unread bytes in z_stream_output_.
  //
  // z_stream_output_  -->
  //
  // [RRRRRRRRRRRRRRRRRRUUUUUUUUUUUUUU000000000000000000]
  //                    ^             ^
  //           next_unread_byte_    z_stream_->next_out
  //
  // R: Read bytes
  // U: Unread bytes
  // 0: garbage bytes where new output will be written
  //
  // Returns the size of [next_unread_byte_, z_stream_->next_out)
  size_t NumUnreadBytes() const;

  // Number of *uncompressed* bytes that have been read from this stream.
  int64 bytes_read_;

  TF_DISALLOW_COPY_AND_ASSIGN(ZlibInputStream);
};

}  // namespace io
}  // namespace tensorflow

#endif  // TENSORFLOW_LIB_IO_ZLIB_INPUTSTREAM_H_
#endif  // THIRD_PARTY_NUCLEUS_VENDOR_ZLIB_INPUTSTREAM_H_
