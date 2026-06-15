/*
 * Copyright 2026 Google LLC.
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
/* Aggregate runtimes and happy summary hashes */

process aggregate {
  /*
  Aggregates runtimes from multiple trials and hashes happy outputs to ensure
  deterministic results.

  It is possible to run without inference, in which case no happy results will
  be output.
  */

  label 'sm'

  container 'google/deepvariant:1.10.0'

  publishDir "${params.output_dir}"

 input:
    path(aggregate_inputs)

  output:
    path("runtimes.tsv"), optional: true
    path("runtimes.md"), optional: true
    path("md5sums.tsv")
    path("md5sums.md")
    path("happy.summary.tsv"), optional: true
    path("happy.summary.md"), optional: true
    path("happy.error_counts.tsv"), optional: true
    path("happy.error_counts.md"), optional: true
    path("happy.extended.tsv"), optional: true

  script:
    """
    aggregate.py
    """

  stub:
    """
    touch runtimes.tsv
    touch runtimes.md
    touch md5sums.tsv
    touch md5sums.md
    touch happy.summary.tsv
    touch happy.summary.md
    touch happy.error_counts.tsv
    touch happy.error_counts.md
    touch happy.extended.tsv
    """
}
