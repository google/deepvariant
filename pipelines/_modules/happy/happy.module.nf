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

process run_happy {
  /*
    Runs hap.py
  */

  tag { "${output_fname}" }

  label 'xl'

  errorStrategy 'retry'
  maxRetries 2
  container 'jmcdani20/hap.py:v0.3.12'
  shell = ['/bin/bash', '-euo', 'pipefail']

  publishDir "${params.output_dir}/${uid}", pattern: "*happy*"

  input:
      tuple val(uid), \
            val(sample), \
            path(input_vcf), \
            path(input_vcf_index), \
            val(dataset_name), \
            path(truth_vcf), \
            path(truth_vcf_index), \
            path(truth_bed), \
            path(reference), \
            path(reference_fai), \
            val(region_str), \
            path(region_file), \
            val(haploid_contigs), \
            path(par_regions_bed)

  output:
      tuple path("${output_fname}.happy.summary.csv"), \
            path("${output_fname}.happy.extended.csv"), \
            path("${output_fname}.happy.vcf.gz"), \
            path("${output_fname}.happy.vcf.gz.tbi")
      tuple path("${output_fname}.happy.summary.csv"), \
            path("${output_fname}.happy.extended.csv"),
            path("${output_fname}.info.tsv"), emit: summary



  script:
    // Specifies regions by chr:start-end or chr.
    location_flag = region_str ? "--location ${region_str.join(",")}" : ""

    if (region_file.getName().startsWith("NO_FILE") ? "" : region_file) {
      target_regions_flag = "--target-regions ${region_file}"
    } else {
      target_regions_flag = ""
    }

    if (truth_bed.getName().startsWith("NO_FILE") ? "" : truth_bed) {
      fp_flag = "--false-positives ${truth_bed}"
    } else {
      fp_flag = ""
    }

    // If dataset_name is specified, use it as part of the output file name.
    output_fname = dataset_name ? "${uid}_${sample}_${dataset_name}" : "${uid}_${sample}"

    gender_flag = ""
    if (haploid_contigs && haploid_contigs.contains("chrY")) {
      gender_flag = "--gender male"
    } else if (haploid_contigs && haploid_contigs.contains("chrX")) {
      gender_flag = "--gender female"
    }

  """
  /opt/hap.py/bin/hap.py \\
      ${truth_vcf} \\
      ${input_vcf} \\
      ${fp_flag} \\
      ${location_flag} \\
      ${target_regions_flag} \\
      --reference ${reference} \\
      --report-prefix "${output_fname}.happy" \\
      --threads ${task.cpus} \\
      --engine=vcfeval \\
      --preprocess-truth \\
      ${gender_flag} \\
      --pass-only \\
      --verbose

  echo -e "${uid}\t${sample}\t${dataset_name}\t${output_fname}.happy.summary.csv" > ${output_fname}.info.tsv
  """

  stub:
    output_fname = dataset_name ? "${uid}_${sample}_${dataset_name}" : "${uid}_${sample}"
    """
    touch ${output_fname}.happy.summary.csv
    touch ${output_fname}.happy.extended.csv
    touch ${output_fname}.happy.vcf.gz
    touch ${output_fname}.happy.vcf.gz.tbi
    touch ${output_fname}.info.tsv
    """

}
