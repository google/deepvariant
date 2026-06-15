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
/*
    DeepVariant module
    =============================================
*/

/*-----------------
POST_PROCESSING PROCESS
------------------*/
process run_post_processing {
  /*
  Runs DeepVariant using the `run_deepvariant` command.
  */

  label 'xl'
  tag { "${uid}" }
  stageInMode { task.executor == 'local' ? 'symlink' : 'copy' }
  container { "${docker_image}" }

  publishDir "${params.output_dir}/${uid}"

  input:
    tuple val(uid), \
          val(sample), \
          path(bam), \
          path(bam_index), \
          val(model_type), \
          path(model_files), \
          path(small_model_file), \
          val(extra_args_flags), \
          val(disable_small_model_flag), \
          val(keep_intermediate_results), \
          val(make_examples_extra_args_flag), \
          val(call_variants_extra_args_flag), \
          val(postprocess_variants_extra_args_flag), \
          val(docker_image), \
          val(region_flag), \
          path(region_file), \
          path(ref), \
          path(ref_fai), \
          val(customized_model_flag), \
          val(customized_small_model_flag), \
          val(num_instances), \
          val(shards_per_instance), \
          path(call_variants_output_files, stageAs: 'intermediate_results/*'), \
          path(make_examples_call_variant_outputs, stageAs: 'intermediate_results/*'), \
          path(gvcf_files, stageAs: 'intermediate_results/*')
  output:
    tuple val(uid), \
          val(sample), \
          path("${uid}_${sample}.deepvariant.vcf.gz"), \
          path("${uid}_${sample}.deepvariant.vcf.gz.tbi"), emit: 'to_happy'
    tuple path("${uid}_${sample}.deepvariant.g.vcf.gz"), \
          path("${uid}_${sample}.deepvariant.g.vcf.gz.tbi")
    tuple path("${uid}_${sample}.md5sum.txt"), emit: 'to_aggregate'

  script:
  """
    run_deepvariant \\
      --steps=postprocess_variants \\
      --model_type="${model_type}" \\
      ${customized_model_flag} \\
      ${customized_small_model_flag} \\
      --ref="${ref}" \\
      --reads=${bam} \\
      --sample_name=${sample} \\
      --output_vcf=${uid}_${sample}.deepvariant.vcf.gz \\
      --output_gvcf=${uid}_${sample}.deepvariant.g.vcf.gz \\
      --num_shards ${shards_per_instance * num_instances} \\
      --postprocess_cpus ${task.cpus} \\
      ${region_flag} \\
      ${extra_args_flags} \\
      ${disable_small_model_flag} \\
      ${make_examples_extra_args_flag} \\
      ${call_variants_extra_args_flag} \\
      ${postprocess_variants_extra_args_flag} \\
      --novcf_stats_report \\
      --intermediate_results_dir=intermediate_results \\
      --logging_dir=. | tee deepvariant.call_variants.log

      md5sum ${uid}_${sample}.deepvariant.vcf.gz > ${uid}_${sample}.md5sum.txt
  """

}
