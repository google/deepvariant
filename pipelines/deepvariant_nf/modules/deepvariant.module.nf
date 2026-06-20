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

include { set_flag; get_make_examples_extra_args_flag; } from '../../_modules/util/util.module.nf'


/*-----------------
DEEPVARIANT PROCESS
------------------*/
process run_deepvariant {
  /*
  Runs DeepVariant using the `run_deepvariant` command.
  */

  label params.label ?: (params.gpu ? 'gpu' : 'xl')
  tag { "${uid}" }
  stageInMode { task.executor == 'local' ? 'symlink' : 'copy' }
  container { "${docker_image}" }

  publishDir "${params.output_dir}/${uid}"

  input:
    tuple val(uid),
          val(sample), \
          path(bam), \
          path(bam_index), \
          val(model_type), \
          path(model_files), \
          path(small_model_file), \
          val(extra_args), \
          val(disable_small_model), \
          val(keep_intermediate_results), \
          val(make_examples_extra_args), \
          val(call_variants_extra_args), \
          val(postprocess_variants_extra_args), \
          val(docker_image), \
          val(region_str), \
          path(region_file), \
          path(ref), \
          path(ref_fai), \
          val(haploid_contigs), \
          path(par_regions_bed)

  output:
    tuple val(uid), \
          val(sample), \
          path("${uid}_${sample}.deepvariant.vcf.gz"), \
          path("${uid}_${sample}.deepvariant.vcf.gz.tbi"), emit: 'to_happy'
    tuple path("${uid}_${sample}.deepvariant.g.vcf.gz"), \
          path("${uid}_${sample}.deepvariant.g.vcf.gz.tbi")
    tuple path("${uid}_${sample}.runtimes.tsv"), path("${uid}_${sample}.md5sum.txt"), emit: 'to_aggregate'
    tuple path("make_examples.log"), \
          path("call_variants.log"), \
          path("postprocess_variants.log"), \
          path("vcf_stats_report.log"), \
          path("deepvariant.log")
    path("intermediate_results/*"), optional: true
    path("${uid}_${sample}.deepvariant.visual_report.html")

  script:
    // Format extra args flags.
    extra_args_flags = extra_args.join(" ")

    disable_small_model_flag = disable_small_model ? "--disable_small_model" : "--nodisable_small_model"

    // Resolve region flag:
    region_bed = region_file.getName().startsWith("NO_FILE") ? "" : region_file
    all_regions="${region_str.join(" ")} ${region_bed}".trim()
    if (all_regions) {
      region_flag = "--regions=\"${all_regions}\""
    } else {
      region_flag = ""
    }

    // Resolve customized model flag:
    if (model_files[0].name.endsWith("savedmodel") || model_files[0].name.endsWith("saved_model.pb")) {
      customized_model_flag = "--customized_model=." // DeepVariant will check for ./saved_model.pb
    } else {
      customized_model_flag = model_files[0].getName().startsWith("NO_FILE") ? "" : "--customized_model=./model.ckpt"
    }

    // Resolve customized small model flag:
    customized_small_model_flag = small_model_file.getName().startsWith("NO_FILE") ? "" : "--customized_small_model=${small_model_file}/model.keras"

    // Intermediate results flag:
    keep_intermediate_results_flag = keep_intermediate_results ? "--intermediate_results_dir=intermediate_results" : ""

    make_examples_extra_args_flag = get_make_examples_extra_args_flag(model_files, make_examples_extra_args)
    call_variants_extra_args_flag = set_flag("call_variants_extra_args", call_variants_extra_args)
    postprocess_variants_extra_args_flag = set_flag("postprocess_variants_extra_args", postprocess_variants_extra_args)
    haploid_contigs_flag = set_flag("haploid_contigs", haploid_contigs)
    par_regions_bed_local = par_regions_bed.getName().startsWith("NO_FILE") ? "" : par_regions_bed
    par_regions_bed_flag = set_flag("par_regions_bed", par_regions_bed_local)

  """
    mv *.data-00000-of-00001 model.ckpt.data-00000-of-00001 2> /dev/null || true
    mv *.index model.ckpt.index 2> /dev/null || true

    run_deepvariant \\
      --model_type="${model_type}" \\
      ${customized_model_flag} \\
      ${customized_small_model_flag} \\
      --ref="${ref}" \\
      --reads=${bam} \\
      --sample_name=${sample} \\
      --output_vcf=${uid}_${sample}.deepvariant.vcf.gz \\
      --output_gvcf=${uid}_${sample}.deepvariant.g.vcf.gz \\
      --num_shards ${task.cpus} \\
      ${region_flag} \\
      ${extra_args_flags} \\
      ${disable_small_model_flag} \\
      ${make_examples_extra_args_flag} \\
      ${call_variants_extra_args_flag} \\
      ${postprocess_variants_extra_args_flag} \\
      ${haploid_contigs_flag} \\
      ${par_regions_bed_flag} \\
      --vcf_stats_report true \\
      ${keep_intermediate_results_flag} \\
      --logging_dir=. | tee deepvariant.log

    collect_runtimes.py "${uid}_${sample}"
    md5sum ${uid}_${sample}.deepvariant.vcf.gz > ${uid}_${sample}.md5sum.txt
  """

}



/*-------------------------
PANGENOME_AWARE_DEEPVARIANT
--------------------------*/

process run_pangenome_aware_deepvariant {
  /*
  A separate process is necessary here because we are using a different
  binary, and several different options.
  */

  label params.label ?: 'xxl_highmem'
  tag { "${uid}" }
  stageInMode { task.executor == 'local' ? 'symlink' : 'copy' }
  errorStrategy 'retry'
  maxRetries { params.debug ? 0 : 2 }

  container { "${docker_image}" }
  containerOptions { "--shm-size ${shared_memory_size_gb}gb" }

  publishDir "${params.output_dir}/${uid}"

  input:
    tuple val(uid),
          val(sample), \
          path(bam), \
          path(bam_index), \
          val(model_type), \
          path(model_files), \
          path(small_model_file), \
          val(extra_args), \
          val(disable_small_model), \
          val(keep_intermediate_results), \
          val(make_examples_extra_args), \
          val(call_variants_extra_args), \
          val(postprocess_variants_extra_args), \
          val(docker_image), \
          val(region_str), \
          path(region_file), \
          path(ref), \
          path(ref_fai), \
          val(haploid_contigs), \
          path(par_regions_bed), \
          val(shared_memory_size_gb), \
          path(pangenome_graph)

    output:
    tuple val(uid), \
          val(sample), \
          path("${uid}_${sample}.deepvariant.vcf.gz"), \
          path("${uid}_${sample}.deepvariant.vcf.gz.tbi"), emit: 'to_happy'
    tuple path("${uid}_${sample}.deepvariant.g.vcf.gz"), \
          path("${uid}_${sample}.deepvariant.g.vcf.gz.tbi")
    tuple path("${uid}_${sample}.runtimes.tsv"), path("${uid}_${sample}.md5sum.txt"), emit: 'to_aggregate'
    tuple path("make_examples_pangenome_aware_dv.log"), \
          path("call_variants.log"), \
          path("postprocess_variants.log"), \
          path("vcf_stats_report.log"), \
          path("deepvariant.log")
    path("intermediate_results/*"), optional: true

  script:

    // Format extra args flags.
    extra_args_flags = extra_args.join(" ")

    // Always disable small model for pangenome for now - we do not package
    // small models with the pangenome currently.
    // disable_small_model_flag = disable_small_model ? "--disable_small_model" : "--nodisable_small_model"
    disable_small_model_flag = "--disable_small_model"

    // Resolve region flag:
    region_bed = region_file.getName().startsWith("NO_FILE") ? "" : region_file
    all_regions="${region_str.join(" ")} ${region_bed}".trim()
    if (all_regions) {
      region_flag = "--regions=\"${all_regions}\""
    } else {
      region_flag = ""
    }

    // Resolve customized model flag:
    if (model_files[0].name.endsWith("savedmodel") || model_files[0].name.endsWith("saved_model.pb")) {
      customized_model_flag = "--customized_model=." // DeepVariant will check for ./saved_model.pb
    } else {
      customized_model_flag = model_files[0].getName().startsWith("NO_FILE") ? "" : "--customized_model=./model.ckpt"
    }

    // Resolve customized small model flag:
    customized_small_model_flag = small_model_file.getName().startsWith("NO_FILE") ? "" : "--customized_small_model=./model.keras"

    // Intermediate results flag:
    keep_intermediate_results_flag = keep_intermediate_results ? "--intermediate_results_dir=intermediate_results" : ""

    make_examples_extra_args_flag = get_make_examples_extra_args_flag(model_files, make_examples_extra_args)
    call_variants_extra_args_flag = set_flag("call_variants_extra_args", call_variants_extra_args)
    postprocess_variants_extra_args_flag = set_flag("postprocess_variants_extra_args", postprocess_variants_extra_args)
    haploid_contigs_flag = set_flag("haploid_contigs", haploid_contigs)
    par_regions_bed_local = par_regions_bed.getName().startsWith("NO_FILE") ? "" : par_regions_bed
    par_regions_bed_flag = set_flag("par_regions_bed", par_regions_bed_local)

    """
    export IS_PANGENOME_AWARE_DEEPVARIANT=1 # Used by get_shape_flags.sh
    mv *.data-00000-of-00001 model.ckpt.data-00000-of-00001 2> /dev/null || true
    mv *.index model.ckpt.index 2> /dev/null || true

    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \\
      --model_type="${model_type}" \\
      ${customized_model_flag} \\
      ${customized_small_model_flag} \\
      --ref="${ref}" \\
      --reads=${bam} \\
      --pangenome=${pangenome_graph} \\
      --output_vcf=${uid}_${sample}.deepvariant.vcf.gz \\
      --output_gvcf=${uid}_${sample}.deepvariant.g.vcf.gz \\
      --num_shards ${task.cpus} \\
      ${region_flag} \\
      ${extra_args_flags} \\
      ${disable_small_model_flag} \\
      ${make_examples_extra_args_flag} \\
      ${call_variants_extra_args_flag} \\
      ${postprocess_variants_extra_args_flag} \\
      ${haploid_contigs_flag} \\
      ${par_regions_bed_flag} \\
      --gbz_shared_memory_size_gb ${shared_memory_size_gb} \\
      --vcf_stats_report true \\
      ${keep_intermediate_results_flag} \\
      --logging_dir=. | tee deepvariant.log

    collect_runtimes.py "${uid}_${sample}"
    md5sum ${uid}_${sample}.deepvariant.vcf.gz > ${uid}_${sample}.md5sum.txt
    """
}
