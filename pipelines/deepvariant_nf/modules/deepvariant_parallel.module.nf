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
include { is_real_file; set_flag; } from '../../_modules/util/util.module.nf'
include { run_make_examples; } from './make_examples.module.nf'
include { run_call_variants; } from './call_variants.module.nf'
include { run_post_processing; } from './post_processing.module.nf'

/**
 * Normalize inputs for DeepVariant parallel processing.
 */
def normalize_inputs(List row) {
    def (uid,
         sample,
         bam,
         bam_index,
         model_type,
         model_files,
         small_model_file,
         extra_args,
         disable_small_model,
         keep_intermediate_results,
         make_examples_extra_args,
         call_variants_extra_args,
         postprocess_variants_extra_args,
         docker_image,
         region_str,
         region_file,
         ref,
         ref_fai) = row

    extra_args_flag = extra_args.join(" ")

    disable_small_model_flag = disable_small_model ? "--disable_small_model" : "--nodisable_small_model"

    def region_bed = is_real_file(region_file[0]) ? region_file[0].name : ""
    def all_regions="${region_str.join(" ")} ${region_bed}".trim()
    region_flag = all_regions ? "--regions=\"${all_regions}\"" : ""

    def is_custom_model = is_real_file(model_files[0])
    if (is_custom_model) {
      // Set shape flags when using a custom model.
      make_examples_extra_args_flag = set_flag("make_examples_extra_args", "\$(get_shape_flags.sh),${row[10]}".strip(","))
    } else {
      make_examples_extra_args_flag = set_flag("make_examples_extra_args", make_examples_extra_args)
    }

    call_variants_extra_args_flag = set_flag("call_variants_extra_args", call_variants_extra_args)
    postprocess_variants_extra_args_flag = set_flag("postprocess_variants_extra_args", postprocess_variants_extra_args)

    customized_model_flag = is_real_file(model_files[0]) ? "--customized_model=./model.ckpt" : ""
    customized_small_model_flag = is_real_file(small_model_file) ? "--customized_small_model=./model.keras" : ""

    [
      uid, \
      sample, \
      bam, \
      bam_index, \
      model_type, \
      model_files, \
      small_model_file, \
      extra_args_flag, \
      disable_small_model_flag, \
      keep_intermediate_results, \
      make_examples_extra_args_flag, \
      call_variants_extra_args_flag, \
      postprocess_variants_extra_args_flag, \
      docker_image, \
      region_flag, \
      region_file, \
      ref, \
      ref_fai, \
      customized_model_flag, \
      customized_small_model_flag
    ]
}

workflow run_deepvariant_parallel {
  take:
    rows_input
    num_instances

  main:
      // Process the input rows to normalize flags and check input files.
      normalized_input = rows_input.map { normalize_inputs(it) }

      // Add instance index and num_instances to each input row.
      // This parameter will be used to split the input into multiple instances.
      def rows_with_instance_index = normalized_input
          .combine(Channel.of(0..<num_instances))
          .map(v -> tuple(*v, num_instances))

      // Run make_examples and call_variants for each record in rows_input.
      run_make_examples(rows_with_instance_index)
      run_call_variants(
        normalized_input.combine(run_make_examples.out.to_call_variants, by: 0))


      // Group the outputs from call_variants by uid and flatten the list to
      // prepare for the post_processing step.
      to_post_processing_grouped = run_call_variants.out.to_post_processing
      .groupTuple(by: 0)
      .map{uid, \
        num_instances, \
        shards_count, \
        intermediate_results, \
        make_examples_call_variant_outputs, \
        gvcf_files ->
        [ uid, \
          num_instances[0], \
          shards_count[0], \
          intermediate_results.flatten(), \
          make_examples_call_variant_outputs.flatten(), \
          gvcf_files.flatten()\
        ]
      }

      run_post_processing(
        normalized_input
        .combine(to_post_processing_grouped, by: 0))

      emit:
        to_happy = run_post_processing.out.to_happy
        to_aggregate = run_post_processing.out.to_aggregate
}
