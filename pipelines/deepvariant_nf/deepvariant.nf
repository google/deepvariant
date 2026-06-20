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
    Run DeepVariant
    ===============
*/
import groovy.yaml.YamlSlurper

// Global Params --> These are not added to sample rows.
params.output_dir=""
params.sample_sheet=""
params.happy_bench_inventory="./case_studies/happy_bench_inventory.yaml"
params.n_trials=1
params.limit=-1
// Do not manually specify params.overrides.
// This is used to indicate non-default values passed in via command line.
params.override_vars=""
params.label=""

// Params
params.uid="" // A unique id for the set of parameters for the given sample.
params.uid_filter="" // If specified, only run the specified uids.
params.sample=""
params.ref=""
params.model="" // Path to model.
params.small_model="" // Path to small model.
params.bam=""
params.regions=""
params.capture_bed=""
params.truth_vcf=""
params.truth_bed=""
params.docker_image=""
params.keep_intermediate_results=false // Passed to run_deepvariant.
params.ignore_model_example_info=false // Use to disable checking for model.example_info.json.
params.run_happy = true
params.run_happy_bench = false
params.num_instances = 1

// run_deepvariant params; For convenience:
params.disable_small_model=""

params.extra_args=""

// Passing Extra Params DeepVariant params.
params.make_examples_extra_args = ""
params.call_variants_extra_args = ""
params.postprocess_variants_extra_args = ""
params.model_type=""
params.haploid_contigs=""
params.par_regions_bed=""
params.shared_memory_size_gb=12

// Nested make_examples flags.
// When set in the YAML under `make_examples_args:`, these are automatically
// merged into make_examples_extra_args. Any valid make_examples flag can be
// used without modifying this file. Example YAML usage:
//   make_examples_args:
//     min_mapping_quality: 0
//     normalize_reads: true
params.make_examples_args=[:]

// Module imports
include { resolve_model_files;
          resolve_small_model_files;
          resolve_bam_index;
          resolve_vcf_index;
          parse_regions;
          save_params;
          create_placeholder_file;
          parse_sample_sheet;
          field_is_set;
          get_field;
          parse_happy_bench_inventory } from '../_modules/util/util.module.nf'
include { run_deepvariant;
          run_pangenome_aware_deepvariant } from './modules/deepvariant.module.nf'
include { run_deepvariant_parallel } from './modules/deepvariant_parallel.module.nf'
include { run_happy; } from '../_modules/happy/happy.module.nf'
include { aggregate; } from '../_modules/aggregate.module.nf'
include { multiqc; } from '../_modules/multiqc.module.nf'

// Constants
DELIM = "\t"

workflow {

  required_fields = ['uid',
                     'sample',
                     'bam',
                     'ref',
                     'docker_image']
  allowed_fields = ['pangenome_graph',
                    'regions',
                    'run_happy',
                    'model_type',
                    'small_model',
                    'extra_args',
                    'disable_small_model',
                    'keep_intermediate_results',
                    'make_examples_args',
                    'make_examples_extra_args',
                    'call_variants_extra_args',
                    'postprocess_variants_extra_args',
                    'haploid_contigs',
                    'par_regions_bed',
                    'capture_bed',
                    'truth_vcf',
                    'truth_bed',
                    'ignore_model_example_info'] + required_fields + params.keySet()
  sample_sheet = parse_sample_sheet(sample_sheet=params.sample_sheet,
                                    output_dir=params.output_dir,
                                    uid_filter=params.uid_filter,
                                    n_trials=params.n_trials,
                                    limit=params.limit,
                                    required_fields=required_fields,
                                    allowed_fields=allowed_fields).map {
                 row ->
                 // Combine regions with capture_bed if specified.
                 row.regions = ((row.regions ?: "") + " " + (row.capture_bed ?: "")).trim()
                 // Resolve bam index.
                 row.bam_index = row.bam_index ?: resolve_bam_index(row.bam);
                 row
                 }

  regions_in = sample_sheet.map { [it.uid] + parse_regions(it.regions) }
  reference_in = sample_sheet.map { [it.uid, it.ref, "${it.ref}.fai"] }

  sample_sheet.map { row ->
    row.model = resolve_model_files(row.model, row.ignore_model_example_info)
    row.small_model = resolve_small_model_files(row.small_model)
    // Format extra_args using `--` prefix for flags if not present.
    def extra_args = []
    if (row.extra_args) {
      extra_args = row.extra_args.split(",").collect { arg -> if (!arg.startsWith("--")) "--${arg}" }
    }

    // Merge nested make_examples_args YAML fields into make_examples_extra_args.
    def me_args = row.make_examples_args
    if (me_args instanceof Map && me_args) {
      def me_flags = me_args.collect { k, v -> "${k}=${v}" }.join(",")
      def existing = row.make_examples_extra_args ?: ""
      row.make_examples_extra_args = [existing, me_flags].findAll { it }.join(",")
    }

    [
      row.uid,
      row.sample,
      row.bam,
      row.bam_index,
      row.model_type,
      row.model,
      row.small_model,
      extra_args,
      row.disable_small_model,
      row.keep_intermediate_results,
      get_field(row, 'make_examples_extra_args'),
      get_field(row, 'call_variants_extra_args'),
      get_field(row, 'postprocess_variants_extra_args'),
      row.docker_image,
      *parse_regions(row.regions),
      row.ref,
      "${row.ref}.fai",
      row.haploid_contigs,
      get_field(row, 'par_regions_bed', "", true),
      get_field(row, 'shared_memory_size_gb', 12),
      row.pangenome_graph
    ]
  }.branch { row ->
    to_deepvariant: !row[-1]
      // Remove the last two elements (shared_memory_size_gb, pangenome_graph).
      return row[0..row.size() - 3]
    to_pangenome: row[-1]
      return row
  }.set { rows_input }

  if (params.num_instances > 1) {
    print("Running in parallel on ${params.num_instances} instances.")
    run_deepvariant_parallel(rows_input.to_deepvariant, params.num_instances)
    dv_out = run_deepvariant_parallel.out
  } else {
    run_deepvariant(rows_input.to_deepvariant)
    dv_out = run_deepvariant.out
  }

  run_pangenome_aware_deepvariant(rows_input.to_pangenome)

  if (params.run_happy_bench) {
    /* This option can run multiple evaluations for a single sample. */
    happy_bench_inventory = parse_happy_bench_inventory(params.happy_bench_inventory)
    truth_in = sample_sheet.map { [it.sample, it.uid] }
                .combine(happy_bench_inventory, by: 0)
                .map { it.swap(0, 1) /* New order [uid, sample, ...] */ }
  } else {
    truth_in = sample_sheet.filter { row ->
      // Check if the truth VCF is not specified or if we have disabled running happy.
      row.truth_vcf && row.run_happy
      }.map { row ->
      [
        row.uid,
        row.sample,
        null, // dataset_name
        row.truth_vcf,
        resolve_vcf_index(row.truth_vcf),
        row.truth_bed
      ]
    }
  }

  happy_in = dv_out.to_happy \
    .concat(run_pangenome_aware_deepvariant.out.to_happy) \
    .combine(truth_in, by: [0, 1] /* Join on uid, sample */)
    .combine(reference_in, by: 0 /* Join on uid */)
    .combine(regions_in, by: 0)

  run_happy(happy_in)
  multiqc(run_happy.out.summary.collect())

  aggregate_in = dv_out.to_aggregate \
                  .concat(run_pangenome_aware_deepvariant.out.to_aggregate) \
                  .concat(run_happy.out.summary) \
                  .collect()
  aggregate(aggregate_in)

}

workflow.onComplete {
  if (workflow.success) {
    if (file(params.output_dir + "/happy.summary.md").exists()) {
      log.info "\nHappy Results: ${params.output_dir}/happy.summary.tsv\n"
      print(file(params.output_dir + "/happy.summary.md").text)
    }

    if (file(params.output_dir + "/happy.error_counts.md").exists()) {
      log.info "\nHappy Error Counts: ${params.output_dir}/happy.error_counts.tsv\n"
      print(file(params.output_dir + "/happy.error_counts.md").text)
    }

    if (file(params.output_dir + "/runtimes.md").exists()) {
      log.info "\nRuntimes: ${params.output_dir}/runtimes.tsv\n"
      print(file(params.output_dir + "/runtimes.md").text)
    }

    if (file(params.output_dir + "/md5sums.md").exists()) {
      log.info "\nmd5Sums: ${params.output_dir}/md5sums.tsv\n"
      print(file(params.output_dir + "/md5sums.md").text)
    }
  }
}
