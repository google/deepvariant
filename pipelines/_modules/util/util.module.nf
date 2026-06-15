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
    Utilities Module
    ================

    This module contains general utilities.
*/
import groovy.json.JsonOutput
import groovy.yaml.YamlSlurper
import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions
import groovy.transform.Memoized

process save_params {

    executor 'local'

    publishDir "${params.output_dir}", mode: 'copy'

    output:
        path('params.json')

    script:
      "echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > params.json"
}

process samtools_index {
  /*
    Args:
      input_bam: Path to input bam file.
      index_type: Type of index to create. Must be 'csi' or 'bai'.
      output_path: Directory to publish index file, relative to the output_dir.
                   The same basename will be used.
  */

  container 'quay.io/biocontainers/samtools:1.15.1--h6899075_1'

  publishDir "${params.output_dir}/${output_path}", mode: 'copy'

  input:
      path(input_bam)
      val(index_type)
      val(output_path)

  output:
      path(output_vcf_idx)

  script:
    if (index_type == 'csi') {
      output_vcf_idx = "${input_bam}.csi"
      index_flag = "-c"
    } else if (index_type == 'bai') {
      output_vcf_idx = "${input_bam}.bai"
      index_flag = "-b"
    } else {
      exit 1, "Unsupported index type: ${index_type}; Must be 'csi' or 'bai'"
    }

  """
  samtools index ${index_flag} -@ ${task.cpus} ${input_bam}
  """

}

process bcftools_index {
  /*
    Args:
      input_bam: Path to input vcf/bcf file.
      index_type: Type of index to create. Must be 'csi' or 'tbi'.
      output_path: Directory to publish index file, relative to the output_dir.
                   The same basename will be used.
  */

  container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

  publishDir "${params.output_dir}/${output_path}", mode: 'copy'

  input:
      path(input_vcf)
      val(index_type)
      val(output_path)

  output:
      path(output_vcf_idx)

  script:
    if (index_type == 'csi') {
      output_vcf_idx = "${input_vcf}.csi"
      index_flag = "--csi"
    } else if (index_type == 'tbi') {
      output_vcf_idx = "${input_vcf}.tbi"
      index_flag = "--tbi"
    } else {
      exit 1, "Unsupported index type: ${index_type}"
    }

  """
  bcftools index ${index_flag} --threads ${task.cpus} ${input_vcf}
  """

}

@Memoized
def resolve_vcf_index(vcf) {
  /*
  Args:
    vcf: Path to bam file.

  Returns:
    An associated index file.
  */
  def idx
  if (file("${vcf}.tbi").exists()) {
      idx = "${vcf}.tbi"
  } else if (file("${vcf}.csi").exists()) {
      idx = "${vcf}.csi"
  } else {
    exit 1, "No index file found for vcf: ${vcf}"
  }
  idx
}



@Memoized
def resolve_bam_index(bam) {
  /*
  Args:
    bam: Path to bam or cram file.

  Returns:
    An associated index file.
  */

  // Return null if called with null input.
  if (bam == null) {
    return null
  }

  def idx

  if (file(bam).getName().startsWith("NO_FILE")) {
    idx = create_placeholder_file("bam_index")
  }
  if (file("${bam}.bai").exists()) {
      idx = "${bam}.bai"
  } else if (bam.endsWith(".bam") && file(bam.replace("bam", "bai")).exists()) {
      idx = bam.replaceAll('\\.bam$', '.bai')
  } else if (file("${bam}.csi").exists()) {
      idx = "${bam}.csi"
  } else if (bam.endsWith(".cram") && file("${bam}.crai").exists()) {
      idx = "${bam}.crai"
  } else {
    if (params.check_for_errors) {
      throw new IllegalArgumentException("No index file found for bam: ${bam}")
    }
  }
  return idx
}

@Memoized
def resolve_fai_index(ref) {
  /*
  Args:
    ref: Path to ref file.

  Returns:
    The fai index file path.
  */
  def fai = file("${ref}.fai")
  if (!fai.exists()) {
    throw new IllegalArgumentException("No .fai index file found for ref: ${ref}")
  }
  return fai.toUriString()
}

@Memoized
def resolve_gzi_index(ref) {
  /*
  Args:
    ref: Path to ref file.

  Returns:
    The gzi index file path if the reference is compressed, otherwise null.
  */
  if (file(ref).getName().endsWith(".gz")) {
    def gzi = file("${ref}.gzi")
    if (!gzi.exists()) {
      throw new IllegalArgumentException("No .gzi index file found for compressed ref: ${ref}")
    }
    return gzi.toUriString()
  }
  return null
}


@Memoized
def resolve_model_files(model_path, ignore_model_example_info_json=false, placeholder_name="model") {
  // Allow for a model to be specified using its checkpoint filename.
  // If revert_to_example_info_json is true, then the example_info.json file will
  // be used instead of the model.example_info.json file.
  //
  // If ignore_model_example_info_json is true, then the model.example_info.json
  // file will not be checked for existence. This is useful for older models.
  if (!model_path) {
    return [create_placeholder_file(placeholder_name)]
  }

  def model_files = [] as List
  def model_file = file(model_path.replaceAll(/\/+$/, ""))
  def parent_dir = ""
  if (model_file.name.endsWith("savedmodel") || model_file.name.endsWith("saved_model.pb")) {
    // Saved Model
    parent_dir = model_file.parent.toUriString()
    model_files.add("${parent_dir}saved_model.pb")
    model_files.add("${parent_dir}variables")
  } else {
    // Checkpoint
    if (model_file.isDirectory()) {
      // If a directory is specified then we expect a single checkpoint
      model_fname = ""
      parent_dir = model_file.toUriString()
    } else if (model_file.name == "checkpoint") {
      // If the file is named 'checkpoint', then we expect a single checkpoint.
      model_fname = ""
      parent_dir = model_file.parent.toUriString()
    } else {
      // If a prefix is provided, we will use that to glob the appropriate
      // checkpoint.
      model_fname = model_file.name
      parent_dir = model_file.parent.toUriString()
    }

    def data_file_list = file("${parent_dir}/${model_fname}*.data-00000-of-00001")
    def index_file_list = file("${parent_dir}/${model_fname}*.index")

    if (params.check_for_errors) {
      if (data_file_list.size() == 0) {
        throw new Exception("Data file not found for model: ${model_file.toUriString()}")
      }
      if (data_file_list.size() > 1) {
        throw new Exception("Multiple data files found for model: ${model_file.toUriString()}. Please specify a more unique model_path.")
      }
      if (index_file_list.size() == 0) {
        throw new Exception("Index file not found for model: ${model_file.toUriString()}")
      }
      if (index_file_list.size() > 1) {
        throw new Exception("Multiple index files found for model: ${model_file.toUriString()}. Please specify a more unique model_path.")
      }
    }

    def data_file = data_file_list.size() > 0 ? data_file_list[0] : null
    def index_file = index_file_list.size() > 0 ? index_file_list[0] : null

    model_files.add(data_file.toUriString())
    model_files.add(index_file.toUriString())
  }

  // Optional Files:
  def meta_file = file("${parent_dir}/*.meta")
  meta_file = meta_file.size() > 0 ? meta_file[0] : null
  def example_info_json = "${parent_dir}/example_info.json"
  def model_example_info_json = file("${parent_dir}/model.example_info.json")
  if (params.check_for_errors) {
    if (!model_example_info_json.exists() && !ignore_model_example_info_json) {
      throw new Exception("Model example info json file not found: ${model_example_info_json.toUriString()}")
    }
  }

  def optional_files = [model_example_info_json, example_info_json, meta_file]
  for (def optional_file : optional_files) {
    if (optional_file == null) continue
    def current_file = file(optional_file)
    if (current_file.exists()) {
      model_files.add(current_file.toUriString())
    }
  }
  return model_files as List
}

@Memoized
def resolve_small_model_files(small_model_path, placeholder_name="small_model") {
  /*
  Allow a custom small model to be specified.
  */
  if (!small_model_path) {
    return create_placeholder_file(placeholder_name)
  }
  if (params.check_for_errors) {
    if (!file(small_model_path + "/model.keras").exists() && !file(small_model_path + "/fingerprint.pb").exists()) {
      throw new Exception("Small model files not found. Checked for newer setup: model.keras and older setup: fingerprint.pb at '${small_model_path}'")
    }
  }
  return file(small_model_path)
}

@Memoized
def create_placeholder_file(name) {
  /*
  Placeholder files are used to enable optional file inputs.
  Because the location needs to be constant, placeholder files are created
  in the `workDir`.
  */
  def placeholder_path = "${workDir.toUriString()}/placeholder_files/NO_FILE_${name}"
  def placeholder_file = file(placeholder_path)
  if (!placeholder_file.exists()) {
    println("Creating placeholder file: ${placeholder_path}")
    file(placeholder_path).getParent().mkdirs()
    placeholder_file.withWriter { writer ->
        writer.write("")
    }
  }
  placeholder_file
}

@Memoized
def parse_regions(String inputString) {
  /*
  This utility function parses regions into strings and files and returns
  a tuple of [region_str, region_file].

  Regions can be specified as:
  - A single string (e.g. chr1)
  - A space delimited list of chroms/regions (e.g. "chr1 chr2:1-100,000")
  - A path to a bed Files (e.g. "gs://path/to/bed")
  - A mixture of regions and bed files (e.g. "chr1 chr2:chr2:1-100,000 gs://path/to/bed")

  * For convenience, regions can contain commas.
  * Only one bed file can be specified.
  */

  // 1. Split the string by spaces to get all individual tokens.
  // .trim() removes leading/trailing whitespace from each token.
  // .findAll { it } filters out any empty strings that might result from multiple spaces.
  inputString = inputString ?: ""
  def allTokens = inputString.split(' ').collect { it.trim() }.findAll { it }

  // 2. Initialize the two sets
  def regionStr = [] as Set // Will contain both chr:start-end and chrX
  def regionFiles = [] as Set

  // 3. Process and classify each token
  allTokens.each { token ->
      // Check for file paths first
      if (token.startsWith('gs://') || token.endsWith('.bed') || token.endsWith('.bed.gz')) {
          regionFiles.add(file(token))
      } else {
          // Remove commas from region strings.
          def sanitizedRegion = token.replaceAll(',', '')
          regionStr.add(sanitizedRegion)
      }
  }

  if (regionFiles.size() > 1) {
    throw new IllegalArgumentException("Only one region bed file is supported.")
  }

  if (regionFiles.size() == 0) {
    regionFiles.add(create_placeholder_file("region"))
  }

  return [regionStr, regionFiles]
}


def field_is_set(item, field) {
  /*
  Checks whether a field exists in a row and is not empty.
  */
  item.containsKey(field) && item[field] != "" && item[field] != null
}

def get_field(row, field, value="", placeholder=false) {
  /*
  Returns value if a field is set.
  Args:
    row - The groovy Map (dictionary) to check.
    field - The field to check.
    value - The default value to return if the field is not set.
    placeholder - If true, then a placeholder file will be returned based on the field name.
  Returns:
    The value of the field if it is set, otherwise the value is returned.
  */
  if (field_is_set(row, field)) {
    return row[field]
  } else {
    if (placeholder) {
      return create_placeholder_file(field)
    }
    return value
  }
}


def parse_sample_sheet(sample_sheet,
                       output_dir,
                       uid_filter = "",
                       n_trials = 1,
                       limit = -1,
                       required_fields = [],
                       allowed_fields = [],
                       delim= "\t" ) {
  /*
  This function is used to:
  - intake sample sheets,
  - override values passed in as params,
  - validate required fields are present.
  - (optionally) validate whether a specified field is allowed.
  - (optionally) filter by uids.
  - (optionally) duplicate sample sheet rows for n_trials.
  - (optionally) limit the number of samples to process.

  This function is long, but it standardizes the way sample sheets are
  specified and handled across pipelines. The intention after parsing a sample
  sheet is to use a `.map` function to pick out fields and construct inputs
  for downstream processes.

  Sample sheets can be provided as a YAML or TSV.
  This function will read in the sample sheet, and override it with any
  manually passed arguments (params.<x>). Then it will check that all
  `required_fields` are present.

  Every row of a sample sheet *must* contain a uid column, and that column
  must be unique.

  Optionally, this function can also create n_trials of each sample by setting
  `params.n_trials`.

  Args:
    sample_sheet - The path to the sample sheet (TSV or YAML).
    output_dir - The directory to write the sample sheet to.
    uid_filter - A comma-separated list of UIDs to filter.
    n_trials - An integer number of trials to duplicate for each sample.
    required_fields - An integer number of samples to limit the sheet to.
    delim - The delimiter to use for the sample sheet (default: "\t").

  Returns:
    A list of sample sheet rows
    (A list of groovy Maps; Similar to a list of python dictionaries).
  */

  def GLOBAL_FIELDS = [
    'uid_filter',
    'run_name',
    'sample_sheet',
    'output_dir',
    'n_trials',
    'debug',
    'limit',
    'override_vars',
  ]

  // Read in YAML or TSV sample sheet.
  def sample_sheet_list
  if (sample_sheet == "") {
    throw new IllegalArgumentException("Must specify a sample sheet --sample_sheet.")
  } else if (sample_sheet.endsWith(".yaml")) {
    print("Reading sample sheet as YAML")
    sample_sheet_in = new YamlSlurper().parse(file(sample_sheet));
    if (sample_sheet_in instanceof Map) {
      sample_sheet_in = [sample_sheet_in]
    }
    sample_sheet_list = Channel.fromList(sample_sheet_in)
  } else {
    print("Reading sample sheet as TSV")
    sample_sheet_list = Channel.fromPath(sample_sheet,
                                          checkIfExists: true)
                                .splitCsv(sep: delim,
                                          header: true)
  }

  // Use a shared set to track which fields have been reported as overridden.
  def reportedOverrides = new HashSet()

  // Read in Sample Sheet.
  sample_sheet_0 = sample_sheet_list
                        .filter { if (uid_filter.toString() != "") {
                                      uid_filter.toString().split(",").contains(it.uid.toString())
                                  } else { true }
                                }
                        .map { row ->
                          /*
                          Override params from the sample sheet if specified.

                          Values are specified in the following order:
                          - 1. If a field is specified on the command-line, it will be used.
                          - 2. Use value if specified in sample sheet.
                          - 3. Finally, if a default value is specified as param.key=<default> in `.nf` file, use that.

                          Nextflow allows params to be specified using
                          params.value=<default>, and params can be overridden
                          by passing in `--my_value=new_value` to the command-line.
                          However, nextflow does not tell you whether a default
                          value has been changed, which makes it difficult to
                          override defaults. Therefore, we pass in an extra
                          field to indicate which fields are explicitly set
                          so they will be overridden when not set in a sample
                          sheet.

                          If additional global parameters are defined outside
                          the standard set, it is possible they will be added to
                          the sample sheet, but they should not interfere with
                          thedownstream processes since process inputs are
                          defined by picking out specific fields.
                          */
                          def override_vars = params.override_vars.split(",") ?: []
                          for (param in params) {
                            field = param.key
                            if (GLOBAL_FIELDS.contains(field)) {
                              // Global fields are not added to the sample sheet.
                              continue
                            }
                            if (override_vars.contains(field)) {
                              // This is an override. Log it only once per field.
                              if (!reportedOverrides.contains(field)) {
                                println("\033[1mOverriding field params.${field}=${row[field]} => ${param.value}\033[0m")
                                reportedOverrides.add(field)
                              }
                              row[field] = param.value
                            } else {
                              if (!row.containsKey(field) && param.value) {
                                // Set the default value if it is not specified.
                                row[field] = param.value
                              }
                            }
                          }
                          row
                        }

  // Check for errors in the sample sheet.
  if (params.check_for_errors) {
    // Check if the sample sheet is empty.
    sample_sheet_0.ifEmpty {
      throw new IllegalArgumentException("No rows found. Check your sample sheet and filters.")
    }

    // Validate sample sheet fields by checking required / allowed fields.
    sample_sheet_0.map { row ->
      if (required_fields.size() > 0) {
        for (field in required_fields) {
          if (!row.containsKey(field)) {
            throw new IllegalArgumentException("Field '${field}' is required but is not set.")
          }
        }
      }
      if (allowed_fields.size() > 0) {
        for (field in row.keySet()) {
          if (!allowed_fields.contains(field)) {
            throw new IllegalArgumentException("Field: '${field}' is set but is not an allowed field.")
          }
        }
      }
      if (row.uid.toString().contains("-trial")) {
        throw new IllegalArgumentException("uid column cannot contain the substring '-trial'.")
      }
    }
  }

  // Duplicate sample sheet rows for n trials.
  sample_sheet = sample_sheet_0.take(limit).map { row ->
      (1..n_trials).toList().collect {
        def rc = row.clone();
        if (it > 1) {
          // Only append after first trial; This way we can leverage caching
          // if we choose to run trials *after* an initial run.
          rc.uid = "${rc.uid}-trial${it}";
        };
        rc }
  }.flatten()

  // Write the sample in YAML format to the output directory (with overrides).
  //def builder = new YamlBuilder()
  println("Writing sample sheet to ${params.output_dir}/sample_sheet.yaml")
  def options = new DumperOptions()
  options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK)
  options.setPrettyFlow(true)
  def yaml = new Yaml(options)
  sample_sheet.collect().toList().map { yaml.dumpAll(it.clone().iterator()) }
              .collectFile(name: "sample_sheet.yaml",
                           storeDir: output_dir,
                           newLine: true,
                           keepHeader: false,
                           sort: 'index')

  if (params.check_for_errors) {
    // Check UIDs are unique, and uid characters are valid.
    // uids can contain `-`, but not `_`.
    uids_ch = sample_sheet.map {
      if (it.uid.toString().contains('_')) {
        throw new IllegalArgumentException("uid value ${it.uid} contains an underscore. The uid column in your sample sheet must not contain underscores. You can use dashes instead.")
      }
      it.uid
    }
    uids_ch.unique().count().join(uids_ch.count()).ifEmpty {
      throw new IllegalArgumentException("The uid column in your sample sheet must be unique.")
    }
  }

  return sample_sheet

}


def parse_happy_bench_inventory(happy_bench_inventory_path) {
  /*
  This function is used to parse the happy_bench_inventory file.
  */
  if (happy_bench_inventory_path == "") {
    throw new IllegalArgumentException("Must specify a happy_bench_inventory file.")
  }
  happy_bench_inventory_in = new YamlSlurper().parse(file(happy_bench_inventory_path));
  // Nest sample sheet in a list if its a single-record file.
  if (happy_bench_inventory_in instanceof Map) {
    happy_bench_inventory_in = [happy_bench_inventory_in]
  }
  happy_bench_inventory_list = Channel.fromList(happy_bench_inventory_in)
                                      .map { row ->
                                            [
                                              row.sample,
                                              row.dataset_name,
                                              row.truth_vcf,
                                              resolve_vcf_index(row.truth_vcf),
                                              row.truth_bed
                                            ]
                                          }
}


def set_flag(flag, value) {
  if (value) {
    return "--${flag}=\"${value}\""
  } else {
    return ""
  }
}

def get_make_examples_extra_args_flag(model_files, make_examples_extra_args) {
  is_custom_model = !model_files[0].getName().startsWith("NO_FILE")
  current_make_examples_extra_args = make_examples_extra_args
  if (is_custom_model) {
    // Set shape flags when using a custom model.
    has_model_example_info = model_files.any { it.getName() == 'model.example_info.json' }
    has_example_info = model_files.any { it.getName() == 'example_info.json' }
    if (!has_model_example_info && has_example_info) {
      // If no model.example_info.json, and an example_info.json file exists,
      // then we can extract the shape flags from it.
      // This is done for backward compatibility with older models.
      current_make_examples_extra_args = "\$(get_shape_flags.sh),${make_examples_extra_args}".strip(",")
    }
  }
  return set_flag("make_examples_extra_args", current_make_examples_extra_args)
}

def is_real_file(input) {
    !input?.getName()?.startsWith("NO_FILE")
}
