# DeepVariant whole exome sequencing (WES) case study

Similar to the [case study on whole genome sequencing data], in this
study we describe applying DeepVariant to a real exome sample using a single
machine.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. We report the runtime with [specific machine type]
for the sake of consistency in reporting run time. This is NOT the fastest or
cheapest configuration. For more scalable execution of DeepVariant see the
[External Solutions] section.

## Running DeepVariant with one command

DeepVariant pipeline consists of 3 steps: `make_examples`, `call_variants`, and
`postprocess_variants`. You can now run DeepVariant with one command using the
`run_deepvariant` script.

### Running on a CPU-only machine

Here is an example command:

```bash
sudo docker run \
  -v "${DATA_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --regions="/input/${CAPTURE_BED}" \
  --output_vcf=/output/HG002.output.vcf.gz \
  --output_gvcf=/output/HG002.output.g.vcf.gz \
  --num_shards=${N_SHARDS}
```

By specifying `--model_type=WES`, you'll be using a model that is best suited
for Illumina Whole Exome Sequencing data.

The script [run_wes_case_study_docker.sh] shows a full example of which data to
download, and run DeepVariant with the data.

Before you run the script, you can read through all sections to understand the
details. Here is a quick way to get the script and run it:

```bash
curl https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wes_case_study_docker.sh | bash
```

### Runtime

See
[this page](deepvariant-details.md#commands-for-requesting-machines-used-in-case-studies)
for the commands used to obtain different machine types on Google Cloud.

Step                               | Hardware | Wall time
---------------------------------- | -------- | ---------
`make_examples`                    | 64 CPUs  | ~ 19m
`call_variants`                    | 64 CPUs  | ~ 2m
`postprocess_variants` (with gVCF) | 1 CPU    | ~ 1m

In this example, `call_variants` does not take much time on 64 CPUs. Running
with GPU might be unnecessary. You can read
[case study on whole genome sequencing data] about the use of GPU. If you want
to use GPU on the exome data in this case study, you can see
[run_wes_case_study_docker_gpu.sh] shows a full example.

## Data description

### BAM file:

`151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam`

Downloaded from
[https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_OsloUniversityHospital_IlluminaExome_bwamem_GRCh37_11252015](https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_OsloUniversityHospital_IlluminaExome_bwamem_GRCh37_11252015)

### FASTA

Same as described in the
[case study for whole genome data](deepvariant-case-study.md#test_data)

### Truth VCF and BED

`HG002_GRCh37_1_22_v4.1_draft_benchmark.*` are from NIST, as part of the
[Genomes in a Bottle project](http://jimb.stanford.edu/giab/). They are
downloaded from
[ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/)

### Capture target BED file

According to the paper "[Extensive sequencing of seven human genomes to
characterize benchmark reference
materials](https://www.nature.com/articles/sdata201625)", the HG002 exome was
generated with Agilent SureSelect. In this case study we'll use the SureSelect
v5 BED (`agilent_sureselect_human_all_exon_v5_b37_targets.bed`) and intersect it
with the GIAB confident regions for evaluation.

## Variant call quality

We used the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This serves as a check
to ensure the three DeepVariant commands ran correctly and produced high-quality
results.

We evaluate against the capture region:

Type  | # TP  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ----- | ---- | ---- | -------- | --------- | ---------
INDEL | 2894  | 126  |  58  | 0.958278 | 0.980602  | 0.969312
SNP   | 38173 | 403  | 130  | 0.989553 | 0.996610  | 0.993069

[specific machine type]: deepvariant-details.md#commands-for-requesting-machines-used-in-case-studies
[install_nvidia_docker.sh]: ../scripts/install_nvidia_docker.sh
[run_wes_case_study_docker.sh]: ../scripts/run_wes_case_study_docker.sh
[run_wes_case_study_docker_gpu.sh]: ../scripts/run_wes_case_study_docker_gpu.sh
[External Solutions]: https://github.com/google/deepvariant#external-solutions
[case study on whole genome sequencing data]: deepvariant-case-study.md
