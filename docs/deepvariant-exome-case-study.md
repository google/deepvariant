# DeepVariant exome case study

Similar to the [case study on whole genome sequencing
data](deepvariant-case-study.md), in this study we describe applying DeepVariant
to a real exome sample using a single machine.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. This might not be the fastest or cheapest
configuration for your needs. For more scalable execution of DeepVariant see the
[Docker-based exome pipeline](https://cloud.google.com/genomics/deepvariant)
created for Google Cloud Platform.

## Update since r0.7 : using docker instead of copying binaries.

Starting from the 0.7 release, we use docker to run the binaries instead of
copying binaries to local machines first. You can still read about the previous
approach in
[the Exome Case Study in r0.6](https://github.com/google/deepvariant/blob/r0.6/docs/deepvariant-exome-case-study.md).

We recognize that there might be some overhead of using docker run. But using
docker makes this case study easier to generalize to different versions of Linux
systems. For example, we have verified that you can use docker to run
DeepVariant on other Linux systems such as CentOS 7.

## Update since v0.7.1: Run the case study with one script

Script:
https://github.com/google/deepvariant/blob/r0.7/tools/run_wes_case_study.sh

Get the script and run everything. This will install and download everything
needed for this case study. And it will run DeepVariant to generate the output
VCFs, and also run the evaluation as well.

Before you run the script, you can read through all sections to understand the
details.

```bash
wget https://raw.githubusercontent.com/google/deepvariant/r0.7/tools/run_wes_case_study.sh -P ${HOME}
chmod +x ${HOME}/run_wes_case_study.sh
${HOME}/run_wes_case_study.sh
```

## (Optional) If on Google Cloud Platform, request a machine with this example command

Any sufficiently capable machine will do. For this case study, we used a 64-core
non-preemptible instance with 128GiB and no GPU.

If you need an example, see
[this section](deepvariant-case-study.md#optional-if-on-google-cloud-platform-request-a-machine-with-this-example-command).

## Description for data

### BAM file:

`151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam`

Downloaded from
[https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_OsloUniversityHospital_IlluminaExome_bwamem_GRCh37_11252015](https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_OsloUniversityHospital_IlluminaExome_bwamem_GRCh37_11252015)

### FASTA

Same as described in the [case study for whole genome
data](deepvariant-case-study.md#test_data)

### Truth VCF and BED

`HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_*`
are from NIST, as part of the [Genomes in a Bottle
project](http://jimb.stanford.edu/giab/). They are downloaded from
[ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/)

### Capture target BED file

According to the paper "[Extensive sequencing of seven human genomes to
characterize benchmark reference
materials](https://www.nature.com/articles/sdata201625)", the HG002 exome was
generated with Agilent SureSelect. In this case study we'll use the SureSelect
v5 BED (`agilent_sureselect_human_all_exon_v5_b37_targets.bed`) and intersect it
with the GIAB confident regions for evaluation.

## Notes for `call_variants`

More discussion can be found in the
[call_variants section in the case study](deepvariant-case-study.md#notes-for-call-variants).

## Notes for `postprocess_variants`

More discussion can be found in the
[postprocess_variants section in the case study](deepvariant-case-study.md#notes-for-postprocess-variants).

## Resources used by each step

Step                               | wall time
---------------------------------- | ---------
`make_examples`                    | 13m 48s
`call_variants`                    | 2m 7s
`postprocess_variants` (no gVCF)   | 0m 14s
`postprocess_variants` (with gVCF) | 1m 19s
total time (single machine)        | ~17m

## Variant call quality

We used the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This serves as a check
to ensure the three DeepVariant commands ran correctly and produced high-quality
results.

We evaluate against the capture region:

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 111  | 51   | 0.957308 | 0.980086  | 0.968563
SNP   | 48   | 15   | 0.998577 | 0.999555  | 0.999066

## Separate models for calling whole genome and exome data

Starting from DeepVariant 0.5.\* and later releases, we recommend a separate
model for calling exome sequencing data. Here is how the exome model is trained:
we used a WGS model as the starting checkpoint (instead of an ImageNet one), and
trained only on examples created from exome data.
