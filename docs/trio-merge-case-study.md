# Best practices for multi-sample variant calling with DeepVariant (WES trio demonstration)

## Overview

This document outlines all the steps and considerations for calling and merging
a trio using DeepVariant and [GLnexus](https://github.com/dnanexus-rnd/GLnexus).
These best practices were developed and evaluated as described in the article
published in _Bioinformatics_:
[Accurate, scalable cohort variant calls using DeepVariant and GLnexus](https://doi.org/10.1093/bioinformatics/btaa1081)
(2021).

The process involves 3 major stages: running DeepVariant to create individual
genome call sets, running GLnexus to merge call sets, and analyzing the merged
call set.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. The steps below were done on a machine with this
[example command to start a machine](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform).

The steps in this document can be extended to merge larger cohorts as well.

See this workflow:

![workflow](images/cohort-workflow.png?raw=true "DeepVariant+GLnexus cohort workflow")

A few things to note before we start:

*   If you are looking for ways to run DeepVariant in larger batches, please
    refer to the
    [third party solutions](https://github.com/google/deepvariant#external-solutions)
    section.
*   It is recommended to use BAM files with original quality scores. In the case
    that BAM files went through recalibration, optional DV flags can be used in
    order to use original scores: `--parse_sam_aux_fields`,
    `--use_original_quality_scores`.
*   DeepVariant optionally allows gVCF output. This option is required for
    further GLnexus analysis in this document.

## Dataset

The Whole Exome Sequencing (WES) dataset we're using is from:

[ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/)

*   HG002_NA24385_son
*   HG003_NA24149_father
*   HG004_NA24143_mother

### Commands for downloading the input BAMs

Just for convenience, we use aria2 to download our data. You can change it to
whatever other tools (wget, curl) that you prefer.

To install aria2, you can run: `sudo apt-get -y install aria2`

```
DIR="${PWD}/trio"
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam -o HG002.bam
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai -o HG002.bai
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam -o HG003.bam
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bai -o HG003.bai
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam -o HG004.bam
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bai -o HG004.bai
```

### Command for downloading the reference file

```
aria2c -c -x10 -s10 -d "${DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.gz
gunzip ${DIR}/hs37d5.fa.gz
aria2c -c -x10 -s10 -d "${DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/hs37d5.fa.fai
```

### Command for downloading the input capture region BED file

```
aria2c -c -x10 -s10 -d "${DIR}" https://storage.googleapis.com/deepvariant/exome-case-study-testdata/agilent_sureselect_human_all_exon_v5_b37_targets.bed
```

### Command for downloading the truth files


HG002:

```
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz -o HG002_truth.vcf.gz
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi -o HG002_truth.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed -o HG002_truth.bed
```

HG003:

```
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz -o HG003_truth.vcf.gz
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi -o HG003_truth.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed -o HG003_truth.bed
```

HG004:

```
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz -o HG004_truth.vcf.gz
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi -o HG004_truth.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed -o HG004_truth.bed
```

(No need to install bcftools and other tools, because they are now installed in
the DeepVariant images.)

## Run DeepVariant on trio to get 3 single sample VCFs

First, install docker if you don't have it yet: `sudo apt-get -y install
docker.io`

With the example command below, it runs DeepVariant on the trio one by one. This
is for demonstration only. If you're running this on a large cohort, running
serially is not the most effective approach.

```
N_SHARDS=$(nproc)  # Or change to the number of cores you want to use
CAPTURE_BED=agilent_sureselect_human_all_exon_v5_b37_targets.bed
VERSION=1.3.0

declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  BAM=${SAMPLE}.bam

  OUTPUT_VCF=${SAMPLE}.vcf.gz
  OUTPUT_GVCF=${SAMPLE}.g.vcf.gz

  time sudo docker run \
    -v "${DIR}":"/data" \
    google/deepvariant:${VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref="/data/hs37d5.fa" \
    --reads="/data/${BAM}" \
    --regions="/data/${CAPTURE_BED}" \
    --output_vcf="/data/${OUTPUT_VCF}" \
    --output_gvcf="/data/${OUTPUT_GVCF}" \
    --num_shards=${N_SHARDS}
done
```

Note: The BAM files should provide unique names for each sample in their `SM`
header tag, which is usually derived from a command-line flag to the read
aligner. If your BAM files don't have unique `SM` tags (and if it's not feasible
to adjust the alignment pipeline), add the `--sample_name=XYZ` flag to
`run_deepvariant` to override the sample name written into the gVCF file header.

## Merge the trio samples using GLnexus

### Run GLnexus to merge 3 gVCFs

And then run GLnexus with this config:

```
sudo docker pull quay.io/mlin/glnexus:v1.2.7

time sudo docker run \
  -v "${DIR}":"/data" \
  quay.io/mlin/glnexus:v1.2.7 \
  /usr/local/bin/glnexus_cli \
  --config DeepVariantWES \
  --bed "/data/${CAPTURE_BED}" \
  /data/HG004.g.vcf.gz /data/HG003.g.vcf.gz /data/HG002.g.vcf.gz \
  | sudo docker run -i google/deepvariant:${VERSION} bcftools view - \
  | sudo docker run -i google/deepvariant:${VERSION} bgzip -c \
  > ${DIR}/deepvariant.cohort.vcf.gz
```

When we ran on this WES trio, it took only about 13 seconds. For more details on
performance, see
[GLnexus performance guide](https://github.com/dnanexus-rnd/GLnexus/wiki/Performance).

For a WGS cohort, we recommend using `--config DeepVariantWGS` instead of
`DeepVariantWES`. Another preset `DeepVariant_unfiltered` is available in
`glnexus:v1.2.7` or later versions for merging DeepVariant gVCFs with no QC
filters or genotype revision (see
[GitHub issue #326](https://github.com/google/deepvariant/issues/326) for a
potential use case). The details of these presets can be found
[here](../deepvariant/cohort_best_practice).

## Annotate the merged VCF with Mendelian discordance information using RTG Tools

Create an SDF template from our reference file:

```
sudo docker run \
  -v "${DIR}":"/data" \
  realtimegenomics/rtg-tools format \
  -o /data/hs37d5.sdf /data/hs37d5.fa
```

Create a PED file `$DIR/trio.ped` that looks like this (with the sample name
of the trio):

```
#PED format pedigree
#
#fam-id/ind-id/pat-id/mat-id: 0=unknown
#sex: 1=male; 2=female; 0=unknown
#phenotype: -9=missing, 0=missing; 1=unaffected; 2=affected
#
#fam-id ind-id pat-id mat-id sex phen
1 Sample_Diag-excap51-HG002-EEogPU Sample_Diag-excap51-HG003-EEogPU Sample_Diag-excap51-HG004-EEogPU 1 0
1 Sample_Diag-excap51-HG003-EEogPU 0 0 1 0
1 Sample_Diag-excap51-HG004-EEogPU 0 0 2 0
```

## Annotate merged VCF with RTG Tools

```
sudo docker run \
  -v "${DIR}":"/data" \
  realtimegenomics/rtg-tools mendelian \
  -i /data/deepvariant.cohort.vcf.gz \
  -o /data/deepvariant.annotated.vcf.gz \
  --pedigree=/data/trio.ped \
  -t /data/hs37d5.sdf \
  | tee ${DIR}/deepvariant.input_rtg_output.txt
```

The output is:

```
Checking: /data/deepvariant.cohort.vcf.gz
Family: [Sample_Diag-excap51-HG003-EEogPU + Sample_Diag-excap51-HG004-EEogPU] -> [Sample_Diag-excap51-HG002-EEogPU]
Concordance Sample_Diag-excap51-HG002-EEogPU: F:56629/57124 (99.13%)  M:57014/57138 (99.78%)  F+M:56384/57014 (98.90%)
Sample Sample_Diag-excap51-HG002-EEogPU has less than 99.0 concordance with both parents. Check for incorrect pedigree or sample mislabelling.
815/57369 (1.42%) records did not conform to expected call ploidy
57271/57369 (99.83%) records were variant in at least 1 family member and checked for Mendelian constraints
210/57271 (0.37%) records had indeterminate consistency status due to incomplete calls
636/57271 (1.11%) records contained a violation of Mendelian constraints
```

From this report, we know that there is a 1.11% Mendelian violation rate, and
0.37% of the records had incomplete calls (with `.`) so RTG couldn't determine
whether there is violation or not.

## Single sample quality metrics

In addition to the cohort quality statistics, for completeness we generate
single-sample quality metrics.

### ti/tv ratio

We run `bcftools stats` on the 3 VCF outputs. Since our DeepVariant run already
constrained to just the capture regions, no need to specify it again here. We
had to pass in the `-f PASS` flag so that only the PASS calls are considered.

```
declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  sudo docker run \
  -v ${DIR}:${DIR} \
  google/deepvariant:${VERSION} \
  bcftools stats -f PASS \
    ${DIR}/${SAMPLE}.vcf.gz \
  > ${DIR}/${SAMPLE}.stats
done
```

| Sample | [3]ts | [4]tv | [5]ts/tv | [6]ts (1st ALT) | [7]tv (1st ALT) | [8]ts/tv (1st ALT) |
| ------ | ----- | ----- | -------- | --------------- | --------------- | ------------------ |
| HG002  | 29819 | 11626 | 2.56     | 29809           | 11606           | 2.57               |
| HG003  | 29686 | 11670 | 2.54     | 29677           | 11649           | 2.55               |
| HG004  | 29920 | 11801 | 2.54     | 29906           | 11783           | 2.54               |

If you want to restrict to the truth BED files, use this command:

```
declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  sudo docker run \
  -v ${DIR}:${DIR} \
  google/deepvariant:${VERSION} \
  bcftools stats -f PASS \
    -T ${DIR}/${SAMPLE}_truth.bed \
    ${DIR}/${SAMPLE}.vcf.gz \
  > ${DIR}/${SAMPLE}.with_truth_bed.stats
done
```

Which resulted in this table:

| Sample | [3]ts | [4]tv | [5]ts/tv | [6]ts (1st ALT) | [7]tv (1st ALT) | [8]ts/tv (1st ALT) |
| ------ | ----- | ----- | -------- | --------------- | --------------- | ------------------ |
| HG002  | 27688 | 10533 | 2.63     | 27680           | 10520           | 2.63               |
| HG003  | 27324 | 10502 | 2.60     | 27319           | 10491           | 2.60               |
| HG004  | 27474 | 10592 | 2.59     | 27467           | 10580           | 2.60               |


### Rtg vcfstats

```
declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  sudo docker run \
  -v "${DIR}":"/data" \
  realtimegenomics/rtg-tools vcfstats \
  /data/${SAMPLE}.vcf.gz \
  > ${DIR}/${SAMPLE}.vcfstats
done
```

which shows the following:

HG002:

```
Location                     : /data/HG002.vcf.gz
Failed Filters               : 14776
Passed Filters               : 45080
SNPs                         : 41414
MNPs                         : 0
Insertions                   : 1866
Deletions                    : 1777
Indels                       : 19
Same as reference            : 4
SNP Transitions/Transversions: 2.56 (41696/16271)
Total Het/Hom ratio          : 1.49 (26956/18120)
SNP Het/Hom ratio            : 1.51 (24883/16531)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.11 (981/885)
Deletion Het/Hom ratio       : 1.52 (1073/704)
Indel Het/Hom ratio          : - (19/0)
Insertion/Deletion ratio     : 1.05 (1866/1777)
Indel/SNP+MNP ratio          : 0.09 (3662/41414)
```

HG003:

```
Location                     : /data/HG003.vcf.gz
Failed Filters               : 15654
Passed Filters               : 44919
SNPs                         : 41321
MNPs                         : 0
Insertions                   : 1850
Deletions                    : 1729
Indels                       : 18
Same as reference            : 1
SNP Transitions/Transversions: 2.52 (41489/16448)
Total Het/Hom ratio          : 1.48 (26775/18143)
SNP Het/Hom ratio            : 1.49 (24728/16593)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.15 (988/862)
Deletion Het/Hom ratio       : 1.51 (1041/688)
Indel Het/Hom ratio          : - (18/0)
Insertion/Deletion ratio     : 1.07 (1850/1729)
Indel/SNP+MNP ratio          : 0.09 (3597/41321)
```

HG004:

```
Location                     : /data/HG004.vcf.gz
Failed Filters               : 15365
Passed Filters               : 45316
SNPs                         : 41686
MNPs                         : 0
Insertions                   : 1860
Deletions                    : 1750
Indels                       : 19
Same as reference            : 1
SNP Transitions/Transversions: 2.55 (41517/16304)
Total Het/Hom ratio          : 1.57 (27656/17659)
SNP Het/Hom ratio            : 1.59 (25573/16113)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.13 (986/874)
Deletion Het/Hom ratio       : 1.60 (1078/672)
Indel Het/Hom ratio          : - (19/0)
Insertion/Deletion ratio     : 1.06 (1860/1750)
Indel/SNP+MNP ratio          : 0.09 (3629/41686)
```

### Run hap.py to calculate the accuracy of DeepVariant generated call sets

```
sudo docker pull jmcdani20/hap.py:v0.3.12

declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  sudo docker run -i \
    -v "${DIR}":"/data" \
    jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    "/data/${SAMPLE}_truth.vcf.gz" \
    "/data/${SAMPLE}.vcf.gz" \
    -f "/data/${SAMPLE}_truth.bed" \
    -T "/data/${CAPTURE_BED}" \
    -r "/data/hs37d5.fa" \
    -o "/data/${SAMPLE}.happy.output" \
    --engine=vcfeval \
    --pass-only > ${DIR}/${SAMPLE}.stdout
done
```

Accuracy F1 scores:

Sample | Indel    | SNP
------ | -------- | --------
HG002  | 0.971084 | 0.993661
HG003  | 0.965795 | 0.993012
HG004  | 0.972641 | 0.993467
