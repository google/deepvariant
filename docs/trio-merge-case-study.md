# Best practices for multi-sample variant calling with DeepVariant (WES trio demonstration)

## Overview

This document outlines all the steps and considerations for calling and merging
a trio using DeepVariant and [GLnexus](https://github.com/dnanexus-rnd/GLnexus).
These best practices were developed and evaluated as described in the bioRxiv
preprint
[Accurate, scalable cohort variant calls using DeepVariant and GLnexus](https://doi.org/10.1101/2020.02.10.942086).

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

There have been newer version of the truth files, including
[v4.1, GRCh37 for HG002](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37),
and
[v4.2, GRCh38 for HG002-4](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/).
In the future we will plan to update this documentation with newer versions.

HG002:

```
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz -o HG002_truth.vcf.gz
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi -o HG002_truth.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed -o HG002_truth.bed
```

HG003:

```
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz -o HG003_truth.vcf.gz
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi -o HG003_truth.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed -o HG003_truth.bed
```

HG004:

```
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh37/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz -o HG004_truth.vcf.gz
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh37/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi -o HG004_truth.vcf.gz.tbi
aria2c -c -x10 -s10 -d "${DIR}" ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh37/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed -o HG004_truth.bed
```

## Install bcftools, samtools, htslib, if you don't already have them

Here are example commands used to install these tools. You should also look at
http://www.htslib.org for official instructions.

```
sudo apt-get install -y build-essential libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev tabix

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -xvf samtools-1.9.tar.bz2
pushd samtools-1.9 && ./configure && make && sudo make install && popd

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -xvf bcftools-1.9.tar.bz2
pushd bcftools-1.9 && ./configure && make && sudo make install && popd

wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -xvf htslib-1.9.tar.bz2
pushd htslib-1.9 && ./configure && make && sudo make install && popd
```

## Run DeepVariant on trio to get 3 single sample VCFs

First, install docker if you don't have it yet: `sudo apt-get -y install
docker.io`

With the example command below, it runs DeepVariant on the trio one by one. This
is for demonstration only. If you're running this on a large cohort, running
serially is not the most effective approach.

```
N_SHARDS=$(nproc)  # Or change to the number of cores you want to use
CAPTURE_BED=agilent_sureselect_human_all_exon_v5_b37_targets.bed
VERSION=1.1.0

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
  | bcftools view - | bgzip -c > ${DIR}/deepvariant.cohort.vcf.gz
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
1 non-pass records were skipped
Concordance Sample_Diag-excap51-HG002-EEogPU: F:59383/59949 (99.06%)  M:59760/59944 (99.69%)  F+M:59053/59806 (98.74%)
Sample Sample_Diag-excap51-HG002-EEogPU has less than 99.0 concordance with both parents. Check for incorrect pedigree or sample mislabelling.
852/60261 (1.41%) records did not conform to expected call ploidy
60107/60261 (99.74%) records were variant in at least 1 family member and checked for Mendelian constraints
242/60107 (0.40%) records had indeterminate consistency status due to incomplete calls
781/60107 (1.30%) records contained a violation of Mendelian constraints
```

From this report, we know that there is a 1.30% Mendelian violation rate, and
0.40% of the records had incomplete calls (with `.`) so RTG couldn't determine
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
  bcftools stats -f PASS \
    ${DIR}/${SAMPLE}.vcf.gz \
  > ${DIR}/${SAMPLE}.stats
done
```

| Sample | [3]ts | [4]tv | [5]ts/tv | [6]ts (1st ALT) | [7]tv (1st ALT) | [8]ts/tv (1st ALT) |
| ------ | ----- | ----- | -------- | --------------- | --------------- | ------------------ |
| HG002  | 30016 | 11709 | 2.56     | 30002           | 11693           | 2.57               |
| HG003  | 29880 | 11747 | 2.54     | 29871           | 11731           | 2.55               |
| HG004  | 30133 | 11860 | 2.54     | 30120           | 11848           | 2.54               |

If you want to restrict to the truth BED files, use this command:

```
declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  bcftools stats -f PASS \
    -T ${DIR}/${SAMPLE}_truth.bed \
    ${DIR}/${SAMPLE}.vcf.gz \
  > ${DIR}/${SAMPLE}.with_truth_bed.stats
done
```

Which resulted in this table:

| Sample | [3]ts | [4]tv | [5]ts/tv | [6]ts (1st ALT) | [7]tv (1st ALT) | [8]ts/tv (1st ALT) |
| ------ | ----- | ----- | -------- | --------------- | --------------- | ------------------ |
| HG002  | 24474 | 9255  | 2.64     | 24469           | 9245            | 2.65               |
| HG003  | 24175 | 9182  | 2.63     | 24172           | 9174            | 2.63               |
| HG004  | 24313 | 9334  | 2.60     | 24306           | 9327            | 2.61               |

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
Failed Filters               : 14405
Passed Filters               : 45447
SNPs                         : 41696
MNPs                         : 0
Insertions                   : 1909
Deletions                    : 1821
Indels                       : 19
Same as reference            : 2
SNP Transitions/Transversions: 2.56 (41882/16353)
Total Het/Hom ratio          : 1.51 (27315/18130)
SNP Het/Hom ratio            : 1.52 (25178/16518)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.10 (1001/908)
Deletion Het/Hom ratio       : 1.59 (1117/704)
Indel Het/Hom ratio          : - (19/0)
Insertion/Deletion ratio     : 1.05 (1909/1821)
Indel/SNP+MNP ratio          : 0.09 (3749/41696)
```

HG003:

```
Location                     : /data/HG003.vcf.gz
Failed Filters               : 15215
Passed Filters               : 45306
SNPs                         : 41600
MNPs                         : 0
Insertions                   : 1902
Deletions                    : 1784
Indels                       : 18
Same as reference            : 2
SNP Transitions/Transversions: 2.52 (41676/16511)
Total Het/Hom ratio          : 1.50 (27156/18148)
SNP Het/Hom ratio            : 1.51 (25031/16569)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.15 (1019/883)
Deletion Het/Hom ratio       : 1.56 (1088/696)
Indel Het/Hom ratio          : - (18/0)
Insertion/Deletion ratio     : 1.07 (1902/1784)
Indel/SNP+MNP ratio          : 0.09 (3704/41600)
```

HG004:

```
Location                     : /data/HG004.vcf.gz
Failed Filters               : 14832
Passed Filters               : 45681
SNPs                         : 41965
MNPs                         : 0
Insertions                   : 1899
Deletions                    : 1796
Indels                       : 18
Same as reference            : 3
SNP Transitions/Transversions: 2.55 (41746/16352)
Total Het/Hom ratio          : 1.58 (28000/17678)
SNP Het/Hom ratio            : 1.60 (25850/16115)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.17 (1023/876)
Deletion Het/Hom ratio       : 1.61 (1109/687)
Indel Het/Hom ratio          : - (18/0)
Insertion/Deletion ratio     : 1.06 (1899/1796)
Indel/SNP+MNP ratio          : 0.09 (3713/41965)
```

### Run hap.py to calculate the accuracy of DeepVariant generated call sets

```
sudo docker pull pkrusche/hap.py

declare -a trio=(HG002 HG003 HG004)
for SAMPLE in "${trio[@]}"
do
  sudo docker run -i \
    -v "${DIR}":"/data" \
    pkrusche/hap.py /opt/hap.py/bin/hap.py \
    "/data/${SAMPLE}_truth.vcf.gz" \
    "/data/${SAMPLE}.vcf.gz" \
    -f "/data/${SAMPLE}_truth.bed" \
    -T "/data/${CAPTURE_BED}" \
    -r "/data/hs37d5.fa" \
    -o "/data/${SAMPLE}.happy.output" \
    --engine=vcfeval > ${DIR}/${SAMPLE}.stdout
done
```

Accuracy F1 scores:

Sample | Indel    | SNP
------ | -------- | --------
HG002  | 0.969993 | 0.999036
HG003  | 0.969311 | 0.998921
HG004  | 0.967794 | 0.999212
