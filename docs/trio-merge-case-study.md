# Advanced Case study: Merge DeepVariant Whole Exome Sequencing (WES) trio call sets

## Overview

This document outlines all the steps and considerations for calling and
merging a trio using DeepVariant and GLnexus. The process involves 3 major
stages: running DeepVariant to create individual genome call sets, running
GLnexus to merge call sets, and analyzing the merged call set.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. The steps below were done on a machine with this
[example command to start a machine](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform).

The steps in this document can be extended to merge larger cohort as well.

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

With the example command below, it runs DeepVariant on the trio one by one.

```
N_SHARDS=64  # Change to your number of cores.
CAPTURE_BED=agilent_sureselect_human_all_exon_v5_b37_targets.bed
VERSION=0.8.0

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

## Single sample quality metrics

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
| HG002  | 30221 | 11820 | 2.56     | 30209           | 11799           | 2.56               |
| HG003  | 30150 | 11921 | 2.53     | 30136           | 11897           | 2.53               |
| HG004  | 30379 | 12031 | 2.53     | 30366           | 12011           | 2.53               |

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
| HG002  | 24472 | 9256  | 2.64     | 24467           | 9247            | 2.65               |
| HG003  | 24172 | 9184  | 2.63     | 24169           | 9176            | 2.63               |
| HG004  | 24313 | 9340  | 2.60     | 24306           | 9333            | 2.60               |

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
Failed Filters               : 9949
Passed Filters               : 45657
SNPs                         : 42007
MNPs                         : 0
Insertions                   : 1870
Deletions                    : 1762
Indels                       : 15
Same as reference            : 3
SNP Transitions/Transversions: 2.56 (42090/16470)
Total Het/Hom ratio          : 1.52 (27556/18098)
SNP Het/Hom ratio            : 1.54 (25477/16530)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.11 (984/886)
Deletion Het/Hom ratio       : 1.58 (1080/682)
Indel Het/Hom ratio          : - (15/0)
Insertion/Deletion ratio     : 1.06 (1870/1762)
Indel/SNP+MNP ratio          : 0.09 (3647/42007)
```

HG003:

```
Location                     : /data/HG003.vcf.gz
Failed Filters               : 10444
Passed Filters               : 45614
SNPs                         : 42029
MNPs                         : 0
Insertions                   : 1842
Deletions                    : 1715
Indels                       : 21
Same as reference            : 7
SNP Transitions/Transversions: 2.51 (41970/16707)
Total Het/Hom ratio          : 1.51 (27450/18157)
SNP Het/Hom ratio            : 1.53 (25407/16622)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.14 (980/862)
Deletion Het/Hom ratio       : 1.55 (1042/673)
Indel Het/Hom ratio          : - (21/0)
Insertion/Deletion ratio     : 1.07 (1842/1715)
Indel/SNP+MNP ratio          : 0.09 (3578/42029)
```

HG004:

```
Location                     : /data/HG004.vcf.gz
Failed Filters               : 10312
Passed Filters               : 46002
SNPs                         : 42373
MNPs                         : 0
Insertions                   : 1859
Deletions                    : 1750
Indels                       : 16
Same as reference            : 4
SNP Transitions/Transversions: 2.54 (41992/16542)
Total Het/Hom ratio          : 1.60 (28335/17663)
SNP Het/Hom ratio            : 1.63 (26231/16142)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 1.15 (994/865)
Deletion Het/Hom ratio       : 1.67 (1094/656)
Indel Het/Hom ratio          : - (16/0)
Insertion/Deletion ratio     : 1.06 (1859/1750)
Indel/SNP+MNP ratio          : 0.09 (3625/42373)
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
HG002  | 0.966611 | 0.999110
HG003  | 0.961020 | 0.998906
HG004  | 0.966357 | 0.999153

## Merge the trio samples using GLnexus

### Download GLnexus

```
wget https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.1.10/glnexus_cli
chmod +x glnexus_cli
```

### Run GLnexus to merge 3 gVCFs

```
time ./glnexus_cli \
 --config DeepVariant \
 --bed ${DIR}/${CAPTURE_BED} \
 ${DIR}/HG004.g.vcf.gz ${DIR}/HG003.g.vcf.gz ${DIR}/HG002.g.vcf.gz \
 | bcftools view - | bgzip -c > ${DIR}/deepvariant.cohort.vcf.gz
```

When we ran on this WES trio, it took only about 13 seconds. However, with WGS
data or larger cohort, you might want to consider using jemalloc to improve
performance. See https://github.com/dnanexus-rnd/GLnexus/wiki/Performance for
more detail.

NOTE: If you need to re-run this, you often need to clean up the GLnexus.DB
directory it generates before you can successfully re-run.

## Annotate the merged VCF with Mendelian discordance information using RTG Tools

Create SDF template from our reference file:

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
4 non-pass records were skipped
Concordance Sample_Diag-excap51-HG002-EEogPU: F:60687/61400 (98.84%)  M:61157/61380 (99.64%)  F+M:60174/61156 (98.39%)
Sample Sample_Diag-excap51-HG002-EEogPU has less than 99.0 concordance with the father (Sample_Diag-excap51-HG003-EEogPU). Check for incorrect pedigree or sample mislabelling.
883/62027 (1.42%) records did not conform to expected call ploidy
61652/62027 (99.40%) records were variant in at least 1 family member and checked for Mendelian constraints
436/61652 (0.71%) records had indeterminate consistency status due to incomplete calls
1033/61652 (1.68%) records contained a violation of Mendelian constraints
```

From this report, we know that there is a 1.68% Mendelian violation rate, and
0.71% of the records had incomplete calls (with `.`) so RTG couldn't determine
whether there is violation or not.
