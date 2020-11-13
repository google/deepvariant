# DeepVariant training data

### WGS models

version | Replicates                             | #examples
------- | -------------------------------------- | -----------
v0.4    | 9 HG001                                | 85,323,867
v0.5    | 9 HG001<br>2 HG005<br>78 HG001 WES<br>1 HG005 WES<sup>[(1)](#vfootnote1)</sup> | 115,975,740
v0.6    | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>4 HG001 PCR+     | 156,571,227
v0.7    | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>4 HG001 PCR+     | 158,571,078
v0.8    | 12 HG001 PCR-free<br>2 HG005 PCR-free<br>4 HG001 PCR+<br>(and, more `dowsample_fraction` since last version)     | 346,505,686
v0.9    | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>2 HG006 PCR-free<br>2 HG007 PCR-free<br>5 HG001 PCR+     | 325,202,093
v0.10   | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>2 HG006 PCR-free<br>2 HG007 PCR-free<br>5 HG001 PCR+     | 339,410,078
v1.0    | 11 HG001<br>2 HG005-HG007<br>2 HG002-HG004<sup>[(7)](#vfootnote7)</sup>     | 317,486,837

### WES models

version | Replicates                  | #examples
------- | --------------------------- | ------------------------------
v0.5    | 78 HG001<br>1 HG005 | 15,714,062
v0.6    | 78 HG001<br>1 HG005<sup>[(2)](#vfootnote2)</sup> | 15,705,449
v0.7    | 78 HG001<br>1 HG005 | 15,704,197
v0.8    | 78 HG001<br>1 HG005<sup>[(3)](#vfootnote3)</sup> | 18,683,247
v0.9    | 81 HG001<br>1 HG005<sup>[(3)](#vfootnote3)[(4)](#vfootnote4)[(5)](#vfootnote5)</sup> | 61,953,965
v0.10   | 81 HG001<br>1 HG005<sup>[(3)](#vfootnote3)[(4)](#vfootnote4)[(5)](#vfootnote5)</sup> | 61,953,965
v1.0    | 32 HG001<br>9 HG002<br>6 HG003<br>6 HG004<br>12 HG005<br>9 HG006<br>9 HG007<sup>[(7)](#vfootnote7)</sup> | 10,716,281

### PACBIO models

version | Replicates                  | #examples
------- | --------------------------- | ------------------------------
v0.8    | 16 HG002 | 160,025,931
v0.9    | 49 HG002 <sup>[(6)](#vfootnote6)</sup> | 357,507,235
v0.10   | 49 HG002, 2 HG003, 2 HG004, 1 HG002 (amplified) <sup>[(6)](#vfootnote6)</sup> | 472,711,858
v1.0    | 1 HG001<br>2 HG002<br>2 HG003<br>2 HG004<br>1 HG005 <sup>[(8)](#vfootnote8)</sup>  | 302,331,948

<a name="vfootnote1">(1)</a>: In v0.5, we experimented with adding whole exome
sequencing data into training data. In v0.6, we took it out because it didn't
improve the WGS accuracy.

<a name="vfootnote2">(2)</a>: The training data are from the same replicates as
v0.5. The number of examples changed because of the update in
[haplotype_labeler](https://github.com/google/deepvariant/tree/r0.6/deepvariant/labeler/haplotype_labeler.py).

<a name="vfootnote3">(3)</a>: In v0.8, we used the
[Platinum Genomes Truthset](https://github.com/Illumina/PlatinumGenomes) to
create more training examples outside the GIAB confident regions.

<a name="vfootnote4">(4)</a>: Previously, we split train/tune by leaving 3 WES
for tuning. Starting from this release, we leave out chr1 and chr20 from
training, and use chr1 for tuning.

<a name="vfootnote5">(5)</a>: Starting from this version, we padded (100bps on
both sides) of the capture BED and used that for generating training examples.
We also added more `downsample_fraction`.

<a name="vfootnote6">(6)</a>: (Before v1.0) PacBio is the only one we currently
uses HG002 in training and tuning.

<a name="vfootnote7">(7)</a>: In v1.0, we train on HG002-HG004 for WGS as well,
but only using examples from the region of NIST truth confident region v4.2
subtracting v3.3.2.

<a name="vfootnote8">(8)</a>: In v1.0, PacBio training data contains training
examples with haplotag sorted images and unsorted images.


## WGS training data:

Here are 6 of WGS BAM files from public sources:

BAM file (`--reads`)                                                                                                                                     | PCR-free? | FASTA file (`--ref`)                                                                                                         | Truth VCF (`--truth_variants`)                                                                          | BED file (`--confident_regions`)
-------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- | --------------------------------
[HG001-NA12878-pFDA.merged.sorted.bam](https://console.cloud.google.com/storage/browser/deepvariant/public-training-data)<sup>[(1)](#myfootnote1)</sup>  | Yes      | [GRCh38_Verily_v1.genome.fa](https://console.cloud.google.com/storage/browser/genomics-public-data/references/GRCh38_Verily) | [NISTv3.3.2/GRCh38](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/) | [NISTv3.3.2/GRCh38](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/)
[NA12878D_HiSeqX_R1.deduplicated.bam](https://console.cloud.google.com/storage/browser/deepvariant/public-training-data)<sup>[(2)](#myfootnote2)</sup>   | No       | [hs37d5.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence)                 | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/) | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/)
[NA12878J_HiSeqX_R1.deduplicated.bam](https://console.cloud.google.com/storage/browser/deepvariant/public-training-data)<sup>[(2)](#myfootnote2)</sup>   | No       | [hs37d5.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence)                 | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/) | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/)
[NA12878-Rep01_S1_L001_001_markdup.bam](https://console.cloud.google.com/storage/browser/deepvariant/public-training-data)<sup>[(2)](#myfootnote2)</sup> | No       | [hs37d5.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence)                 | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/) | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/)
N3C9-2plex1-L1-171212B-NA12878-1_S1_L001_001_markdup.bam<sup>[(3)](#myfootnote3)</sup>                                                                   | Yes      | [hs37d5.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence)                 | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/) | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/)
NexteraFlex-2plex1-L1-NA12878-1_S1_L001_001_markdup.bam<sup>[(4)](#myfootnote4)</sup>                                                                    | No       | [hs37d5.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence)                 | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/) | [NISTv3.3.2/GRCh37](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/)

<a name="myfootnote1">(1)</a>: FASTQ files from
[Precision FDA Truth Challenge](https://precision.fda.gov/challenges/truth).

<a name="myfootnote2">(2)</a>: BAM files provided by
[DNAnexus](https://www.dnanexus.com/).

<a name="myfootnote3">(3)</a>: FASTQ files from
[BaseSpace public data](https://basespace.illumina.com/datacentral): `NovaSeq S1
Xp: TruSeq Nano 350 (Replicates of
NA12878)/Samples/N3C9_2plex1_L1_171212B_NA12878-1/Files/N3C9-2plex1-L1-171212B-NA12878-1_S1_L001_R1_001.fastq.gz`
and `N3C9-2plex1-L1-171212B-NA12878-1_S1_L001_R2_001.fastq.gz`

<a name="myfootnote4">(4)</a>: FASTQ files from
[BaseSpace public data](https://basespace.illumina.com/datacentral): `NovaSeq S1
Xp: Nextera DNA Flex (Replicates of
NA12878)/Samples/NexteraFlex_2plex1_L1_NA12878-1/Files/NexteraFlex-2plex1-L1-NA12878-1_S1_L001_R1_001.fastq.gz`
and `NexteraFlex-2plex1-L1-NA12878-1_S1_L001_R2_001.fastq.gz`

We generated our own BAM files using BWA-MEM to map the reads to the reference,
and sorts the output. We also mark duplicated reads.
