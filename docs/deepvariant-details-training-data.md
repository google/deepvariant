# DeepVariant training data

## WGS training data v0.7

We used: 10 HG001 PCR-free, 2 HG005 PCR-free, 4 HG001 PCR+ for training.

Among these 16 BAM files, 6 of them are from public sources:

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
