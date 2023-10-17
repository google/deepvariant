# VG Giraffe + DeepVariant case study
---

This is an example to run `vg giraffe`, so we can go from FASTQs --> BAM.

For simplicity and consistency, we run the following with a
[Google Cloud instance with 64 cores](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform).

I added more disks because 300G is not enough for the example below. I changed
it to `--boot-disk-size "1000"`.

## Install softwares that will be used later

```bash
sudo apt update -y
sudo apt-get -y install aria2 docker.io samtools
```

## Download input FASTQ files

```bash
DATA_DIR=${PWD}/data
mkdir -p ${DATA_DIR}
gcloud storage cp gs://brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R?.fastq.gz ${DATA_DIR}/
```

## Download VG files

```bash
gcloud storage cp gs://hprc-pangenomes/hprc-v1.0-mc-chm13-minaf.0.1.gbz ${DATA_DIR}/
```

I used `aria2c` to download these files. You can use other approaches as well.

```bash
aria2c -c -x10 -s10 -d "${DATA_DIR}" https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/filtered/hprc-v1.0-mc-chm13-minaf.0.1.min
aria2c -c -x10 -s10 -d "${DATA_DIR}" https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/filtered/hprc-v1.0-mc-chm13-minaf.0.1.dist
```

## Run `vg giraffe` with one command to get from FASTQs to BAM.

```
wget -P ${DATA_DIR} https://storage.googleapis.com/hprc-pangenomes/GRCh38.path_list.txt
SAMP=HG003
READS1=${DATA_DIR}/HG003.novaseq.pcr-free.35x.R1.fastq.gz
READS2=${DATA_DIR}/HG003.novaseq.pcr-free.35x.R2.fastq.gz

VG_DOCKER_VERSION=v1.47.0
sudo docker pull quay.io/vgteam/vg:${VG_DOCKER_VERSION}

time sudo docker run --rm -v ${DATA_DIR}:${DATA_DIR}  \
  quay.io/vgteam/vg:${VG_DOCKER_VERSION} \
  vg giraffe --progress \
  --read-group "ID:1 LB:lib1 SM:${SAMP} PL:illumina PU:unit1" \
  --sample "${SAMP}" \
  -o BAM --ref-paths ${DATA_DIR}/GRCh38.path_list.txt \
  -P -L 3000 \
  -f $READS1 -f $READS2 \
  -Z ${DATA_DIR}/hprc-v1.0-mc-chm13-minaf.0.1.gbz \
  -d ${DATA_DIR}/hprc-v1.0-mc-chm13-minaf.0.1.dist \
  -m ${DATA_DIR}/hprc-v1.0-mc-chm13-minaf.0.1.min \
  -t $(nproc) > reads.unsorted.bam
```

NOTE: No need to sort this yet, because we'll need to sort it in the next step.

## Runtime

On my machine, the last few lines of the log showed:

```
Mapped 838385300 reads across 64 threads in 17842.7 seconds with 2.03384 additional single-threaded seconds.
Mapping speed: 734.18 reads per second per thread
Used 1.13562e+06 CPU-seconds (including output).
Achieved 738.262 reads per CPU-second (including output)
Memory footprint: 60.2227 GB

real    300m17.121s
user    0m47.536s
sys     3m56.936s
```

File size:

```
$ ls -lh reads.unsorted.bam
-rw-rw-r-- 1 user user Oct 16 22:35 reads.unsorted.bam
```

## Fix the BAM's contig name

```bash
INBAM=reads.unsorted.bam
BAM=reads.sorted.chrfixed.bam

## Download the reference fasta, index, and dict (although for this you only need the .dict).
wget https://storage.googleapis.com/hprc-pangenomes/hg38.fa
wget https://storage.googleapis.com/hprc-pangenomes/hg38.fa.fai
wget https://storage.googleapis.com/hprc-pangenomes/hg38.dict

# Prepare the new header.
samtools view -H $INBAM | grep ^@HD > new_header.sam
grep ^@SQ hg38.dict | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
samtools view -H $INBAM  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam
## write new BAM, removing the "GRCh38." prefixes
time cat <(cat new_header.sam) <(samtools view $INBAM) | sed -e "s/GRCh38.//g" | samtools sort --threads 10 -m 2G -O BAM > ${BAM}
# Index the BAM.
samtools index -@$(nproc) ${BAM}
```

The step with `time` above took:

```
real    81m40.953s
user    170m42.641s
sys     27m33.029s
```

File size:

```
$ ls -lh reads.sorted.chrfixed.bam
-rw-rw-r-- 1 user user 40G Oct 17 02:23 reads.sorted.chrfixed.bam
```

## Run DeepVariant With `min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true`

Get the same reference we used for
[DeepVariant Case Study](deepvariant-case-study.md)

```bash
FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > ${DATA_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx ${DATA_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

And then, run DeepVariant.

(If you want to test on one smaller chromosome first, you can add
`--regions chr20` like what we did in
[DeepVariant Case Study](deepvariant-case-study.md).)

```bash
BIN_VERSION="1.6.0"

sudo docker pull google/deepvariant:"${BIN_VERSION}"

time sudo docker run --rm \
  -v "${DATA_DIR}":"${DATA_DIR}" \
  -v "${PWD}:${PWD}" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=${DATA_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --reads=${PWD}/${BAM} \
  --output_vcf=${PWD}/min_mapping_quality-keep_legacy_allele_counter_behavior-normalize_reads-vg.vcf.gz \
  --output_gvcf=${PWD}/min_mapping_quality-keep_legacy_allele_counter_behavior-normalize_reads-vg.g.vcf.gz \
  --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true" \
  --num_shards=$(nproc)
```

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 116m20.443s
call_variants                    | 173m55.861s
postprocess_variants (with gVCF) | 28m18.185s


### Run hap.py

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38

curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

```bash
mkdir -p happy

sudo docker pull jmcdani20/hap.py:v0.3.12

sudo docker run --rm \
  -v "${DATA_DIR}":"${DATA_DIR}" \
  -v "${PWD}:${PWD}" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  ${PWD}/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  ${PWD}/min_mapping_quality-keep_legacy_allele_counter_behavior-normalize_reads-vg.vcf.gz \
  -f ${PWD}/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r ${DATA_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -o ${PWD}/happy/happy.output \
  --engine=vcfeval \
  --pass-only
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL       504501    502108      2393       948864      1453     423968    903    342       0.995257          0.997232        0.446816         0.996243                     NaN                     NaN                   1.489759                   1.965681
INDEL   PASS       504501    502108      2393       948864      1453     423968    903    342       0.995257          0.997232        0.446816         0.996243                     NaN                     NaN                   1.489759                   1.965681
  SNP    ALL      3327496   3314384     13112      3725023      4167     404651   1751    300       0.996059          0.998745        0.108630         0.997400                2.102576                2.031762                   1.535137                   1.492075
  SNP   PASS      3327496   3314384     13112      3725023      4167     404651   1751    300       0.996059          0.998745        0.108630         0.997400                2.102576                2.031762                   1.535137                   1.492075
```

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502108   | 2393     | 1453     | 0.995257      | 0.997232         | 0.996243        |
| SNP   | 3314384  | 13112    | 4167     | 0.996059      | 0.998745         | 0.9974          |

This can be compared with
https://github.com/google/deepvariant/blob/r1.6/docs/metrics.md#accuracy.

Which shows that `vg giraffe` improves F1:

- Indel F1: 0.995998 --> 0.996243
- SNP F1: 0.996237 --> 0.9974
