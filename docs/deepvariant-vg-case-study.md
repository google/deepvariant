# Using graph genomes: VG Giraffe + DeepVariant case study
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

Get binaries `vg` 1.51.0 and `kmc`:

```bash
wget https://github.com/refresh-bio/KMC/releases/download/v3.2.2/KMC3.2.2.linux.x64.tar.gz
tar zxf KMC3.2.2.linux.x64.tar.gz bin/kmc
mv bin/kmc ${DATA_DIR}/
wget https://github.com/vgteam/vg/releases/download/v1.51.0/vg -O ${DATA_DIR}/vg
chmod +x ${DATA_DIR}/vg ${DATA_DIR}/kmc
```

Get the graph (.gbz) and haplotype index (.hapl).
I used `aria2c` to download these files. You can use other approaches as well.

```bash
aria2c -c -x10 -s10 -d "${DATA_DIR}" https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
aria2c -c -x10 -s10 -d "${DATA_DIR}" https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.hapl
```

## Run `vg giraffe` with one command to get from FASTQs to BAM.

Put the paths name into a file named HG003.fq.paths:

```bash
cat > HG003.fq.paths <<- EOM
${DATA_DIR}/HG003.novaseq.pcr-free.35x.R1.fastq.gz
${DATA_DIR}/HG003.novaseq.pcr-free.35x.R2.fastq.gz
EOM
```

Run `kmc`` on this file. I used -t$(nproc) to use all cores, and $TMPDIR for a
scratch directory:

```bash
TMPDIR=$(mktemp -d)
time ${DATA_DIR}/kmc -k29 -m128 -okff -t$(nproc) @HG003.fq.paths ${DATA_DIR}/HG003.fq $TMPDIR
```

Output on the terminal:

```
*******************************************************************************
Stage 1: 55%
Stage 1: 100%
Stage 2: 100%
1st stage: 944.485s
2nd stage: 506.81s
Total    : 1451.29s
Tmp size : 110071MB

Stats:
   No. of k-mers below min. threshold :   5828096589
   No. of k-mers above max. threshold :            0
   No. of unique k-mers               :   8581831809
   No. of unique counted k-mers       :   2753735220
   Total no. of k-mers                : 103092565745
   Total no. of reads                 :    838385300
   Total no. of super-k-mers          :   9929565346

real    24m11.431s
user    142m37.817s
sys     8m14.566s
```

Run `giraffe`` on the graph, haplotype index, kmers and reads:

```bash
${DATA_DIR}/vg paths \
  -x ${DATA_DIR}/hprc-v1.1-mc-grch38.gbz \
  -L -Q GRCh38 > ${DATA_DIR}/GRCh38.path_list.txt
```

```bash
time ${DATA_DIR}/vg giraffe --progress \
 --read-group "ID:1 LB:lib1 SM:HG003 PL:illumina PU:unit1" \
  --sample "HG003" \
  -o BAM --ref-paths ${DATA_DIR}/GRCh38.path_list.txt \
  -P -L 3000 \
  -f ${DATA_DIR}/HG003.novaseq.pcr-free.35x.R1.fastq.gz \
  -f ${DATA_DIR}/HG003.novaseq.pcr-free.35x.R2.fastq.gz \
  -Z ${DATA_DIR}/hprc-v1.1-mc-grch38.gbz \
  --kff-name ${DATA_DIR}/HG003.fq.kff \
  --haplotype-name ${DATA_DIR}/hprc-v1.1-mc-grch38.hapl \
  -t $(nproc) > reads.unsorted.bam
```

NOTE: No need to sort this yet, because we'll need to sort it in the next step.

## Runtime

On my machine, the last few lines of the log showed:

```
Mapped 838385300 reads across 64 threads in 14093.4 seconds with 3.25431 additional single-threaded seconds.
Mapping speed: 929.496 reads per second per thread
Used 896175 CPU-seconds (including output).
Achieved 935.515 reads per CPU-second (including output)
Memory footprint: 61.0703 GB

real    283m10.368s
user    15260m35.845s
sys     214m57.882s
```

File size:

```
$ ls -lh reads.unsorted.bam
-rw-rw-r-- 1 pichuan pichuan 69G Nov  1 23:56 reads.unsorted.bam
```

Then, clean up contig names, and sort:

```bash
INBAM=reads.unsorted.bam
BAM=reads.sorted.chrfixed.bam
time samtools view -h $INBAM | sed -e "s/GRCh38#0#//g" | samtools sort --threads 10 -m 2G -O BAM > ${BAM}
# Index the BAM.
samtools index -@$(nproc) ${BAM}
```

The step with `time` above took:

```
real    73m19.172s
user    178m59.088s
sys     24m36.986s
```


File size:

```
$ ls -lh reads.sorted.chrfixed.bam
-rw-rw-r-- 1 pichuan pichuan 40G Nov  2 02:09 reads.sorted.chrfixed.bam
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
make_examples                    | 116m37.385s
call_variants                    | 214m37.055s
postprocess_variants (with gVCF) | 30m59.968s


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
INDEL    ALL       504501    502199      2302       960061      1526     434935    906    371       0.995437          0.997094        0.453029         0.996265                     NaN                     NaN                   1.489759                   1.952023
INDEL   PASS       504501    502199      2302       960061      1526     434935    906    371       0.995437          0.997094        0.453029         0.996265                     NaN                     NaN                   1.489759                   1.952023
  SNP    ALL      3327496   3316515     10981      3858659      5550     534709   2104    475       0.996700          0.998330        0.138574         0.997514                2.102576                1.970783                   1.535137                   1.436586
  SNP   PASS      3327496   3316515     10981      3858659      5550     534709   2104    475       0.996700          0.998330        0.138574         0.997514                2.102576                1.970783                   1.535137                   1.436586
```

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502199   | 2302     | 1526     | 0.995437      | 0.997094         | 0.996265        |
| SNP   | 3316515  | 10981    | 5550     | 0.9967        | 0.99833          | 0.997514        |

This can be compared with
https://github.com/google/deepvariant/blob/r1.6/docs/metrics.md#accuracy.

Which shows that `vg giraffe` improves F1:

- Indel F1: 0.995998 --> 0.996265
- SNP F1: 0.996237 --> 0.997514
