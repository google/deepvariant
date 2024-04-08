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

Get binaries `vg` 1.55.0 and `kmc`:

```bash
wget https://github.com/refresh-bio/KMC/releases/download/v3.2.2/KMC3.2.2.linux.x64.tar.gz
tar zxf KMC3.2.2.linux.x64.tar.gz bin/kmc
mv bin/kmc ${DATA_DIR}/
wget https://github.com/vgteam/vg/releases/download/v1.55.0/vg -O ${DATA_DIR}/vg
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
Stage 1: 100%
Stage 2: 100%
1st stage: 1080.93s
2nd stage: 383.703s
Total    : 1464.63s
Tmp size : 110071MB

Stats:
   No. of k-mers below min. threshold :   5828096589
   No. of k-mers above max. threshold :            0
   No. of unique k-mers               :   8581831809
   No. of unique counted k-mers       :   2753735220
   Total no. of k-mers                : 103092565745
   Total no. of reads                 :    838385300
   Total no. of super-k-mers          :   9929565346

real    24m24.961s
user    157m44.462s
sys     10m9.054s
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
Mapped 838385300 reads across 64 threads in 16094.2 seconds with 2.27046 additional single-threaded seconds.
Mapping speed: 813.942 reads per second per thread
Used 1.02948e+06 CPU-seconds (including output).
Achieved 814.375 reads per CPU-second (including output)
Memory footprint: 60.8593 GB

real    322m57.058s
user    17482m34.387s
sys     257m57.409s
```

File size:

```
$ ls -lh reads.unsorted.bam
-rw-rw-r-- 1 pichuan pichuan 69G Apr  4 08:39 reads.unsorted.bam
```

Then, clean up contig names, and sort:

```bash
INBAM=reads.unsorted.bam
BAM=HG003.novaseq.pcr-free.35x.vg-1.55.0.bam
time samtools view -h $INBAM | sed -e "s/GRCh38#0#//g" | samtools sort --threads 10 -m 2G -O BAM > ${BAM}
# Index the BAM.
samtools index -@$(nproc) ${BAM}
```

The step with `time` above took:

```
real    77m14.107s
user    184m18.458s
sys     28m51.487s
```


File size:

```
$ ls -lh HG003.novaseq.pcr-free.35x.vg-1.55.0.bam
-rw-rw-r-- 1 pichuan pichuan 39G Apr  4 17:06 HG003.novaseq.pcr-free.35x.vg-1.55.0.bam
```

This file can also be found at:

`gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam`

## Run DeepVariant With `min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true`

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
BIN_VERSION="1.6.1"

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
  --make_examples_extra_args="min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true" \
  --num_shards=$(nproc)
```

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 101m31.676s
call_variants                    | 215m33.631s
postprocess_variants (with gVCF) | 24m44.242s


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
INDEL    ALL       504501    502283      2218       958181      1471     433079    913    351       0.995604          0.997199        0.451980         0.996400                     NaN                     NaN                   1.489759                   1.954212
INDEL   PASS       504501    502283      2218       958181      1471     433079    913    351       0.995604          0.997199        0.451980         0.996400                     NaN                     NaN                   1.489759                   1.954212
  SNP    ALL      3327496   3316374     11122      3820052      4177     497662   1686    344       0.996658          0.998743        0.130276         0.997699                2.102576                1.991054                   1.535137                   1.457635
  SNP   PASS      3327496   3316374     11122      3820052      4177     497662   1686    344       0.996658          0.998743        0.130276         0.997699                2.102576                1.991054                   1.535137                   1.457635
```

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502283   | 2218     | 1471     | 0.995604      | 0.997199         | 0.9964          |
| SNP   | 3316374  | 11122    | 4177     | 0.996658      | 0.998743         | 0.997699        |

This can be compared with
https://github.com/google/deepvariant/blob/r1.6.1/docs/metrics.md#accuracy.

Which shows that `vg giraffe` improves F1:

- Indel F1: 0.995998 --> 0.9964
- SNP F1: 0.996237 --> 0.997699
