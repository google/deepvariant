# DeepTrio whole genome sequencing case study

In this case study, we describe applying DeepTrio to a real WGS trio. Then we
assess the quality of the DeepTrio variant calls with `hap.py`. In addition we
evaluate a mendelian violation rate for a merged VCF.

To make it faster to run over this case study, we run only on chromosome 20.


## Prepare environment

### Tools

[Docker](https://docs.docker.com/get-docker/) will be used to run DeepTrio and
[hap.py](https://github.com/illumina/hap.py),

### Download Reference

We will be using GRCh38 for this case study.

```bash
mkdir -p reference

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids

curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai
```

### Download Genome in a Bottle Benchmarks

We will benchmark our variant calls against v4.2.1 of the Genome in a Bottle
small variant benchmarks for HG002, HG003, and HG004 trio.

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio

curl ${FTPDIR}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

curl ${FTPDIR}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

curl ${FTPDIR}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

### Download HG002, HG003, and HG004 BAM files

We'll use HG002, HG003, HG004 Illumina WGS reads publicly available from the
[PrecisionFDA Truth v2 Challenge](https://precision.fda.gov/challenges/10).

```bash
mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata

curl ${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam > input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
curl ${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai > input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai

curl ${HTTPDIR}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam > input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
curl ${HTTPDIR}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai > input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai

curl ${HTTPDIR}/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam > input/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
curl ${HTTPDIR}/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai > input/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai
```

## Running DeepTrio with one command

DeepTrio pipeline consists of 4 steps: `make_examples`, `call_variants`,
`postprocess_variants` and `GLnexus merge`. It is possible to run DeepTrio with
one command using the `run_deepvariant` script. GLnexus is run as a separate
command.

### Running on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="1.7.0"

sudo docker pull google/deepvariant:deeptrio-"${BIN_VERSION}"

time sudo docker run \
  -v "${PWD}/input":"/input"   \
  -v "${PWD}/output":"/output"  \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:deeptrio-"${BIN_VERSION}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type WGS \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads_child /input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
  --reads_parent1 /input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
  --reads_parent2 /input/HG004.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
  --output_vcf_child /output/HG002.output.vcf.gz \
  --output_vcf_parent1 /output/HG003.output.vcf.gz \
  --output_vcf_parent2 /output/HG004.output.vcf.gz \
  --sample_name_child 'HG002' \
  --sample_name_parent1 'HG003' \
  --sample_name_parent2 'HG004' \
  --num_shards $(nproc)  \
  --regions chr20 \
  --intermediate_results_dir /output/intermediate_results_dir \
  --output_gvcf_child /output/HG002.g.vcf.gz \
  --output_gvcf_parent1 /output/HG003.g.vcf.gz \
  --output_gvcf_parent2 /output/HG004.g.vcf.gz
```

By specifying `--model_type WGS`, you'll be using a model that is best suited
for Illumina Whole Genome Sequencing data.

NOTE: If you want to run each of the steps separately, add `--dry_run=true`
to the command above to figure out what flags you need in each step. Based on
the different model types, different flags are needed in the `make_examples`
step.

`--intermediate_results_dir` flag is optional. By specifying it, the
intermediate outputs of `make_examples` and `call_variants` stages can be found
in the directory. After the command, you can find these files in the directory:

```
call_variants_output_child.tfrecord.gz
call_variants_output_parent1.tfrecord.gz
call_variants_output_parent2.tfrecord.gz

gvcf_child.tfrecord-?????-of-?????.gz
gvcf_parent1.tfrecord-?????-of-?????.gz
gvcf_parent2.tfrecord-?????-of-?????.gz

make_examples_child.tfrecord-?????-of-?????.gz
make_examples_parent1.tfrecord-?????-of-?????.gz
make_examples_parent2.tfrecord-?????-of-?????.gz
```

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md).

## Merge VCFs using GLnexus

At this step we take all 3 VCFs generated in the previous step and merge them
using GLnexus.

```bash
# bcftools and bgzip are now included in our docker images.
# You can also install them separately.
sudo docker run \
  -v "${PWD}/output":"/output" \
  quay.io/mlin/glnexus:v1.2.7 \
  /usr/local/bin/glnexus_cli \
  --config DeepVariant_unfiltered \
  /output/HG002.g.vcf.gz \
  /output/HG003.g.vcf.gz \
  /output/HG004.g.vcf.gz \
  | sudo docker run -i google/deepvariant:deeptrio-"${BIN_VERSION}" \
    bcftools view - \
  | sudo docker run -i google/deepvariant:deeptrio-"${BIN_VERSION}" \
    bgzip -c > output/HG002_trio_merged.vcf.gz
```

After completion of GLnexus command we should have a new merged VCF file in the
output directory.

```
HG002_trio_merged.vcf.gz
```

## Benchmark on chr20

### Calculate mendelian violation rate

```bash
sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/reference":"/reference" \
  realtimegenomics/rtg-tools format \
  -o /reference/GRCh38_no_alt_analysis_set.sdf "/reference/GRCh38_no_alt_analysis_set.fasta"

FILE="reference/trio.ped"
cat <<EOM >$FILE
#PED format pedigree
#
#fam-id/ind-id/pat-id/mat-id: 0=unknown
#sex: 1=male; 2=female; 0=unknown
#phenotype: -9=missing, 0=missing; 1=unaffected; 2=affected
#
#fam-id ind-id pat-id mat-id sex phen
1 HG002 HG003 HG004 1 0
1 HG003 0 0 1 0
1 HG004 0 0 2 0
EOM

sudo docker run \
-v "${PWD}/input":"/input" \
-v "${PWD}/reference":"/reference" \
-v "${PWD}/output":"/output" \
realtimegenomics/rtg-tools mendelian \
-i "/output/HG002_trio_merged.vcf.gz" \
-o "/output/HG002_trio_annotated.output.vcf.gz" \
--pedigree=/reference/trio.ped \
-t /reference/GRCh38_no_alt_analysis_set.sdf \
| tee output/deepvariant.input_rtg_output.txt
```

As a result we should get the following output:

```bash
Checking: /output/HG002_trio_merged.vcf.gz
Family: [HG003 + HG004] -> [HG002]
95 non-pass records were skipped
Concordance HG002: F:137908/139703 (98.72%)  M:137988/139909 (98.63%)  F+M:134596/137968 (97.56%)
Sample HG002 has less than 99.0 concordance with both parents. Check for incorrect pedigree or sample mislabelling.
0/146013 (0.00%) records did not conform to expected call ploidy
143704/146013 (98.42%) records were variant in at least 1 family member and checked for Mendelian constraints
5066/143704 (3.53%) records had indeterminate consistency status due to incomplete calls
3886/143704 (2.70%) records contained a violation of Mendelian constraints
```

### Perform analysis with hap.py against 4.2.1 truth set

```bash
mkdir -p happy

sudo docker pull jmcdani20/hap.py:v0.3.12

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG002.output.vcf.gz \
  -f /benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/HG002.output \
  --engine=vcfeval \
  --pass-only \
  -l chr20

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG003.output.vcf.gz \
  -f /benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/HG003.output \
  --engine=vcfeval \
  --pass-only \
  -l chr20

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /output/HG004.output.vcf.gz \
  -f /benchmark/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/HG004.output \
  --engine=vcfeval \
  --pass-only \
  -l chr20
```

```
Benchmarking Summary for HG002:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        11256     11208        48        21239        13       9586      7      4       0.995736          0.998884        0.451340         0.997308                     NaN                     NaN                   1.561710                   2.047281
INDEL   PASS        11256     11208        48        21239        13       9586      7      4       0.995736          0.998884        0.451340         0.997308                     NaN                     NaN                   1.561710                   2.047281
  SNP    ALL        71333     71087       246        88976        42      17795      5      4       0.996551          0.999410        0.199998         0.997979                2.314904                2.029984                   1.715978                   1.716560
  SNP   PASS        71333     71087       246        88976        42      17795      5      4       0.996551          0.999410        0.199998         0.997979                2.314904                2.029984                   1.715978                   1.716560

Benchmarking Summary for HG003:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10584        44        21028        20       9969     13      6       0.995860          0.998192        0.474082         0.997024                     NaN                     NaN                   1.748961                   2.197401
INDEL   PASS        10628     10584        44        21028        20       9969     13      6       0.995860          0.998192        0.474082         0.997024                     NaN                     NaN                   1.748961                   2.197401
  SNP    ALL        70166     69975       191        85299        55      15231     15      4       0.997278          0.999215        0.178560         0.998246                2.296566                2.064978                   1.883951                   1.845348
  SNP   PASS        70166     69975       191        85299        55      15231     15      4       0.997278          0.999215        0.178560         0.998246                2.296566                2.064978                   1.883951                   1.845348

Benchmarking Summary for HG004:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        11000     10945        55        21426        27       9969     22      4       0.995000          0.997643        0.465276         0.996320                     NaN                     NaN                   1.792709                   2.279678
INDEL   PASS        11000     10945        55        21426        27       9969     22      4       0.995000          0.997643        0.465276         0.996320                     NaN                     NaN                   1.792709                   2.279678
  SNP    ALL        71659     71446       213        86406        52      14858      9      4       0.997028          0.999273        0.171956         0.998149                2.310073                2.064306                   1.878340                   1.735500
  SNP   PASS        71659     71446       213        86406        52      14858      9      4       0.997028          0.999273        0.171956         0.998149                2.310073                2.064306                   1.878340                   1.735500
```
