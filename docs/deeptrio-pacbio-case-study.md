# Using DeepTrio for small variant calling from the trio sequenced with PacBio HiFi

In this case study, we describe applying [DeepTrio](deeptrio-details.md) to a
real PacBio WGS trio. Then we assess the quality of the DeepTrio variant calls
with `hap.py`. In addition we evaluate a Mendelian violation rate for a merged
VCF.

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

We'll use HG002, HG003, HG004 PacBio HiFi WGS reads publicly available from the
[PrecisionFDA Truth v2 Challenge](https://precision.fda.gov/challenges/10).
These reads have been aligned to the GRCh38_no_alt_analysis reference using
[pbmm2](https://github.com/PacificBiosciences/pbmm2).

```bash
mkdir -p input
HTTPDIR=https://storage.googleapis.com/deepvariant/pacbio-case-study-testdata

curl ${HTTPDIR}/HG002.pfda_challenge.grch38.phased.chr20.bam > input/HG002.pfda_challenge.grch38.phased.chr20.bam
curl ${HTTPDIR}/HG002.pfda_challenge.grch38.phased.chr20.bam.bai > input/HG002.pfda_challenge.grch38.phased.chr20.bam.bai

curl ${HTTPDIR}/HG003.pfda_challenge.grch38.phased.chr20.bam > input/HG003.pfda_challenge.grch38.phased.chr20.bam
curl ${HTTPDIR}/HG003.pfda_challenge.grch38.phased.chr20.bam.bai > input/HG003.pfda_challenge.grch38.phased.chr20.bam.bai

curl ${HTTPDIR}/HG004.pfda_challenge.grch38.phased.chr20.bam > input/HG004.pfda_challenge.grch38.phased.chr20.bam
curl ${HTTPDIR}/HG004.pfda_challenge.grch38.phased.chr20.bam.bai > input/HG004.pfda_challenge.grch38.phased.chr20.bam.bai
```

## Running DeepTrio with one command

DeepTrio pipeline consists of 4 steps: `make_examples`, `call_variants`,
`postprocess_variants` and `GLnexus merge`. It is possible to run the first
three steps with one command using the `run_deeptrio` script. GLnexus
is run as a separate command.

### Running on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION="1.4.0"

sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:deeptrio-"${BIN_VERSION}"

time sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:deeptrio-"${BIN_VERSION}" \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type PACBIO \
  --ref /reference/GRCh38_no_alt_analysis_set.fasta \
  --reads_child /input/HG002.pfda_challenge.grch38.phased.chr20.bam \
  --reads_parent1 /input/HG003.pfda_challenge.grch38.phased.chr20.bam \
  --reads_parent2 /input/HG004.pfda_challenge.grch38.phased.chr20.bam \
  --output_vcf_child /output/HG002.output.vcf.gz \
  --output_vcf_parent1 /output/HG003.output.vcf.gz \
  --output_vcf_parent2 /output/HG004.output.vcf.gz \
  --sample_name_child 'HG002' \
  --sample_name_parent1 'HG003' \
  --sample_name_parent2 'HG004' \
  --num_shards $(nproc) \
  --intermediate_results_dir /output/intermediate_results_dir \
  --output_gvcf_child /output/HG002.g.vcf.gz \
  --output_gvcf_parent1 /output/HG003.g.vcf.gz \
  --output_gvcf_parent2 /output/HG004.g.vcf.gz \
  --regions chr20 \
  --use_hp_information
```

The `--use_hp_information` arg makes use of a phased reads, thus allowing a
further improvement of the accuracy. You can use tools like whatshap to phase.

By specifying `--model_type PACBIO`, you'll be using a model that is best suited
for PacBio HiFi Whole Genome Sequencing data.

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
[Quick Start](deepvariant-quick-start.md) or
[DeepVariant PacBio case study](deepvariant-pacbio-model-case-study.md).

## Merge VCFs using GLnexus

At this step we take all 3 VCFs generated in the previous step and merge them
using GLnexus.

```bash
sudo docker pull quay.io/mlin/glnexus:v1.2.7

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

### Calculate Mendelian Violation rate

```bash
sudo docker pull realtimegenomics/rtg-tools

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
175 non-pass records were skipped
Concordance HG002: F:152203/153489 (99.16%)  M:153447/154331 (99.43%)  F+M:148795/151183 (98.42%)
Sample HG002 has less than 99.0 concordance with both parents. Check for incorrect pedigree or sample mislabelling.
0/165025 (0.00%) records did not conform to expected call ploidy
159539/165025 (96.68%) records were variant in at least 1 family member and checked for Mendelian constraints
7726/159539 (4.84%) records had indeterminate consistency status due to incomplete calls
2678/159539 (1.68%) records contained a violation of Mendelian constraints
```

### Benchmark variant calls against 4.2.1 truth set with hap.py

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
INDEL    ALL        11256     11233        23        23105        62      11332     15     43       0.997957          0.994734        0.490457         0.996343                     NaN                     NaN                   1.561710                   2.316632
INDEL   PASS        11256     11233        23        23105        62      11332     15     43       0.997957          0.994734        0.490457         0.996343                     NaN                     NaN                   1.561710                   2.316632
  SNP    ALL        71333     71272        61       100735        22      29369     15      7       0.999145          0.999692        0.291547         0.999418                2.314904                1.932428                   1.715978                   2.093416
  SNP   PASS        71333     71272        61       100735        22      29369     15      7       0.999145          0.999692        0.291547         0.999418                2.314904                1.932428                   1.715978                   2.093416

Benchmarking Summary for HG003:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10595        33        23047        54      11901     22     30       0.996895          0.995155        0.516380         0.996024                     NaN                     NaN                   1.748961                   2.557273
INDEL   PASS        10628     10595        33        23047        54      11901     22     30       0.996895          0.995155        0.516380         0.996024                     NaN                     NaN                   1.748961                   2.557273
  SNP    ALL        70166     70144        22        99173        19      28954      9      2       0.999686          0.999729        0.291954         0.999708                2.296566                1.923292                   1.883951                   2.164234
  SNP   PASS        70166     70144        22        99173        19      28954      9      2       0.999686          0.999729        0.291954         0.999708                2.296566                1.923292                   1.883951                   2.164234

Benchmarking Summary for HG004:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        11000     10967        33        23397        53      11867     19     33       0.997000          0.995403        0.507202         0.996201                     NaN                     NaN                   1.792709                   2.553655
INDEL   PASS        11000     10967        33        23397        53      11867     19     33       0.997000          0.995403        0.507202         0.996201                     NaN                     NaN                   1.792709                   2.553655
  SNP    ALL        71659     71594        65       100530        38      28814      9     11       0.999093          0.999470        0.286621         0.999281                2.310073                1.949802                   1.878340                   2.053544
  SNP   PASS        71659     71594        65       100530        38      28814      9     11       0.999093          0.999470        0.286621         0.999281                2.310073                1.949802                   1.878340                   2.053544
```
