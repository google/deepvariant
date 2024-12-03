# DeepVariant with PacBio HiFi data

In this case study we describe applying DeepVariant to PacBio HiFi reads to call
variants. We will call small variants from a publicly available whole genome
HiFi dataset from PacBio.

### Updated dataset in release 1.8.0

In release 1.8.0, we have updated the PacBio test data from HG003 Sequel-II to
latest Revio with SPRQ chemistry data to showcase performance on the updated
platform and chemistry. The full bam data is available [here](https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/GIAB_trio/HG003/analysis/GRCh38.m84039_241002_000337_s3.hifi_reads.bc2020.bam).

The dataset used in this case-study has following attributes:

```bash
Sample: HG003
Region: Chr20
Chemistry: REVIO SPRQ
Coverage: 32x
```

## Prepare environment

In this case-study, we will use [Docker](https://docs.docker.com/get-docker/) to
run DeepVariant for variant calling and
[hap.py](https://github.com/illumina/hap.py) for benchmarking.

If you want to run on GPU machines, or use `Singularity` instead of `Docker`,
please follow [Quick Start](deepvariant-quick-start.md) documentation.

### Create input and output directory structures and download inputs

```bash
BASE="${HOME}/pacbio-case-study"

# Set up input and output directory data
INPUT_DIR="${BASE}/input/data"
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Download reference to input directory
FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > ${INPUT_DIR}/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > ${INPUT_DIR}/GRCh38_no_alt_analysis_set.fasta.fai

HTTPDIR=https://storage.googleapis.com/deepvariant/pacbio-case-study-testdata
curl ${HTTPDIR}/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam > ${INPUT_DIR}/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam
curl ${HTTPDIR}/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam.bai > ${INPUT_DIR}/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam.bai

# Set up input variables
REF="GRCh38_no_alt_analysis_set.fasta"
BAM="HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam"
THREADS=$(nproc)
REGION="chr20"

# Set up output variable
OUTPUT_VCF="HG003_PACBIO_SPRQ_GRCh38.chr20.output.vcf.gz"
OUTPUT_GVCF="HG003_PACBIO_SPRQ_GRCh38.chr20.output.g.vcf.gz"
INTERMEDIATE_DIRECTORY="intermediate_results_dir"

mkdir -p "${OUTPUT_DIR}/${INTERMEDIATE_DIRECTORY}"
```

## Run DeepVariant

We will run DeepVariant from docker using the `run_deepvariant` script.

```bash
BIN_VERSION="1.8.0"

sudo docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref "${INPUT_DIR}/${REF}" \
  --reads "${INPUT_DIR}/${BAM}" \
  --output_vcf "${OUTPUT_DIR}/${OUTPUT_VCF}" \
  --output_gvcf "${OUTPUT_DIR}/${OUTPUT_GVCF}" \
  --num_shards "${THREADS}" \
  --regions "${REGION}" \
  --intermediate_results_dir "${OUTPUT_DIR}/${INTERMEDIATE_DIRECTORY}"
```

By specifying `--model_type PACBIO`, you'll be using a model that is best
suited for PacBio data.

NOTE: If you want to run each of the steps separately, add `--dry_run=true` to
the command above to figure out what flags you need in each step. Based on the
different model types, different flags are needed in the `make_examples` step.

`--intermediate_results_dir` flag is optional. By specifying it, the
intermediate outputs of `make_examples` and `call_variants` stages can be found
in the directory. After the command, you can find these files in the directory:

```
call_variants_output.tfrecord.gz
gvcf.tfrecord-?????-of-?????.gz
make_examples.tfrecord-?????-of-?????.gz
```

## Benchmark HG003 chr20 output from DeepVariant

We will use Genome-in-a-Bottle (GIAB) dataset to evaluate the performance of
DeepVariant.

### Download Genome in a Bottle Benchmarks

We will benchmark our variant calls against v4.2.1 of the Genome in a Bottle
small variant benchmarks for HG003.

```bash
FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38

curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > ${INPUT_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > ${INPUT_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > ${INPUT_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
```

```bash
sudo docker pull jmcdani20/hap.py:v0.3.12

sudo docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  "${INPUT_DIR}/${TRUTH_VCF}" \
  "${OUTPUT_DIR}/${OUTPUT_VCF}" \
  -f "${INPUT_DIR}/${TRUTH_BED}" \
  -r "${INPUT_DIR}/${REF}" \
  -o "${OUTPUT_DIR}/hg003.pacbio.chr20.happy.output" \
  --engine=vcfeval \
  --pass-only \
  -l "${REGION}"
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628     10543        85        22403        74      11375     40     29       0.992002          0.993290        0.507744         0.992646                     NaN                     NaN                   1.748961                   2.138647
INDEL   PASS        10628     10543        85        22403        74      11375     40     29       0.992002          0.993290        0.507744         0.992646                     NaN                     NaN                   1.748961                   2.138647
  SNP    ALL        70166     70101        65       105602        71      35342     12     12       0.999074          0.998989        0.334672         0.999032                2.296566                1.713281                   1.883951                   1.503192
  SNP   PASS        70166     70101        65       105602        71      35342     12     12       0.999074          0.998989        0.334672         0.999032                2.296566                1.713281                   1.883951                   1.503192
```
