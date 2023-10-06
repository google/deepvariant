# Calling variants in non-autosomal contigs

For details about the support for haploid contigs, please read
[DeepVariant haploid support](deepvariant-haploid-support.md).

In this case study, we describe how to call variants in non-autosomal regions
like X, Y chromosomes. Then we assess the quality of the DeepVariant variant
calls with `hap.py`.

The dataset used in this case-study has following attributes:

```bash
Sample: HG002
Region: ChrX, ChrY
Platform: PacBio
Sample Karyotype: X, Y
```

## Prepare environment

In this case study, we will use [Docker](https://docs.docker.com/get-docker/) to
run DeepVariant for variant calling and
[hap.py](https://github.com/illumina/hap.py) for benchmarking.

If you want to run on GPU machines, or use `Singularity` instead of `Docker`,
please follow [Quick Start](deepvariant-quick-start.md) documentation.

### Create input and output directory structures and download inputs

```bash
BASE="${HOME}/XY-walkthrough"

# Set up input and output directory data
INPUT_DIR="${BASE}/input"
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/data"

# Download reference to input directory
FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > ${INPUT_DIR}/GRCh38_no_alt_analysis_set.fasta
curl ${FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > ${INPUT_DIR}/GRCh38_no_alt_analysis_set.fasta.fai

# Download bam file to input directory
HTTPDIR=https://storage.googleapis.com/deepvariant/xy-case-study-testdata
curl ${HTTPDIR}/HG002.pfda_challenge.grch38.chrXY.bam > ${INPUT_DIR}/HG002.pfda_challenge.grch38.chrXY.bam
curl ${HTTPDIR}/HG002.pfda_challenge.grch38.chrXY.bam.bai > ${INPUT_DIR}/HG002.pfda_challenge.grch38.chrXY.bam.bai

HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata
curl ${HTTPDIR}/GRCh38_PAR.bed > ${INPUT_DIR}/GRCh38_PAR.bed

# Set up input variables
REF="GRCh38_no_alt_analysis_set.fasta"
BAM="HG002.pfda_challenge.grch38.chrXY.bam"
THREADS=$(nproc)
REGION="chrX chrY"
HAPLOID_CONTIGS="chrX,chrY"
PAR_BED="GRCh38_PAR.bed"

# Set up output variable
OUTPUT_VCF="HG002_pacbio_hifi.chrXY.output.vcf.gz"
OUTPUT_GVCF="HG002_pacbio_hifi.chrXY.output.g.vcf.gz"
INTERMEDIATE_DIRECTORY="intermediate_results_dir"

mkdir -p "${OUTPUT_DIR}/${INTERMEDIATE_DIRECTORY}"
```

## Run DeepVariant

We will run DeepVariant from docker using the `run_deepvariant` script.

```bash
BIN_VERSION="1.6.0"

docker pull google/deepvariant:"${BIN_VERSION}"

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
--haploid_contigs "${HAPLOID_CONTIGS}" \
--par_regions_bed "${INPUT_DIR}/${PAR_BED}" \
--regions "${REGION}" \
--intermediate_results_dir "${OUTPUT_DIR}/${INTERMEDIATE_DIRECTORY}"
```

## Benchmark X, Y outputs from DeepVariant

We will use Genome-in-a-Bottle (GIAB) dataset to evaluate the performance of
DeepVariant.

### Download Genome in a Bottle Benchmarks

We will benchmark our variant calls against v1.0 of the Genome in a Bottle
small variant benchmarks for HG002_chrXY.

```bash
FTPDIR=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/chrXY_v1.0/GRCh38/SmallVariant

curl ${FTPDIR}/HG002_GRCh38_chrXY_smallvar_v1.0.bed > ${INPUT_DIR}/HG002_GRCh38_chrXY_smallvar_v1.0.bed
curl ${FTPDIR}/HG002_GRCh38_chrXY_smallvar_v1.0.vcf.gz > ${INPUT_DIR}/HG002_GRCh38_chrXY_smallvar_v1.0.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_chrXY_smallvar_v1.0.vcf.gz.tbi > ${INPUT_DIR}/HG002_GRCh38_chrXY_smallvar_v1.0.vcf.gz.tbi

TRUTH_VCF="HG002_GRCh38_chrXY_smallvar_v1.0.vcf.gz"
TRUTH_BED="HG002_GRCh38_chrXY_smallvar_v1.0.bed"
```

```bash
sudo docker pull jmcdani20/hap.py:v0.3.12

REGION="chrX,chrY"
sudo docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  -v "${PWD}/happy:/happy" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  "${INPUT_DIR}/${TRUTH_VCF}" \
  "${OUTPUT_DIR}/${OUTPUT_VCF}" \
  -f "${INPUT_DIR}/${TRUTH_BED}" \
  -r "${INPUT_DIR}/${REF}" \
  -o "${OUTPUT_DIR}/hg002.chrXY.happy.output" \
  --engine=vcfeval \
  --pass-only \
  -l "${REGION}"
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        24273     23511       762        30833       517       6378     35    367       0.968607          0.978859        0.206856         0.973706                     NaN                     NaN                   1.559454                   0.071354
INDEL   PASS        24273     23511       762        30833       517       6378     35    367       0.968607          0.978859        0.206856         0.973706                     NaN                     NaN                   1.559454                   0.071354
  SNP    ALL        87443     86921       522       116830       551      29511     13    112       0.994030          0.993690        0.252598         0.993860                1.937122                1.649898                   1.825434                   0.052139
  SNP   PASS        87443     86921       522       116830       551      29511     13    112       0.994030          0.993690        0.252598         0.993860                1.937122                1.649898                   1.825434                   0.052139
```
