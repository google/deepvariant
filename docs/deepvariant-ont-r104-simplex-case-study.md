# DeepVariant with Oxford Nanopore R10.4.1 Simplex reads

In this case study, we describe applying DeepVariant to Oxford Nanopore R10.4.1
simplex reads. Then we assess the quality of the DeepVariant variant calls with
`hap.py`.

To make it faster to go over this case study, we run only on chromosome 20.

The dataset used in this case-study has following attributes:

```bash
Sample: HG003
Region: Chr20
Chemistry: ONT R10.4.1
Coverage: 80x
```

**Model note:**

*   The model is trained with Guppy 6+ "SUP" Simplex and Dorado v0.1.1 Duplex
    reads.

*   The model is trained on both Ultra-long and sheared reads with varying read
    N50 and coverage.

## Prepare environment

In this case-study, we will use [Docker](https://docs.docker.com/get-docker/) to
run DeepVariant for variant calling and
[hap.py](https://github.com/illumina/hap.py) for benchmarking.

If you want to run on GPU machines, or use `Singularity` instead of `Docker`,
please follow [Quick Start](deepvariant-quick-start.md) documentation.

### Create input and output directory structures and download inputs

```bash
BASE="${HOME}/ont-case-study"

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

# Download HG003 Ultra-long chr20 bam file to input directory
HTTPDIR=https://storage.googleapis.com/deepvariant/ont-case-study-testdata
curl ${HTTPDIR}/HG003_R104_sup_merged.80x.chr20.bam > ${INPUT_DIR}/HG003_R104_sup_merged.80x.chr20.bam
curl ${HTTPDIR}/HG003_R104_sup_merged.80x.chr20.bam.bai > ${INPUT_DIR}/HG003_R104_sup_merged.80x.chr20.bam.bai

# Set up input variables
REF="GRCh38_no_alt_analysis_set.fasta"
BAM="HG003_R104_sup_merged.80x.chr20.bam"
THREADS=$(nproc)
REGION="chr20"

# Set up output variable
OUTPUT_VCF="HG003_UL_R1041_Guppy6_sup_2_GRCh38.chr20.output.vcf.gz"
OUTPUT_GVCF="HG003_UL_R1041_Guppy6_sup_2_GRCh38.output.g.vcf.gz"
INTERMEDIATE_DIRECTORY="intermediate_results_dir"

mkdir -p "${OUTPUT_DIR}/${INTERMEDIATE_DIRECTORY}"
```

## Run DeepVariant

We will run DeepVariant from docker using the `run_deepvariant` script.

```bash
BIN_VERSION="1.6.0"

sudo docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type ONT_R104 \
  --ref "${INPUT_DIR}/${REF}" \
  --reads "${INPUT_DIR}/${BAM}" \
  --output_vcf "${OUTPUT_DIR}/${OUTPUT_VCF}" \
  --output_gvcf "${OUTPUT_DIR}/${OUTPUT_GVCF}" \
  --num_shards "${THREADS}" \
  --regions "${REGION}" \
  --intermediate_results_dir "${OUTPUT_DIR}/${INTERMEDIATE_DIRECTORY}"
```

By specifying `--model_type ONT_R104`, you'll be using a model that is best
suited for Oxford Nanopore R10.4.1 chemistry Simplex and Duplex reads.

NOTE: If you want to run each of the steps separately, add `--dry_run=true`
to the command above to figure out what flags you need in each step. Based on
the different model types, different flags are needed in the `make_examples`
step.

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
  -o "${OUTPUT_DIR}/hg003.ul.r104.ont.chr20.happy.output" \
  --engine=vcfeval \
  --pass-only \
  -l "${REGION}"
```

Output:

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL        10628      9165      1463        18460       934       8010    432    339       0.862345          0.910622        0.433911         0.885826                     NaN                     NaN                   1.748961                   2.129371
INDEL   PASS        10628      9165      1463        18460       934       8010    432    339       0.862345          0.910622        0.433911         0.885826                     NaN                     NaN                   1.748961                   2.129371
  SNP    ALL        70166     69925       241        91027       151      20935     82     29       0.996565          0.997846        0.229987         0.997205                2.296566                1.944646                   1.883951                   1.843331
  SNP   PASS        70166     69925       241        91027       151      20935     82     29       0.996565          0.997846        0.229987         0.997205                2.296566                1.944646                   1.883951                   1.843331
```

## Acknowledgement

**For providing analysis results and expertise, we are thankful to:**

*   Karen Miga, Brandy McNulty, Jean Monlong, Benedict Paten from UC Santa Cruz
    Genomics Institute, University of California, Santa Cruz, CA.
*   Miten Jain from Department of Bioengineering, Department of Physics,
    Northeastern University, Boston, MA.
