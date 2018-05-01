# DeepVariant exome case study

Similar to the [case study on whole genome sequencing
data](deepvariant-case-study.md), in this study we describe applying DeepVariant
to a real exome sample using a single machine.

NOTE: This case study demonstrates an example of how to run DeepVariant
end-to-end on one machine. This might not be the fastest or cheapest
configuration for your needs. For more scalable execution of DeepVariant see the
[Docker-based exome pipeline](https://cloud.google.com/genomics/deepvariant)
created for Google Cloud Platform.

## Request a machine

Any sufficiently capable machine will do. For this case study, we used a 64-core
non-preemptible instance with 128GiB and no GPU.

## Preliminaries

Set a number of shell variables, to make what follows easier to read.

```bash
BASE="${HOME}/exome-case-study"
BUCKET="gs://deepvariant"
BIN_VERSION="0.6.1"
MODEL_VERSION="0.6.0"
MODEL_CL="191676894"

# Note that we don't specify the CL number for the binary, only the bin version.
BIN_BUCKET="${BUCKET}/binaries/DeepVariant/${BIN_VERSION}/DeepVariant-${BIN_VERSION}+cl-*"
MODEL_BUCKET="${BUCKET}/models/DeepVariant/${MODEL_VERSION}/DeepVariant-inception_v3-${MODEL_VERSION}+cl-${MODEL_CL}.data-wes_standard"
DATA_BUCKET="${BUCKET}/exome-case-study-testdata"

INPUT_DIR="${BASE}/input"
BIN_DIR="${INPUT_DIR}/bin"
MODELS_DIR="${INPUT_DIR}/models"
MODEL="${MODELS_DIR}/model.ckpt"
DATA_DIR="${INPUT_DIR}/data"
REF="${DATA_DIR}/hs37d5.fa.gz"
BAM="${DATA_DIR}/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
TRUTH_VCF="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/HG002.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="${OUTPUT_DIR}/HG002.output.vcf.gz"
OUTPUT_GVCF="${OUTPUT_DIR}/HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

CAPTURE_BED="${DATA_DIR}/agilent_sureselect_human_all_exon_v5_b37_targets.bed"
```

## Create local directory structure

```bash
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BIN_DIR}"
mkdir -p "${DATA_DIR}"
mkdir -p "${MODELS_DIR}"
mkdir -p "${LOG_DIR}"
```

## Download extra packages

There are some extra programs we will need.

We are going to use [GNU Parallel](https://www.gnu.org/software/parallel/) to
run `make_examples`. We are going to install `samtools` and `docker.io` to help
do some analysis at the end.

```bash
sudo apt-get -y install parallel
sudo apt-get -y install samtools
sudo apt-get -y install docker.io
```

## Download binaries, models, and test data

### Binaries

Copy our binaries from the cloud bucket.

```bash
time gsutil -m cp -r "${BIN_BUCKET}/*" "${BIN_DIR}"
chmod a+x "${BIN_DIR}"/*
```

This step should be very fast - it took us about 6 seconds when we tested.

Now, we need to install all prerequisites on the machine. Run this command:

```bash
cd "${BIN_DIR}"; time bash run-prereq.sh; cd -
```

In our test run it took about 1 min.

### Models

Copy the model files to your local disk.

```bash
time gsutil -m cp -r "${MODEL_BUCKET}/*" "${MODELS_DIR}"
```

This step should be really fast. It took us about 5 seconds.

### Test data

Copy the input files you need to your local disk from our gs:// bucket.

The original source of these files are:

##### BAM file:

`151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam`

Downloaded from
[https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_OsloUniversityHospital_IlluminaExome_bwamem_GRCh37_11252015](https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_OsloUniversityHospital_IlluminaExome_bwamem_GRCh37_11252015)

##### FASTA

Same as described in the [case study for whole genome
data](deepvariant-case-study.md#test_data)

##### Truth VCF and BED

`HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_*`
are from NIST, as part of the [Genomes in a Bottle
project](http://jimb.stanford.edu/giab/). They are downloaded from
[ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/)

##### Capture target BED file

According to the paper "[Extensive sequencing of seven human genomes to
characterize benchmark reference
materials](https://www.nature.com/articles/sdata201625)", the HG002 exome was
generated with Agilent SureSelect. In this case study we'll use the SureSelect
v5 BED (`agilent_sureselect_human_all_exon_v5_b37_targets.bed`) and intersect it
with the GIAB confident regions for evaluation.

##### Copy the data

You can simply run the command below to get all the data you need for this case
study.

```bash
time gsutil -m cp -r "${DATA_BUCKET}/*" "${DATA_DIR}"
```

It took us a few minuntes to copy the files.

## Run `make_examples`

In this step, we used the `--regions` flag to constrain the regions we processed
to the capture region BED file:

```bash
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    python "${BIN_DIR}"/make_examples.zip \
      --mode calling \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${EXAMPLES}" \
      --regions "${CAPTURE_BED}" \
      --gvcf "${GVCF_TFRECORDS}" \
      --task {}
) >"${LOG_DIR}/make_examples.log" 2>&1
```

Timing information is included in [a later
section](#resources_used_by_each_step).

## Run `call_variants`

There are different ways to run `call_variants`. In this case study, we ran just
one `call_variants` job. Here's the command that we used:

```bash
( time python "${BIN_DIR}"/call_variants.zip \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --examples "${EXAMPLES}" \
    --checkpoint "${MODEL}" \
    --batch_size 32
) >"${LOG_DIR}/call_variants.log" 2>&1
```

More discussion can be found in the [call_variants section in the case
study](deepvariant-case-study.md#run_call_variants).

## Run `postprocess_variants`

```bash
( time python "${BIN_DIR}"/postprocess_variants.zip \
    --ref "${REF}" \
    --infile "${CALL_VARIANTS_OUTPUT}" \
    --outfile "${OUTPUT_VCF}" \
    --nonvariant_site_tfrecord_path "${GVCF_TFRECORDS}" \
    --gvcf_outfile "${OUTPUT_GVCF}"
) >"${LOG_DIR}/postprocess_variants.withGVCF.log" 2>&1
```

The last two flags are optional, only if gVCF outputs are desired. More
discussion can be found in the [postprocess_variants section in the case
study](deepvariant-case-study.md#run_postprocess_variants).

## Resources used by each step

Step                               | wall time
---------------------------------- | ---------
`make_examples`                    | 65m 45s
`call_variants`                    | 6m 21s
`postprocess_variants` (no gVCF)   | 0m 13s
`postprocess_variants` (with gVCF) | 1m 24s
total time (single machine)        | ~ 1h 13m

## Variant call quality

Here we use the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This serves as a check
to ensure the three DeepVariant commands ran correctly and produced high-quality
results.

To set up:

```bash
UNCOMPRESSED_REF="${OUTPUT_DIR}/hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory and index it.
zcat <"${REF}" >"${UNCOMPRESSED_REF}"
samtools faidx "${UNCOMPRESSED_REF}"

sudo docker pull pkrusche/hap.py
```

We evaluate against the capture region:

```bash
sudo docker run -it \
-v "${DATA_DIR}:${DATA_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${OUTPUT_VCF}" \
  -f "${TRUTH_BED}" \
  -T "${CAPTURE_BED}" \
  -r "${UNCOMPRESSED_REF}" \
  -o "${OUTPUT_DIR}/happy.output" \
  --engine=vcfeval
```

Here are the results:

Type  | # FN | # FP | Recall   | Precision | F1\_Score
----- | ---- | ---- | -------- | --------- | ---------
INDEL | 127  | 35   | 0.951154 | 0.986193  | 0.968357
SNP   | 43   | 26   | 0.998725 | 0.999229  | 0.998977

## Separate models for calling whole genome and exome data

In the DeepVariant 0.5.\* release, we recommend a separate model for calling
exome sequencing data. Here is how the exome model is trained: we used a WGS
model as the starting checkpoint (instead of an ImageNet one), and trained only
on examples created from exome data.
