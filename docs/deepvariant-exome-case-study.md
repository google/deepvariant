# DeepVariant exome case study

Similar to the [case study on whole genome sequencing
data](deepvariant-case-study.md), in this study we describe applying
DeepVariant to a real exome sample.

## Request a machine

For this case study, we use a 64-core DeepVariant non-preemptible instance in
the "us-west1-b" zone with no GPU. From our local command line, we do:

```bash
gcloud beta compute instances create "${USER}-deepvariant-exome-casestudy"  \
--scopes "compute-rw,storage-full,cloud-platform" \
--image-family "ubuntu-1604-lts" --image-project "ubuntu-os-cloud" \
--machine-type "custom-64-131072" \
--boot-disk-size "50" --boot-disk-type "pd-ssd" \
--boot-disk-device-name "deepvariant-exome-casestudy" \
--zone "us-west1-b"
```

The `custom-64-131072` machine type gives you 64 vCPU, 128.0 GiB.

Then connect to your instance via SSH:

```bash
gcloud compute ssh --zone "us-west1-b" "${USER}-deepvariant-exome-casestudy"
```

## Preliminaries

Set a number of shell variables, to make what follows easier.

```bash
BASE="${HOME}/exome-case-study"
BUCKET="gs://deepvariant"
BIN_VERSION="0.4.0"
MODEL_VERSION="0.4.0"
MODEL_CL="174375304"

# Note that we don't specify the CL number for the binary, only the bin version.
BIN_BUCKET="${BUCKET}/binaries/DeepVariant/${BIN_VERSION}/DeepVariant-${BIN_VERSION}+cl-*"
MODEL_BUCKET="${BUCKET}/models/DeepVariant/${MODEL_VERSION}/DeepVariant-inception_v3-${MODEL_VERSION}+cl-${MODEL_CL}.data-wgs_standard"
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
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="${OUTPUT_DIR}/HG002.output.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

REFSEQ_BED="${DATA_DIR}/refseq.coding_exons.b37.bed"
EXTENDED_REFSEQ_BED="${DATA_DIR}/refseq.coding_exons.b37.extended50.bed"
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
run `make_examples`.
We are going to install `samtools` and `docker.io` to help do some
analysis at the end.

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

##### RefSeq target BED file

We prepared two different BED files:

1.  `refseq.coding_exons.b37.bed`

    *   Download hg19 refseq.coding_exones.hg19.bed from
        [http://genome.ucsc.edu/cgi-bin/hgTables](http://genome.ucsc.edu/cgi-bin/hgTables)
    *   Edit file to remove genes on unk/alt genes
    *   Then replace-regex '^chr' => '' throughout the file.
    *   Make sure the chromosomes are sorted the same way as `hs37d5.fa`.
    *   Save as refseq.coding_exons.b37.bed

    Coverage is 33,889,421 bases.

1.  `refseq.coding_exons.b37.extended50.bed`

    Coverage is 53,393,967 bases.

##### Copy the data

You can simply run the command below to get all the data you need for this case
study.

```bash
time gsutil -m cp -r "${DATA_BUCKET}/*" "${DATA_DIR}"
```

It took us a few minuntes to copy the files.

## Run `make_examples`

In this step, we used the `--regions` flag to constrain the regions
we processed to the extended RefSeq BED file:

```bash
( time seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog "${LOG_DIR}/log" --res "${LOG_DIR}" \
    python "${BIN_DIR}"/make_examples.zip \
      --mode calling \
      --ref "${REF}" \
      --reads "${BAM}" \
      --examples "${EXAMPLES}" \
      --regions "${EXTENDED_REFSEQ_BED}" \
      --task {}
) >"${LOG_DIR}/make_examples.log" 2>&1
```

Timing information is included in [a later
section](#resources_used_by_each_step).

## Run `call_variants`

Follow the same instructions (reuse the same commands) in the [call_variants
section in the case study](deepvariant-case-study.md#run_call_variants).

## Run `postprocess_variants`

Follow the same instructions (reuse the same commands) in the
[postprocess_variants section in the case
study](deepvariant-case-study.md#run_postprocess_variants).

## Resources used by each step

Step                   | wall time
---------------------- | ---------
`make_examples`        | 64m 7s
`call_variants`        | 5m 43s
`postprocess_variants` | 0m 18s
total time             | ~ 1h 10m

## Variant call quality

Here we use the `hap.py` ([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting vcf file. This
serves as a check to ensure the three DeepVariant commands ran correctly and
produced high-quality results.

To set up:

```bash
UNCOMPRESSED_REF="${OUTPUT_DIR}/hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory and index it.
zcat <"${REF}" >"${UNCOMPRESSED_REF}"
samtools faidx "${UNCOMPRESSED_REF}"

sudo docker pull pkrusche/hap.py
```

First, we evaluate against just the RefSeq region:

```bash
sudo docker run -it \
-v "${DATA_DIR}:${DATA_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${OUTPUT_VCF}" \
  -f "${TRUTH_BED}" \
  -T "${REFSEQ_BED}" \
  -r "${UNCOMPRESSED_REF}" \
  -o "${OUTPUT_DIR}/happy.output"
```

Then, we also evaluate against the extended RefSeq region:

```bash
sudo docker run -it \
-v "${DATA_DIR}:${DATA_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${OUTPUT_VCF}" \
  -f "${TRUTH_BED}" \
  -T "${EXTENDED_REFSEQ_BED}" \
  -r "${UNCOMPRESSED_REF}" \
  -o "${OUTPUT_DIR}/extended.happy.output"
```

Putting the quality results in one table, we have:

BED                 | Type  | Recall   | Precision | F1_Score
------------------- | ----- | -------- | --------- | --------
REFSEQ_BED          | INDEL | 0.969194 | 0.971496  | 0.970344
REFSEQ_BED          | SNP   | 0.995194 | 0.997902  | 0.996546
EXTENDED_REFSEQ_BED | INDEL | 0.912479 | 0.921393  | 0.916914
EXTENDED_REFSEQ_BED | SNP   | 0.994076 | 0.997693  | 0.995881

## Limitations and Future Work

The current released model is trained on whole genome sequencing data. Based on
the evaluation on RefSeq, the F1 scores are reasonable even though we didn't
include exome data in training.

However, from our experience in the [PrecisionFDA Hidden Treasures
Challenge](https://precision.fda.gov/challenges/1/view/results) we know that we
can train a better model for calling exome data if we include exome in our
training data.

We are actively working on extending the training set, so we can release a
better model for calling on exome data in the future.
