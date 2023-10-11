# DeepVariant quick start

This is an explanation of how to use DeepVariant.

## Background

To get started, you'll need the DeepVariant programs (and some packages they
depend on), some test data, and of course a place to run them.

We've provided a Docker image, and some test data in a bucket on Google Cloud
Storage. The instructions below show how to download the data through the
corresponding public URLs from these data.

This setup requires a machine with the AVX instruction set. To see if your
machine meets this requirement, you can check the `/proc/cpuinfo` file, which
lists this information under "flags". If you do not have the necessary
instructions, see the next section for more information on how to build your own
Docker image.

### Use Docker to run DeepVariant in one command.

Starting from the 0.8 release, we introduced one convenient command that will
run through all 3 steps that are required to go from a BAM file to the VCF/gVCF
output files. You can still read about the r0.7 approach in
[Quick Start in r0.7].

If you want to compile the DeepVariant binaries for yourself, we also have a
[Dockerfile] that you can use to build your own Docker image. You can read the
[docker build] documentation on how to build.

## Get Docker image, models, and test data

### Get Docker image

```bash
BIN_VERSION="1.6.0"

sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:"${BIN_VERSION}"
```

### Download test data

Before you start running, you need to have the following input files:

1.  A reference genome in [FASTA] format and its corresponding index file
    (.fai).

1.  An aligned reads file in [BAM] format and its corresponding index file
    (.bai). You get this by aligning the reads from a sequencing instrument,
    using an aligner like [BWA] for example.

We've prepared a small test data bundle for use in this quick start guide that
can be downloaded to your instance from the public URLs.

Download the test bundle:

```bash
INPUT_DIR="${PWD}/quickstart-testdata"
DATA_HTTP_DIR="https://storage.googleapis.com/deepvariant/quickstart-testdata"

mkdir -p ${INPUT_DIR}
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/NA12878_S1.chr20.10_10p1mb.bam
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/NA12878_S1.chr20.10_10p1mb.bam.bai
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/test_nist.b37_chr20_100kbp_at_10mb.bed
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.fai
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.gz
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.gz.fai
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.gz.gzi
```

This should create a subdirectory in the current directory containing the actual
data files:

```bash
ls -1 ${INPUT_DIR}
```

outputting:

```
NA12878_S1.chr20.10_10p1mb.bam
NA12878_S1.chr20.10_10p1mb.bam.bai
test_nist.b37_chr20_100kbp_at_10mb.bed
test_nist.b37_chr20_100kbp_at_10mb.vcf.gz
test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi
ucsc.hg19.chr20.unittest.fasta
ucsc.hg19.chr20.unittest.fasta.fai
ucsc.hg19.chr20.unittest.fasta.gz
ucsc.hg19.chr20.unittest.fasta.gz.fai
ucsc.hg19.chr20.unittest.fasta.gz.gzi
```

### Model location (optional)

Starting from r0.8, we put the model files inside the released Docker images.
So there is no need to download model files anymore. If you want to find the
model files of all releases, you can find them in our bucket on the Google Cloud
Storage. You can view them in the browser:
https://console.cloud.google.com/storage/browser/deepvariant/models/DeepVariant

## Run DeepVariant with one command

DeepVariant consists of 3 main binaries: `make_examples`, `call_variants`, and
`postprocess_variants`. To make it easier to run, we create one entrypoint that
can be directly run as a docker command. If you want to see the details, you can
read through [run_deepvariant.py].

```bash
OUTPUT_DIR="${PWD}/quickstart-output"
mkdir -p "${OUTPUT_DIR}"
```

You can run everything with the following command:

```bash
sudo docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}":"/output" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/input/ucsc.hg19.chr20.unittest.fasta \
  --reads=/input/NA12878_S1.chr20.10_10p1mb.bam \
  --regions "chr20:10,000,000-10,010,000" \
  --output_vcf=/output/output.vcf.gz \
  --output_gvcf=/output/output.g.vcf.gz \
  --intermediate_results_dir /output/intermediate_results_dir \
  --num_shards=1
```

NOTE: If you want to look at all the commands being run, you can add
`--dry_run=true` to the command above, which will print out all the commands
but not execute them.

This will generate 5 files and 1 directory in `${OUTPUT_DIR}`:

```bash
ls -1 ${OUTPUT_DIR}
```

outputting:

```
intermediate_results_dir
output.g.vcf.gz
output.g.vcf.gz.tbi
output.vcf.gz
output.vcf.gz.tbi
output.visual_report.html
```

The directory "intermediate_results_dir" exists because
`--intermediate_results_dir /output/intermediate_results_dir` is specified. This
directory contains the intermediate output of make_examples and call_variants
steps.

For more information about `output.visual_report.html`, see the
[VCF stats report documentation](deepvariant-vcf-stats-report.md).

## Notes on GPU image

If you are using GPUs, you can pull the GPU version, and make sure you run with
`--gpus 1`. `call_variants` is the only step that uses the GPU, and can only use
one at a time. `make_examples` and `postprocess_variants` do not run on GPU.

For an example to install GPU driver and docker, see [install_nvidia_docker.sh].



```
sudo docker run --gpus 1 \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  ...
```

## Notes on Singularity

### CPU version

```
# Pull the image.
singularity pull docker://google/deepvariant:"${BIN_VERSION}"

# Run DeepVariant.
singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
  docker://google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \ **Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
  --ref="${INPUT_DIR}"/ucsc.hg19.chr20.unittest.fasta \
  --reads="${INPUT_DIR}"/NA12878_S1.chr20.10_10p1mb.bam \
  --regions "chr20:10,000,000-10,010,000" \
  --output_vcf="${OUTPUT_DIR}"/output.vcf.gz \
  --output_gvcf="${OUTPUT_DIR}"/output.g.vcf.gz \
  --intermediate_results_dir "${OUTPUT_DIR}/intermediate_results_dir" \ **Optional.
  --num_shards=1 \ **How many cores the `make_examples` step uses. Change it to the number of CPU cores you have.**
```

### GPU version

```
# Pull the image.
singularity pull docker://google/deepvariant:"${BIN_VERSION}-gpu"

# Run DeepVariant.
# Using "--nv" and "${BIN_VERSION}-gpu" is important.
singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ \
  docker://google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  ...
```

## Evaluating the results

Here we use the `hap.py`
([https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py))
program from Illumina to evaluate the resulting 10 kilobase vcf file. This
serves as a quick check to ensure the three DeepVariant commands ran correctly.

```bash
sudo docker pull jmcdani20/hap.py:v0.3.12
sudo docker run -it \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /input/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz \
  /output/output.vcf.gz \
  -f "/input/test_nist.b37_chr20_100kbp_at_10mb.bed" \
  -r "/input/ucsc.hg19.chr20.unittest.fasta" \
  -o "/output/happy.output" \
  --engine=vcfeval \
  --pass-only \
  -l chr20:10000000-10010000
```

You should see output similar to the following.

```
Benchmarking Summary:
Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  FP.al  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
INDEL    ALL            4         4         0           13         0          9      0      0            1.0               1.0        0.692308              1.0                     NaN                     NaN                   0.333333                   1.000000
INDEL   PASS            4         4         0           13         0          9      0      0            1.0               1.0        0.692308              1.0                     NaN                     NaN                   0.333333                   1.000000
  SNP    ALL           44        44         0           60         0         16      0      0            1.0               1.0        0.266667              1.0                     1.2                1.307692                   0.333333                   0.363636
  SNP   PASS           44        44         0           60         0         16      0      0            1.0               1.0        0.266667              1.0                     1.2                1.307692                   0.333333                   0.363636
```

[BAM]: http://genome.sph.umich.edu/wiki/BAM
[BWA]: https://academic.oup.com/bioinformatics/article/25/14/1754/225615/Fast-and-accurate-short-read-alignment-with
[docker build]: https://docs.docker.com/engine/reference/commandline/build/
[Dockerfile]: https://github.com/google/deepvariant/blob/r1.5/Dockerfile
[FASTA]: https://en.wikipedia.org/wiki/FASTA_format
[Quick Start in r0.7]: https://github.com/google/deepvariant/blob/r0.7/docs/deepvariant-quick-start.md
[VCF]: https://samtools.github.io/hts-specs/VCFv4.3.pdf
[run_deepvariant.py]: ../scripts/run_deepvariant.py
[install_nvidia_docker.sh]: ../scripts/install_nvidia_docker.sh
