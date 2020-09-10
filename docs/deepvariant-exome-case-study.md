# DeepVariant whole exome sequencing (WES) case study

Similar to the [case study on whole genome sequencing data], in this
study we describe applying DeepVariant to a real exome sample using a single
machine.

## Prepare environment

### Tools

[Docker](https://docs.docker.com/get-docker/) will be used to run DeepVariant
and [hap.py](https://github.com/illumina/hap.py),

### Download Reference

We will be using GRCh37 for this case study.

```bash
mkdir -p reference

HTTPDIR=https://storage.googleapis.com/deepvariant/exome-case-study-testdata

curl ${HTTPDIR}/hs37d5.fa.gz | gunzip > reference/hs37d5.fa
curl ${HTTPDIR}/hs37d5.fa.fai > reference/hs37d5.fa.fai
```

### Download Genome in a Bottle Benchmarks

We will benchmark our variant calls against v4.1 of the Genome in a Bottle small
variant benchmarks for HG002.

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37

curl ${FTPDIR}/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed > benchmark/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed
curl ${FTPDIR}/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz > benchmark/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi
```

### Download HG002 BAM

```bash
mkdir -p input

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome
BAM=151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup

curl ${FTPDIR}/${BAM}.bam > input/${BAM}.bam
curl ${FTPDIR}/${BAM}.bai > input/${BAM}.bai
```

### Download capture target BED file

According to the paper "[Extensive sequencing of seven human genomes to
characterize benchmark reference
materials](https://www.nature.com/articles/sdata201625)", the HG002 exome was
generated with Agilent SureSelect. In this case study we'll use the SureSelect
v5 BED (`agilent_sureselect_human_all_exon_v5_b37_targets.bed`). For evaluation,
`hap.py` will intersect this BED with the GIAB confident regions.

```bash
HTTPDIR=https://storage.googleapis.com/deepvariant/exome-case-study-testdata

curl ${HTTPDIR}/agilent_sureselect_human_all_exon_v5_b37_targets.bed > input/agilent_sureselect_human_all_exon_v5_b37_targets.bed
```


## Running on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION=1.0.0

sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref /reference/hs37d5.fa \
  --reads /input/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam \
  --regions /input/agilent_sureselect_human_all_exon_v5_b37_targets.bed \
  --output_vcf /output/HG002.output.vcf.gz \
  --output_gvcf /output/HG002.output.g.vcf.gz \
  --num_shards $(nproc) \
  --intermediate_results_dir /output/intermediate_results_dir
```

By specifying `--model_type WES`, you'll be using a model that is best suited
for Illumina Whole Exome Sequencing data.

`--intermediate_results_dir` flag is optional. By specifying it, the
intermediate outputs of `make_examples` and `call_variants` stages can be found
in the directory. After the command, you can find these files in the directory:

```
call_variants_output.tfrecord.gz
gvcf.tfrecord-?????-of-?????.gz
make_examples.tfrecord-?????-of-?????.gz
```

For running on GPU machines, or using Singularity instead of Docker, see
[Quick Start](deepvariant-quick-start.md).

## Benchmark on all chromosomes

```bash
mkdir -p happy

sudo docker pull pkrusche/hap.py

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
  /benchmark/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz \
  /output/HG002.output.vcf.gz \
  -f /benchmark/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed \
  -T /input/agilent_sureselect_human_all_exon_v5_b37_targets.bed \
  -r /reference/hs37d5.fa \
  -o /happy/happy.output \
  --engine=vcfeval
```

Output:

```
Benchmarking Summary:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL         3020      2896       124         3718        81        703     35       0.958940          0.973134        0.189080         0.965985                     NaN                     NaN                   1.325453                   1.368421
 INDEL   PASS         3020      2896       124         3718        81        703     35       0.958940          0.973134        0.189080         0.965985                     NaN                     NaN                   1.325453                   1.368421
   SNP    ALL        38576     38180       396        41157       116       2816     27       0.989735          0.996975        0.068421         0.993341                2.624273                2.568779                   1.548701                   1.563630
   SNP   PASS        38576     38180       396        41157       116       2816     27       0.989735          0.996975        0.068421         0.993341                2.624273                2.568779                   1.548701                   1.563630
```

[case study on whole genome sequencing data]: deepvariant-case-study.md
