# Using DeepTrio for small variant calling from the trio sequenced with PacBio HiFi

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

We will benchmark our variant calls against v4.2 of the Genome in a Bottle small
variant benchmarks for HG002, HG003, and HG004 trio.

```bash
mkdir -p benchmark

FTPDIR=ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020

curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi

curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG003_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi

curl ${FTPDIR}/HG004_GRCh38_1_22_v4.2_benchmark.bed > benchmark/HG004_GRCh38_1_22_v4.2_benchmark.bed
curl ${FTPDIR}/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz > benchmark/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz
curl ${FTPDIR}/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi > benchmark/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi
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
`postprocess_variants` and `GLNexus merge`. It is possible to run the first
three steps with one command using the `run_deeptrio` script. GLNexus
is run as a separate command.

### Running on a CPU-only machine

```bash
mkdir -p output
mkdir -p output/intermediate_results_dir

BIN_VERSION=1.0.1rc

sudo docker pull gcr.io/deepvariant-docker/deeptrio:"${BIN_VERSION}"

time sudo docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  gcr.io/deepvariant-docker/deeptrio:"${BIN_VERSION}" \
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
  --make_examples_extra_args "sort_by_haplotypes=true,parse_sam_aux_fields=true"
```

Note that there are extra arguments added: "sort_by_haplotypes" and
"parse_sam_aux_fields". These extra arguments make use of a phased reads, thus
allowing a further improvement of the accuracy. In order to use this feature
input BAM files have to be phased. For the detailed description on how to do
that please see
[DeepVariant PacBio case study](deepvariant-pacbio-model-case-study.md).

By specifying `--model_type PACBIO`, you'll be using a model that is best suited
for PacBio HiFi Whole Genome Sequencing data.

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

## Merge VCFs using GLNexus

At this step we take all 3 VCFs generated in the previous step and merge them
using GLNexus.

```bash
# BCFTools are required:
sudo apt-get -y install bcftools
sudo apt-get -y install tabix

sudo docker pull quay.io/mlin/glnexus:v1.2.7

sudo docker run \
  -v "${PWD}/output":"/output" \
  quay.io/mlin/glnexus:v1.2.7 \
  /usr/local/bin/glnexus_cli \
  --config DeepVariantWGS \
  /output/HG002.g.vcf.gz \
  /output/HG003.g.vcf.gz \
  /output/HG004.g.vcf.gz \
  | bcftools view - | bgzip -c > output/HG002_trio_merged.vcf.gz
```

After completion of GLNexus command we should have a new merged VCF file in the
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
45 non-pass records were skipped
Concordance HG002: F:143838/144142 (99.79%)  M:143867/144177 (99.78%)  F+M:142108/142800 (99.52%)
0/148910 (0.00%) records did not conform to expected call ploidy
146501/148910 (98.38%) records were variant in at least 1 family member and checked for Mendelian constraints
3317/146501 (2.26%) records had indeterminate consistency status due to incomplete calls
748/146501 (0.51%) records contained a violation of Mendelian constraints
```

### Benchmark variant calls against 4.2 truth set with hap.py

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
  /benchmark/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz \
  /output/HG002.output.vcf.gz \
  -f /benchmark/HG002_GRCh38_1_22_v4.2_benchmark.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval \
  -l chr20

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
  /benchmark/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz \
  /output/HG003.output.vcf.gz \
  -f /benchmark/HG003_GRCh38_1_22_v4.2_benchmark.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval \
  -l chr20

sudo docker run \
  -v "${PWD}/benchmark":"/benchmark" \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/happy:/happy" \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
  /benchmark/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz \
  /output/HG004.output.vcf.gz \
  -f /benchmark/HG004_GRCh38_1_22_v4.2_benchmark.bed \
  -r /reference/GRCh38_no_alt_analysis_set.fasta \
  -o /happy/happy.output \
  --engine=vcfeval \
  -l chr20
```

```
Benchmarking Summary for HG002:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL        11256     11214        42        22325        73      10571     24       0.996269          0.993789        0.473505         0.995027                     NaN                     NaN                   1.561710                   2.282047
 INDEL   PASS        11256     11214        42        22325        73      10571     24       0.996269          0.993789        0.473505         0.995027                     NaN                     NaN                   1.561710                   2.282047
   SNP    ALL        71333     71273        60        94297        20      22929     16       0.999159          0.999720        0.243157         0.999439                2.314904                2.023814                   1.715978                   2.001115
   SNP   PASS        71333     71273        60        94297        20      22929     16       0.999159          0.999720        0.243157         0.999439                2.314904                2.023814                   1.715978                   2.001115

Benchmarking Summary for HG003:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL        10634     10584        50        22122        72      10974     34       0.995298          0.993541        0.496067         0.994419                     NaN                     NaN                   1.749861                   2.498041
 INDEL   PASS        10634     10584        50        22122        72      10974     34       0.995298          0.993541        0.496067         0.994419                     NaN                     NaN                   1.749861                   2.498041
   SNP    ALL        70209     70186        23        94386        29      24118     15       0.999672          0.999587        0.255525         0.999630                2.297347                1.977766                   1.884533                   2.173085
   SNP   PASS        70209     70186        23        94386        29      24118     15       0.999672          0.999587        0.255525         0.999630                2.297347                1.977766                   1.884533                   2.173085

Benchmarking Summary for HG004:
  Type Filter  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  QUERY.UNK  FP.gt  METRIC.Recall  METRIC.Precision  METRIC.Frac_NA  METRIC.F1_Score  TRUTH.TOTAL.TiTv_ratio  QUERY.TOTAL.TiTv_ratio  TRUTH.TOTAL.het_hom_ratio  QUERY.TOTAL.het_hom_ratio
 INDEL    ALL        11036     10984        52        22518        76      10950     31       0.995288          0.993430        0.486278         0.994358                     NaN                     NaN                   1.791542                   2.507959
 INDEL   PASS        11036     10984        52        22518        76      10950     31       0.995288          0.993430        0.486278         0.994358                     NaN                     NaN                   1.791542                   2.507959
   SNP    ALL        71933     71861        72        95483        32      23506     11       0.998999          0.999555        0.246180         0.999277                2.309582                1.995672                   1.878938                   2.019747
   SNP   PASS        71933     71861        72        95483        32      23506     11       0.998999          0.999555        0.246180         0.999277                2.309582                1.995672                   1.878938                   2.019747
```
