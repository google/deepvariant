# Runtime and accuracy metrics for all release models

## Setup

The runtime and accuracy reported in this page are generated using
`n2-standard-96` GCP instances which has the following configuration:

```bash
GCP instance type: n2-standard-96
CPUs: 96-core (vCPU)
Memory: 384GiB
GPUs: 0
```

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    |  45m55.29s
call_variants                    |  16m34.10s
postprocess_variants (with gVCF) |  6m49.36s
vcf_stats_report (optional)      |  5m30.57s (optional)
total                            |  83m6.88s (1h23m6.88s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501527   | 2974     | 1262     | 0.994105      | 0.997591         | 0.995845        |
| SNP   | 3306720  | 20776    | 4900     | 0.993756      | 0.998521         | 0.996133        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/WGS/deepvariant.output.visual_report.html)

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 2m58.39s
call_variants                    | 0m34.54s
postprocess_variants (with gVCF) | 0m39.29s
vcf_stats_report (optional)      | 0m5.06s (optional)
total                            | 5m3.19s

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1024     | 27       | 8        | 0.97431       | 0.992417         | 0.98328         |
| SNP   | 24983    | 296      | 60       | 0.988291      | 0.997604         | 0.992926        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/WES/deepvariant.output.visual_report.html)

## PacBio (HiFi)

### Updated dataset

We have updated the PacBio test data from HG003 Sequel-II to
latest Revio with SPRQ chemistry data to showcase performance on the updated
platform and chemistry. The numbers reported here are generated using the bam
that can be found in:

```bash
gs://deepvariant/pacbio-case-study-testdata/HG003.SPRQ.pacbio.GRCh38.nov2024.bam
```

Which is also available through [here](https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/GIAB_trio/HG003/analysis/GRCh38.m84039_241002_000337_s3.hifi_reads.bc2020.bam).

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | 31m41.58s
call_variants                    | 12m42.97s
postprocess_variants (with gVCF) | 4m58.54s
vcf_stats_report (optional)      | 5m46.53s (optional)
total                            | 64m3.68s (1h04m3.68s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501343   | 3158     | 3167     | 0.99374       | 0.993951         | 0.993846        |
| SNP   | 3321731  | 5764     | 4061     | 0.998268      | 0.998781         | 0.998524        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | 52m15.03s
call_variants                    | 16m42.26s
postprocess_variants (with gVCF) | 5m49.54s
vcf_stats_report (optional)      | 6m22.82s (optional)
total                            | 89m32.02s (1h29m32.02s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 450841   | 53660    | 40739    | 0.893637      | 0.91953          | 0.906399        |
| SNP   | 3319365  | 8130     | 2964     | 0.997557      | 0.999108         | 0.998332        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/ONT_R104/deepvariant.output.visual_report.html)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 84m23.11s
call_variants                    | 103m22.48s
postprocess_variants (with gVCF) | 5m31.07s
vcf_stats_report (optional)      | 8m20.23s (optional)
total                            | 209m0.92s (3h29m0.92s)

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503160   | 1341     | 2243     | 0.997342      | 0.99577          | 0.996555        |
| SNP   | 3323907  | 3588     | 1981     | 0.998922      | 0.999405         | 0.999163        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/HYBRID/deepvariant.output.visual_report.html)

## Inspect outputs that produced the metrics above

The DeepVariant VCFs, gVCFs, and hap.py evaluation outputs are available at:

```
gs://deepvariant/case-study-outputs
```

You can also inspect them in a web browser here:
https://42basepairs.com/browse/gs/deepvariant/case-study-outputs

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.9/scripts/inference_deepvariant.sh

# WGS
bash inference_deepvariant.sh --model_preset WGS

# WES
bash inference_deepvariant.sh --model_preset WES

# PacBio
bash inference_deepvariant.sh --model_preset PACBIO

# ONT_R104
bash inference_deepvariant.sh --model_preset ONT_R104

# Hybrid
bash inference_deepvariant.sh --model_preset HYBRID_PACBIO_ILLUMINA
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant. The runtime numbers reported above are the average of 5 runs each.
The accuracy metrics come from the hap.py summary.csv output file.
The runs are deterministic so all 5 runs produced the same output.
