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
make_examples                    |  47m4.92s
call_variants                    |  15m56.52s
postprocess_variants (with gVCF) |  7m0.99s
vcf_stats_report (optional)      |  5m17.67s (optional)
total                            |  83m57.12s (1h23m57.12s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501594   | 2907     | 1190     | 0.994238      | 0.997729         | 0.99598         |
| SNP   | 3306720  | 20776    | 4880     | 0.993756      | 0.998527         | 0.996136        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/WGS/deepvariant.output.visual_report.html)

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 3m0.33s
call_variants                    | 0m33.72s
postprocess_variants (with gVCF) | 0m39.24s
vcf_stats_report (optional)      | 0m5.10s (optional)
total                            | 5m7.71s

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
make_examples                    | 33m46.75s
call_variants                    | 11m38.86s
postprocess_variants (with gVCF) | 5m12.45s
vcf_stats_report (optional)      | 5m34.81s (optional)
total                            | 65m27.90s (1h05m27.90s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501455   | 3046     | 2986     | 0.993962      | 0.994296         | 0.994129        |
| SNP   | 3321751  | 5744     | 4032     | 0.998274      | 0.998789         | 0.998532        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | 46m29.14s
call_variants                    | 53m48.26s
postprocess_variants (with gVCF) | 11m25.74s
vcf_stats_report (optional)      | 7m22.90s (optional)
total                            | 127m34.97s (2h07m34.97s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 461818   | 42683    | 31344    | 0.915396      | 0.938385         | 0.926748        |
| SNP   | 3321289  | 6206     | 5476     | 0.998135      | 0.998355         | 0.998245        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/ONT_R104/deepvariant.output.visual_report.html)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 60m4.06s
call_variants                    | 62m23.86s
postprocess_variants (with gVCF) | 4m10.56s
vcf_stats_report (optional)      | 5m16.31s (optional)
total                            | 162m45.17s (2h42m45.17s)

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503264   | 1237     | 2052     | 0.997548      | 0.996129         | 0.996838        |
| SNP   | 3324021  | 3474     | 1856     | 0.998956      | 0.999442         | 0.999199        |

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
