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
make_examples                    |  45m6.46s
call_variants                    |  16m25.54s
postprocess_variants (with gVCF) |  6m46.32s
vcf_stats_report (optional)      |  5m14.59s (optional)
total                            |  82m5.75s (1h22m5.75s)

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
make_examples                    | 2m58.54s
call_variants                    | 0m33.96s
postprocess_variants (with gVCF) | 0m39.10s
vcf_stats_report (optional)      | 0m6.71s (optional)
total                            | 5m1.55s

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
make_examples                    | 36m25.99s
call_variants                    | 11m27.85s
postprocess_variants (with gVCF) | 4m44.98s
vcf_stats_report (optional)      | 5m24.03s (optional)
total                            | 67m39.74s (1h07m39.74s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501410   | 3091     | 2969     | 0.993873      | 0.994327         | 0.9941          |
| SNP   | 3321461  | 6034     | 3855     | 0.998187      | 0.998842         | 0.998514        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.9.0/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | 47m44.53s
call_variants                    | 58m42.49s
postprocess_variants (with gVCF) | 11m11.15s
vcf_stats_report (optional)      | 7m30.22s (optional)
total                            | 133m33.89s (2h13m33.89s)

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
make_examples                    | 63m29.95s
call_variants                    | 65m14.15s
postprocess_variants (with gVCF) | 4m3.75s
vcf_stats_report (optional)      | 5m16.99s (optional)
total                            | 168m42.48s (2h48m42.48s)

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
