# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

## Setup

The runtime and accuracy reported in this page are generated using
`n2-standard-96` GCP instances which has the following configuration:

```bash
GCP instance type: n2-standard-96
CPUs: 96-core (vCPU)
Memory: 384GiB
GPUs: 0
```

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | -----------------
make_examples                    | 172m53.87s
call_variants: HG002             | 269m26.55s
call_variants: HG003             | 268m2.29s
call_variants: HG004             | 270m22.72s
postprocess_variants (parallel)  | 34m12.36s; 35m4.75s; 35m8.14s
vcf_stats_report(optional):HG002 | 6m36.58s
vcf_stats_report(optional):HG003 | 6m39.92s
vcf_stats_report(optional):HG003 | 6m40.64s
total                            | 1028m3.08s (17h08m3.08s)

### Accuracy

We report hap.py results on HG002/HG003/HG004 trio (chr20, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 11208    | 48       | 13       | 0.995736      | 0.998884         | 0.997308        |
| SNP   | 71088    | 245      | 41       | 0.996565      | 0.999424         | 0.997993        |


#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10578    | 50       | 24       | 0.995295      | 0.99783          | 0.996561        |
| SNP   | 69977    | 189      | 64       | 0.997306      | 0.999087         | 0.998196        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10949    | 51       | 23       | 0.995364      | 0.997993         | 0.996676        |
| SNP   | 71445    | 214      | 48       | 0.997014      | 0.999329         | 0.99817         |

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/WGS/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/WGS/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/WGS/HG004.output.visual_report.html)

## PacBio (HiFi)

Read haplotagging in DeepTrio PacBio is on by default. You no longer
need to run DeepVariant->WhatsHap->DeepTrio, and can just run DeepTrio once.

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | 16m48.88s+288m15.08s
call_variants: HG002             | 279m5.76s
call_variants: HG003             | 274m47.90s
call_variants: HG004             | 283m37.89s
postprocess_variants (parallel)  | 44m12.28s; 51m39.02s; 51m52.66s
vcf_stats_report(optional):HG002 | 6m49.94s
vcf_stats_report(optional):HG003 | 6m53.24s
vcf_stats_report(optional):HG003 | 7m19.57s
total                            | 1206m35.85s (20h6m35.85s)

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/PACBIO/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/PACBIO/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/PACBIO/HG004.output.visual_report.html)

### Accuracy

We report hap.py results on HG002/HG003/HG004 trio (chr20, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 11213    | 43       | 84       | 0.99618       | 0.992863         | 0.994519        |
| SNP   | 71305    | 28       | 21       | 0.999607      | 0.999706         | 0.999657        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10577    | 51       | 77       | 0.995201      | 0.993089         | 0.994144        |
| SNP   | 70143    | 23       | 35       | 0.999672      | 0.999502         | 0.999587        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10954    | 46       | 70       | 0.995818      | 0.993931         | 0.994874        |
| SNP   | 71617    | 42       | 22       | 0.999414      | 0.999693         | 0.999554        |

## Whole Exome Sequencing (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | --------------
make_examples                    | 7m11.47s
call_variants: HG002             | 3m49.25s
call_variants: HG003             | 3m53.32s
call_variants: HG004             | 3m52.68s
postprocess_variants (parallel)  | 0m40.52s; 0m42.09s; 0m42.30s
vcf_stats_report(optional):HG002 | 0m5.65s
vcf_stats_report(optional):HG003 | 0m5.69s
vcf_stats_report(optional):HG003 | 0m7.15s
total                            | 20m6.26s

### Accuracy

We report hap.py results on HG002/HG003/HG004 trio (chr20, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 34       | 0        | 0        | 1.0           | 1.0              | 1.0             |
| SNP   | 670      | 2        | 0        | 0.997024      | 1.0              | 0.99851         |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 29       | 0        | 0        | 1.0           | 1.0              | 1.0             |
| SNP   | 683      | 2        | 1        | 0.99708       | 0.998538         | 0.997809        |


#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 32       | 1        | 1        | 0.969697      | 0.969697         | 0.969697        |
| SNP   | 676      | 3        | 0        | 0.995582      | 1.0              | 0.997786        |

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/WES/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/WES/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.8.0/WES/HG004.output.visual_report.html)

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
For bigger datasets (WGS and PACBIO), we used bigger disk size (900G).
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.8/scripts/inference_deeptrio.sh

# WGS
bash inference_deeptrio.sh --model_preset WGS

# WES
bash inference_deeptrio.sh --model_preset WES

# PacBio
bash inference_deeptrio.sh --model_preset PACBIO

```

Runtime metrics are taken from the resulting log after each stage of
DeepTrio. The runtime numbers reported above are the average of 5 runs each.
The accuracy metrics come from the hap.py summary.csv output file.
The runs are deterministic so all 5 runs produced the same output.

[CPU instance with 96 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
