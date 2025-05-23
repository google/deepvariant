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
make_examples                    | 165m53.25s
call_variants: HG002             | 20m35.09s
call_variants: HG003             | 22m9.17s
call_variants: HG004             | 21m49.53s
postprocess_variants (parallel)  | 8m12.84s; 8m48.05s; 8m50.40s
vcf_stats_report(optional):HG002 | 6m23.53s
vcf_stats_report(optional):HG003 | 6m26.04s
vcf_stats_report(optional):HG003 | 6m43.78s
total                            | 251m13.52s (4h11m13.52s)

### Accuracy

We report hap.py results on HG002/HG003/HG004 trio (chr20, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 11208    | 48       | 13       | 0.995736      | 0.998885         | 0.997308        |
| SNP   | 71087    | 246      | 37       | 0.996551      | 0.99948          | 0.998014        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10576    | 52       | 23       | 0.995107      | 0.997919         | 0.996511        |
| SNP   | 69976    | 190      | 55       | 0.997292      | 0.999215         | 0.998253        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10950    | 50       | 24       | 0.995455      | 0.997906         | 0.996679        |
| SNP   | 71445    | 214      | 49       | 0.997014      | 0.999315         | 0.998163        |

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/WGS/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/WGS/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/WGS/HG004.output.visual_report.html)

## PacBio (HiFi)

Read haplotagging in DeepTrio PacBio is on by default. You no longer
need to run DeepVariant->WhatsHap->DeepTrio, and can just run DeepTrio once.

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | 17m8.34s+189m10.68s
call_variants: HG002             | 30m23.02s
call_variants: HG003             | 37m13.83s
call_variants: HG004             | 36m59.66s
postprocess_variants (parallel)  | 7m59.91s; 8m23.51s; 8m30.73s
vcf_stats_report(optional):HG002 | 6m52.84s
vcf_stats_report(optional):HG003 | 6m58.93s
vcf_stats_report(optional):HG003 | 7m6.09s
total                            | 331m0.71s (5h31m0.71s)

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/PACBIO/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/PACBIO/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/PACBIO/HG004.output.visual_report.html)

### Accuracy

We report hap.py results on HG002/HG003/HG004 trio (chr20, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 11214    | 42       | 78       | 0.996269      | 0.993369         | 0.994817        |
| SNP   | 71310    | 23       | 18       | 0.999678      | 0.999748         | 0.999713        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10575    | 53       | 71       | 0.995013      | 0.993623         | 0.994317        |
| SNP   | 70149    | 17       | 32       | 0.999758      | 0.999544         | 0.999651        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10952    | 48       | 61       | 0.995636      | 0.994709         | 0.995172        |
| SNP   | 71622    | 37       | 18       | 0.999484      | 0.999749         | 0.999616        |

## Whole Exome Sequencing (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | --------------
make_examples                    | 7m14.00s
call_variants: HG002             | 2m15.50s
call_variants: HG003             | 2m16.04s
call_variants: HG004             | 2m17.15s
postprocess_variants (parallel)  | 0m9.25s; 0m9.44s; 0m9.57s
vcf_stats_report(optional):HG002 | 0m5.60s
vcf_stats_report(optional):HG003 | 0m5.69s
vcf_stats_report(optional):HG003 | 0m5.84s
total                            | 14m23.30s

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
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/WES/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/WES/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.9.0/WES/HG004.output.visual_report.html)

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
For bigger datasets (WGS and PACBIO), we used bigger disk size (900G).
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.9/scripts/inference_deeptrio.sh

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
