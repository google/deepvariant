# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | -----------------
make_examples                    | 358m41.88s
call_variants: HG002             | 371m20.17s
call_variants: HG003             | 375m57.41s
call_variants: HG004             | 374m18.63s
postprocess_variants (parallel)  | 47m21.31s; 49m4.97s; 49m24.06s
vcf_stats_report(optional):HG002 | 9m32.26s
vcf_stats_report(optional):HG003 | 9m43.99s
vcf_stats_report(optional):HG003 | 10m1.37s
total                            | 1540m27.85s (25h40m27.85s)

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
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WGS/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WGS/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WGS/HG004.output.visual_report.html)

## PacBio (HiFi)

In v1.7.0, we introduced read haplotagging in DeepTrio PacBio. You no longer
need to run DeepVariant->WhatsHap->DeepTrio, and can just run DeepTrio once.

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | 49m35.25s+605m53.88s
call_variants: HG002             | 375m17.22s
call_variants: HG003             | 375m29.55s
call_variants: HG004             | 382m30.66s
postprocess_variants (parallel)  | 59m38.61s; 67m56.74s; 69m15.34s
vcf_stats_report(optional):HG002 | 9m53.77s
vcf_stats_report(optional):HG003 | 10m22.43s
vcf_stats_report(optional):HG003 | 10m32.42s
total                            | 1871m44.75s (31h11m44.75s)

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/PACBIO/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/PACBIO/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/PACBIO/HG004.output.visual_report.html)

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
make_examples                    | 14m30.21s
call_variants: HG002             | 5m29.36s
call_variants: HG003             | 5m33.41s
call_variants: HG004             | 5m32.59s
postprocess_variants (parallel)  | 0m53.02s; 0m53.46s; 0m53.64s
vcf_stats_report(optional):HG002 | 0m7.96s
vcf_stats_report(optional):HG003 | 0m8.04s
vcf_stats_report(optional):HG003 | 0m8.15s
total                            | 32m24.97s

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
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WES/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WES/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WES/HG004.output.visual_report.html)

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
For bigger datasets (WGS and PACBIO), we used bigger disk size (900G).
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.7/scripts/inference_deeptrio.sh

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

[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
