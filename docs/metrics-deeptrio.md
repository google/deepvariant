# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -----------------
make_examples                    | ~379m
call_variants for HG002          | ~385m
call_variants for HG003          | ~380m
call_variants for HG004          | ~377m
postprocess_variants (parallel)  | ~50m
total                            | ~1571m = ~26.18 hours

### Accuracy

We report hap.py results on HG002/HG003/HG004 trio (chr20, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 11207    | 49       | 14       | 0.995647      | 0.998798         | 0.99722         |
| SNP   | 71088    | 245      | 41       | 0.996565      | 0.999424         | 0.997993        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10578    | 50       | 24       | 0.995295      | 0.99783          | 0.996561        |
| SNP   | 69975    | 191      | 65       | 0.997278      | 0.999072         | 0.998174        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10948    | 52       | 24       | 0.995273      | 0.997905         | 0.996587        |
| SNP   | 71447    | 212      | 50       | 0.997042      | 0.999301         | 0.99817         |

* See VCF stats report (for all chromosomes)
  - [HG002](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WGS/HG002.output.visual_report.html)
  - [HG003](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WGS/HG003.output.visual_report.html)
  - [HG004](https://storage.googleapis.com/deepvariant/visual_reports/DeepTrio/1.7.0/WGS/HG004.output.visual_report.html)

## PacBio (HiFi)

In v1.7.0, we introduced read haplotagging in DeepTrio PacBio. You no longer
need to run DeepVariant->WhatsHap->DeepTrio, and can just run DeepTrio once.

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | ~46+824m
call_variants for HG002          | ~373m
call_variants for HG003          | ~372m
call_variants for HG004          | ~376m
postprocess_variants (parallel)  | ~66m
total                            | ~2057m = ~34.28 hours

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
| INDEL | 11214    | 42       | 84       | 0.996269      | 0.992864         | 0.994563        |
| SNP   | 71306    | 27       | 18       | 0.999621      | 0.999748         | 0.999685        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10581    | 47       | 76       | 0.995578      | 0.993183         | 0.994379        |
| SNP   | 70143    | 23       | 35       | 0.999672      | 0.999502         | 0.999587        |
#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 10953    | 47       | 76       | 0.995727      | 0.993411         | 0.994568        |
| SNP   | 71615    | 44       | 25       | 0.999386      | 0.999651         | 0.999519        |

## Whole Exome Sequencing (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | --------------
make_examples                    | ~15m
call_variants for HG002          | ~5m
call_variants for HG003          | ~5m
call_variants for HG004          | ~5m
postprocess_variants (parallel)  | ~1m
total                            | ~31m

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
