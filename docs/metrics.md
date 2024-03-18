# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | ~103m
call_variants                    | ~196m
postprocess_variants (with gVCF) | ~27m
total                            | ~326m = ~5.43 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501683   | 2818     | 1265     | 0.994414      | 0.997586         | 0.995998        |
| SNP   | 3306788  | 20708    | 4274     | 0.993777      | 0.99871          | 0.996237        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.1/WGS/deepvariant.output.visual_report.html)

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | ~7m
call_variants                    | ~1m
postprocess_variants (with gVCF) | ~1m
total                            | ~9m

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1022     | 29       | 13       | 0.972407      | 0.987713         | 0.98            |
| SNP   | 24987    | 292      | 59       | 0.988449      | 0.997645         | 0.993025        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.1/WES/deepvariant.output.visual_report.html)

## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | ~149m
call_variants                    | ~217m
postprocess_variants (with gVCF) | ~33m
total                            | ~399m = ~6.65 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501516   | 2985     | 2745     | 0.994083      | 0.994773         | 0.994428        |
| SNP   | 3324302  | 3193     | 1502     | 0.99904       | 0.999549         | 0.999295        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.1/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | ~329m
call_variants                    | ~281m
postprocess_variants (with gVCF) | ~34m
total                            | ~644m = ~10.73 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 441658   | 62843    | 41301    | 0.875435      | 0.917411         | 0.895932        |
| SNP   | 3314131  | 13364    | 8115     | 0.995984      | 0.997558         | 0.99677         |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.1/ONT_R104/deepvariant.output.visual_report.html)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | ~172m
call_variants                    | ~211m
postprocess_variants (with gVCF) | ~24m
total                            | ~407m = ~6.78 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503014   | 1487     | 2767     | 0.997053      | 0.994781         | 0.995916        |
| SNP   | 3323624  | 3871     | 2273     | 0.998837      | 0.999317         | 0.999077        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.1/HYBRID/deepvariant.output.visual_report.html)

## Inspect outputs that produced the metrics above

The DeepVariant VCFs, gVCFs, and hap.py evaluation outputs are available at:

```
gs://deepvariant/case-study-outputs
```

You can also inspect them in a web browser here:
https://42basepairs.com/browse/gs/deepvariant/case-study-outputs

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.6.1/scripts/inference_deepvariant.sh

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

[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform

