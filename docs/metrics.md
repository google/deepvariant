# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | ~98m
call_variants                    | ~207m
postprocess_variants (with gVCF) | ~27m
total                            | ~332m = ~5.53 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501680   | 2821     | 1276     | 0.994408      | 0.997565         | 0.995984        |
| SNP   | 3306792  | 20704    | 4270     | 0.993778      | 0.998711         | 0.996238        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.0/WGS/deepvariant.output.visual_report.html)

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
| INDEL | 1021     | 30       | 13       | 0.971456      | 0.987689         | 0.979505        |
| SNP   | 24981    | 298      | 62       | 0.988212      | 0.997524         | 0.992846        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.0/WES/deepvariant.output.visual_report.html)

## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | ~132m
call_variants                    | ~242m
postprocess_variants (with gVCF) | ~34m
total                            | ~408m = ~6.8 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501502   | 2999     | 2781     | 0.994056      | 0.994704         | 0.99438         |
| SNP   | 3324311  | 3184     | 1497     | 0.999043      | 0.99955          | 0.999297        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.0/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | ~251m
call_variants                    | ~295m
postprocess_variants (with gVCF) | ~35m
total                            | ~581m = ~9.68 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 438790   | 65711    | 40704    | 0.869751      | 0.918025         | 0.893236        |
| SNP   | 3314049  | 13446    | 8147     | 0.995959      | 0.997549         | 0.996753        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.0/ONT_R104/deepvariant.output.visual_report.html)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | ~147m
call_variants                    | ~209m
postprocess_variants (with gVCF) | ~22m
total                            | ~378m = ~6.3 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503007   | 1494     | 2765     | 0.997039      | 0.994785         | 0.995911        |
| SNP   | 3323624  | 3871     | 2277     | 0.998837      | 0.999316         | 0.999076        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.6.0/HYBRID/deepvariant.output.visual_report.html)

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
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.6/scripts/inference_deepvariant.sh

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

